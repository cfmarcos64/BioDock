# -*- coding: utf-8 -*-
"""
BioDock Studio: Protein Parser, Inspector & Docking Prep
v9 — sesión PyMOL (.pse), docking con Vina integrado (paquete Python o ejecutable), RMSD de redocking,
caja con ligando externo, bolsillo y contorno corregidos.
"""

import base64
import hashlib
import html as html_lib
import importlib.util
import io
import json
import math
import os
import re
import shutil
import subprocess
import sys
import tempfile
import time
import urllib.request
import zipfile

import py3Dmol
import streamlit as st
from Bio import PDB
from Bio.SeqUtils import ProtParam
from stmol import showmol

# Librerías opcionales para la preparación de docking
try:
    from rdkit import Chem
    from rdkit.Chem import AllChem
    HAS_RDKIT = True
except Exception:
    HAS_RDKIT = False
try:
    from openbabel import pybel
    HAS_OB = True
except Exception:
    HAS_OB = False
try:
    from meeko import MoleculePreparation, PDBQTWriterLegacy
    HAS_MEEKO = True
except Exception:
    HAS_MEEKO = False
try:
    from meeko import PDBQTMolecule, RDKitMolCreate
    HAS_MEEKO_READ = True
except Exception:
    HAS_MEEKO_READ = False
HAS_PLIP = HAS_OB and importlib.util.find_spec("plip") is not None
HAS_PDB2PQR = all(importlib.util.find_spec(m) is not None for m in ("pdb2pqr", "propka"))

# --- CONFIGURACIÓN DE PÁGINA STREAMLIT ---
st.set_page_config(
    page_title="BioDock Studio | Prep & Inspector",
    page_icon="🧬",
    layout="wide",
)

# --- TEMA NEGRO (CSS) ---
st.markdown(
    """
<style>
:root { --bg:#000000; --panel:#0a0c10; --border:#1e2430; --accent:#19e3c5;
        --accent2:#4d7cff; --text:#e8eaed; --muted:#8b93a1; }

html, body, .stApp, [data-testid="stAppViewContainer"], [data-testid="stHeader"],
[data-testid="stToolbar"], [data-testid="stSidebar"] { background-color: var(--bg) !important; color: var(--text); }

.block-container { padding-top: 2rem; max-width: 1500px; }

h1 { background: linear-gradient(90deg, #19e3c5, #4d7cff, #b06bff);
     -webkit-background-clip: text; -webkit-text-fill-color: transparent; font-weight: 800; }
h2, h3, h4, h5 { color: #ffffff !important; }
p, li, label, span, .stMarkdown, [data-testid="stCaptionContainer"] { color: var(--text); }
[data-testid="stCaptionContainer"] * { color: var(--muted) !important; }
hr { border-color: var(--border) !important; }

/* Contenedores con borde */
[data-testid="stVerticalBlockBorderWrapper"] { border-color: var(--border) !important; border-radius: 14px; }

/* Métricas */
[data-testid="stMetric"] { background: var(--panel); border: 1px solid var(--border);
                           border-radius: 12px; padding: 12px 16px; }
[data-testid="stMetricLabel"] * { color: var(--muted) !important; }
[data-testid="stMetricValue"] { color: var(--accent) !important; }

/* Botones */
.stButton > button, .stDownloadButton > button {
    border-radius: 10px; border: 1px solid #2a3242; background: #0f131b; color: var(--text); }
.stButton > button:hover, .stDownloadButton > button:hover {
    border-color: var(--accent); color: var(--accent); }
button[kind="primary"], button[data-testid="stBaseButton-primary"] {
    background: linear-gradient(90deg, #19e3c5, #4d7cff) !important; color: #000 !important;
    border: none !important; font-weight: 700; }
button[kind="primary"]:hover, button[data-testid="stBaseButton-primary"]:hover { filter: brightness(1.15); color:#000 !important; }

/* Inputs */
div[data-baseweb="select"] > div, div[data-baseweb="input"] > div, textarea,
[data-testid="stFileUploaderDropzone"] { background-color: var(--panel) !important;
    border-color: #222836 !important; color: var(--text) !important; }
div[data-baseweb="popover"] ul, div[data-baseweb="popover"] > div { background-color: var(--panel) !important; }
span[data-baseweb="tag"] { background-color: #0f2f36 !important; color: var(--accent) !important; }

/* Pestañas */
button[data-baseweb="tab"] { color: var(--muted); background: transparent; }
button[data-baseweb="tab"][aria-selected="true"] { color: var(--accent); }
div[data-baseweb="tab-highlight"] { background-color: var(--accent) !important; }
div[data-baseweb="tab-border"] { background-color: var(--border) !important; }

/* Código / secuencia */
[data-testid="stCode"] pre, [data-testid="stCode"] code { background: var(--panel) !important; color: #b8f5e8 !important; }
[data-testid="stAlert"] { background: var(--panel); border: 1px solid var(--border); }
</style>
""",
    unsafe_allow_html=True,
)

AA_3TO1 = {
    'ALA': 'A', 'CYS': 'C', 'ASP': 'D', 'GLU': 'E', 'PHE': 'F',
    'GLY': 'G', 'HIS': 'H', 'ILE': 'I', 'LYS': 'K', 'LEU': 'L',
    'MET': 'M', 'ASN': 'N', 'PRO': 'P', 'GLN': 'Q', 'ARG': 'R',
    'SER': 'S', 'THR': 'T', 'VAL': 'V', 'TRP': 'W', 'TYR': 'Y',
    'MSE': 'M',
}
WATER_NAMES = {"HOH", "WAT", "DOD", "H2O"}
COMMON_IONS = {"NA", "CL", "K", "MG", "ZN", "CA", "FE", "CU", "MN", "SO4", "PO4", "EDTA"}


# --- FUNCIONES AUXILIARES ---

def fetch_pdb_from_web(pdb_code):
    if not pdb_code or len(pdb_code.strip()) != 4:
        return None
    url = f"https://files.rcsb.org/download/{pdb_code.strip().upper()}.pdb"
    try:
        req = urllib.request.Request(url, headers={'User-Agent': 'Mozilla/5.0'})
        with urllib.request.urlopen(req, timeout=20) as response:
            return response.read().decode('utf-8', errors='ignore')
    except Exception:
        return None


def is_water(res):
    return res.get_id()[0] == 'W' or res.get_resname().strip() in WATER_NAMES


def is_amino_acid(res):
    return res.get_resname().strip() in AA_3TO1 and not is_water(res)


def is_ion(res):
    return res.get_resname().strip() in COMMON_IONS and not is_amino_acid(res)


def is_ligand(res):
    return (not is_water(res)) and (not is_ion(res)) and (not is_amino_acid(res)) \
        and res.get_id()[0].startswith('H_')


# --- ESTRUCTURA SECUNDARIA ---

def parse_ss_records(raw_pdb):
    segments = []
    for line in raw_pdb.splitlines():
        try:
            if line.startswith("HELIX"):
                segments.append(("H", line[19], line[15:18].strip(), int(line[21:25]), line[25].strip(),
                                 line[27:30].strip(), int(line[33:37]), line[37].strip()))
            elif line.startswith("SHEET"):
                segments.append(("E", line[21], line[17:20].strip(), int(line[22:26]), line[26].strip(),
                                 line[28:31].strip(), int(line[33:37]), line[37].strip()))
        except (ValueError, IndexError):
            continue
    return segments


def compute_ss_segments(model):
    segments = []
    for chain in model:
        for pp in PDB.PPBuilder().build_peptides(chain):
            labels = []
            for phi, psi in pp.get_phi_psi_list():
                if phi is None or psi is None:
                    labels.append('C')
                    continue
                phi, psi = math.degrees(phi), math.degrees(psi)
                if -160 < phi < -20 and -120 < psi < 50:
                    labels.append('H')
                elif -180 < phi < -40 and (psi > 90 or psi < -150):
                    labels.append('E')
                else:
                    labels.append('C')

            i = 0
            while i < len(labels):
                j = i
                while j + 1 < len(labels) and labels[j + 1] == labels[i]:
                    j += 1
                min_len = 4 if labels[i] == 'H' else 3
                if labels[i] != 'C' and (j - i + 1) >= min_len:
                    r1, r2 = pp[i], pp[j]
                    segments.append((
                        labels[i], chain.id,
                        r1.get_resname().strip(), r1.get_id()[1], r1.get_id()[2].strip(),
                        r2.get_resname().strip(), r2.get_id()[1], r2.get_id()[2].strip(),
                    ))
                i = j + 1
    return segments


def segments_to_pdb_records(segments):
    lines = []
    h_n = s_n = 0
    for kind, ch, n1, s, ic1, n2, e, ic2 in segments:
        if kind == 'H':
            h_n += 1
            lines.append(
                f"HELIX  {h_n:3d} {h_n:3d} {n1:>3s} {ch:1s} {s:4d}{ic1:1s} "
                f"{n2:>3s} {ch:1s} {e:4d}{ic2:1s} 1"
            )
        else:
            s_n += 1
            lines.append(
                f"SHEET  {s_n:3d}   A 1 {n1:>3s} {ch:1s}{s:4d}{ic1:1s} "
                f"{n2:>3s} {ch:1s}{e:4d}{ic2:1s}  0"
            )
    return "\n".join(lines) + ("\n" if lines else "")


def residue_ss(segments, chain_id, resnum):
    for kind, ch, _, s, _, _, e, _ in segments:
        if ch == chain_id and s <= resnum <= e:
            return kind
    return 'C'


# --- PROCESADO Y FILTRADO ---

@st.cache_data(show_spinner=False)
def get_chain_summary(pdb_text):
    """Devuelve {cadena: nº de aminoácidos} para todas las cadenas del primer modelo."""
    structure = PDB.PDBParser(QUIET=True).get_structure("s", io.StringIO(pdb_text))
    model = next(iter(structure))
    return {chain.id: sum(1 for r in chain if is_amino_acid(r)) for chain in model}


@st.cache_data(show_spinner=False)
def process_pdb_structure(pdb_text, selected_chains, remove_water, remove_ions):
    selected = set(selected_chains)
    structure = PDB.PDBParser(QUIET=True).get_structure("protein_structure", io.StringIO(pdb_text))
    model = next(iter(structure))
    model_id = model.id

    segments = parse_ss_records(pdb_text)
    ss_source = "registros del PDB"
    if not segments:
        segments = compute_ss_segments(model)
        ss_source = "calculada (phi/psi)"
    ss_header = segments_to_pdb_records([s for s in segments if s[1] in selected])

    chains_data = {}
    ligands = []
    water_count = ion_count = 0

    for chain in model:
        if chain.id not in selected:
            continue
        residues = []
        for res in chain:
            if is_water(res):
                water_count += 1
            elif is_amino_acid(res):
                num = res.get_id()[1]
                residues.append((AA_3TO1[res.get_resname().strip()], num,
                                 residue_ss(segments, chain.id, num)))
            elif is_ion(res):
                ion_count += 1
            elif is_ligand(res):
                ligands.append((res.get_resname().strip(), res.get_id()[1], chain.id))
        if residues:
            chains_data[chain.id] = residues

    class ProteinSelect(PDB.Select):
        def accept_model(self, m):
            return m.id == model_id

        def accept_chain(self, chain):
            return chain.id in selected

        def accept_residue(self, residue):
            if is_amino_acid(residue):
                return True
            if not remove_water and is_water(residue):
                return True
            if not remove_ions and is_ion(residue):
                return True
            return False

    class LigandSelect(PDB.Select):
        def accept_model(self, m):
            return m.id == model_id

        def accept_chain(self, chain):
            return chain.id in selected

        def accept_residue(self, residue):
            return is_ligand(residue)

    writer = PDB.PDBIO()
    writer.set_structure(structure)

    buf = io.StringIO()
    writer.save(buf, ProteinSelect())
    cleaned_pdb_str = ss_header + buf.getvalue()

    buf = io.StringIO()
    writer.save(buf, LigandSelect())
    ligand_pdb_str = buf.getvalue()

    return {
        "cleaned": cleaned_pdb_str,
        "ligand": ligand_pdb_str,
        "chains": chains_data,
        "water": water_count,
        "ions": ion_count,
        "ligands": ligands,
        "ss_source": ss_source,
    }


def calculate_protein_properties(sequence_str):
    if not sequence_str or len(sequence_str) < 3:
        return None
    try:
        clean_seq = re.sub(r'[^ACDEFGHIKLMNPQRSTVWY]', 'A', sequence_str)
        analyser = ProtParam.ProteinAnalysis(clean_seq)
        return {
            "length": len(clean_seq),
            "mw_kDa": round(analyser.molecular_weight() / 1000.0, 2),
            "pi": round(analyser.isoelectric_point(), 2),
            "aromaticity": round(analyser.aromaticity() * 100, 1),
            "instability": round(analyser.instability_index(), 1),
        }
    except Exception:
        return None


def format_chain_block(residues, width=60):
    out = []
    for i in range(0, len(residues), width):
        part = residues[i:i + width]
        seq = "".join(r[0] for r in part)
        ss = "".join({'H': 'H', 'E': 'E'}.get(r[2], '-') for r in part)
        out.append(f"{part[0][1]:>5}  {seq}  {part[-1][1]}")
        out.append(f"       {ss}")
        out.append("")
    return "\n".join(out)


# --- VISOR 3D ---

STYLE_OPTIONS = [
    "Cartoon (flechas)",
    "Cartoon ovalado",
    "Tubo (trace)",
    "Esqueleto (Backbone)",
    "Palos (Sticks)",
    "Bolas y palos",
    "Líneas (Wireframe)",
    "Esferas (van der Waals)",
    "Superficie + Cartoon",
]

COLOR_OPTIONS = [
    "Arcoíris (N→C)", "Por cadena", "Estructura secundaria", "Propiedad del residuo",
    "Factor B", "Cian", "Blanco", "Naranja", "Magenta",
]
SOLID_COLORS = {"Cian": "#19e3c5", "Blanco": "#f2f2f2", "Naranja": "#ff9f43", "Magenta": "#ff4fd8"}

BACKBONE_ATOMS = ['N', 'CA', 'C', 'O']
SS_COLORS = {'coil': '#6b7280', 'h': '#ff4b5c', 's': '#ffd166'}
PROPERTY_GROUPS = [
    (['ALA', 'VAL', 'LEU', 'ILE', 'MET', 'PHE', 'TRP', 'PRO', 'MSE'], '#ff9f43'),  # hidrofóbicos
    (['SER', 'THR', 'ASN', 'GLN', 'CYS', 'TYR', 'GLY'], '#19e3c5'),                # polares
    (['LYS', 'ARG', 'HIS'], '#4d7cff'),                                             # básicos (+)
    (['ASP', 'GLU'], '#ff4b5c'),                                                    # ácidos (−)
]


def style_dict(style_type, c):
    """Devuelve el estilo 3Dmol para un tipo de visualización y una especificación de color c."""
    if style_type == "Cartoon (flechas)":
        return {'cartoon': dict(c, arrows=True, thickness=0.35)}
    if style_type == "Cartoon ovalado":
        return {'cartoon': dict(c, style='oval', arrows=True, thickness=0.4)}
    if style_type == "Tubo (trace)":
        return {'cartoon': dict(c, style='trace', thickness=0.6)}
    if style_type == "Esqueleto (Backbone)":
        return {'stick': dict(c, radius=0.25)}
    if style_type == "Palos (Sticks)":
        return {'stick': dict(c, radius=0.14)}
    if style_type == "Bolas y palos":
        return {'stick': dict(c, radius=0.12), 'sphere': dict(c, scale=0.25)}
    if style_type == "Líneas (Wireframe)":
        return {'line': dict(c, linewidth=1.5)}
    if style_type == "Esferas (van der Waals)":
        return {'sphere': dict(c, scale=1.0)}
    return {'cartoon': dict(c, arrows=True, thickness=0.3)}  # Superficie + Cartoon


def simple_color(mode, rmin, rmax):
    """Especificación de color para los modos que no dependen de selecciones."""
    if mode == "Por cadena":
        return {'colorscheme': 'chain'}
    if mode == "Factor B":
        return {'colorscheme': {'prop': 'b', 'gradient': 'roygb', 'min': 80, 'max': 10}}  # alto B = rojo
    if mode in SOLID_COLORS:
        return {'color': SOLID_COLORS[mode]}
    return {'colorscheme': {'prop': 'resi', 'gradient': 'roygb', 'min': rmin, 'max': rmax}}


def build_3dmol_view(pdb_str, ligand_str, style_type, color_mode, show_protein, show_ligand,
                     rmin, rmax, surf_opacity, show_pocket, outline, spin, inter=None):
    view = py3Dmol.view(width=700, height=560)
    view.setBackgroundColor('#000000')
    if outline:
        view.setViewStyle({'style': 'outline', 'color': '#ffffff', 'width': 0.10})   # silueta blanca (el negro no se ve sobre fondo negro)

    model_idx = 0
    prot_idx = None
    lig_idx = None

    if show_protein and pdb_str:
        view.addModel(pdb_str, "pdb")
        prot_idx = model_idx
        model_idx += 1

        base = {'model': prot_idx}
        if style_type == "Esqueleto (Backbone)":
            view.setStyle(base, {})                 # ocultar todo y pintar solo el esqueleto
            base = dict(base, atom=BACKBONE_ATOMS)

        if color_mode == "Estructura secundaria":
            view.setStyle(base, style_dict(style_type, {'color': SS_COLORS['coil']}))
            view.setStyle(dict(base, ss='h'), style_dict(style_type, {'color': SS_COLORS['h']}))
            view.setStyle(dict(base, ss='s'), style_dict(style_type, {'color': SS_COLORS['s']}))
        elif color_mode == "Propiedad del residuo":
            view.setStyle(base, style_dict(style_type, {'color': '#9aa0aa'}))
            for resns, col in PROPERTY_GROUPS:
                view.setStyle(dict(base, resn=resns), style_dict(style_type, {'color': col}))
        else:
            view.setStyle(base, style_dict(style_type, simple_color(color_mode, rmin, rmax)))

        # Aguas e iones conservados -> esferas pequeñas
        view.addStyle({'model': prot_idx, 'resn': sorted(WATER_NAMES)},
                      {'sphere': {'radius': 0.3, 'color': '#ff4b5c'}})
        view.addStyle({'model': prot_idx, 'resn': sorted(COMMON_IONS)},
                      {'sphere': {'radius': 0.6, 'colorscheme': 'Jmol'}})

        if style_type == "Superficie + Cartoon":
            if color_mode == "Por cadena":
                surf = {'colorscheme': 'chain'}
            elif color_mode in SOLID_COLORS:
                surf = {'color': SOLID_COLORS[color_mode]}
            elif color_mode == "Arcoíris (N→C)":
                surf = {'colorscheme': {'prop': 'resi', 'gradient': 'roygb', 'min': rmin, 'max': rmax}}
            elif color_mode == "Factor B":
                surf = {'colorscheme': {'prop': 'b', 'gradient': 'roygb', 'min': 80, 'max': 10}}
            else:
                surf = {'color': '#b8bcc6'}
            view.addSurface(py3Dmol.SES, dict(surf, opacity=surf_opacity), {'model': prot_idx})

    if show_ligand and ligand_str and ligand_str.strip():
        view.addModel(ligand_str, "pdb")
        lig_idx = model_idx
        view.setStyle({'model': lig_idx},
                      {'stick': {'colorscheme': 'greenCarbon', 'radius': 0.25},
                       'sphere': {'colorscheme': 'greenCarbon', 'scale': 0.25}})

    if show_pocket and prot_idx is not None:   # lista de residuos calculada con pocket_residues()
        by_chain = {}
        for r in show_pocket:
            by_chain.setdefault(r["chain"], []).append(r["resnum"])
        for ch, nums in by_chain.items():
            view.addStyle({'model': prot_idx, 'chain': ch, 'resi': sorted(nums)},
                          {'stick': {'colorscheme': 'cyanCarbon', 'radius': 0.16}})
        if len(show_pocket) <= 40:
            for r in show_pocket:
                x, y, z = r["xyz"]
                view.addLabel(f"{r['resname']}{r['resnum']}",
                              {'position': {'x': x, 'y': y, 'z': z}, 'fontSize': 11, 'fontColor': 'white',
                               'backgroundColor': '#000000', 'backgroundOpacity': 0.5,
                               'borderThickness': 0, 'inFront': True, 'showBackground': True})

    if inter and inter.get("items"):
        draw_interactions(view, prot_idx, inter["items"], inter.get("labels", True))

    if show_pocket and prot_idx is not None:
        sels = [{'model': prot_idx, 'chain': ch, 'resi': sorted(nums)} for ch, nums in by_chain.items()]
        view.zoomTo(sels[0] if len(sels) == 1 else {'or': sels})     # centra la vista en el bolsillo
    else:
        view.zoomTo()
    if spin:
        view.spin(True)
    return view



# --- PREPARACIÓN PARA DOCKING (hidrógenos, cargas, SDF, PDBQT) ---

@st.cache_data(show_spinner=False)
def extract_single_ligand(pdb_text, resname, resnum, chain_id):
    """PDB con un único ligando (residuo concreto)."""
    structure = PDB.PDBParser(QUIET=True).get_structure("l", io.StringIO(pdb_text))
    model = next(iter(structure))

    class Sel(PDB.Select):
        def accept_model(self, m):
            return m.id == model.id

        def accept_chain(self, c):
            return c.id == chain_id

        def accept_residue(self, r):
            return is_ligand(r) and r.get_resname().strip() == resname and r.get_id()[1] == resnum

    buf = io.StringIO()
    w = PDB.PDBIO()
    w.set_structure(structure)
    w.save(buf, Sel())
    return buf.getvalue()


def fetch_ligand_smiles(resname):
    """SMILES de referencia del ligando (Chemical Component Dictionary de RCSB)."""
    url = f"https://data.rcsb.org/rest/v1/core/chemcomp/{resname.strip().upper()}"
    try:
        req = urllib.request.Request(url, headers={'User-Agent': 'Mozilla/5.0'})
        with urllib.request.urlopen(req, timeout=15) as r:
            data = json.loads(r.read().decode('utf-8'))
    except Exception:
        return ""
    d = data.get("rcsb_chem_comp_descriptor") or {}
    smi = d.get("smiles_stereo") or d.get("smiles")
    if not smi:
        for x in data.get("pdbx_chem_comp_descriptor") or []:
            if "SMILES" in str(x.get("type", "")):
                smi = x.get("descriptor")
                break
    return smi or ""


def pdb_coords(pdb_str):
    """Coordenadas de átomos pesados de un bloque PDB."""
    pts = []
    for l in pdb_str.splitlines():
        if l.startswith(("ATOM", "HETATM")):
            el = (l[76:78].strip() or l[12:16].strip()[:1]).upper()
            if el != "H":
                pts.append((float(l[30:38]), float(l[38:46]), float(l[46:54])))
    return pts


def _finalize_ligand(obmol, name, ph, info):
    """Cargas Gasteiger + SDF / PDB / PDBQT (Meeko si está disponible, si no Open Babel)."""
    obmol.calccharges("gasteiger")
    obmol.title = name
    net_charge = sum(a.formalcharge for a in obmol.atoms)

    sdf = obmol.write("sdf")
    pdb_h = obmol.write("pdb")

    pdbqt, engine = None, "Open Babel"
    if HAS_MEEKO and HAS_RDKIT:
        try:
            rd = Chem.MolFromMolBlock(sdf, removeHs=False)
            setups = MoleculePreparation().prepare(rd)
            txt, ok, err = PDBQTWriterLegacy.write_string(setups[0])
            if not ok:
                raise RuntimeError(err)
            pdbqt, engine = txt, "Meeko"
        except Exception as e:
            info.append(f"⚠️ Meeko falló ({str(e)[:100]}); PDBQT generado con Open Babel.")
    if pdbqt is None:
        pdbqt = obmol.write("pdbqt")

    tors = re.search(r"TORSDOF\s+(\d+)", pdbqt)
    return {
        "name": name, "sdf": sdf, "pdb": pdb_h, "pdbqt": pdbqt, "engine": engine, "info": info,
        "net_charge": net_charge, "n_atoms": len(obmol.atoms),
        "torsions": int(tors.group(1)) if tors else None,
    }


def prepare_ligand(lig_pdb, smiles, ph, name):
    """Ligando cristalográfico -> órdenes de enlace (plantilla SMILES) -> H a pH -> SDF/PDBQT."""
    info = []
    obmol = None
    if HAS_RDKIT and smiles.strip():
        try:
            raw = Chem.MolFromPDBBlock(lig_pdb, sanitize=False, removeHs=False)
            tpl = Chem.MolFromSmiles(smiles.strip())
            fixed = AllChem.AssignBondOrdersFromTemplate(tpl, raw)
            Chem.SanitizeMol(fixed)
            obmol = pybel.readstring("mol", Chem.MolToMolBlock(fixed))
            info.append("✅ Órdenes de enlace asignados con la plantilla SMILES.")
        except Exception as e:
            info.append(f"⚠️ No se pudo aplicar la plantilla SMILES ({str(e)[:120]}). "
                        "Se usa la percepción de enlaces de Open Babel.")
    if obmol is None:
        obmol = pybel.readstring("pdb", lig_pdb)
        info.append("⚠️ Enlaces y órdenes inferidos por geometría: revisa el resultado "
                    "(un SMILES correcto da un resultado mucho más fiable).")

    obmol.OBMol.DeleteHydrogens()
    obmol.OBMol.AddHydrogens(False, True, float(ph))   # todos los H, estados de ionización a este pH
    return _finalize_ligand(obmol, name, ph, info)


def safe_name(text, default="lig"):
    t = re.sub(r"[^A-Za-z0-9_.-]", "_", (text or "").strip())[:40].strip("._")
    return t or default


def read_external_ligands(filename, data, limit=500):
    """Lee ligandos de SDF/MOL/MOL2/PDB/SMI con Open Babel. Devuelve [(nombre, molécula pybel)]."""
    ext = filename.rsplit(".", 1)[-1].lower() if "." in filename else "smi"
    fmt = {"sdf": "sdf", "sd": "sdf", "mol": "mol", "mol2": "mol2", "pdb": "pdb", "ent": "pdb",
           "smi": "smi", "smiles": "smi", "txt": "smi"}.get(ext)
    if fmt is None:
        raise ValueError(f"Formato no soportado: .{ext}")
    base = safe_name(filename.rsplit(".", 1)[0], "lig")
    out, used = [], set()
    with tempfile.TemporaryDirectory() as td:
        path = os.path.join(td, f"in.{fmt}")
        with open(path, "wb") as f:
            f.write(data)
        for k, m in enumerate(pybel.readfile(fmt, path)):
            if k >= limit:
                break
            name = safe_name(m.title, f"{base}_{k + 1}")
            if name in used:
                name = f"{name}_{k + 1}"
            used.add(name)
            out.append((name, m))
    return out


def prepare_external_ligand(obmol, name, ph):
    """Ligando externo -> fragmento principal -> H a pH -> 3D si hace falta -> SDF/PDBQT."""
    info = []
    try:
        obmol.OBMol.StripSalts(0)           # conserva el fragmento mayor (quita contraiones)
    except Exception:
        pass
    dim = obmol.OBMol.GetDimension()
    obmol.OBMol.DeleteHydrogens()
    obmol.OBMol.AddHydrogens(False, True, float(ph))
    if dim == 3:
        info.append("✅ Se conservan las coordenadas 3D de la entrada; H añadidos a pH %s." % ph)
    else:
        built = False
        if HAS_RDKIT:
            try:
                rd = Chem.MolFromMolBlock(obmol.write("mol"), removeHs=False)
                params = AllChem.ETKDGv3()
                params.randomSeed = 42
                if AllChem.EmbedMolecule(rd, params) != 0:
                    params.useRandomCoords = True
                    if AllChem.EmbedMolecule(rd, params) != 0:
                        raise RuntimeError("no se pudo generar la conformación")
                try:
                    AllChem.MMFFOptimizeMolecule(rd, maxIters=1000)
                except Exception:
                    AllChem.UFFOptimizeMolecule(rd, maxIters=1000)
                obmol = pybel.readstring("mol", Chem.MolToMolBlock(rd))
                built = True
                info.append("✅ Conformación 3D generada con ETKDG + MMFF (sin estereoquímica definida "
                            "si el SMILES no la incluye).")
            except Exception as e:
                info.append(f"⚠️ RDKit no pudo generar el 3D ({str(e)[:80]}); se usa Open Babel.")
        if not built:
            obmol.make3D("mmff94", 500)
            info.append("✅ Conformación 3D generada con Open Babel (MMFF94).")
    return _finalize_ligand(obmol, name, ph, info)


def zip_files(files):
    """{nombre: texto} -> bytes de un ZIP."""
    buf = io.BytesIO()
    with zipfile.ZipFile(buf, "w", zipfile.ZIP_DEFLATED) as z:
        for fname, content in files.items():
            z.writestr(fname, content)
    return buf.getvalue()


def prepare_receptor(pdb_str, ph, method):
    """Proteína -> H a pH (PDB2PQR+PROPKA u Open Babel) -> cargas Gasteiger -> PDBQT rígido."""
    info = []
    # Sin aguas para docking
    pdb_str = "\n".join(l for l in pdb_str.splitlines()
                        if not (l.startswith("HETATM") and l[17:20].strip() in WATER_NAMES))
    # Solo átomos: HELIX/SHEET y demás cabeceras acaban en el PDBQT y Vina las rechaza
    pdb_str = "\n".join(l for l in pdb_str.splitlines() if l.startswith(("ATOM", "HETATM", "TER", "END")))
    add_h = True
    protonated = pdb_str

    if method.startswith("PDB2PQR"):
        with tempfile.TemporaryDirectory() as td:
            pin, pout, pqr = (os.path.join(td, n) for n in ("in.pdb", "out.pdb", "out.pqr"))
            with open(pin, "w") as f:
                f.write(pdb_str)
            cmd = [sys.executable, "-m", "pdb2pqr", "--ff", "AMBER", "--with-ph", str(ph),
                   "--titration-state-method", "propka", "--drop-water", "--keep-chain",
                   "--pdb-output", pout, pin, pqr]
            proc = subprocess.run(cmd, capture_output=True, text=True, timeout=900)
            if proc.returncode != 0 or not os.path.exists(pout):
                raise RuntimeError("PDB2PQR falló:\n" + (proc.stderr or proc.stdout)[-600:])
            with open(pout) as f:
                protonated = f.read()
        add_h = False
        info.append(f"✅ Protonación con PDB2PQR + PROPKA a pH {ph}.")
    else:
        info.append(f"✅ Protonación con Open Babel a pH {ph} (modelo simple de pKa; "
                    "para His, Cys o sitios activos críticos es preferible PDB2PQR+PROPKA).")

    mol = pybel.readstring("pdb", protonated)
    if add_h:
        mol.OBMol.AddHydrogens(False, True, float(ph))
    n_h = sum(1 for a in mol.atoms if a.atomicnum == 1)
    mol.calccharges("gasteiger")
    pdbqt = mol.write("pdbqt", opt={"r": None})   # receptor rígido, sin ROOT/BRANCH
    pdbqt = "\n".join(l for l in pdbqt.splitlines() if l.startswith(("ATOM", "HETATM", "REMARK"))) + "\n"
    pdb_h = mol.write("pdb") if add_h else protonated
    n_qt = sum(1 for l in pdbqt.splitlines() if l.startswith(("ATOM", "HETATM")))
    return {"pdbqt": pdbqt, "pdb": pdb_h, "n_h": n_h, "n_atoms_pdbqt": n_qt, "info": info}



# --- INTERACCIONES LIGANDO–PROTEÍNA ---

INTER_COLORS = {
    "Puente de hidrógeno": "#4d9bff",
    "Contacto hidrofóbico": "#c8c8c8",
    "Puente salino": "#ffd166",
    "Apilamiento π": "#ff4fd8",
    "Interacción π–catión": "#ff9f43",
    "Enlace de halógeno": "#19e3c5",
    "Puente de agua": "#7fd1ff",
    "Coordinación metálica": "#c084fc",
    "Contacto polar (aprox.)": "#4d9bff",
    "Contacto hidrofóbico (aprox.)": "#c8c8c8",
}
INTER_EMOJI = {
    "Puente de hidrógeno": "🔵", "Contacto hidrofóbico": "⚪", "Puente salino": "🟡",
    "Apilamiento π": "🟣", "Interacción π–catión": "🟠", "Enlace de halógeno": "🟢",
    "Puente de agua": "💧", "Coordinación metálica": "🟪",
    "Contacto polar (aprox.)": "🔵", "Contacto hidrofóbico (aprox.)": "⚪",
}


def _xyz(o):
    """Coordenadas (x, y, z) de un átomo, anillo o grupo cargado de PLIP."""
    for attr in ("coords", "center"):
        v = getattr(o, attr, None)
        if v is not None:
            try:
                return tuple(float(x) for x in v)
            except TypeError:
                pass
    for attr in ("o", "x", "oxy", "m", "atom"):
        sub = getattr(o, attr, None)
        if sub is not None:
            r = _xyz(sub)
            if r:
                return r
    return None


def _dist(a, b):
    return math.sqrt(sum((a[i] - b[i]) ** 2 for i in range(3)))


@st.cache_data(show_spinner=False)
def build_complex_pdb(raw_pdb, chains, resname, resnum, chain_id):
    """Complejo para el análisis: proteína (cadenas elegidas) + aguas + iones + UN ligando (con CONECT)."""
    selected = set(chains)
    structure = PDB.PDBParser(QUIET=True).get_structure("c", io.StringIO(raw_pdb))
    model = next(iter(structure))

    def is_target(r):
        return (is_ligand(r) and r.get_resname().strip() == resname
                and r.get_id()[1] == resnum and r.get_parent().id == chain_id)

    lig_serials = set()
    for r in model[chain_id]:
        if is_target(r):
            lig_serials = {a.serial_number for a in r}

    class Sel(PDB.Select):
        def accept_model(self, m):
            return m.id == model.id

        def accept_chain(self, c):
            return c.id in selected or c.id == chain_id

        def accept_residue(self, r):
            if is_target(r):
                return True
            if r.get_parent().id not in selected:
                return False
            return is_amino_acid(r) or is_water(r) or is_ion(r)

    buf = io.StringIO()
    w = PDB.PDBIO()
    w.set_structure(structure)
    w.save(buf, Sel(), preserve_atom_numbering=True)
    allowed = selected | {chain_id}
    # Los TER originales ocupan un nº de serie: hay que conservarlos para que los CONECT apunten bien
    raw_ter = {int(l[6:11]): l for l in raw_pdb.splitlines()
               if l.startswith("TER") and l[6:11].strip().isdigit() and len(l) > 21 and l[21] in allowed}
    atom_lines, kept = [], set()
    for l in buf.getvalue().splitlines():
        if l.startswith(("ATOM", "HETATM")) and l[6:11].strip().isdigit():
            atom_lines.append(l)
            serial = int(l[6:11])
            kept.add(serial)
            if serial + 1 in raw_ter:
                atom_lines.append(raw_ter[serial + 1])

    # PLIP necesita los registros de conectividad originales (enlaces covalentes ligando–proteína, SS, etc.)
    extra = []
    for l in raw_pdb.splitlines():
        if l.startswith("CONECT"):
            nums = [int(l[i:i + 5]) for i in range(6, min(len(l), 31), 5) if l[i:i + 5].strip().isdigit()]
            if nums and all(n in kept for n in nums):
                extra.append(l)
        elif l.startswith("LINK") and len(l) > 51 and l[21] in allowed and l[51] in allowed:
            extra.append(l)
        elif l.startswith("SSBOND") and len(l) > 29 and l[15] in allowed and l[29] in allowed:
            extra.append(l)
    header = [l for l in raw_pdb.splitlines() if l.startswith("HEADER")][:1]
    order = {"HEADER": 0, "LINK": 1, "SSBOND": 1}
    extra.sort(key=lambda l: order.get(l.split()[0], 2))
    head_ext = [l for l in extra if not l.startswith("CONECT")]
    conect = [l for l in extra if l.startswith("CONECT")]
    return "\n".join(header + head_ext + atom_lines + conect + ["END"]) + "\n"


_PYBEL_PATCHED = False


def _patch_pybel_for_plip():
    """PLIP pide la InChIKey del ligando solo para sus informes. Algunas instalaciones de
    Open Babel no incluyen el formato InChI y PLIP se detiene; aquí ese paso se vuelve opcional."""
    global _PYBEL_PATCHED
    if _PYBEL_PATCHED:
        return
    original_write = pybel.Molecule.write

    def safe_write(self, format="smi", filename=None, overwrite=False, opt=None):
        if format in ("inchikey", "inchi"):
            try:
                return original_write(self, format, filename, overwrite, opt)
            except Exception:
                return ""
        return original_write(self, format, filename, overwrite, opt)

    pybel.Molecule.write = safe_write
    _PYBEL_PATCHED = True


def _plip_interactions(complex_pdb, key):
    _patch_pybel_for_plip()
    from plip.structure.preparation import PDBComplex
    with tempfile.TemporaryDirectory() as td:
        path = os.path.join(td, "complex.pdb")
        with open(path, "w") as f:
            f.write(complex_pdb)
        mol = PDBComplex()
        mol.load_pdb(path)
        mol.analyze()
    sets = mol.interaction_sets
    iset = sets.get(key) or (next(iter(sets.values())) if sets else None)
    if iset is None:
        return []

    items = []

    def add(tipo, i, segs, label, dist=None, lig=None):
        segs = [(a, b) for a, b in segs if a and b]
        if not segs:
            return
        d = float(dist) if dist is not None else _dist(*segs[0])
        items.append({"tipo": tipo, "resname": str(i.restype), "resnum": int(i.resnr),
                      "chain": str(i.reschain), "dist": round(d, 2),
                      "segs": segs, "lig": lig or segs[0][1],
                      "label": label or tuple((segs[0][0][k] + segs[0][1][k]) / 2 for k in range(3))})

    for i in iset.hbonds_pdon:       # proteína = dador, ligando = aceptor
        add("Puente de hidrógeno", i, [(_xyz(i.d), _xyz(i.a))], _xyz(i.d),
            getattr(i, "distance_ad", None), _xyz(i.a))
    for i in iset.hbonds_ldon:       # proteína = aceptor, ligando = dador
        add("Puente de hidrógeno", i, [(_xyz(i.d), _xyz(i.a))], _xyz(i.a),
            getattr(i, "distance_ad", None), _xyz(i.d))
    for i in iset.hydrophobic_contacts:
        add("Contacto hidrofóbico", i, [(_xyz(i.bsatom), _xyz(i.ligatom))], _xyz(i.bsatom),
            i.distance, _xyz(i.ligatom))
    for i in list(iset.saltbridge_lneg) + list(iset.saltbridge_pneg):
        pos, neg = _xyz(i.positive), _xyz(i.negative)
        add("Puente salino", i, [(pos, neg)], pos if i.protispos else neg, i.distance,
            neg if i.protispos else pos)
    for i in iset.pistacking:
        add("Apilamiento π", i, [(_xyz(i.proteinring), _xyz(i.ligandring))], _xyz(i.proteinring),
            i.distance, _xyz(i.ligandring))
    for i in list(iset.pication_laro) + list(iset.pication_paro):
        ring, chg = _xyz(i.ring), _xyz(i.charge)
        add("Interacción π–catión", i, [(ring, chg)], chg if i.protcharged else ring, i.distance,
            ring if i.protcharged else chg)
    for i in iset.halogen_bonds:
        add("Enlace de halógeno", i, [(_xyz(i.don), _xyz(i.acc))], _xyz(i.acc), i.distance, _xyz(i.don))
    for i in iset.water_bridges:
        d, a, w = _xyz(i.d), _xyz(i.a), _xyz(i.water)
        add("Puente de agua", i, [(d, w), (w, a)], d if i.protisdon else a,
            getattr(i, "distance_aw", None), a if i.protisdon else d)
    for i in iset.metal_complexes:
        add("Coordinación metálica", i, [(_xyz(i.metal), _xyz(i.target))], None, i.distance, _xyz(i.target))
    return items


def _geometric_interactions(complex_pdb, resname, resnum, chain_id):
    """Respaldo sin PLIP: contactos polares (N/O ≤ 3.5 Å) e hidrofóbicos (C–C ≤ 4.0 Å). Aproximado."""
    structure = PDB.PDBParser(QUIET=True).get_structure("c", io.StringIO(complex_pdb))
    model = next(iter(structure))
    lig_atoms, prot_atoms = [], []
    for ch in model:
        for r in ch:
            if is_ligand(r) and r.get_resname().strip() == resname and r.get_id()[1] == resnum and ch.id == chain_id:
                lig_atoms += [a for a in r if a.element != "H"]
            elif is_amino_acid(r):
                prot_atoms += [(a, r, ch.id) for a in r if a.element != "H"]
    if not lig_atoms or not prot_atoms:
        return []
    ns = PDB.NeighborSearch([a for a, _, _ in prot_atoms])
    info = {id(a): (r, c) for a, r, c in prot_atoms}
    items, seen_hphob = [], set()
    for la in lig_atoms:
        for pa in ns.search(la.coord, 4.0):
            d = float(la - pa)
            r, c = info[id(pa)]
            el_l, el_p = la.element, pa.element
            if el_l in ("N", "O") and el_p in ("N", "O") and d <= 3.5:
                tipo = "Contacto polar (aprox.)"
            elif el_l == "C" and el_p == "C" and d <= 4.0 and id(pa) not in seen_hphob:
                tipo = "Contacto hidrofóbico (aprox.)"
                seen_hphob.add(id(pa))
            else:
                continue
            p1, p2 = tuple(float(x) for x in pa.coord), tuple(float(x) for x in la.coord)
            items.append({"tipo": tipo, "resname": r.get_resname().strip(), "resnum": int(r.get_id()[1]),
                          "chain": c, "dist": round(d, 2), "segs": [(p1, p2)], "label": p1, "lig": p2})
    return items


_ANALYSIS_VERSION = 3   # súbelo si cambias la lógica de análisis: invalida la caché


@st.cache_data(show_spinner=False)
def _plip_cached(raw_pdb, chains, resname, resnum, chain_id, _version=_ANALYSIS_VERSION):
    """Solo se cachean los éxitos: si PLIP lanza una excepción, Streamlit no guarda nada."""
    cpx = build_complex_pdb(raw_pdb, chains, resname, resnum, chain_id)
    return _plip_interactions(cpx, f"{resname}:{chain_id}:{resnum}")


@st.cache_data(show_spinner=False)
def _geometric_cached(raw_pdb, chains, resname, resnum, chain_id, _version=_ANALYSIS_VERSION):
    cpx = build_complex_pdb(raw_pdb, chains, resname, resnum, chain_id)
    return _geometric_interactions(cpx, resname, resnum, chain_id)


def analyze_interactions(raw_pdb, chains, resname, resnum, chain_id):
    """Devuelve (interacciones, método, aviso). Sin caché propia: los errores nunca se memorizan."""
    if HAS_PLIP:
        try:
            return _plip_cached(raw_pdb, chains, resname, resnum, chain_id), "PLIP", ""
        except Exception as e:
            fallback_msg = (f"PLIP falló ({type(e).__name__}: {str(e)[:120]}); "
                            "se usa el análisis geométrico aproximado.")
    else:
        fallback_msg = "PLIP no está instalado (`pip install plip`); se usa el análisis geométrico aproximado."
    try:
        return (_geometric_cached(raw_pdb, chains, resname, resnum, chain_id),
                "Geométrico (aproximado)", fallback_msg)
    except Exception as e:
        return [], "", f"No se pudieron calcular las interacciones: {e}"


def draw_interactions(view, prot_idx, items, labels):
    """Líneas discontinuas por interacción, residuos implicados en sticks y etiquetas."""
    for it in items:
        col = INTER_COLORS.get(it["tipo"], "#ffffff")
        for p1, p2 in it["segs"]:
            view.addCylinder({
                'start': {'x': p1[0], 'y': p1[1], 'z': p1[2]},
                'end': {'x': p2[0], 'y': p2[1], 'z': p2[2]},
                'radius': 0.06, 'color': col, 'dashed': True,
                'dashLength': 0.25, 'gapLength': 0.15, 'fromCap': 1, 'toCap': 1,
            })

    residues = {}
    for it in items:
        residues.setdefault((it["chain"], it["resnum"]), it)
    if prot_idx is not None:
        by_chain = {}
        for (ch, num) in residues:
            by_chain.setdefault(ch, []).append(num)
        for ch, nums in by_chain.items():
            view.addStyle({'model': prot_idx, 'chain': ch, 'resi': sorted(nums)},
                          {'stick': {'colorscheme': 'cyanCarbon', 'radius': 0.18}})
    if labels:
        for (ch, num), it in residues.items():
            x, y, z = it["label"]
            view.addLabel(f"{it['resname']}{num}",
                          {'position': {'x': x, 'y': y, 'z': z}, 'fontSize': 12, 'fontColor': 'white',
                           'backgroundColor': '#000000', 'backgroundOpacity': 0.55,
                           'borderThickness': 0, 'inFront': True, 'showBackground': True})



# --- MAPA 2D DE INTERACCIONES ---

_DASH = {
    "Puente de hidrógeno": "7,4", "Contacto polar (aprox.)": "7,4",
    "Contacto hidrofóbico": "2,4", "Contacto hidrofóbico (aprox.)": "2,4",
    "Puente salino": "9,3", "Apilamiento π": "10,3,2,3", "Interacción π–catión": "10,3,2,3",
    "Enlace de halógeno": "6,3", "Puente de agua": "3,3", "Coordinación metálica": "0",
}


def _res_color(resname):
    for resns, col in PROPERTY_GROUPS:
        if resname in resns:
            return col
    return "#9aa0aa"


def heavy_atoms(mol):
    """Molécula sin H. Con sanitización (recalcula los H implícitos de los aromáticos); si falla, sin ella."""
    try:
        return Chem.RemoveHs(mol)
    except Exception:
        h = Chem.RemoveHs(mol, sanitize=False)
        h.UpdatePropertyCache(strict=False)
        return h


def ligand_mol_for_map(lig_pdb, smiles=""):
    """Molécula RDKit del ligando (con órdenes de enlace si el SMILES encaja) y coordenadas 3D."""
    mol = None
    if smiles and smiles.strip():
        try:
            raw = Chem.MolFromPDBBlock(lig_pdb, sanitize=False, removeHs=True)
            mol = AllChem.AssignBondOrdersFromTemplate(Chem.MolFromSmiles(smiles.strip()), raw)
        except Exception:
            mol = None
    if mol is None:
        mol = Chem.MolFromPDBBlock(lig_pdb, sanitize=False, removeHs=True)
        if mol is not None:
            try:
                mol.UpdatePropertyCache(strict=False)
                Chem.SanitizeMol(mol, catchErrors=True)
            except Exception:
                pass
    return mol


def interaction_map_svg(mol, items, title="", dark=True, width=900, height=720):
    """Diagrama 2D estilo LigPlot: ligando dibujado con RDKit y residuos en corona con líneas por tipo."""
    from rdkit.Chem import rdDepictor
    from rdkit.Chem.Draw import rdMolDraw2D

    heavy = heavy_atoms(mol)
    if heavy.GetNumConformers() == 0:
        raise ValueError("El ligando no tiene coordenadas 3D.")
    conf = heavy.GetConformer()
    n = heavy.GetNumAtoms()
    pos3d = [(conf.GetAtomPosition(i).x, conf.GetAtomPosition(i).y, conf.GetAtomPosition(i).z)
             for i in range(n)]

    m2 = Chem.Mol(heavy)
    try:
        m2.UpdatePropertyCache(strict=False)
        Chem.FastFindRings(m2)
    except Exception:
        pass
    rdDepictor.Compute2DCoords(m2)

    d = rdMolDraw2D.MolDraw2DSVG(width, height)
    d.drawOptions().padding = 0.30
    if dark:
        try:
            rdMolDraw2D.SetDarkMode(d)
        except Exception:
            try:
                d.drawOptions().setBackgroundColour((0, 0, 0, 1))
            except Exception:
                pass
    d.DrawMolecule(m2)
    px = []
    for i in range(n):
        pt = d.GetDrawCoords(i)
        px.append((pt.x, pt.y))
    d.FinishDrawing()
    svg = d.GetDrawingText().replace("encoding='iso-8859-1'", "encoding='UTF-8'")

    fg = "#e8eaed" if dark else "#111111"
    bg_text = "#000000" if dark else "#ffffff"
    cx_l = sum(x for x, _ in px) / n
    cy_l = sum(y for _, y in px) / n

    def nearest(pt):
        return min(range(n), key=lambda k: sum((pt[a] - pos3d[k][a]) ** 2 for a in range(3)))

    nodes = {}
    for it in items:
        key = (it["chain"], it["resnum"], it["resname"])
        nodes.setdefault(key, {"edges": []})["edges"].append((it["tipo"], nearest(it["lig"]), it["dist"]))

    for key, nd in nodes.items():
        mx = sum(px[a][0] for _, a, _ in nd["edges"]) / len(nd["edges"])
        my = sum(px[a][1] for _, a, _ in nd["edges"]) / len(nd["edges"])
        nd["ang"] = math.atan2(my - cy_l, mx - cx_l)

    keys = sorted(nodes, key=lambda k: nodes[k]["ang"])
    N = len(keys)
    angs = [nodes[k]["ang"] for k in keys]
    if N > 1:
        gmin = min(0.42, 2 * math.pi / N * 0.9)
        for _ in range(300):
            for i in range(N):
                j = (i + 1) % N
                gap = angs[j] - angs[i] if j != 0 else angs[0] + 2 * math.pi - angs[-1]
                if gap < gmin:
                    delta = (gmin - gap) / 2
                    angs[i] -= delta
                    angs[j] += delta

    W2, H2 = width / 2, height / 2
    rx, ry = width / 2 - 70, height / 2 - 52
    bw, bh = 98, 34
    multi_chain = len({k[0] for k in keys}) > 1
    out = [f'<g font-family="Helvetica, Arial, sans-serif">']
    present = []
    labels_svg = []
    nodes_svg = []

    for ang, key in zip(angs, keys):
        ncx, ncy = W2 + rx * math.cos(ang), H2 + ry * math.sin(ang)
        color = _res_color(key[2])
        for e_i, (tipo, a, dist) in enumerate(nodes[key]["edges"]):
            if tipo not in present:
                present.append(tipo)
            ax, ay = px[a]
            vx, vy = ncx - ax, ncy - ay
            sx = (bw / 2) / abs(vx) if abs(vx) > 1e-6 else float("inf")
            sy = (bh / 2) / abs(vy) if abs(vy) > 1e-6 else float("inf")
            sc = min(sx, sy, 1.0)
            ex, ey = ncx - vx * sc, ncy - vy * sc
            lc = INTER_COLORS.get(tipo, "#ffffff")
            dash = _DASH.get(tipo, "6,4")
            out.append(f'<line x1="{ax:.1f}" y1="{ay:.1f}" x2="{ex:.1f}" y2="{ey:.1f}" '
                       f'stroke="{lc}" stroke-width="2.2" stroke-dasharray="{dash}" stroke-linecap="round"/>')
            t = 0.42 + 0.12 * e_i
            lx, ly = ax + (ex - ax) * t, ay + (ey - ay) * t
            labels_svg.append(   # halo (trazo) + texto: se ve igual en navegadores y en visores SVG simples
                f'<text x="{lx:.1f}" y="{ly - 4:.1f}" font-size="11" fill="{bg_text}" stroke="{bg_text}" '
                f'stroke-width="3.5" text-anchor="middle">{dist:.2f}</text>'
                f'<text x="{lx:.1f}" y="{ly - 4:.1f}" font-size="11" fill="{lc}" text-anchor="middle">'
                f'{dist:.2f}</text>')
        label = f"{key[2]}{key[1]}" + (f" ({key[0]})" if multi_chain else "")
        nodes_svg.append(
            f'<rect x="{ncx - bw / 2:.1f}" y="{ncy - bh / 2:.1f}" width="{bw}" height="{bh}" rx="9" '
            f'fill="{color}" fill-opacity="0.22" stroke="{color}" stroke-width="2"/>'
            f'<text x="{ncx:.1f}" y="{ncy + 5:.1f}" font-size="14" font-weight="bold" fill="{fg}" '
            f'text-anchor="middle">{html_lib.escape(label)}</text>')

    out += nodes_svg + labels_svg

    # Título y leyenda
    if title:
        out.append(f'<text x="14" y="24" font-size="16" font-weight="bold" fill="{fg}">'
                   f'{html_lib.escape(title)}</text>')
    y0 = height - 14 - 18 * (len(present) - 1)
    for k, tipo in enumerate(present):
        yy = y0 + 18 * k
        lc = INTER_COLORS.get(tipo, "#ffffff")
        out.append(f'<line x1="14" y1="{yy - 4}" x2="48" y2="{yy - 4}" stroke="{lc}" stroke-width="2.2" '
                   f'stroke-dasharray="{_DASH.get(tipo, "6,4")}"/>'
                   f'<text x="56" y="{yy}" font-size="12" fill="{fg}">{html_lib.escape(tipo)}</text>')
    leg = [("Hidrofóbico", "#ff9f43"), ("Polar", "#19e3c5"), ("Básico (+)", "#4d7cff"), ("Ácido (−)", "#ff4b5c")]
    for k, (txt, col) in enumerate(leg):
        yy = height - 14 - 18 * (len(leg) - 1 - k)
        out.append(f'<rect x="{width - 130}" y="{yy - 11}" width="14" height="14" rx="3" fill="{col}" '
                   f'fill-opacity="0.35" stroke="{col}"/><text x="{width - 110}" y="{yy}" font-size="12" '
                   f'fill="{fg}">{txt}</text>')
    out.append("</g>")
    return svg.replace("</svg>", "".join(out) + "</svg>")


def svg_img_html(svg, max_width="100%"):
    b64 = base64.b64encode(svg.encode("utf-8")).decode()
    return f'<img src="data:image/svg+xml;base64,{b64}" style="width:{max_width};border-radius:12px;border:1px solid #1e2430">'


def ensure_smiles(pdb_name, lig_label, resname):
    """SMILES del ligando cristalográfico (RCSB), guardado en sesión; sin red queda vacío."""
    k = f"smiles_{pdb_name}_{lig_label}"
    if k not in st.session_state:
        st.session_state[k] = fetch_ligand_smiles(resname)
    return st.session_state[k]


# --- RESULTADOS DE DOCKING: LECTURA DE POSES Y ANÁLISIS ---

_SCORE_RE = re.compile(r"REMARK VINA RESULT:\s+(-?\d+\.?\d*)\s+(-?\d+\.?\d*)\s+(-?\d+\.?\d*)")
_SCORE_KEYS = ("minimizedAffinity", "minimized_affinity", "vina_score", "Vina_score", "score",
               "Score", "docking_score", "r_i_docking_score", "affinity")


def parse_vina_scores(text):
    parts = re.split(r"^MODEL", text, flags=re.M)
    chunks = parts[1:] if len(parts) > 1 else [text]
    out = []
    for c in chunks:
        m = _SCORE_RE.search(c)
        out.append(tuple(float(x) for x in m.groups()) if m else (None, None, None))
    return out


def _mol_block(m):
    return Chem.MolToMolBlock(m, kekulize=False)


def detect_pose_format(filename, text):
    ext = filename.rsplit(".", 1)[-1].lower() if "." in filename else ""
    if ext in ("pdbqt", "sdf", "pdb"):
        return ext
    if "$$$$" in text or "M  END" in text:
        return "sdf"
    if "ROOT" in text and "ENDROOT" in text:
        return "pdbqt"
    return "pdb"


@st.cache_data(show_spinner=False)
def load_poses(filename, text, smiles=""):
    """Devuelve (poses, método, avisos). Cada pose: dict(molblock, score, rmsd_lb, rmsd_ub)."""
    fmt = detect_pose_format(filename, text)
    warnings, mols, method = [], [], ""
    scores = []

    if fmt == "pdbqt":
        scores = parse_vina_scores(text)
        if HAS_MEEKO_READ and "REMARK SMILES" in text:   # solo los PDBQT preparados con Meeko lo traen
            try:
                pm = PDBQTMolecule(text, is_dlg=False, skip_typing=True)
                base = RDKitMolCreate.from_pdbqt_mol(pm)[0]
                for conf in base.GetConformers():
                    m = Chem.Mol(base)
                    m.RemoveAllConformers()
                    m.AddConformer(Chem.Conformer(conf), assignId=True)
                    mols.append(m)
                method = "Meeko (órdenes de enlace y H tomados del SMILES guardado en el PDBQT)"
            except Exception as e:
                warnings.append(f"Meeko no pudo leer el PDBQT ({str(e)[:90]}); se usa Open Babel.")
                mols = []
    if fmt == "sdf":
        sup = Chem.SDMolSupplier()
        sup.SetData(text, removeHs=False)
        for m in sup:
            if m is None:
                continue
            sc = None
            props = m.GetPropsAsDict()
            for k in _SCORE_KEYS:
                if k in props:
                    try:
                        sc = float(props[k])
                        break
                    except (TypeError, ValueError):
                        pass
            scores.append((sc, None, None))
            mols.append(m)
        method = "SDF"
    if not mols and fmt in ("pdbqt", "pdb"):
        with tempfile.TemporaryDirectory() as td:
            path = os.path.join(td, f"p.{fmt}")
            with open(path, "w") as f:
                f.write(text)
            obs = list(pybel.readfile(fmt, path))
        for m in obs:
            rd = Chem.MolFromMolBlock(m.write("mol"), sanitize=False, removeHs=False)
            if rd is None:
                continue
            if smiles.strip():
                try:
                    heavy = Chem.RemoveHs(rd, sanitize=False)
                    rd = AllChem.AssignBondOrdersFromTemplate(Chem.MolFromSmiles(smiles.strip()), heavy)
                except Exception as e:
                    warnings.append(f"No se pudo aplicar el SMILES de plantilla ({str(e)[:80]}).")
            mols.append(rd)
        method = "Open Babel" + (" + plantilla SMILES" if smiles.strip() else
                                 " (enlaces inferidos; aporta un SMILES para mayor fiabilidad)")
        if fmt == "pdb":
            scores = [(None, None, None)] * len(mols)
        elif not scores:
            scores = [(None, None, None)] * len(mols)

    poses = []
    for k, m in enumerate(mols):
        sc = scores[k] if k < len(scores) else (None, None, None)
        poses.append({"molblock": _mol_block(m), "score": sc[0], "rmsd_lb": sc[1], "rmsd_ub": sc[2]})
    return poses, method, warnings


def pose_mol(molblock):
    m = Chem.MolFromMolBlock(molblock, removeHs=False, sanitize=False)
    m.UpdatePropertyCache(strict=False)
    Chem.SanitizeMol(m, catchErrors=True)
    return m


def build_pose_complex(protein_pdb, mol, chain="Z", resname="LIG", resnum=1):
    """Proteína + pose como HETATM con CONECT (numeración sin colisiones). Devuelve (pdb, clave PLIP)."""
    plines = [l for l in protein_pdb.splitlines() if l.startswith(("ATOM", "HETATM", "TER"))]
    sers = [int(l[6:11]) for l in plines if l[6:11].strip().isdigit()]
    last = max(sers) if sers else 0
    ter = f"TER   {last + 1:5d}"
    serial = last + 2

    heavy = heavy_atoms(mol)
    block = Chem.MolToPDBBlock(heavy)
    ser_map, lig_lines = {}, []
    for l in block.splitlines():
        if l.startswith(("ATOM", "HETATM")):
            old = int(l[6:11])
            ser_map[old] = serial
            lig_lines.append("HETATM" + f"{serial:5d}" + l[11:17] + f"{resname:>3s}" + " " + chain
                             + f"{resnum:4d}" + l[26:])
            serial += 1
    for l in block.splitlines():
        if l.startswith("CONECT"):
            nums = [int(l[i:i + 5]) for i in range(6, len(l.rstrip()), 5) if l[i:i + 5].strip().isdigit()]
            if nums and all(x in ser_map for x in nums):
                lig_lines.append("CONECT" + "".join(f"{ser_map[x]:5d}" for x in nums))
    text = "\n".join(plines + [ter] + lig_lines + ["END"]) + "\n"
    return text, f"{resname}:{chain}:{resnum}"


@st.cache_data(show_spinner=False)
def _plip_complex_cached(complex_pdb, key, _version=_ANALYSIS_VERSION):
    return _plip_interactions(complex_pdb, key)


def analyze_complex(complex_pdb, resname, resnum, chain_id):
    """(interacciones, método, aviso) sobre un complejo ya montado. Los errores no se cachean."""
    key = f"{resname}:{chain_id}:{resnum}"
    if HAS_PLIP:
        try:
            return _plip_complex_cached(complex_pdb, key), "PLIP", ""
        except Exception as e:
            msg = (f"PLIP falló ({type(e).__name__}: {str(e)[:120]}); "
                   "se usa el análisis geométrico aproximado.")
    else:
        msg = "PLIP no está instalado (`pip install plip`); se usa el análisis geométrico aproximado."
    try:
        return _geometric_interactions(complex_pdb, resname, resnum, chain_id), "Geométrico (aproximado)", msg
    except Exception as e:
        return [], "", f"No se pudieron calcular las interacciones: {e}"


def build_pose_view(protein_pdb, mol, items, labels=True, ref_ligand_pdb=""):
    view = py3Dmol.view(width=700, height=560)
    view.setBackgroundColor('#000000')
    view.addModel(protein_pdb, "pdb")
    view.setStyle({'model': 0}, {'cartoon': {'color': '#9aa0aa', 'opacity': 0.85}})
    heavy = heavy_atoms(mol)
    view.addModel(Chem.MolToMolBlock(heavy, kekulize=False), "sdf")
    view.setStyle({'model': 1}, {'stick': {'colorscheme': 'magentaCarbon', 'radius': 0.25}})
    if ref_ligand_pdb.strip():
        view.addModel(ref_ligand_pdb, "pdb")
        view.setStyle({'model': 2}, {'stick': {'colorscheme': 'greenCarbon', 'radius': 0.12}})
    if items:
        draw_interactions(view, 0, items, labels)
    view.zoomTo({'model': 1})
    return view


TYPE_CODE = {
    "Puente de hidrógeno": "HB", "Contacto hidrofóbico": "HP", "Puente salino": "SB",
    "Apilamiento π": "πS", "Interacción π–catión": "πC", "Enlace de halógeno": "HX",
    "Puente de agua": "WB", "Coordinación metálica": "MC",
    "Contacto polar (aprox.)": "HB~", "Contacto hidrofóbico (aprox.)": "HP~",
}


def Counter_types(items):
    """{código: nº} por tipo de interacción, para la tabla resumen."""
    out = {}
    for it in items:
        c = TYPE_CODE.get(it["tipo"], "?")
        out[c] = out.get(c, 0) + 1
    return out


def interaction_matrix(per_pose):
    """per_pose: [(etiqueta, items)] -> filas residuo × pose con códigos de interacción."""
    rows = {}
    for label, items in per_pose:
        for it in items:
            k = (it["chain"], it["resnum"], it["resname"])
            rows.setdefault(k, {}).setdefault(label, set()).add(TYPE_CODE.get(it["tipo"], "?"))
    table = []
    for (ch, num, rn), d in rows.items():
        row = {"Residuo": f"{rn}{num}", "Cadena": ch, "Nº poses": len(d)}
        for label, _ in per_pose:
            row[label] = ", ".join(sorted(d.get(label, [])))
        table.append((num, row))
    table.sort(key=lambda t: (-t[1]["Nº poses"], t[0]))
    return [r for _, r in table]



def render_ligand_results(results, pdb_name, prefix, kind="ext"):
    """Tabla, métricas y descargas para uno o varios ligandos preparados."""
    ok = [r for r in results if r.get("ok")]
    bad = [r for r in results if not r.get("ok")]
    r = ok[0] if ok else None
    if len(results) > 1:
        st.dataframe([{"Ligando": x["name"], "Átomos": x.get("n_atoms"), "Carga neta": x.get("net_charge"),
                       "Torsiones": x.get("torsions"),
                       "Estado": "✅" if x.get("ok") else "❌ " + x.get("error", "")} for x in results],
                     hide_index=True, use_container_width=True)
        if ok:
            st.download_button("📦 Todos los PDBQT (.zip)",
                               zip_files({f"{x['name']}.pdbqt": x["pdbqt"] for x in ok}),
                               file_name=f"{pdb_name}_ligandos_pdbqt.zip", mime="application/zip",
                               type="primary", use_container_width=True, key=f"{prefix}_zip_qt")
            st.download_button("📦 Todos los SDF (.zip)",
                               zip_files({f"{x['name']}.sdf": x["sdf"] for x in ok}),
                               file_name=f"{pdb_name}_ligandos_sdf.zip", mime="application/zip",
                               use_container_width=True, key=f"{prefix}_zip_sdf")
            chosen = st.selectbox("Ver / descargar uno:", [x["name"] for x in ok], key=f"{prefix}_sel")
            r = next(x for x in ok if x["name"] == chosen)
    elif bad:
        st.error(f"No se pudo preparar «{bad[0]['name']}»: {bad[0].get('error', '')}")
    if r:
        st.session_state["last_prepared"] = {"name": r["name"], "sdf": r["sdf"], "kind": kind}
        for msg in r["info"]:
            st.write(msg)
        a, b, c = st.columns(3)
        a.metric("Átomos (con H)", r["n_atoms"])
        b.metric("Carga neta", f"{r['net_charge']:+d}")
        c.metric("Torsiones", r["torsions"] if r["torsions"] is not None else "—")
        st.caption(f"PDBQT generado con {r['engine']} · cargas Gasteiger")
        st.download_button("📥 Ligando (.sdf, con H a pH)", r["sdf"].encode(), file_name=f"{r['name']}.sdf",
                           type="primary" if len(results) == 1 else "secondary", use_container_width=True,
                           key=f"{prefix}_sdf")
        st.download_button("📥 Ligando (.pdbqt)", r["pdbqt"].encode(), file_name=f"{r['name']}.pdbqt",
                           use_container_width=True, key=f"{prefix}_qt")
        st.download_button("📥 Ligando (.pdb con H)", r["pdb"].encode(), file_name=f"{r['name']}_H.pdb",
                           use_container_width=True, key=f"{prefix}_pdb")



# --- BOLSILLO, SESIÓN PYMOL, RMSD Y VINA ---

HAS_PYMOL = importlib.util.find_spec("pymol") is not None
HAS_VINA = importlib.util.find_spec("vina") is not None


@st.cache_data(show_spinner=False)
def pocket_residues(protein_pdb, ligand_pdb, dist=5.0):
    """Residuos con algún átomo pesado a ≤ dist Å de algún átomo del ligando."""
    if not ligand_pdb or not ligand_pdb.strip():
        return []
    prot = PDB.PDBParser(QUIET=True).get_structure("p", io.StringIO(protein_pdb))
    lig = PDB.PDBParser(QUIET=True).get_structure("l", io.StringIO(ligand_pdb))
    p_atoms = [a for a in prot.get_atoms() if a.element != "H" and is_amino_acid(a.get_parent())]
    if not p_atoms:
        return []
    ns = PDB.NeighborSearch(p_atoms)
    found = {}
    for la in lig.get_atoms():
        if la.element == "H":
            continue
        for pa in ns.search(la.coord, float(dist)):
            r = pa.get_parent()
            found[(r.get_parent().id, r.get_id()[1], r.get_resname().strip())] = r
    out = []
    for (ch, num, name), r in sorted(found.items(), key=lambda kv: (kv[0][0], kv[0][1])):
        atom = r["CA"] if "CA" in r else next(iter(r))
        out.append({"chain": ch, "resnum": int(num), "resname": name,
                    "xyz": tuple(float(x) for x in atom.coord)})
    return out


_PLIP_PSE_SCRIPT = r"""
import sys, os
outdir, path, key = sys.argv[1], sys.argv[2], sys.argv[3]
from openbabel import pybel
_orig = pybel.Molecule.write
def _safe(self, format="smi", filename=None, overwrite=False, opt=None):
    if format in ("inchikey", "inchi"):
        try:
            return _orig(self, format, filename, overwrite, opt)
        except Exception:
            return ""
    return _orig(self, format, filename, overwrite, opt)
pybel.Molecule.write = _safe
from plip.basic import config
config.PYMOL = True
config.SILENT = True
config.OUTPATH = outdir.rstrip("/") + "/"
config.BASEPATH = config.OUTPATH
from plip.structure.preparation import PDBComplex
from plip.basic.remote import VisualizerData
from plip.visualization.visualize import visualize_in_pymol
mol = PDBComplex()
mol.output_path = outdir
mol.load_pdb(path)
mol.analyze()
done = 0
for site in sorted(mol.interaction_sets):
    if site == key and len(mol.interaction_sets[site].interacting_res) > 0:
        visualize_in_pymol(VisualizerData(mol, site))
        done += 1
print("PSE_DONE", done)
"""


def build_pml(items, lig_sel, complex_name="complejo.pdb", session_name="interacciones.pse"):
    """Script de PyMOL que reproduce la visualización (respaldo cuando PyMOL no está instalado)."""
    py_colors = {"Puente de hidrógeno": "marine", "Contacto hidrofóbico": "gray70", "Puente salino": "yellow",
                 "Apilamiento π": "magenta", "Interacción π–catión": "orange", "Enlace de halógeno": "cyan",
                 "Puente de agua": "lightblue", "Coordinación metálica": "violet",
                 "Contacto polar (aprox.)": "marine", "Contacto hidrofóbico (aprox.)": "gray70"}
    lines = [f"load {complex_name}, complejo", "hide everything", "bg_color black",
             "show cartoon, polymer", "color gray70, polymer", f"select ligando, {lig_sel}",
             "show sticks, ligando", "util.cbag ligando", "set stick_radius, 0.2, ligando"]
    resi_sel = sorted({(it["chain"], it["resnum"]) for it in items})
    if resi_sel:
        parts = [f"(chain {c} and resi {n})" for c, n in resi_sel]
        lines += [f"select sitio, polymer and ({' or '.join(parts)})", "show sticks, sitio",
                  "util.cbac sitio", "label sitio and name CA, resn+resi"]
    for k, it in enumerate(items):
        col = py_colors.get(it["tipo"], "white")
        gname = "int_" + "".join(ch for ch in it["tipo"] if ch.isalnum())[:14]
        for j, (p1, p2) in enumerate(it["segs"]):
            a, b = f"ps_{k}_{j}a", f"ps_{k}_{j}b"
            lines += [f"pseudoatom {a}, pos=[{p1[0]:.3f},{p1[1]:.3f},{p1[2]:.3f}]",
                      f"pseudoatom {b}, pos=[{p2[0]:.3f},{p2[1]:.3f},{p2[2]:.3f}]",
                      f"distance {gname}, {a}, {b}", f"color {col}, {gname}",
                      f"set dash_color, {col}, {gname}", f"hide everything, {a} or {b}"]
    lines += ["set dash_gap, 0.3", "set dash_radius, 0.07", "set label_color, white", "set label_size, 14",
              "zoom ligando, 8", f"save {session_name}"]
    return "\n".join(lines) + "\n"


def export_pymol_session(complex_pdb, key, items, lig_sel):
    """(bytes, nombre de archivo, descripción). .pse real de PLIP si PyMOL está disponible; si no, ZIP con .pml."""
    if HAS_PYMOL and HAS_PLIP:
        try:
            with tempfile.TemporaryDirectory() as td:
                path = os.path.join(td, "complejo.pdb")
                with open(path, "w") as f:
                    f.write(complex_pdb)
                proc = subprocess.run([sys.executable, "-c", _PLIP_PSE_SCRIPT, td, path, key],
                                      capture_output=True, text=True, timeout=600, cwd=td)
                pses = [f for f in os.listdir(td) if f.endswith(".pse")]
                if proc.returncode == 0 and pses:
                    with open(os.path.join(td, pses[0]), "rb") as f:
                        return f.read(), pses[0], "Sesión de PyMOL generada por PLIP (.pse)."
                err = (proc.stderr or proc.stdout)[-300:]
        except Exception as e:
            err = str(e)[:300]
        note = f"PLIP/PyMOL no pudo crear el .pse ({err.strip()[-120:]}); se entrega un script .pml. "
    else:
        note = ("PyMOL no está instalado (`pip install pymol-open-source`, o `conda install -c conda-forge "
                "pymol-open-source`); se entrega un script .pml. ")
    pml = build_pml(items, lig_sel)
    zbytes = zip_files({"complejo.pdb": complex_pdb, "interacciones.pml": pml,
                        "LEEME.txt": "Abre PyMOL en esta carpeta y ejecuta:  pymol interacciones.pml\n"
                                     "Se crea interacciones.pse con el complejo y las interacciones.\n"})
    return zbytes, "interacciones_pymol.zip", note + "Ejecuta `pymol interacciones.pml` para obtener el .pse."


def _skeleton(mol):
    """Átomos pesados y enlaces simples (sin órdenes): permite comparar grafos de distinta procedencia."""
    import numpy as np
    heavy = heavy_atoms(mol)
    rw = Chem.RWMol()
    for a in heavy.GetAtoms():
        rw.AddAtom(Chem.Atom(a.GetAtomicNum()))
    for b in heavy.GetBonds():
        rw.AddBond(b.GetBeginAtomIdx(), b.GetEndAtomIdx(), Chem.BondType.SINGLE)
    m = rw.GetMol()
    m.UpdatePropertyCache(strict=False)
    Chem.FastFindRings(m)
    conf = heavy.GetConformer()
    xyz = np.array([[conf.GetAtomPosition(i).x, conf.GetAtomPosition(i).y, conf.GetAtomPosition(i).z]
                    for i in range(heavy.GetNumAtoms())])
    return m, xyz


def ligand_rmsd(pose, ref):
    """RMSD de átomos pesados SIN alinear (misma referencia espacial), con simetría molecular.
    Devuelve (rmsd | None, nota)."""
    import numpy as np
    pm, pxyz = _skeleton(pose)
    rm, rxyz = _skeleton(ref)
    if pm.GetNumAtoms() != rm.GetNumAtoms():
        return None, f"distinto nº de átomos pesados (pose {pm.GetNumAtoms()}, referencia {rm.GetNumAtoms()})"
    matches = pm.GetSubstructMatches(rm, uniquify=False, maxMatches=20000, useChirality=False)
    if matches:
        best = min(float(np.sqrt(((pxyz[list(m)] - rxyz) ** 2).sum(axis=1).mean())) for m in matches)
        return best, "exacto (con simetría)"
    try:   # respaldo: emparejamiento óptimo por elemento (no respeta el grafo; es una cota inferior)
        from scipy.optimize import linear_sum_assignment
        p_el = [a.GetAtomicNum() for a in pm.GetAtoms()]
        r_el = [a.GetAtomicNum() for a in rm.GetAtoms()]
        tot = 0.0
        for el in set(r_el):
            ri = [i for i, e in enumerate(r_el) if e == el]
            pi = [i for i, e in enumerate(p_el) if e == el]
            if len(ri) != len(pi):
                return None, "composición distinta entre pose y referencia"
            cost = ((rxyz[ri][:, None, :] - pxyz[pi][None, :, :]) ** 2).sum(axis=2)
            r_, c_ = linear_sum_assignment(cost)
            tot += cost[r_, c_].sum()
        return float(np.sqrt(tot / len(r_el))), "aproximado (por elemento; los grafos no coinciden)"
    except Exception:
        return None, "los grafos de la pose y la referencia no coinciden"


def find_vina_exe(custom=""):
    """Localiza el ejecutable oficial de Vina: ruta indicada, junto al script/directorio actual o en el PATH."""
    cands = []
    if custom and custom.strip():
        cands.append(custom.strip().strip('"'))
    dirs = [os.getcwd()]
    try:
        dirs.append(os.path.dirname(os.path.abspath(__file__)))
    except NameError:
        pass
    skip = (".py", ".pdbqt", ".txt", ".zip", ".md", ".log", ".pdb", ".sdf")
    for d in dirs:
        try:
            for f in sorted(os.listdir(d)):
                fl = f.lower()
                if fl.startswith("vina") and not fl.startswith("vina_split") and not fl.endswith(skip) \
                        and os.path.isfile(os.path.join(d, f)):
                    cands.append(os.path.join(d, f))
        except OSError:
            pass
    cands += ["vina", "vina.exe"]
    for c in cands:
        path = c if os.path.isfile(c) else shutil.which(c)
        if not path:
            continue
        try:
            r = subprocess.run([path, "--version"], capture_output=True, text=True, timeout=20)
            if "vina" in (r.stdout + r.stderr).lower():
                return path
        except Exception:
            continue
    return None


def _run_vina_exe(exe, receptor_pdbqt, ligand_pdbqt, center, size, exhaustiveness, n_poses, energy_range,
                  seed, cpu, scoring):
    with tempfile.TemporaryDirectory() as td:
        rp, lp, op = (os.path.join(td, n) for n in ("receptor.pdbqt", "ligand.pdbqt", "out.pdbqt"))
        with open(rp, "w") as f:
            f.write(receptor_pdbqt)
        with open(lp, "w") as f:
            f.write(ligand_pdbqt)
        cmd = [exe, "--receptor", rp, "--ligand", lp, "--out", op,
               "--center_x", f"{center[0]:.3f}", "--center_y", f"{center[1]:.3f}", "--center_z", f"{center[2]:.3f}",
               "--size_x", f"{size[0]:.1f}", "--size_y", f"{size[1]:.1f}", "--size_z", f"{size[2]:.1f}",
               "--exhaustiveness", str(int(exhaustiveness)), "--num_modes", str(max(int(n_poses), 1)),
               "--energy_range", str(float(energy_range)), "--seed", str(int(seed)),
               "--cpu", str(int(cpu)), "--scoring", scoring]
        proc = subprocess.run(cmd, capture_output=True, text=True, timeout=7200)
        if proc.returncode != 0 or not os.path.exists(op):
            raise RuntimeError((proc.stderr or proc.stdout or "Vina terminó sin resultados").strip()[-300:])
        with open(op) as f:
            pdbqt = f.read()
    energies = [[a, b, c] for a, b, c in
                (tuple(float(x) for x in m.groups()) for m in _SCORE_RE.finditer(pdbqt))]
    return pdbqt, energies


def run_vina(receptor_pdbqt, ligand_pdbqt, center, size, exhaustiveness, n_poses, energy_range,
             seed, cpu, scoring, exe=None):
    """Docking con Vina: paquete Python (si no se fuerza el ejecutable) o ejecutable oficial.
    Devuelve (pdbqt de poses, energías por pose; la primera columna es la afinidad)."""
    if exe or not HAS_VINA:
        if not exe:
            raise RuntimeError("No hay motor de Vina disponible (ni paquete Python ni ejecutable).")
        return _run_vina_exe(exe, receptor_pdbqt, ligand_pdbqt, center, size, exhaustiveness, n_poses,
                             energy_range, seed, cpu, scoring)
    from vina import Vina
    with tempfile.TemporaryDirectory() as td:
        rp = os.path.join(td, "receptor.pdbqt")
        with open(rp, "w") as f:
            f.write(receptor_pdbqt)
        v = Vina(sf_name=scoring, seed=int(seed), cpu=int(cpu), verbosity=0)
        v.set_receptor(rp)
        v.set_ligand_from_string(ligand_pdbqt)
        v.compute_vina_maps(center=[float(c) for c in center], box_size=[float(x) for x in size])
        v.dock(exhaustiveness=int(exhaustiveness), n_poses=max(int(n_poses), 1))
        pdbqt = v.poses(n_poses=int(n_poses), energy_range=float(energy_range))
        energies = v.energies(n_poses=int(n_poses), energy_range=float(energy_range))
    return pdbqt, [[round(float(x), 3) for x in row] for row in energies]


# --- INTERFAZ STREAMLIT ---

st.title("🧬 BioDock Studio: Inspector PDB & Preparación de Docking")
st.markdown("Herramienta bioinformática avanzada para visualización y preparación de complejos macromoleculares.")
st.markdown("---")

st.markdown("### 📥 1. Cargar Estructura Proteica")
col_input1, col_input2 = st.columns([1, 2])

with col_input1:
    input_method = st.radio("Método de entrada:",
                            ["Código PDB (Ej. 1HSG, 4AKE)", "Subir archivo .pdb local"],
                            key="input_method")

with col_input2:
    if input_method == "Código PDB (Ej. 1HSG, 4AKE)":
        c_code, c_btn = st.columns([2, 1])
        with c_code:
            pdb_code = st.text_input("Código PDB:", value="1HSG").strip()
        with c_btn:
            st.markdown("<br>", unsafe_allow_html=True)
            btn_fetch = st.button("Buscar en RCSB", type="primary")

        if btn_fetch and pdb_code:
            with st.spinner("Descargando..."):
                data = fetch_pdb_from_web(pdb_code)
                if data:
                    st.session_state["raw_pdb"] = data
                    st.session_state["pdb_name"] = pdb_code.upper()
                    st.success("¡Estructura cargada con éxito!")
                else:
                    st.error("Error al descargar el PDB.")
    else:
        uploaded_file = st.file_uploader("Sube tu archivo .pdb", type=["pdb"])
        if uploaded_file is not None:
            # getvalue() (no read()): read() devuelve vacío en las recargas siguientes
            st.session_state["raw_pdb"] = uploaded_file.getvalue().decode('utf-8', errors='ignore')
            st.session_state["pdb_name"] = uploaded_file.name.rsplit('.', 1)[0]

if st.session_state.get("raw_pdb"):
    raw_pdb = st.session_state["raw_pdb"]
    pdb_name = st.session_state.get("pdb_name", "Protein")
    chain_summary = get_chain_summary(raw_pdb)
    all_chains = list(chain_summary.keys())

    st.markdown("---")
    tab1, tab2, tab3, tab4, tab5 = st.tabs(["👁️ Visor 3D & Preparación", "📊 Propiedades Fisicoquímicas & Secuencia",
                                            "🧪 Preparación Docking (Vina)", "🚀 Docking (Vina)",
                                            "📈 Resultados de Docking"])

    with tab1:
        col_vis, col_panel = st.columns([1.5, 1])

        with col_panel:
            st.markdown(f"### ⚙️ Panel de Limpieza: `{pdb_name}`")

            with st.container(border=True):
                st.markdown("**Opciones de Filtrado & Limpieza** (se aplican al instante)")
                selected_chains = st.multiselect(
                    "Cadenas a conservar:", all_chains, default=all_chains,
                    format_func=lambda c: f"{c} ({chain_summary[c]} aa)" if chain_summary[c] else f"{c} (sin proteína)",
                    key=f"chains_{pdb_name}",   # clave por estructura: reinicia al cargar otra
                )
                remove_water = st.checkbox("Eliminar moléculas de agua (HOH)", value=True)
                remove_ions = st.checkbox("Eliminar iones y sales (Na, Cl, Mg...)", value=True)

            if not selected_chains:
                st.warning("Selecciona al menos una cadena.")
                st.stop()

            result = process_pdb_structure(raw_pdb, tuple(selected_chains), remove_water, remove_ions)
            total_res = sum(len(v) for v in result["chains"].values())
            unique_ligands = {f"{n} {num}:{ch}" for n, num, ch in result["ligands"]}

            with st.container(border=True):
                st.markdown("**Resumen Activo:**")
                st.write(f"• **Cadenas activas:** `{', '.join(result['chains'].keys()) or 'Ninguna'}`")
                st.write(f"• **Residuos de aminoácido:** `{total_res}`")
                st.write(f"• **Aguas:** `{result['water']}` detectadas → "
                         f"{'eliminadas' if remove_water else 'conservadas'}")
                st.write(f"• **Iones:** `{result['ions']}` detectados → "
                         f"{'eliminados' if remove_ions else 'conservados'}")
                st.write(f"• **Ligandos orgánicos:** `{len(unique_ligands)}` "
                         f"{('(' + ', '.join(sorted(unique_ligands)) + ')') if unique_ligands else ''}")
                st.caption(f"Estructura secundaria: {result['ss_source']}")

            st.markdown("#### 💾 Exportación para Docking")
            st.download_button(
                label="📥 Descargar Proteína Limpia (.pdb)",
                data=result["cleaned"].encode('utf-8'),
                file_name=f"{pdb_name}_clean.pdb",
                mime="text/plain",
                type="primary",
                use_container_width=True,
            )
            if result["ligand"].strip():
                st.download_button(
                    label="📥 Descargar Ligando Extraído (.pdb)",
                    data=result["ligand"].encode('utf-8'),
                    file_name=f"{pdb_name}_ligand.pdb",
                    mime="text/plain",
                    use_container_width=True,
                )
                st.caption("Para docking usa la pestaña 🧪: ahí obtienes SDF/PDBQT con H y cargas.")

        with col_vis:
            st.markdown("### 👁️ Visor 3D Interactivo")

            c_style, c_col = st.columns(2)
            with c_style:
                style_opt = st.selectbox("Estilo:", STYLE_OPTIONS)
            with c_col:
                color_opt = st.selectbox("Color:", COLOR_OPTIONS)

            c1, c2, c3, c4 = st.columns(4)
            with c1:
                show_protein = st.checkbox("Proteína", value=True)
            with c2:
                show_ligand = st.checkbox("Ligando", value=True)
            with c3:
                outline = st.checkbox("Contorno", value=True,
                                      help="Dibuja una línea blanca alrededor de las formas: separa las "
                                           "zonas que se solapan y da relieve a la estructura.")
            with c4:
                spin = st.checkbox("Girar", value=False)

            pk1, pk2 = st.columns([2, 1])
            with pk1:
                pocket_on = st.checkbox("Resaltar residuos del bolsillo (cercanos al ligando)", value=False)
            pocket_list = None
            if pocket_on:
                with pk2:
                    pocket_dist = st.slider("Distancia (Å)", 3.0, 8.0, 5.0, 0.5, key="pocket_dist")
                pk_opts = sorted({f"{n} {num}:{ch}" for n, num, ch in result["ligands"]})
                if pk_opts:
                    with pk1:
                        pk_lig = st.selectbox("Ligando de referencia del bolsillo:", pk_opts,
                                              key=f"pocket_lig_{pdb_name}")
                    pk_res, pk_rest = pk_lig.split(" ", 1)
                    pk_num, pk_chain = pk_rest.split(":")
                    pk_pdb = extract_single_ligand(raw_pdb, pk_res, int(pk_num), pk_chain)
                    pocket_list = pocket_residues(result["cleaned"], pk_pdb, pocket_dist)
                    st.caption(f"{len(pocket_list)} residuos a ≤ {pocket_dist:.1f} Å de {pk_lig}: " +
                               ", ".join(f"{r['resname']}{r['resnum']}" for r in pocket_list[:25]) +
                               ("…" if len(pocket_list) > 25 else ""))
                else:
                    st.caption("No hay ligando orgánico en las cadenas seleccionadas: no hay bolsillo que resaltar.")
            surf_opacity = 0.6
            if style_opt == "Superficie + Cartoon":
                surf_opacity = st.slider("Opacidad de la superficie", 0.1, 1.0, 0.6, 0.05)

            if color_opt == "Estructura secundaria":
                st.caption("🔴 Hélice α · 🟡 Lámina β · ⚫ Lazo")
            elif color_opt == "Propiedad del residuo":
                st.caption("🟠 Hidrofóbico · 🟢 Polar · 🔵 Básico (+) · 🔴 Ácido (−)")
            elif color_opt == "Factor B":
                st.caption("🔴 B alto (móvil) → 🔵 B bajo (rígido)")

            # ---- Interacciones ligando–proteína ----
            inter_data = None
            lig_inter_opts = sorted({f"{n} {num}:{ch}" for n, num, ch in result["ligands"]})
            if lig_inter_opts:
                show_inter = st.checkbox("Mostrar interacciones ligando–proteína", value=False)
                if show_inter:
                    ci1, ci2 = st.columns([1, 2])
                    with ci1:
                        sel_lig = st.selectbox("Ligando:", lig_inter_opts, key=f"inter_lig_{pdb_name}")
                    i_res, i_rest = sel_lig.split(" ", 1)
                    i_num, i_chain = i_rest.split(":")
                    with st.spinner("Analizando interacciones..."):
                        i_items, i_method, i_err = analyze_interactions(
                            raw_pdb, tuple(selected_chains), i_res, int(i_num), i_chain)
                    if i_err:
                        st.warning(i_err)
                    if i_items:
                        found = sorted({it["tipo"] for it in i_items})
                        with ci2:
                            chosen = st.multiselect("Tipos:", found, default=found,
                                                    key=f"inter_types_{pdb_name}_{sel_lig}")
                        show_labels = st.checkbox("Etiquetas de residuos", value=True)
                        inter_data = {"items": [it for it in i_items if it["tipo"] in chosen],
                                      "labels": show_labels}
                        st.caption(f"Método: {i_method} · " +
                                   " ".join(f"{INTER_EMOJI.get(t, '▫️')} {t}" for t in chosen))
                    elif not i_err:
                        st.info("No se detectaron interacciones para este ligando.")

            all_nums = [r[1] for v in result["chains"].values() for r in v] or [1, 2]
            view_obj = build_3dmol_view(
                result["cleaned"], result["ligand"], style_opt, color_opt,
                show_protein, show_ligand, min(all_nums), max(all_nums),
                surf_opacity, pocket_list, outline, spin, inter_data,
            )
            showmol(view_obj, height=560, width=700)

            if inter_data and inter_data["items"]:
                st.markdown("**Interacciones detectadas**")
                st.dataframe(
                    [{"Tipo": it["tipo"], "Residuo": f"{it['resname']}{it['resnum']}",
                      "Cadena": it["chain"], "Distancia (Å)": it["dist"]}
                     for it in sorted(inter_data["items"], key=lambda t: (t["chain"], t["resnum"], t["tipo"]))],
                    use_container_width=True, hide_index=True)

                with st.expander("🗺️ Mapa 2D de interacciones", expanded=False):
                    dark_map = st.checkbox("Fondo oscuro", value=True, key=f"map_dark_{pdb_name}")
                    smi_map = ensure_smiles(pdb_name, sel_lig, i_res)
                    lig_pdb_map = extract_single_ligand(raw_pdb, i_res, int(i_num), i_chain)
                    mol_map = ligand_mol_for_map(lig_pdb_map, smi_map)
                    if mol_map is None:
                        st.warning("No se pudo construir la molécula del ligando para el mapa 2D.")
                    else:
                        try:
                            svg_map = interaction_map_svg(mol_map, inter_data["items"],
                                                          title=f"{pdb_name} · {sel_lig}", dark=dark_map)
                            st.markdown(svg_img_html(svg_map), unsafe_allow_html=True)
                            st.download_button("📥 Mapa 2D (.svg)", svg_map.encode("utf-8"),
                                               file_name=f"{pdb_name}_{i_res}_mapa2D.svg", mime="image/svg+xml",
                                               key=f"dl_map_{pdb_name}")
                            if not smi_map.strip():
                                st.caption("Sin SMILES de referencia los enlaces se infieren por geometría "
                                           "(anillos sin dobles enlaces). Con red, la app lo descarga de RCSB.")
                        except Exception as e:
                            st.warning(f"No se pudo generar el mapa 2D: {e}")
                    cpx_dl = build_complex_pdb(raw_pdb, tuple(selected_chains), i_res, int(i_num), i_chain)
                    st.download_button("📥 Complejo proteína–ligando (.pdb) para LigPlot+ / PyMOL",
                                       cpx_dl.encode("utf-8"), file_name=f"{pdb_name}_{i_res}_complejo.pdb",
                                       key=f"dl_cpx_{pdb_name}")
                    st.caption("LigPlot+ es una aplicación Java externa con licencia propia y no puede "
                               "integrarse en la app; este PDB sirve como entrada.")

                    st.markdown("**Sesión de PyMOL (.pse) con las interacciones de PLIP**")
                    pse_key = f"pse_{pdb_name}_{sel_lig}"
                    if st.button("🧩 Generar sesión de PyMOL (.pse)", key=f"btn_{pse_key}"):
                        with st.spinner("Generando la sesión con PLIP + PyMOL..."):
                            st.session_state[pse_key] = export_pymol_session(
                                cpx_dl, f"{i_res}:{i_chain}:{int(i_num)}", i_items,
                                f"resn {i_res} and chain {i_chain} and resi {int(i_num)}")
                    if pse_key in st.session_state:
                        pse_data, pse_name, pse_desc = st.session_state[pse_key]
                        st.download_button(f"📥 Descargar {pse_name}", pse_data, file_name=pse_name,
                                           mime="application/octet-stream", type="primary", key=f"dl_{pse_key}")
                        st.caption(pse_desc)

    with tab2:
        st.markdown("### 📊 Propiedades Fisicoquímicas y Secuencia Primaria")
        chains = result["chains"]

        if chains:
            chain_id = st.selectbox("Seleccionar Cadena para análisis:", list(chains.keys()),
                                    key=f"tab2_chain_{pdb_name}")
            residues = chains[chain_id]
            sequence_str = "".join(r[0] for r in residues)
            props = calculate_protein_properties(sequence_str)

            if props:
                n_h = sum(1 for r in residues if r[2] == 'H')
                n_e = sum(1 for r in residues if r[2] == 'E')

                m1, m2, m3, m4 = st.columns(4)
                m1.metric("Longitud", f"{props['length']} aa")
                m2.metric("Peso Molecular", f"{props['mw_kDa']} kDa")
                m3.metric("Punto Isoeléctrico", f"{props['pi']}")
                m4.metric("Inestabilidad", f"{props['instability']}",
                          delta="Estable" if props['instability'] < 40 else "Inestable")

                m5, m6, m7 = st.columns(3)
                m5.metric("Aromaticidad", f"{props['aromaticity']} %")
                m6.metric("% Hélice α", f"{100 * n_h / len(residues):.1f} %")
                m7.metric("% Lámina β", f"{100 * n_e / len(residues):.1f} %")

            st.markdown(f"**Secuencia primaria — cadena {chain_id}** (H = hélice, E = lámina)")
            st.code(format_chain_block(residues), language=None)

            fasta = "".join(
                f">{pdb_name}_{cid}\n{''.join(r[0] for r in res)}\n" for cid, res in chains.items()
            )
            st.download_button("📥 Descargar FASTA (cadenas seleccionadas)", data=fasta,
                               file_name=f"{pdb_name}.fasta", mime="text/plain")
        else:
            st.warning("No se encontraron residuos de aminoácidos en las cadenas seleccionadas.")

    with tab3:
        st.markdown("### 🧪 Preparación para AutoDock Vina")
        st.caption("Añade hidrógenos según el pH, calcula cargas parciales (Gasteiger) y exporta "
                   "SDF / PDBQT. Las aguas se eliminan siempre en esta etapa.")

        lig_options = sorted({f"{n} {num}:{ch}" for n, num, ch in result["ligands"]})

        missing = []
        if not HAS_OB:
            missing.append("openbabel-wheel")
        if not HAS_RDKIT:
            missing.append("rdkit")
        if not HAS_MEEKO:
            missing.append("meeko gemmi")
        if not HAS_PDB2PQR:
            missing.append("pdb2pqr")
        if missing:
            st.info("Librerías opcionales no instaladas: `pip install " + " ".join(missing) + "`. "
                    "Open Babel es imprescindible; RDKit mejora el ligando (órdenes de enlace), "
                    "Meeko genera el mejor PDBQT de ligando y PDB2PQR+PROPKA la mejor protonación de proteína.")
        if not HAS_OB:
            st.error("Esta pestaña necesita Open Babel (`pip install openbabel-wheel`).")
        else:
            ph = st.slider("pH de protonación", 4.0, 10.0, 7.4, 0.1)
            col_lig, col_rec = st.columns(2)

            # ---------------- LIGANDO ----------------
            with col_lig:
                st.markdown("#### 💊 Ligando → SDF / PDBQT")
                sources = (["Ligando cristalográfico del PDB"] if lig_options else []) + \
                          ["Archivo externo (SDF, MOL2, PDB, SMI)", "SMILES pegados"]
                source = st.radio("Origen del ligando:", sources, key="lig_source")

                if source.startswith("Ligando cristalográfico"):
                    lig_label = st.selectbox("Ligando:", lig_options, key=f"lig_{pdb_name}")
                    resname, rest = lig_label.split(" ", 1)
                    resnum, chain_id = rest.split(":")
                    resnum = int(resnum)

                    smi_key = f"smiles_{pdb_name}_{lig_label}"
                    ensure_smiles(pdb_name, lig_label, resname)

                    def _refetch(k=smi_key, rn=resname):
                        st.session_state[k] = fetch_ligand_smiles(rn)

                    smiles = st.text_input("SMILES de referencia (para los órdenes de enlace):", key=smi_key)
                    st.button("🔎 Reintentar búsqueda en RCSB", on_click=_refetch, key=f"refetch_{smi_key}")
                    if not smiles.strip():
                        st.caption("Sin SMILES se infieren los enlaces por geometría (menos fiable). "
                                   "Pega uno si lo conoces.")

                    lig_sig = (pdb_name, lig_label, ph, smiles.strip())
                    if st.button("⚗️ Preparar ligando", type="primary", use_container_width=True,
                                 key="btn_lig"):
                        try:
                            with st.spinner("Preparando ligando..."):
                                lig_pdb = extract_single_ligand(raw_pdb, resname, resnum, chain_id)
                                prep = prepare_ligand(lig_pdb, smiles, ph, resname)
                                prep["ok"] = True
                            st.session_state["lig_prep"] = (lig_sig, prep)
                            st.session_state["last_ligand_name"] = resname
                        except Exception as e:
                            st.session_state.pop("lig_prep", None)
                            st.error(f"No se pudo preparar el ligando: {e}")

                    stored = st.session_state.get("lig_prep")
                    if stored and stored[0] == lig_sig:
                        render_ligand_results([stored[1]], pdb_name, "cry", kind="cry")
                    elif stored:
                        st.caption("Los ajustes han cambiado: pulsa de nuevo «Preparar ligando».")

                else:
                    ext_data = None
                    if source.startswith("Archivo"):
                        up_lig = st.file_uploader("Ligando(s) externo(s):",
                                                  type=["sdf", "sd", "mol", "mol2", "pdb", "smi", "txt"],
                                                  key="ext_lig_file")
                        if up_lig is not None:
                            ext_data = (up_lig.name, up_lig.getvalue())
                        st.caption("Un SDF/MOL2/SMI puede contener varios ligandos (se preparan todos).")
                    else:
                        txt = st.text_area("SMILES (uno por línea; nombre opcional tras un espacio):",
                                           value="CC(=O)Oc1ccccc1C(=O)O aspirina", height=110, key="ext_smiles")
                        if txt.strip():
                            ext_data = ("pegado.smi", txt.encode("utf-8"))

                    if ext_data:
                        ext_sig = (hashlib.md5(ext_data[1]).hexdigest(), ext_data[0], ph)
                        if st.button("⚗️ Preparar ligando(s)", type="primary", use_container_width=True,
                                     key="btn_ext"):
                            try:
                                mols = read_external_ligands(*ext_data)
                                if not mols:
                                    raise ValueError("No se encontró ninguna molécula en la entrada.")
                                results = []
                                bar = st.progress(0.0)
                                for k, (nm, om) in enumerate(mols):
                                    try:
                                        r_ = prepare_external_ligand(om, nm, ph)
                                        r_["ok"] = True
                                    except Exception as e:
                                        r_ = {"name": nm, "ok": False, "error": str(e)[:150]}
                                    results.append(r_)
                                    bar.progress((k + 1) / len(mols))
                                bar.empty()
                                st.session_state["ext_prep"] = (ext_sig, results)
                                good = [x for x in results if x.get("ok")]
                                if good:
                                    st.session_state["last_ligand_name"] = good[0]["name"]
                            except Exception as e:
                                st.session_state.pop("ext_prep", None)
                                st.error(f"No se pudieron leer los ligandos: {e}")

                        stored = st.session_state.get("ext_prep")
                        if stored and stored[0] == ext_sig:
                            render_ligand_results(stored[1], pdb_name, "ext", kind="ext")
                        elif stored:
                            st.caption("Los ajustes han cambiado: pulsa de nuevo «Preparar ligando(s)».")

            # ---------------- RECEPTOR ----------------
            with col_rec:
                st.markdown("#### 🧬 Receptor → PDBQT")
                methods = (["PDB2PQR + PROPKA (recomendado)"] if HAS_PDB2PQR else []) + \
                          ["Open Babel (rápido)"]
                method = st.radio("Método de protonación:", methods, key="rec_method")
                rec_sig = (pdb_name, tuple(selected_chains), remove_ions, ph, method)

                if st.button("⚗️ Preparar receptor", type="primary", use_container_width=True,
                             key="btn_rec"):
                    try:
                        with st.spinner("Protonando y calculando cargas (puede tardar en proteínas grandes)..."):
                            rec = prepare_receptor(result["cleaned"], ph, method)
                        st.session_state["rec_prep"] = (rec_sig, rec)
                    except Exception as e:
                        st.session_state.pop("rec_prep", None)
                        st.error(f"No se pudo preparar el receptor: {e}")

                stored = st.session_state.get("rec_prep")
                if stored and stored[0] == rec_sig:
                    rec = stored[1]
                    for msg in rec["info"]:
                        st.write(msg)
                    a, b = st.columns(2)
                    a.metric("Átomos H", rec["n_h"])
                    b.metric("Átomos en PDBQT", rec["n_atoms_pdbqt"])
                    st.caption("PDBQT rígido: los H apolares se fusionan, como espera Vina.")
                    st.download_button("📥 Receptor (.pdbqt)", rec["pdbqt"].encode(),
                                       file_name=f"{pdb_name}_receptor.pdbqt", type="primary",
                                       use_container_width=True)
                    st.download_button("📥 Receptor protonado (.pdb)", rec["pdb"].encode(),
                                       file_name=f"{pdb_name}_receptor_H.pdb",
                                       use_container_width=True)
                elif stored:
                    st.caption("Los ajustes han cambiado: pulsa de nuevo «Preparar receptor».")

            # ---------------- CAJA DE DOCKING ----------------
            st.markdown("---")
            st.markdown("#### 📦 Caja de docking")

            box_opts = lig_options + ["Proteína completa (docking ciego)"]
            box_ref = st.selectbox("Centrar la caja en:", box_opts, key=f"boxref_{pdb_name}")
            if box_ref in lig_options:
                b_res, b_rest = box_ref.split(" ", 1)
                b_num, b_chain = b_rest.split(":")
                ref_pdb = extract_single_ligand(raw_pdb, b_res, int(b_num), b_chain)
                ref_pts = pdb_coords(ref_pdb)
                ref_name = b_res
                st.caption("Caja centrada en el ligando cristalográfico elegido (sitio activo conocido). "
                           "Ajusta margen, centro y tamaño; el visor se actualiza al instante.")
            else:
                ref_pdb = ""
                ref_pts = pdb_coords(result["cleaned"])
                ref_name = "blind"
                st.caption("La caja cubre toda la proteína (docking ciego). "
                           "Puedes reducirla y moverla a mano.")
            dock_name = st.session_state.get("last_ligand_name", ref_name)

            if ref_pts:
                pad = st.slider("Margen alrededor del ligando (Å)", 4.0, 20.0, 8.0, 0.5)
                mins = [min(p[i] for p in ref_pts) for i in range(3)]
                maxs = [max(p[i] for p in ref_pts) for i in range(3)]
                d_center = [(mins[i] + maxs[i]) / 2 for i in range(3)]
                d_size = [max(20.0, (maxs[i] - mins[i]) + 2 * pad) for i in range(3)]

                lp = st.session_state.get("last_prepared")
                lp_mol, lp_diam = None, 0.0
                if lp and HAS_RDKIT:
                    try:
                        import numpy as np
                        lp_mol = Chem.RemoveHs(Chem.MolFromMolBlock(lp["sdf"], removeHs=False))
                        _xyz_lp = np.array(lp_mol.GetConformer().GetPositions())
                        lp_diam = float(np.linalg.norm(_xyz_lp[:, None, :] - _xyz_lp[None, :, :], axis=2).max())
                    except Exception:
                        lp_mol = None

                col_ctrl, col_box = st.columns([1, 2])
                with col_ctrl:
                    kp = f"box_{pdb_name}_{safe_name(ref_name)}_{pad}"
                    st.markdown("**Centro (Å)**")
                    cx, cy, cz = (st.number_input(ax, value=float(round(d_center[i], 2)), step=0.5,
                                                  format="%.2f", key=f"{kp}_c{i}")
                                  for i, ax in enumerate(("x", "y", "z")))
                    st.markdown("**Tamaño (Å)**")
                    sx, sy, sz = (st.number_input(ax, value=float(round(d_size[i], 1)), min_value=5.0,
                                                  max_value=126.0, step=1.0, format="%.1f",
                                                  key=f"{kp}_s{i}")
                                  for i, ax in enumerate(("x", "y", "z")))
                    show_inside = st.checkbox("Resaltar residuos dentro de la caja", value=True)

                    volume = sx * sy * sz
                    st.caption(f"Volumen: {volume:,.0f} Å³")
                    if volume > 27000:
                        st.warning("Vina avisa con cajas > 27 000 Å³: la búsqueda pierde eficacia. "
                                   "Reduce la caja o sube `exhaustiveness`.")

                    show_prepared = False
                    if lp_mol is not None:
                        origen = "cristalográfico" if lp["kind"] == "cry" else "externo, centrado en la caja"
                        show_prepared = st.checkbox(f"Mostrar ligando preparado: {lp['name']} ({origen})",
                                                    value=True, key=f"show_prep_{pdb_name}")
                        st.caption(f"Tamaño del ligando ≈ {lp_diam:.1f} Å (distancia máxima entre átomos).")
                        if lp_diam > min(sx, sy, sz):
                            st.warning("El ligando es más largo que la dimensión menor de la caja: "
                                       "puede no caber en algunas orientaciones. Agranda la caja.")

                inside = {}
                if show_inside:
                    lo = (cx - sx / 2, cy - sy / 2, cz - sz / 2)
                    hi = (cx + sx / 2, cy + sy / 2, cz + sz / 2)
                    for l in result["cleaned"].splitlines():
                        if l.startswith("ATOM"):
                            try:
                                x, y, z = float(l[30:38]), float(l[38:46]), float(l[46:54])
                            except ValueError:
                                continue
                            if lo[0] <= x <= hi[0] and lo[1] <= y <= hi[1] and lo[2] <= z <= hi[2]:
                                inside.setdefault(l[21], set()).add(int(l[22:26]))

                with col_box:
                    bv = py3Dmol.view(width=700, height=520)
                    bv.setBackgroundColor('#000000')
                    bv.addModel(result["cleaned"], "pdb")
                    bv.setStyle({'model': 0}, {'cartoon': {'color': '#9aa0aa', 'opacity': 0.85}})
                    for ch, resis in inside.items():
                        bv.addStyle({'model': 0, 'chain': ch, 'resi': sorted(resis)},
                                    {'stick': {'colorscheme': 'cyanCarbon', 'radius': 0.15}})
                    if ref_pdb:
                        bv.addModel(ref_pdb, "pdb")
                        bv.setStyle({'model': 1}, {'stick': {'colorscheme': 'greenCarbon', 'radius': 0.28}})
                    center_d = {'x': cx, 'y': cy, 'z': cz}
                    dims = {'w': sx, 'h': sy, 'd': sz}
                    bv.addBox({'center': center_d, 'dimensions': dims, 'color': '#19e3c5', 'opacity': 0.10})
                    bv.addBox({'center': center_d, 'dimensions': dims, 'color': '#19e3c5', 'wireframe': True})
                    if show_prepared and lp["kind"] == "ext":
                        from rdkit.Geometry import Point3D
                        m_show = Chem.Mol(lp_mol)
                        conf_s = m_show.GetConformer()
                        pos_s = conf_s.GetPositions()
                        cen_s = pos_s.mean(axis=0)
                        for i_a in range(m_show.GetNumAtoms()):
                            q = pos_s[i_a] - cen_s
                            conf_s.SetAtomPosition(i_a, Point3D(float(q[0] + cx), float(q[1] + cy),
                                                                float(q[2] + cz)))
                        bv.addModel(Chem.MolToMolBlock(m_show), "sdf")
                        bv.setStyle({'model': 2 if ref_pdb else 1},
                                    {'stick': {'colorscheme': 'magentaCarbon', 'radius': 0.25}})
                    bv.zoomTo()
                    showmol(bv, height=520, width=700)
                    n_in = sum(len(v) for v in inside.values())
                    st.caption("🟢 Ligando cristalográfico · 🟣 Ligando externo (centrado) · "
                               "🔷 Caja y residuos dentro de ella"
                               + (f" ({n_in} residuos)" if show_inside else ""))

                cfg = (f"receptor = {pdb_name}_receptor.pdbqt\n"
                       f"ligand = {dock_name}.pdbqt\n"
                       f"out = {dock_name}_out.pdbqt\n\n"
                       f"center_x = {cx:.3f}\ncenter_y = {cy:.3f}\ncenter_z = {cz:.3f}\n"
                       f"size_x = {sx:.1f}\nsize_y = {sy:.1f}\nsize_z = {sz:.1f}\n\n"
                       f"exhaustiveness = 8\nnum_modes = 9\n")
                st.session_state["box_params"] = ((cx, cy, cz), (sx, sy, sz))
                st.code(cfg, language="ini")
                st.download_button("📥 Descargar config de Vina (.txt)", cfg.encode(),
                                   file_name=f"{pdb_name}_{dock_name}_vina_config.txt")
                st.caption("Receptor y ligando cristalográfico comparten las coordenadas del PDB original. "
                           "Un ligando externo se acopla con esta misma caja. Con varios ligandos usa "
                           "`vina --batch` (sin `ligand`/`out`) y `--dir`.")
            else:
                st.caption("No se pudieron calcular coordenadas para la caja.")

    with tab4:
        st.markdown("### 🚀 Docking con AutoDock Vina")
        st.caption("Ejecuta Vina dentro de la app con el receptor, el/los ligando(s) y la caja de la pestaña "
                   "anterior. Los resultados pasan directamente a la pestaña 📈.")
        exe_path_in = ""
        if not HAS_VINA:
            exe_path_in = st.text_input(
                "Ruta del ejecutable de Vina (opcional)", key="vina_exe_path",
                help="El paquete `vina` de pip solo tiene ruedas para Linux. En Windows/macOS descarga el "
                     "ejecutable oficial de https://github.com/ccsb-scripps/AutoDock-Vina/releases "
                     "(p. ej. vina_1.2.7_win.exe) y guárdalo junto a este script o indica aquí su ruta.")
        vina_exe = None if HAS_VINA else find_vina_exe(exe_path_in)
        if not HAS_VINA and not vina_exe:
            st.error("No hay motor de Vina. Opciones: (1) `pip install vina` (solo Linux, Python ≤ 3.12), "
                     "(2) descargar el ejecutable oficial de Vina para tu sistema "
                     "(https://github.com/ccsb-scripps/AutoDock-Vina/releases; `vina_1.2.7_win.exe` en Windows) "
                     "y dejarlo junto a este script, o escribir su ruta arriba. En macOS/Linux recuerda `chmod +x`.")
        else:
            st.caption("Motor de Vina: " + ("paquete Python `vina`" if HAS_VINA else f"ejecutable `{vina_exe}`"))
            rec_stored = st.session_state.get("rec_prep")
            box_p = st.session_state.get("box_params")
            avail = {}
            cry_prep = st.session_state.get("lig_prep")
            if cry_prep and cry_prep[1].get("pdbqt"):
                avail[f"{cry_prep[1]['name']} (cristalográfico)"] = {
                    "name": cry_prep[1]["name"], "pdbqt": cry_prep[1]["pdbqt"], "ref_label": cry_prep[0][1]}
            ext_prep = st.session_state.get("ext_prep")
            if ext_prep:
                for r_ in ext_prep[1]:
                    if r_.get("ok"):
                        avail[f"{r_['name']} (externo)"] = {"name": r_["name"], "pdbqt": r_["pdbqt"],
                                                           "ref_label": None}

            ok_rec = bool(rec_stored)
            st.write(("✅" if ok_rec else "❌") + " **Receptor** preparado (pestaña 🧪)"
                     + (f" — {rec_stored[1]['n_atoms_pdbqt']} átomos en el PDBQT" if ok_rec else ""))
            st.write(("✅" if avail else "❌") + f" **Ligando(s)** preparado(s): {len(avail)}")
            st.write(("✅" if box_p else "❌") + " **Caja de docking** definida"
                     + (f" — centro ({box_p[0][0]:.1f}, {box_p[0][1]:.1f}, {box_p[0][2]:.1f}), "
                        f"tamaño ({box_p[1][0]:.0f}, {box_p[1][1]:.0f}, {box_p[1][2]:.0f}) Å" if box_p else ""))
            if ok_rec and (rec_stored[0][0] != pdb_name or rec_stored[0][1] != tuple(selected_chains)):
                st.warning("El receptor preparado corresponde a otra estructura o a otras cadenas: "
                           "vuelve a prepararlo en la pestaña 🧪.")

            if ok_rec and avail and box_p:
                try:
                    max_cpu = max(len(os.sched_getaffinity(0)), 1)    # núcleos realmente disponibles (Linux/nube)
                except AttributeError:
                    max_cpu = max(os.cpu_count() or 1, 1)
                d1, d2, d3 = st.columns(3)
                with d1:
                    scoring = st.selectbox("Función de puntuación", ["vina", "vinardo"], key="dk_sf")
                with d2:
                    exh = st.slider("Exhaustividad", 1, 64, 8, key="dk_exh",
                                    help="Más alta = búsqueda más profunda y más lenta.")
                with d3:
                    npose = st.slider("Nº de poses", 1, 20, 9, key="dk_np")
                d4, d5, d6 = st.columns(3)
                with d4:
                    erange = st.slider("Rango de energía (kcal/mol)", 1.0, 10.0, 3.0, 0.5, key="dk_er")
                with d5:
                    seed = st.number_input("Semilla aleatoria", 0, 999999, 42, key="dk_seed")
                with d6:
                    if max_cpu > 1:
                        cpu = st.slider("Núcleos de CPU", 1, max_cpu, min(4, max_cpu), key="dk_cpu")
                    else:
                        cpu = 1
                        st.caption("Núcleos de CPU: 1 (es lo único disponible en este equipo)")
                chosen = st.multiselect("Ligando(s) a acoplar:", list(avail), default=list(avail), key="dk_lig")
                st.caption("Referencia: con exhaustividad 8 y 4 núcleos, el redocking de un ligando de 6 torsiones "
                           "en una proteína de ~5000 átomos tardó ≈ 25 s.")

                if st.button("🚀 Ejecutar docking", type="primary", use_container_width=True, key="dk_run",
                             disabled=not chosen):
                    out_items = {}
                    bar = st.progress(0.0, text="Preparando...")
                    for k, label in enumerate(chosen):
                        d_ = avail[label]
                        bar.progress(k / len(chosen), text=f"Acoplando {d_['name']} ({k + 1}/{len(chosen)})...")
                        t0 = time.time()
                        try:
                            pq, en = run_vina(rec_stored[1]["pdbqt"], d_["pdbqt"], box_p[0], box_p[1],
                                              exh, npose, erange, seed, cpu, scoring, exe=vina_exe)
                            out_items[d_["name"]] = {"ok": True, "pdbqt": pq, "energies": en,
                                                     "seconds": time.time() - t0, "ref_label": d_["ref_label"]}
                        except Exception as e:
                            out_items[d_["name"]] = {"ok": False, "error": str(e)[:300]}
                    bar.progress(1.0, text="Docking completado")
                    st.session_state["dock_results"] = {
                        "items": out_items, "pdb_name": pdb_name,
                        "params": f"{scoring} · exhaustividad {exh} · {npose} poses · semilla {seed}"}

            dock_store = st.session_state.get("dock_results")
            if dock_store:
                st.markdown("---")
                st.markdown("#### Resultados")
                st.caption(f"Parámetros: {dock_store['params']}")
                rows, rmsd_by_lig = [], {}
                for nm, it_ in dock_store["items"].items():
                    if not it_.get("ok"):
                        rows.append({"Ligando": nm, "Mejor afinidad (kcal/mol)": None, "Poses": 0,
                                     "Tiempo (s)": None, "RMSD mejor pose (Å)": None,
                                     "Estado": "❌ " + it_.get("error", "")})
                        continue
                    best_rmsd = best_any = None
                    if it_["ref_label"] and HAS_RDKIT:
                        try:
                            rr_, rrest_ = it_["ref_label"].split(" ", 1)
                            rn2, rc2 = rrest_.split(":")
                            ref_m = ligand_mol_for_map(extract_single_ligand(raw_pdb, rr_, int(rn2), rc2),
                                                       ensure_smiles(pdb_name, it_["ref_label"], rr_))
                            poses_d, _, _ = load_poses("out.pdbqt", it_["pdbqt"], "")
                            vals = [ligand_rmsd(pose_mol(p_["molblock"]), ref_m)[0] for p_ in poses_d]
                            rmsd_by_lig[nm] = vals
                            best_rmsd = vals[0] if vals and vals[0] is not None else None
                            best_any = min([v for v in vals if v is not None], default=None)
                        except Exception:
                            pass
                    rows.append({"Ligando": nm,
                                 "Mejor afinidad (kcal/mol)": it_["energies"][0][0] if it_["energies"] else None,
                                 "Poses": len(it_["energies"]), "Tiempo (s)": round(it_["seconds"], 1),
                                 "RMSD mejor pose (Å)": None if best_rmsd is None else round(best_rmsd, 2),
                                 "Mejor RMSD de todas (Å)": None if best_any is None else round(best_any, 2),
                                 "Estado": "✅" if best_rmsd is None or best_rmsd <= 2.0 else "⚠️ > 2 Å"})
                st.dataframe(rows, hide_index=True, use_container_width=True)
                if rmsd_by_lig:
                    st.caption("RMSD frente al ligando cristalográfico (redocking, átomos pesados, sin alinear, "
                               "con simetría). ≤ 2 Å se considera una reproducción correcta de la pose experimental.")
                good = {nm: it_ for nm, it_ in dock_store["items"].items() if it_.get("ok")}
                if good:
                    pick = st.selectbox("Ver poses de:", list(good), key="dk_pick")
                    it_ = good[pick]
                    vals = rmsd_by_lig.get(pick)
                    st.dataframe([{"Pose": k + 1, "Afinidad (kcal/mol)": e[0],
                                   **({"RMSD vs cristal (Å)": None if vals[k] is None else round(vals[k], 2),
                                       "≤ 2 Å": "✅" if vals[k] is not None and vals[k] <= 2.0 else "❌"}
                                      if vals and k < len(vals) else {})}
                                  for k, e in enumerate(it_["energies"])],
                                 hide_index=True, use_container_width=True)
                    q1, q2 = st.columns(2)
                    with q1:
                        st.download_button("📥 Poses de este ligando (.pdbqt)", it_["pdbqt"].encode(),
                                           file_name=f"{pick}_out.pdbqt", type="primary",
                                           use_container_width=True, key="dk_dl1")
                    with q2:
                        st.download_button("📦 Todos los resultados (.zip)",
                                           zip_files({f"{n}_out.pdbqt": v["pdbqt"] for n, v in good.items()}),
                                           file_name=f"{pdb_name}_docking.zip", mime="application/zip",
                                           use_container_width=True, key="dk_dl2")
                    st.info("Analiza las interacciones de estas poses en la pestaña 📈 «Resultados de Docking» "
                            "(elige «Docking ejecutado en la app»).")

    with tab5:
        st.markdown("### 📈 Análisis de resultados de docking")
        st.caption("Carga las poses (Vina `.pdbqt`, SDF o PDB multimodelo) y analiza sus interacciones "
                   "con el receptor: visor 3D, mapa 2D y resumen de todas las poses.")
        if not (HAS_OB and HAS_RDKIT):
            st.error("Esta pestaña necesita Open Babel y RDKit: `pip install openbabel-wheel rdkit` "
                     "(y `meeko gemmi` para leer bien los PDBQT de Vina).")
        else:
            rec_choices = ["Proteína limpia del visor"]
            prepared_rec = st.session_state.get("rec_prep")
            if prepared_rec:
                rec_choices.append("Receptor protonado preparado (con H)")
            rec_choices.append("Subir otro PDB")
            rec_src = st.radio("Receptor:", rec_choices, horizontal=True, key="res_rec_src")
            protein_pdb = None
            if rec_src == "Subir otro PDB":
                up_rec = st.file_uploader("Receptor (.pdb)", type=["pdb"], key="res_rec_file")
                if up_rec is not None:
                    protein_pdb = up_rec.getvalue().decode("utf-8", errors="ignore")
            elif rec_src.startswith("Receptor protonado"):
                protein_pdb = prepared_rec[1]["pdb"]
            else:
                protein_pdb = result["cleaned"]
            st.caption("El receptor debe estar en el mismo sistema de coordenadas que usaste en el docking.")

            dock_store5 = st.session_state.get("dock_results")
            good5 = {n: v for n, v in (dock_store5["items"] if dock_store5 else {}).items() if v.get("ok")}
            pose_src = "Subir o pegar archivo"
            if good5:
                pose_src = st.radio("Origen de las poses:", ["Docking ejecutado en la app", "Subir o pegar archivo"],
                                    horizontal=True, key="res_pose_src")
            pose_name, pose_text, smi_tpl = None, "", ""
            if pose_src == "Docking ejecutado en la app":
                app_lig = st.selectbox("Ligando acoplado:", list(good5), key="res_app_lig")
                pose_name, pose_text = "out.pdbqt", good5[app_lig]["pdbqt"]
            else:
                c_up, c_paste = st.columns(2)
                with c_up:
                    up_pose = st.file_uploader("Resultado del docking:", type=["pdbqt", "sdf", "pdb"],
                                               key="res_pose_file")
                with c_paste:
                    pasted = st.text_area("…o pega aquí el contenido", height=130, key="res_pose_paste")
                smi_tpl = st.text_input("SMILES del ligando (opcional; restaura enlaces si el archivo no los trae)",
                                        key="res_smiles")
                if up_pose is not None:
                    pose_name, pose_text = up_pose.name, up_pose.getvalue().decode("utf-8", errors="ignore")
                elif pasted.strip():
                    pose_name, pose_text = "pegado", pasted

            if protein_pdb and pose_text:
                try:
                    poses, pose_method, pose_warn = load_poses(pose_name, pose_text, smi_tpl.strip())
                except Exception as e:
                    poses, pose_method, pose_warn = [], "", []
                    st.error(f"No se pudo leer el archivo de poses: {e}")

                if poses:
                    st.caption(f"{len(poses)} poses · lectura con: {pose_method}")
                    for w in pose_warn:
                        st.warning(w)
                    ref_choice = st.selectbox("Ligando cristalográfico de referencia (superponer y calcular RMSD):",
                                              ["Ninguno"] + lig_options, key="res_ref_lig")
                    rmsds, rm_note, ref_lig_pdb = [None] * len(poses), "", ""
                    if ref_choice != "Ninguno":
                        rr, rrest = ref_choice.split(" ", 1)
                        rn_, rc_ = rrest.split(":")
                        ref_lig_pdb = extract_single_ligand(raw_pdb, rr, int(rn_), rc_)
                        ref_mol_res = ligand_mol_for_map(ref_lig_pdb, ensure_smiles(pdb_name, ref_choice, rr))
                        if ref_mol_res is not None:
                            for k_, p_ in enumerate(poses):
                                try:
                                    rmsds[k_], rm_note = ligand_rmsd(pose_mol(p_["molblock"]), ref_mol_res)
                                except Exception as e:
                                    rm_note = str(e)[:80]
                    rows5 = []
                    for k, p in enumerate(poses):
                        row = {"Pose": k + 1, "Afinidad (kcal/mol)": p["score"],
                               "RMSD l.b.": p["rmsd_lb"], "RMSD u.b.": p["rmsd_ub"]}
                        if ref_choice != "Ninguno":
                            row["RMSD vs cristal (Å)"] = None if rmsds[k] is None else round(rmsds[k], 2)
                            row["≤ 2 Å"] = "✅" if rmsds[k] is not None and rmsds[k] <= 2.0 else "❌"
                        rows5.append(row)
                    st.dataframe(rows5, hide_index=True, use_container_width=True)
                    if ref_choice != "Ninguno":
                        vals_ok = [(k, v) for k, v in enumerate(rmsds) if v is not None]
                        if vals_ok:
                            top = rmsds[0]
                            bk, bv_ = min(vals_ok, key=lambda t: t[1])
                            if top is not None:
                                verdict = ("✅ redocking exitoso" if top <= 2.0 else
                                           "⚠️ la pose de mejor afinidad no reproduce la experimental")
                                st.success(f"Pose 1 (mejor afinidad): RMSD = {top:.2f} Å — {verdict}. "
                                           f"Mejor RMSD entre todas: {bv_:.2f} Å (pose {bk + 1}). "
                                           f"Método: {rm_note}.")
                            else:
                                st.info(f"Mejor RMSD entre todas: {bv_:.2f} Å (pose {bk + 1}). "
                                        f"Método: {rm_note}.")
                        else:
                            st.warning(f"No se pudo calcular el RMSD: {rm_note or 'sin correspondencia entre átomos'}. "
                                       "Comprueba que la referencia sea el mismo ligando y que el receptor "
                                       "esté en el mismo sistema de coordenadas.")

                    def _pose_label(k):
                        sc = poses[k]["score"]
                        return f"Pose {k + 1}" + (f" · {sc:.2f} kcal/mol" if sc is not None else "")

                    r1, r2 = st.columns([2, 1])
                    with r1:
                        pose_idx = st.selectbox("Pose a analizar:", list(range(len(poses))),
                                                format_func=_pose_label, key="res_pose_idx")
                    with r2:
                        res_labels = st.checkbox("Etiquetas de residuos", value=True, key="res_labels")

                    pmol = pose_mol(poses[pose_idx]["molblock"])
                    cpx, _key = build_pose_complex(protein_pdb, pmol)
                    with st.spinner("Analizando interacciones..."):
                        p_items, p_method, p_warn = analyze_complex(cpx, "LIG", 1, "Z")
                    if p_warn:
                        st.warning(p_warn)

                    p_found = sorted({it["tipo"] for it in p_items})
                    p_chosen = st.multiselect("Tipos de interacción:", p_found, default=p_found,
                                              key=f"res_types_{pose_idx}_{len(p_items)}")
                    p_filtered = [it for it in p_items if it["tipo"] in p_chosen]
                    st.caption(f"Método: {p_method} · " +
                               " ".join(f"{INTER_EMOJI.get(t, '▫️')} {t}" for t in p_chosen)
                               + "  · 🟣 pose" + ("  · 🟢 cristalográfico" if ref_choice != "Ninguno" else ""))

                    pv = build_pose_view(protein_pdb, pmol, p_filtered, res_labels, ref_lig_pdb)
                    showmol(pv, height=560, width=700)

                    if p_filtered:
                        st.dataframe([{"Tipo": it["tipo"], "Residuo": f"{it['resname']}{it['resnum']}",
                                       "Cadena": it["chain"], "Distancia (Å)": it["dist"]}
                                      for it in sorted(p_filtered,
                                                       key=lambda t: (t["chain"], t["resnum"], t["tipo"]))],
                                     use_container_width=True, hide_index=True)
                    elif not p_items:
                        st.info("No se detectaron interacciones para esta pose.")

                    with st.expander("🗺️ Mapa 2D de interacciones", expanded=False):
                        dark_p = st.checkbox("Fondo oscuro", value=True, key="res_map_dark")
                        if p_filtered:
                            try:
                                svg_p = interaction_map_svg(pmol, p_filtered, title=_pose_label(pose_idx),
                                                            dark=dark_p)
                                st.markdown(svg_img_html(svg_p), unsafe_allow_html=True)
                                st.download_button("📥 Mapa 2D (.svg)", svg_p.encode("utf-8"),
                                                   file_name=f"pose_{pose_idx + 1}_mapa2D.svg",
                                                   mime="image/svg+xml", key="res_dl_svg")
                            except Exception as e:
                                st.warning(f"No se pudo generar el mapa 2D: {e}")
                        else:
                            st.caption("No hay interacciones que dibujar.")

                    d1, d2 = st.columns(2)
                    with d1:
                        st.download_button("📥 Pose (.sdf)", (Chem.MolToMolBlock(pmol) + "$$$$\n").encode(),
                                           file_name=f"pose_{pose_idx + 1}.sdf", use_container_width=True,
                                           key="res_dl_sdf")
                    with d2:
                        st.download_button("📥 Complejo receptor + pose (.pdb)", cpx.encode(),
                                           file_name=f"complejo_pose_{pose_idx + 1}.pdb",
                                           use_container_width=True, key="res_dl_cpx")

                    pse5 = f"res_pse_{pose_idx}_{hashlib.md5(cpx.encode()).hexdigest()[:8]}"
                    if st.button("🧩 Generar sesión de PyMOL (.pse) de esta pose", key=f"btn_{pse5}"):
                        with st.spinner("Generando la sesión con PLIP + PyMOL..."):
                            st.session_state[pse5] = export_pymol_session(
                                cpx, "LIG:Z:1", p_items, "resn LIG and chain Z and resi 1")
                    if pse5 in st.session_state:
                        pse_d, pse_n, pse_ds = st.session_state[pse5]
                        st.download_button(f"📥 Descargar {pse_n}", pse_d, file_name=f"pose_{pose_idx + 1}_{pse_n}",
                                           mime="application/octet-stream", type="primary", key=f"dl_{pse5}")
                        st.caption(pse_ds)

                    # ---- resumen de todas las poses ----
                    with st.expander("📋 Resumen de todas las poses"):
                        n_max = len(poses)
                        n_an = st.slider("Poses a analizar", 1, n_max, min(n_max, 9),
                                         key="res_n_an") if n_max > 1 else 1
                        sum_sig = (hashlib.md5(pose_text.encode()).hexdigest(),
                                   hashlib.md5(protein_pdb.encode()).hexdigest(), n_an)
                        if st.button("Analizar poses", key="res_btn_all"):
                            per_pose, rows = [], []
                            bar = st.progress(0.0, text="Analizando poses...")
                            for k in range(n_an):
                                m_k = pose_mol(poses[k]["molblock"])
                                c_k, _ = build_pose_complex(protein_pdb, m_k)
                                it_k, meth_k, _w = analyze_complex(c_k, "LIG", 1, "Z")
                                lab = f"Pose {k + 1}"
                                per_pose.append((lab, it_k))
                                cnt = Counter_types(it_k)
                                row_s = {"Pose": k + 1, "Afinidad (kcal/mol)": poses[k]["score"],
                                         "Nº interacciones": len(it_k), **cnt}
                                if ref_choice != "Ninguno" and rmsds[k] is not None:
                                    row_s["RMSD vs cristal (Å)"] = round(rmsds[k], 2)
                                rows.append(row_s)
                                bar.progress((k + 1) / n_an, text=f"Pose {k + 1}/{n_an}")
                            bar.empty()
                            st.session_state["res_summary"] = (sum_sig, per_pose, rows)
                        stored_sum = st.session_state.get("res_summary")
                        if stored_sum and stored_sum[0] == sum_sig:
                            _, per_pose, rows = stored_sum
                            st.markdown("**Interacciones por pose**")
                            st.dataframe(rows, hide_index=True, use_container_width=True)
                            matrix = interaction_matrix(per_pose)
                            if matrix:
                                st.markdown("**Residuos × poses** (HB puente H · HP hidrofóbico · SB puente salino · "
                                            "πS apilamiento π · πC π–catión · HX halógeno · WB agua · MC metal; "
                                            "`~` = aproximado)")
                                st.dataframe(matrix, hide_index=True, use_container_width=True)
                                st.bar_chart({f"{r['Residuo']}:{r['Cadena']}": r["Nº poses"] for r in matrix[:15]})
                            else:
                                st.info("No se detectaron interacciones en las poses analizadas.")
                        elif stored_sum:
                            st.caption("Los ajustes han cambiado: pulsa de nuevo «Analizar poses».")
                else:
                    st.warning("No se encontraron poses en el archivo.")
            else:
                st.info("Sube el archivo de resultados (por ejemplo `out.pdbqt` de Vina) o pega su contenido.")
