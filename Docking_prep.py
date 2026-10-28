# -*- coding: utf-8 -*-
"""
Created on Wed Oct 28 13:15:12 2026

@author: Carlos Fernandez
"""

"""
BioDock Studio: Protein Parser, Inspector & Docking Prep
v11 — v10 + exportación a GROMACS (topologías proteína/ligando, .mdp, script del flujo completo y
construcción/validación opcional del sistema si GROMACS está instalado).
"""

import base64
import hashlib
import html as html_lib
import importlib.util
import glob
import io
import json
import math
import os
import platform
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
        --accent2:#4d7cff; --text:#e8eaed; --muted:#c0c7d0; }

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
    

def extract_molecule_name(pdb_text, default_code="PROTEIN"):
    """
    Busca en el texto del PDB las líneas COMPND y extrae el valor del campo MOLECULE.
    Ejemplo: COMPND   2 MOLECULE: HIV-1 PROTEASE; -> 'HIV-1 PROTEASE'
    """
    if not pdb_text:
        return default_code.upper()
    
    # Buscar el patrón MOLECULE: <nombre>; dentro de los registros COMPND
    match = re.search(r"COMPND(?:\s+\d+)?\s+MOLECULE:\s*([^;\n]+);?", pdb_text, re.IGNORECASE)
    if match:
        # Limpiar espacios extra y comillas si las hubiera
        name = match.group(1).strip().strip('"\'')
        return name
    
    return default_code.upper()


def get_structure_heterogens(pdb_text):
    """Escanea el PDB y devuelve una lista ordenada de nombres de residuos no estándar (HETATM que no son agua)."""
    structure = PDB.PDBParser(QUIET=True).get_structure("het", io.StringIO(pdb_text))
    hets = set()
    for model in structure:
        for chain in model:
            for res in chain:
                if not is_water(res) and not is_amino_acid(res):
                    hets.add(res.get_resname().strip())
    return sorted(list(hets))


def is_water(res):
    return res.get_id()[0] == 'W' or res.get_resname().strip() in WATER_NAMES


def is_amino_acid(res):
    return res.get_resname().strip() in AA_3TO1 and not is_water(res)


def is_ion(res):
    return res.get_resname().strip() in COMMON_IONS and not is_amino_acid(res)


def is_ligand(res, cofactors=()):
    resname = res.get_resname().strip()
    if resname in cofactors:
        return False
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
def process_pdb_structure(pdb_text, selected_chains, remove_water, remove_ions, retained_heterogens=()):
    selected = set(selected_chains)
    retained_set = set(retained_heterogens)
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
            resname = res.get_resname().strip()
            if is_water(res):
                water_count += 1
            elif is_amino_acid(res):
                num = res.get_id()[1]
                residues.append((AA_3TO1[resname], num,
                                 residue_ss(segments, chain.id, num)))
            elif is_ion(res):
                ion_count += 1
            
            # Si el residuo está en la lista de retenidos, se queda en la proteína
            if resname in retained_set:
                pass
            elif is_ligand(res, retained_set):
                ligands.append((resname, res.get_id()[1], chain.id))
        if residues:
            chains_data[chain.id] = residues

    class ProteinSelect(PDB.Select):
        def accept_model(self, m):
            return m.id == model_id

        def accept_chain(self, chain):
            return chain.id in selected

        def accept_residue(self, residue):
            resname = residue.get_resname().strip()
            if is_amino_acid(residue):
                return True
            if resname in retained_set:
                return True  # <-- Conserva cofactores e iones específicos aunque se filtren los demás
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
            resname = residue.get_resname().strip()
            if resname in retained_set:
                return False  # <-- Evita que los retenidos se extraigan como ligandos móviles
            return is_ligand(residue, retained_set)

    
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
                     rmin, rmax, surf_opacity, show_pocket, outline, spin, inter=None, cofactors=()):
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
        
        # ESTILO PARA LOS COFACTORES CONSERVADOS (Sticks) ---
        if cofactors:
            view.addStyle({'model': prot_idx, 'resn': list(cofactors)},
                          {'stick': {'colorscheme': 'yellowCarbon', 'radius': 0.18},
                           'sphere': {'colorscheme': 'yellowCarbon', 'scale': 0.25}})

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

    cl = parse_cluster_remarks(text) if fmt == "pdbqt" else []
    poses = []
    for k, m in enumerate(mols):
        sc = scores[k] if k < len(scores) else (None, None, None)
        cn, rn = cl[k] if k < len(cl) else (None, None)
        poses.append({"molblock": _mol_block(m), "score": sc[0], "rmsd_lb": sc[1], "rmsd_ub": sc[2],
                      "cluster_n": cn, "runs": rn})
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


def build_pose_view(protein_pdb, mol, items, labels=True, ref_ligand_pdb="", cofactors=()):
    view = py3Dmol.view(width=700, height=560)
    view.setBackgroundColor('#000000')
    view.addModel(protein_pdb, "pdb")
    view.setStyle({'model': 0}, {'cartoon': {'color': '#9aa0aa', 'opacity': 0.85}})
    
    # PINTAR COFACTORES EN STICKS EN LA VISTA DE POSES
    if cofactors:
        view.addStyle({'model': 0, 'resn': list(cofactors)},
                      {'stick': {'colorscheme': 'yellowCarbon', 'radius': 0.20}})
        
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


def _vina_importable():
    try:
        from vina import Vina  # noqa: F401
        return True
    except Exception:
        return False


HAS_VINA = _vina_importable()


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
        gname = {"Puente de hidrógeno": "HBonds", "Contacto hidrofóbico": "Hydrophobic",
                 "Puente salino": "SaltBridges", "Apilamiento π": "PiStacking",
                 "Interacción π–catión": "PiCation", "Enlace de halógeno": "HalogenBonds",
                 "Puente de agua": "WaterBridges", "Coordinación metálica": "MetalComplexes",
                 "Contacto polar (aprox.)": "PolarContacts",
                 "Contacto hidrofóbico (aprox.)": "HydrophobicApprox"}.get(it["tipo"], "Interactions")
        for j, (p1, p2) in enumerate(it["segs"]):
            a, b = f"ps_{k}_{j}a", f"ps_{k}_{j}b"
            lines += [f"pseudoatom {a}, pos=[{p1[0]:.3f},{p1[1]:.3f},{p1[2]:.3f}]",
                      f"pseudoatom {b}, pos=[{p2[0]:.3f},{p2[1]:.3f},{p2[2]:.3f}]",
                      f"distance {gname}, {a}, {b}", f"color {col}, {gname}",
                      f"set dash_color, {col}, {gname}", f"hide everything, {a} or {b}"]
    lines += ["set dash_gap, 0.3", "set dash_radius, 0.07", "set label_color, white", "set label_size, 14",
              "zoom ligando, 8", f"save {session_name}"]
    return "\n".join(lines) + "\n"


_PYMOL_CHECK = {}


def pymol_module_status():
    """(ok, mensaje): ¿se puede importar realmente el motor de PyMOL? (find_spec no basta: en Windows
    el paquete puede estar instalado y fallar al cargar sus DLL)."""
    if "r" in _PYMOL_CHECK:
        return _PYMOL_CHECK["r"]
    if importlib.util.find_spec("pymol") is None:
        res = (False, "el módulo `pymol` no está instalado (`pip install pymol-open-source`).")
    else:
        try:
            r = subprocess.run([sys.executable, "-c", "import pymol._cmd"], capture_output=True,
                               text=True, timeout=90)
            if r.returncode == 0:
                res = (True, "")
            else:
                lines = (r.stderr or r.stdout).strip().splitlines()
                err = lines[-1][:220] if lines else "error desconocido"
                hint = ""
                if "DLL" in err:
                    hint = (" En Windows suele faltar «Microsoft Visual C++ Redistributable» (x64, 2015-2022) "
                            "o hay una incompatibilidad con tu versión de Python; la vía más fiable es un entorno "
                            "conda-forge (`conda install -c conda-forge pymol-open-source`).")
                res = (False, f"el módulo `pymol` está instalado pero no carga ({err}).{hint}")
        except Exception as e:
            res = (False, f"no se pudo comprobar PyMOL ({str(e)[:120]}).")
    _PYMOL_CHECK["r"] = res
    return res


def find_pymol_exe(custom=""):
    """Ejecutable de PyMOL instalado (para crear el .pse con un script cuando el módulo Python no carga)."""
    cands = []
    if custom and custom.strip():
        cands.append(custom.strip().strip('"'))
    cands += ["pymol", "pymol.exe", "PyMOLWin.exe"]
    if os.name == "nt":
        for base in (os.environ.get("ProgramFiles"), os.environ.get("ProgramFiles(x86)"),
                     os.environ.get("LOCALAPPDATA"), os.environ.get("ProgramData")):
            if base:
                cands += glob.glob(os.path.join(base, "PyMOL*", "PyMOLWin.exe"))
                cands += glob.glob(os.path.join(base, "*", "PyMOL*", "PyMOLWin.exe"))
    for c in cands:
        path = c if os.path.isfile(c) else shutil.which(c)
        if path:
            return path
    return None


def pymol_ui(prefix):
    """Estado de PyMOL en la interfaz; devuelve la ruta opcional del ejecutable."""
    ok, msg = pymol_module_status()
    if ok:
        st.caption("PyMOL: ✅ módulo Python operativo (la sesión .pse la crea PLIP).")
        return ""
    st.caption("PyMOL: ⚠️ " + msg)
    custom = st.text_input("Ruta del ejecutable de PyMOL (opcional)", key=f"{prefix}_pymol_exe",
                           help="Si tienes PyMOL instalado (p. ej. PyMOLWin.exe), la app lo usa para crear el .pse "
                                "a partir de un script con las interacciones (colores y estilo propios de la app).")
    exe = find_pymol_exe(custom)
    st.caption(f"Ejecutable de PyMOL: {exe}" if exe else
               "No se encontró ningún ejecutable de PyMOL: se entregará un ZIP con el script .pml.")
    return custom


def export_pymol_session(complex_pdb, key, items, lig_sel, pymol_exe=""):
    """(bytes, nombre, descripción). 1) .pse real de PLIP (módulo pymol); 2) .pse creado con el ejecutable de
    PyMOL a partir de un script; 3) ZIP con complejo + .pml."""
    notes = []
    mod_ok, mod_msg = pymol_module_status()
    if mod_ok and HAS_PLIP:
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
                notes.append("PLIP/PyMOL no pudo crear el .pse: " + (proc.stderr or proc.stdout).strip()[-150:])
        except Exception as e:
            notes.append(f"PLIP/PyMOL falló: {str(e)[:150]}")
    elif not mod_ok:
        notes.append(mod_msg[0].upper() + mod_msg[1:])

    pml = build_pml(items, lig_sel)
    exe = find_pymol_exe(pymol_exe)
    if exe:
        try:
            with tempfile.TemporaryDirectory() as td:
                with open(os.path.join(td, "complejo.pdb"), "w") as f:
                    f.write(complex_pdb)
                with open(os.path.join(td, "interacciones.pml"), "w") as f:
                    f.write(pml)
                subprocess.run([exe, "-cq", "interacciones.pml"], cwd=td, capture_output=True, text=True,
                               timeout=300)
                pse_path = os.path.join(td, "interacciones.pse")
                if os.path.exists(pse_path):
                    with open(pse_path, "rb") as f:
                        return (f.read(), "interacciones.pse",
                                " ".join(notes) + " Sesión creada con tu PyMOL a partir de un script con las "
                                "interacciones de PLIP (colores y estilo propios de la app, no los de PLIP).")
                notes.append("El ejecutable de PyMOL no generó el .pse.")
        except Exception as e:
            notes.append(f"Fallo al ejecutar PyMOL: {str(e)[:120]}")
    elif mod_ok is False:
        notes.append("No se encontró un ejecutable de PyMOL.")

    zbytes = zip_files({"complejo.pdb": complex_pdb, "interacciones.pml": pml,
                        "LEEME.txt": "Abre PyMOL en esta carpeta y ejecuta:  pymol interacciones.pml\n"
                                     "Se crea interacciones.pse con el complejo y las interacciones.\n"})
    return (zbytes, "interacciones_pymol.zip",
            " ".join(notes) + " Se entrega un ZIP: ejecuta `pymol interacciones.pml` para obtener el .pse.")


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


VINA_VERSION = "1.2.7"


def _vina_cache_dir():
    return os.path.join(tempfile.gettempdir(), "biodock_vina")


def vina_asset_name():
    """Nombre del ejecutable oficial de Vina para este sistema (None si no hay)."""
    system, mach = platform.system(), platform.machine().lower()
    arm = mach in ("arm64", "aarch64")
    if system == "Windows":
        return f"vina_{VINA_VERSION}_win.exe"
    if system == "Darwin":
        return f"vina_{VINA_VERSION}_mac_{'aarch64' if arm else 'x86_64'}"
    if system == "Linux":
        return f"vina_{VINA_VERSION}_linux_{'aarch64' if arm else 'x86_64'}"
    return None


def download_vina_exe():
    """Descarga el ejecutable oficial de Vina (~4 MB) desde la release de GitHub a una carpeta temporal."""
    name = vina_asset_name()
    if not name:
        raise RuntimeError(f"Sistema no soportado: {platform.system()} {platform.machine()}")
    d = _vina_cache_dir()
    os.makedirs(d, exist_ok=True)
    dest = os.path.join(d, name)
    url = f"https://github.com/ccsb-scripps/AutoDock-Vina/releases/download/v{VINA_VERSION}/{name}"
    req = urllib.request.Request(url, headers={"User-Agent": "Mozilla/5.0"})
    tmp = dest + ".part"
    with urllib.request.urlopen(req, timeout=120) as r, open(tmp, "wb") as f:
        shutil.copyfileobj(r, f)
    os.replace(tmp, dest)
    if os.name != "nt":
        os.chmod(dest, 0o755)
    out = subprocess.run([dest, "--version"], capture_output=True, text=True, timeout=30)
    if "vina" not in (out.stdout + out.stderr).lower():
        raise RuntimeError("El archivo descargado no es un ejecutable de Vina válido.")
    return dest


def _ensure_exec(path):
    """Garantiza permiso de ejecución (en la nube un archivo subido por la web de GitHub no lo trae)."""
    if os.name == "nt" or os.access(path, os.X_OK):
        return path
    try:
        os.chmod(path, 0o755)
        if os.access(path, os.X_OK):
            return path
    except OSError:
        pass
    try:
        d = _vina_cache_dir()
        os.makedirs(d, exist_ok=True)
        dest = os.path.join(d, os.path.basename(path))
        shutil.copy(path, dest)
        os.chmod(dest, 0o755)
        return dest
    except Exception:
        return None


def find_vina_exe(custom=""):
    """Localiza el ejecutable oficial de Vina: ruta indicada, junto al script, carpeta de descarga o PATH."""
    cands = []
    if custom and custom.strip():
        cands.append(custom.strip().strip('"'))
    dirs = [os.getcwd(), _vina_cache_dir()]
    try:
        dirs.append(os.path.dirname(os.path.abspath(__file__)))
    except NameError:
        pass
    skip = (".py", ".pdbqt", ".txt", ".zip", ".md", ".log", ".pdb", ".sdf", ".part", ".toml")
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
        path = _ensure_exec(path)
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



# --- Ki ESTIMADA, RUNS INDEPENDIENTES Y CLÚSTERES (estilo AutoDock4) ---

R_KCAL = 0.0019872041    # kcal/(mol·K)


def ki_from_affinity(dg, temp=298.15):
    """Constante de inhibición estimada (M): Ki = exp(ΔG / RT), con ΔG en kcal/mol."""
    return math.exp(float(dg) / (R_KCAL * temp))


def format_ki(dg, temp=298.15):
    if dg is None:
        return "—"
    ki = ki_from_affinity(dg, temp)
    if ki >= 1.0:
        return "> 1 M"
    for unit, f in (("mM", 1e-3), ("µM", 1e-6), ("nM", 1e-9)):
        if ki >= f:
            return f"{ki / f:.3g} {unit}"
    return f"{ki / 1e-12:.3g} pM"


_CLUSTER_RE = re.compile(r"REMARK CLUSTER\s+(\d+)\s+(\d+)")
_AD_EL = {"A": "C", "NA": "N", "NS": "N", "OA": "O", "OS": "O", "SA": "S", "HD": "H", "HS": "H"}


def parse_cluster_remarks(text):
    parts = re.split(r"^MODEL", text, flags=re.M)
    chunks = parts[1:] if len(parts) > 1 else [text]
    out = []
    for c in chunks:
        m = _CLUSTER_RE.search(c)
        out.append((int(m.group(1)), int(m.group(2))) if m else (None, None))
    return out


def split_pdbqt_models(text):
    models = re.findall(r"MODEL\s+\d+\s*\n(.*?)ENDMDL", text, flags=re.S)
    return models if models else [text]


def _model_heavy(block):
    """[(elemento, x, y, z)] de los átomos pesados de un bloque PDBQT."""
    atoms = []
    for l in block.splitlines():
        if l.startswith(("ATOM", "HETATM")):
            t = l.split()[-1]
            el = _AD_EL.get(t, t.capitalize())
            if el.startswith("H") and len(el) == 1:
                continue
            atoms.append((el, float(l[30:38]), float(l[38:46]), float(l[46:54])))
    return atoms


def rmsd_bounds(a, b):
    """(l.b., u.b.) de la pose `a` respecto a la referencia `b`. u.b.: átomo con átomo (mismo orden);
    l.b.: cada átomo con el átomo más cercano del mismo elemento. Ambos en Å, solo átomos pesados."""
    import numpy as np
    if not a or len(a) != len(b):
        return None, None
    A = np.array([x[1:] for x in a])
    B = np.array([x[1:] for x in b])
    ub = float(np.sqrt(((A - B) ** 2).sum(axis=1).mean()))
    tot = 0.0
    for i, (el, *_rest) in enumerate(a):
        idx = [j for j, y in enumerate(b) if y[0] == el] or list(range(len(b)))
        tot += float(((B[idx] - A[i]) ** 2).sum(axis=1).min())
    return float(np.sqrt(tot / len(a))), ub


_VINA_WORKER = r"""
import sys
rec, lig, out = sys.argv[1:4]
cx, cy, cz, sx, sy, sz = [float(x) for x in sys.argv[4:10]]
seed, scoring = int(sys.argv[10]), sys.argv[11]
from vina import Vina
v = Vina(sf_name=scoring, seed=seed, cpu=1, verbosity=0)
v.set_receptor(rec)
v.set_ligand_from_file(lig)
v.compute_vina_maps(center=[cx, cy, cz], box_size=[sx, sy, sz])
v.dock(exhaustiveness=1, n_poses=1)
v.write_poses(out, n_poses=1, overwrite=True)
"""


def run_vina_independent(receptor_pdbqt, ligand_pdbqt, center, size, n_runs, seed, cpu, scoring,
                         exe=None, tol=2.0, progress=None):
    """N runs independientes (exhaustividad 1 y semillas distintas), agrupadas en clústeres por RMSD como
    en AutoDock4. Devuelve (pdbqt con el representante de cada clúster, energías [[E, l.b., u.b.]],
    resumen de clústeres, energías de todas las runs)."""
    from concurrent.futures import ThreadPoolExecutor, as_completed
    n_runs = max(int(n_runs), 1)
    with tempfile.TemporaryDirectory() as td:
        rp, lp = os.path.join(td, "receptor.pdbqt"), os.path.join(td, "ligand.pdbqt")
        with open(rp, "w") as f:
            f.write(receptor_pdbqt)
        with open(lp, "w") as f:
            f.write(ligand_pdbqt)

        def one(i):
            out = os.path.join(td, f"run_{i}.pdbqt")
            sd = int(seed) + i
            if exe:
                cmd = [exe, "--receptor", rp, "--ligand", lp, "--out", out,
                       "--center_x", f"{center[0]:.3f}", "--center_y", f"{center[1]:.3f}",
                       "--center_z", f"{center[2]:.3f}", "--size_x", f"{size[0]:.1f}",
                       "--size_y", f"{size[1]:.1f}", "--size_z", f"{size[2]:.1f}",
                       "--exhaustiveness", "1", "--num_modes", "1", "--seed", str(sd), "--cpu", "1",
                       "--scoring", scoring]
            else:
                cmd = [sys.executable, "-c", _VINA_WORKER, rp, lp, out] + \
                      [f"{x:.3f}" for x in list(center) + list(size)] + [str(sd), scoring]
            proc = subprocess.run(cmd, capture_output=True, text=True, timeout=3600)
            if proc.returncode != 0 or not os.path.exists(out):
                raise RuntimeError((proc.stderr or proc.stdout or "Vina terminó sin resultados").strip()[-300:])
            with open(out) as f:
                return f.read()

        texts = []
        with ThreadPoolExecutor(max_workers=max(1, min(int(cpu), n_runs))) as pool:
            futs = [pool.submit(one, i) for i in range(n_runs)]
            for k, fut in enumerate(as_completed(futs), 1):
                texts.append(fut.result())
                if progress:
                    progress(k, n_runs)

    runs = []
    for t in texts:
        m = _SCORE_RE.search(t)
        blk = split_pdbqt_models(t)[0]
        runs.append({"e": float(m.group(1)) if m else 0.0, "block": blk, "atoms": _model_heavy(blk)})
    order = sorted(range(len(runs)), key=lambda i: runs[i]["e"])
    clusters = []
    for i in order:
        for c in clusters:
            _lb, ub = rmsd_bounds(runs[i]["atoms"], runs[c["leader"]]["atoms"])
            if ub is not None and ub <= tol:
                c["members"].append(i)
                break
        else:
            clusters.append({"leader": i, "members": [i]})

    best_atoms = runs[clusters[0]["leader"]]["atoms"]
    models, energies, summary = [], [], []
    for r, c in enumerate(clusters, 1):
        ld = runs[c["leader"]]
        lb, ub = rmsd_bounds(ld["atoms"], best_atoms)
        lb, ub = lb or 0.0, ub or 0.0
        head = (f"REMARK VINA RESULT: {ld['e']:10.3f} {lb:10.3f} {ub:10.3f}\n"
                f"REMARK CLUSTER {len(c['members'])} {n_runs}\n")
        body = re.sub(r"REMARK VINA RESULT:.*\n", head, ld["block"], count=1)
        models.append(f"MODEL {r}\n{body.rstrip()}\nENDMDL\n")
        energies.append([ld["e"], round(lb, 3), round(ub, 3)])
        mem = [runs[j]["e"] for j in c["members"]]
        summary.append({"rank": r, "best": ld["e"], "mean": sum(mem) / len(mem), "n": len(mem)})
    return "".join(models), energies, summary, sorted(r_["e"] for r_ in runs)


def pose_rows(poses, rmsds=None):
    """Filas de la tabla de poses: afinidad, Ki estimada, clúster/runs y RMSD."""
    rows = []
    for k, p in enumerate(poses):
        sc = p["score"]
        cn, rn = p.get("cluster_n"), p.get("runs")
        row = {"Pose": k + 1, "Afinidad (kcal/mol)": sc, "Ki estimada": format_ki(sc),
               "Nº conf. en clúster": "—" if cn is None else str(cn),
               "Runs totales": "—" if rn is None else str(rn),
               "RMSD l.b. (vs pose 1)": p["rmsd_lb"], "RMSD u.b. (vs pose 1)": p["rmsd_ub"]}
        if rmsds is not None:
            v = rmsds[k] if k < len(rmsds) else None
            row["RMSD vs cristal (Å)"] = None if v is None else round(v, 2)
            row["≤ 2 Å"] = "✅" if v is not None and v <= 2.0 else "❌"
        rows.append(row)
    return rows



# --- EXPORTACIÓN A GROMACS ---

# Este código viaja tal cual dentro del paquete (prepare_system.py): la app y el script del usuario usan la misma lógica.
_GMX_HELPERS = r"""
import re
import sys


def amberize_pdb(pdb_text):
    '''Renombra residuos según los H presentes (salida de PDB2PQR/PROPKA) con la nomenclatura AMBER de GROMACS.'''
    res_atoms = {}
    for l in pdb_text.splitlines():
        if l.startswith(("ATOM", "HETATM")):
            key = (l[21], l[22:27])
            res_atoms.setdefault(key, (l[17:20].strip(), set()))[1].add(l[12:16].strip())
    new_name, changes = {}, {}
    for key, (rn, atoms) in res_atoms.items():
        nn = rn
        if rn == "HIS":
            hd, he = "HD1" in atoms, "HE2" in atoms
            nn = "HIP" if hd and he else "HID" if hd else "HIE" if he else "HIS"
        elif rn == "ASP" and "HD2" in atoms:
            nn = "ASH"
        elif rn == "GLU" and "HE2" in atoms:
            nn = "GLH"
        elif rn == "LYS" and {"HZ1", "HZ2"} <= atoms and "HZ3" not in atoms:
            nn = "LYN"
        if nn != rn:
            new_name[key] = nn
            changes[nn] = changes.get(nn, 0) + 1
    out = []
    for l in pdb_text.splitlines():
        if l.startswith(("ATOM", "HETATM")):
            key = (l[21], l[22:27])
            if key in new_name:
                l = l[:17] + f"{new_name[key]:>3s}" + l[20:]
        out.append(l)
    return "\n".join(out) + "\n", changes


def parse_gro(text):
    lines = text.rstrip("\n").splitlines()
    n = int(lines[1].strip())
    return lines[0], lines[2:2 + n], lines[2 + n]


def merge_gro(protein_gro, ligand_gro, lig_resname="LIG"):
    '''Une las coordenadas del ligando a las de la proteína (caja de la proteína) y renumera los átomos.'''
    t1, a1, box1 = parse_gro(protein_gro)
    _t2, a2, _box2 = parse_gro(ligand_gro)
    last_res = int(a1[-1][:5]) if a1 else 0
    out, serial = [], 0
    for l in a1:
        serial += 1
        out.append(l[:15] + f"{serial % 100000:5d}" + l[20:])
    for l in a2:
        serial += 1
        out.append(f"{(last_res + 1) % 100000:5d}{lig_resname:<5s}{l[10:15]}{serial % 100000:5d}{l[20:]}")
    return f"{t1}\n{len(out):5d}\n" + "\n".join(out) + f"\n{box1}\n"


def normalize_ligand_itp(itp_text, lig_resname="LIG", posre_name="posre_LIG.itp"):
    '''Residuo del ligando -> LIG y restricciones de posición bajo -DPOSRES (idempotente).'''
    out, in_atoms = [], False
    for l in itp_text.splitlines():
        s = l.strip()
        if s.startswith("["):
            in_atoms = s.replace(" ", "").lower().startswith("[atoms]")
        elif in_atoms and s and not s.startswith(";"):
            p = l.split()
            if len(p) >= 8 and p[0].isdigit():
                l = re.sub(r"^(\s*\d+\s+\S+\s+\d+\s+)(\S+)", lambda m: m.group(1) + lig_resname, l, count=1)
        out.append(l)
    text = "\n".join(out)
    if f'#include "{posre_name}"' not in text:
        text += ("\n\n; Restricciones de posición del ligando (átomos pesados)\n#ifdef POSRES\n"
                 f'#include "{posre_name}"\n#endif\n')
    return text if text.endswith("\n") else text + "\n"


def patch_topology(top_text, itp_name="LIG_GMX.itp", mol_name="LIG"):
    '''Incluye el itp del ligando tras el campo de fuerzas y añade el ligando a [ molecules ] (idempotente).'''
    text = top_text
    if f'#include "{itp_name}"' not in text:
        out, inserted = [], False
        for l in text.splitlines():
            out.append(l)
            if not inserted and "forcefield.itp" in l and l.strip().startswith("#include"):
                out += ["", "; Topología del ligando (GAFF/GAFF2, generada con acpype)", f'#include "{itp_name}"']
                inserted = True
        if not inserted:
            raise ValueError("No se encontró el #include del campo de fuerzas en topol.top")
        text = "\n".join(out) + "\n"
    if not re.search(rf"^\s*{mol_name}\s+1\s*$", text, flags=re.M):
        text = text.rstrip("\n") + f"\n{mol_name:<18s}1\n"
    return text


def main(argv):
    if len(argv) >= 6 and argv[0] == "merge":
        _, pro_gro, lig_gro, top, lig_itp, out_gro = argv[:6]
        pg, lg = open(pro_gro).read(), open(lig_gro).read()
        top_t, itp_t = open(top).read(), open(lig_itp).read()
        open(out_gro, "w").write(merge_gro(pg, lg))
        open(lig_itp, "w").write(normalize_ligand_itp(itp_t))
        open(top, "w").write(patch_topology(top_t))
        print("Complejo escrito en", out_gro)
    else:
        print("Uso: python3 prepare_system.py merge proteina.gro ligando.gro topol.top LIG_GMX.itp complejo.gro")


if __name__ == "__main__":
    main(sys.argv[1:])
"""

_GMX_NS = {"__name__": "gmx_helpers"}
exec(_GMX_HELPERS, _GMX_NS)
amberize_pdb = _GMX_NS["amberize_pdb"]
merge_gro = _GMX_NS["merge_gro"]
normalize_ligand_itp = _GMX_NS["normalize_ligand_itp"]
patch_topology = _GMX_NS["patch_topology"]

GMX_AMBER_FFS = ["amber99sb-ildn", "amber99sb", "amber03", "amber99", "amber96", "amber94"]
GMX_WATERS = {"tip3p": "spc216.gro", "spce": "spc216.gro", "tip4p": "tip4p.gro", "tip4pew": "tip4p.gro"}


def mdp_files(temp=300.0, pres=1.0, dt=0.002, em_steps=5000, nvt_ps=100.0, npt_ps=100.0, prod_ns=10.0,
              out_ps=10.0, grp="Protein_LIG", rc=1.0):
    """Parámetros .mdp de todas las fases (dt en ps; tiempos en ps/ns)."""
    common = f"""; Interacciones no enlazantes (PME + Verlet)
cutoff-scheme            = Verlet
nstlist                  = 20
coulombtype              = PME
rcoulomb                 = {rc:g}
pme_order                = 4
fourierspacing           = 0.16
vdwtype                  = Cut-off
vdw-modifier             = Potential-shift
rvdw                     = {rc:g}
DispCorr                 = EnerPres
pbc                      = xyz
"""
    constr = """; Restricciones
constraints              = h-bonds
constraint_algorithm     = lincs
lincs_iter               = 1
lincs_order              = 4
"""
    tc = f"""; Termostato
tcoupl                   = V-rescale
tc-grps                  = {grp} Water_and_ions
tau_t                    = 0.1     0.1
ref_t                    = {temp:g}   {temp:g}
"""

    def steps(ps):
        return max(int(round(ps / dt)), 1)

    return {
        "ions.mdp": f"""; Solo para generar el .tpr que usa genion
integrator               = steep
emtol                    = 1000.0
emstep                   = 0.01
nsteps                   = 50000
{common}""",
        "em.mdp": f"""; Minimización de energía (máxima pendiente)
integrator               = steep
emtol                    = 1000.0
emstep                   = 0.01
nsteps                   = {int(em_steps)}
nstenergy                = 500
nstlog                   = 500
{common}""",
        "nvt.mdp": f"""; Equilibración NVT con restricciones de posición (proteína + ligando)
define                   = -DPOSRES
integrator               = md
dt                       = {dt:g}
nsteps                   = {steps(nvt_ps)}
nstenergy                = 500
nstlog                   = 500
continuation             = no
{constr}{common}{tc}pcoupl                   = no
gen_vel                  = yes
gen_temp                 = {temp:g}
gen_seed                 = -1
""",
        "npt.mdp": f"""; Equilibración NPT con restricciones de posición (C-rescale, requiere GROMACS >= 2021)
define                   = -DPOSRES
integrator               = md
dt                       = {dt:g}
nsteps                   = {steps(npt_ps)}
nstenergy                = 500
nstlog                   = 500
continuation             = yes
{constr}{common}{tc}pcoupl                   = C-rescale
pcoupltype               = isotropic
tau_p                    = 2.0
ref_p                    = {pres:g}
compressibility          = 4.5e-5
refcoord_scaling         = com
gen_vel                  = no
""",
        "md.mdp": f"""; Producción (sin restricciones)
integrator               = md
dt                       = {dt:g}
nsteps                   = {steps(prod_ns * 1000.0)}
nstenergy                = 5000
nstlog                   = 5000
nstxout-compressed       = {steps(out_ps)}
compressed-x-grps        = System
continuation             = yes
{constr}{common}{tc}pcoupl                   = Parrinello-Rahman
pcoupltype               = isotropic
tau_p                    = 2.0
ref_p                    = {pres:g}
compressibility          = 4.5e-5
gen_vel                  = no
""",
    }


def find_gmx():
    for name in ("gmx", "gmx_mpi", "gmx.exe"):
        path = shutil.which(name)
        if path:
            try:
                r = subprocess.run([path, "--version"], capture_output=True, text=True, timeout=30)
                m = re.search(r"GROMACS version:\s*(\S+)", r.stdout + r.stderr)
                if m:
                    return path, m.group(1)
            except Exception:
                continue
    return None, None


def gmx_tools_status():
    gmx, ver = find_gmx()
    return {"gmx": gmx, "gmx_version": ver, "acpype": shutil.which("acpype"),
            "antechamber": shutil.which("antechamber"), "sqm": shutil.which("sqm")}


def mol_to_mol2_sdf(mol):
    """Molécula RDKit con H -> (mol2, sdf, carga neta) mediante Open Babel."""
    molblock = Chem.MolToMolBlock(mol)
    ob = pybel.readstring("mol", molblock)
    net = int(sum(a.formalcharge for a in ob.atoms))
    return ob.write("mol2"), ob.write("sdf"), net, len(ob.atoms)


def ligand_gromacs_params(mol2_text, net_charge, gaff="gaff2", charges="bcc"):
    """Parametriza el ligando con acpype (GAFF/GAFF2). Devuelve {archivo: texto} listo para GROMACS."""
    if not shutil.which("acpype") or not shutil.which("antechamber"):
        raise RuntimeError("Hacen falta acpype y AmberTools (antechamber/sqm).")
    with tempfile.TemporaryDirectory() as td:
        with open(os.path.join(td, "LIG.mol2"), "w") as f:
            f.write(mol2_text)
        proc = subprocess.run(["acpype", "-i", "LIG.mol2", "-b", "LIG", "-n", str(int(net_charge)),
                               "-a", gaff, "-c", charges], cwd=td, capture_output=True, text=True, timeout=3600)
        out_dir = os.path.join(td, "LIG.acpype")
        need = ("LIG_GMX.itp", "LIG_GMX.gro", "posre_LIG.itp")
        if proc.returncode != 0 or not all(os.path.exists(os.path.join(out_dir, n)) for n in need):
            raise RuntimeError("acpype falló: " + (proc.stdout + proc.stderr).strip()[-400:])
        files = {n: open(os.path.join(out_dir, n)).read() for n in need}
    files["LIG_GMX.itp"] = normalize_ligand_itp(files["LIG_GMX.itp"])
    return files


_GMX_SCRIPT = r"""#!/usr/bin/env bash
# Flujo GROMACS generado por BioDock Studio. Es reanudable: omite los pasos cuyo resultado ya existe.
#   Uso:  bash run_gromacs.sh            (o: GMX=gmx_mpi MDRUN_OPTS="-nb gpu -pme gpu" bash run_gromacs.sh)
set -euo pipefail
GMX="${GMX:-gmx}"
MDRUN_OPTS="${MDRUN_OPTS:-}"
export GMX_MAXBACKUP=-1

@@LIGAND_PARAMS@@
echo "==> 1/8 Topología de la proteína (@@FF@@ + @@WATER@@)"
[ -f protein.gro ] || $GMX pdb2gmx -f protein_amber.pdb -o protein.gro -p topol.top -i posre.itp -ff @@FF@@ -water @@WATER@@ -ignh

echo "==> 2/8 Complejo proteína + ligando"
[ -f complex.gro ] || @@MERGE@@

echo "==> 3/8 Caja (@@BOX@@, @@DIST@@ nm) y solvatación"
[ -f newbox.gro ] || $GMX editconf -f complex.gro -o newbox.gro -c -d @@DIST@@ -bt @@BOX@@
[ -f solv.gro ] || $GMX solvate -cp newbox.gro -cs @@SOLVENT@@ -p topol.top -o solv.gro

echo "==> 4/8 Iones: neutralización + @@CONC@@ M NaCl"
if [ ! -f solv_ions.gro ]; then
  # -maxwarn 1: el único aviso esperado es «sistema con carga neta», que genion corrige a continuación
  $GMX grompp -f mdp/ions.mdp -c solv.gro -p topol.top -o ions.tpr -maxwarn 1
  echo SOL | $GMX genion -s ions.tpr -o solv_ions.gro -p topol.top -pname NA -nname CL -neutral -conc @@CONC@@
fi

echo "==> 5/8 Minimización de energía"
[ -f em.gro ] || { $GMX grompp -f mdp/em.mdp -c solv_ions.gro -p topol.top -o em.tpr; $GMX mdrun -v -deffnm em $MDRUN_OPTS; }

echo "==> 6/8 Grupos de índice"
[ -f index.ndx ] || @@INDEX@@

echo "==> 7/8 Equilibración NVT y NPT (con restricciones de posición)"
[ -f nvt.gro ] || { $GMX grompp -f mdp/nvt.mdp -c em.gro -r em.gro -p topol.top -n index.ndx -o nvt.tpr; $GMX mdrun -v -deffnm nvt $MDRUN_OPTS; }
[ -f npt.gro ] || { $GMX grompp -f mdp/npt.mdp -c nvt.gro -r nvt.gro -t nvt.cpt -p topol.top -n index.ndx -o npt.tpr; $GMX mdrun -v -deffnm npt $MDRUN_OPTS; }

echo "==> 8/8 Producción (@@PROD@@ ns)"
[ -f md.gro ] || { $GMX grompp -f mdp/md.mdp -c npt.gro -t npt.cpt -p topol.top -n index.ndx -o md.tpr; $GMX mdrun -v -deffnm md $MDRUN_OPTS; }
echo "Hecho. Trayectoria: md.xtc  ·  Energías: md.edr  ·  Estructura final: md.gro"
"""


def build_gmx_script(p, has_ligand, net_charge):
    ligand_params = ""
    merge = "cp protein.gro complex.gro"
    index = "echo q | $GMX make_ndx -f em.gro -o index.ndx"
    if has_ligand:
        ligand_params = (
            'if [ ! -f LIG_GMX.itp ]; then\n'
            f'  echo "==> 0/8 Parametrizando el ligando ({p["gaff"]}, cargas {p["charges"]}, carga neta {net_charge:+d})"\n'
            '  command -v acpype >/dev/null || { echo "Falta acpype/AmberTools: conda install -c conda-forge ambertools '
            '&& pip install acpype"; exit 1; }\n'
            f'  acpype -i ligand_pose.mol2 -b LIG -n {int(net_charge)} -a {p["gaff"]} -c {p["charges"]}\n'
            '  cp LIG.acpype/LIG_GMX.itp LIG.acpype/LIG_GMX.gro LIG.acpype/posre_LIG.itp .\n'
            'fi\n')
        merge = "python3 prepare_system.py merge protein.gro LIG_GMX.gro topol.top LIG_GMX.itp complex.gro"
        index = "printf '\"Protein\" | \"LIG\"\\nq\\n' | $GMX make_ndx -f em.gro -o index.ndx"
    return (_GMX_SCRIPT.replace("@@LIGAND_PARAMS@@", ligand_params).replace("@@MERGE@@", merge)
            .replace("@@INDEX@@", index).replace("@@FF@@", p["ff"]).replace("@@WATER@@", p["water"])
            .replace("@@BOX@@", p["box"]).replace("@@DIST@@", f'{p["dist"]:g}')
            .replace("@@SOLVENT@@", GMX_WATERS[p["water"]]).replace("@@CONC@@", f'{p["conc"]:g}')
            .replace("@@PROD@@", f'{p["prod_ns"]:g}'))


def build_gmx_readme(p, has_ligand, net_charge, changes, n_prot_atoms, built=None):
    lig = (f"- **Ligando**: {p['gaff'].upper()} con cargas {'AM1-BCC' if p['charges'] == 'bcc' else 'Gasteiger'} "
           f"(acpype), carga neta {net_charge:+d}.\n" if has_ligand else "- Sistema sin ligando (apo).\n")
    ren = ", ".join(f"{v}×{k}" for k, v in sorted(changes.items())) or "ninguno"
    built_txt = ""
    if built:
        built_txt = ("\n## Sistema ya construido\n"
                     f"Esta carpeta incluye el sistema solvatado y neutralizado (`solv_ions.gro`, `topol.top`, "
                     f"`index.ndx`, `em.tpr`): {built}\nEl script omite los pasos ya hechos.\n")
    return f"""# Entrada para GROMACS — generada por BioDock Studio

## Contenido
- `protein_amber.pdb`: proteína con H a pH (PDB2PQR/PROPKA); residuos renombrados según su estado de protonación
  ({ren}) para que `pdb2gmx` los respete.
{"- `ligand_pose.mol2` / `.sdf`: ligando con H en su pose (coordenadas del complejo)." + chr(10) if has_ligand else ""}{"- `LIG_GMX.itp`, `LIG_GMX.gro`, `posre_LIG.itp`: topología del ligando (acpype)." + chr(10) if has_ligand else ""}- `mdp/`: ions, em, nvt, npt, md.
- `prepare_system.py`: une ligando y proteína (solo biblioteca estándar de Python).
- `run_gromacs.sh`: todo el flujo, reanudable.
{built_txt}
## Modelo
- **Proteína**: {p['ff']} · **agua**: {p['water']} · caja {p['box']} a {p['dist']:g} nm · NaCl {p['conc']:g} M (+ neutralización).
{lig}- {p['temp']:g} K, {p['pres']:g} bar, paso {p['dt'] * 1000:g} fs, restricciones h-bonds (LINCS), PME, cutoffs 1,0 nm.
- Fases: minimización ({p['em_steps']} pasos) → NVT {p['nvt_ps']:g} ps → NPT {p['npt_ps']:g} ps (C-rescale) →
  producción {p['prod_ns']:g} ns (Parrinello-Rahman), trayectoria cada {p['out_ps']:g} ps.

## Cómo usarlo
```bash
bash run_gromacs.sh
# con GPU:   MDRUN_OPTS="-nb gpu -pme gpu" bash run_gromacs.sh
```
Requisitos: GROMACS ≥ 2021 (C-rescale){", acpype y AmberTools (solo si falta `LIG_GMX.itp`)" if has_ligand else ""}.

## Antes de producir resultados
- Revisa el sistema tras la minimización y las equilibraciones (temperatura, presión, densidad, RMSD).
- Átomos de la proteína en el PDB de entrada: {n_prot_atoms}. Comprueba que no haya huecos de cadena ni residuos
  no estándar (`pdb2gmx` avisará).
- Las cargas AM1-BCC/GAFF2 son una aproximación estándar para ligandos orgánicos; metales, cofactores o ligandos
  covalentes necesitan una parametrización específica.
- Una pose de docking no es un estado de equilibrio: la minimización y la equilibración con restricciones sirven
  para relajarla antes de la producción.
"""


def build_gmx_package(p, protein_pdb, ligand=None, lig_files=None):
    """Devuelve ({ruta: contenido}, resumen). ligand = (mol2, sdf, carga_neta) o None; lig_files = salida de acpype."""
    prot_lines = [l for l in protein_pdb.splitlines() if l.startswith(("ATOM", "HETATM", "TER"))]
    prot_txt, changes = amberize_pdb("\n".join(prot_lines) + "\nEND\n")
    n_atoms = sum(1 for l in prot_lines if l.startswith(("ATOM", "HETATM")))
    has_lig = ligand is not None
    net = ligand[2] if has_lig else 0
    grp = "Protein_LIG" if has_lig else "Protein"
    files = {"protein_amber.pdb": prot_txt, "prepare_system.py": _GMX_HELPERS.strip() + "\n"}
    for name, txt in mdp_files(p["temp"], p["pres"], p["dt"], p["em_steps"], p["nvt_ps"], p["npt_ps"],
                               p["prod_ns"], p["out_ps"], grp).items():
        files[f"mdp/{name}"] = txt
    if has_lig:
        files["ligand_pose.mol2"], files["ligand_pose.sdf"] = ligand[0], ligand[1]
        if lig_files:
            files.update(lig_files)
    files["run_gromacs.sh"] = build_gmx_script(p, has_lig, net)
    return files, {"changes": changes, "n_prot_atoms": n_atoms, "has_ligand": has_lig, "net": net}


def run_gmx_build(files, p, gmx, has_ligand, progress=None):
    """Construye y valida el sistema con GROMACS: pdb2gmx → complejo → caja → solvatación → iones → grompp (EM y
    fases de dinámica). Devuelve (ok, resumen, archivos construidos, registro)."""
    log = []
    env = dict(os.environ, GMX_MAXBACKUP="-1")
    with tempfile.TemporaryDirectory() as td:
        for rel, txt in files.items():
            path = os.path.join(td, rel)
            os.makedirs(os.path.dirname(path), exist_ok=True)
            with open(path, "w") as f:
                f.write(txt)

        def run(args, stdin=None, label=""):
            if progress:
                progress(label)
            r = subprocess.run([gmx] + args, cwd=td, input=stdin, capture_output=True, text=True,
                               timeout=3600, env=env)
            log.append(f"$ gmx {' '.join(args)}\n{(r.stdout + r.stderr).strip()[-1200:]}")
            if r.returncode != 0:
                tail = (r.stderr or r.stdout).strip()
                fatal = re.search(r"Fatal error:\n(.*?)\n\nFor more", tail, flags=re.S)
                raise RuntimeError(f"{label} falló: " + (fatal.group(1).strip() if fatal else tail[-400:]))
            return r

        try:
            run(["pdb2gmx", "-f", "protein_amber.pdb", "-o", "protein.gro", "-p", "topol.top", "-i", "posre.itp",
                 "-ff", p["ff"], "-water", p["water"], "-ignh"], label="pdb2gmx (topología de la proteína)")
            if has_ligand:
                pg, lg = (open(os.path.join(td, n)).read() for n in ("protein.gro", "LIG_GMX.gro"))
                itp = open(os.path.join(td, "LIG_GMX.itp")).read()
                top = open(os.path.join(td, "topol.top")).read()
                open(os.path.join(td, "complex.gro"), "w").write(merge_gro(pg, lg))
                open(os.path.join(td, "LIG_GMX.itp"), "w").write(normalize_ligand_itp(itp))
                open(os.path.join(td, "topol.top"), "w").write(patch_topology(top))
            else:
                shutil.copy(os.path.join(td, "protein.gro"), os.path.join(td, "complex.gro"))
            run(["editconf", "-f", "complex.gro", "-o", "newbox.gro", "-c", "-d", f'{p["dist"]:g}',
                 "-bt", p["box"]], label="editconf (caja)")
            run(["solvate", "-cp", "newbox.gro", "-cs", GMX_WATERS[p["water"]], "-p", "topol.top",
                 "-o", "solv.gro"], label="solvate (agua)")
            run(["grompp", "-f", "mdp/ions.mdp", "-c", "solv.gro", "-p", "topol.top", "-o", "ions.tpr",
                 "-maxwarn", "1"], label="grompp (iones)")
            run(["genion", "-s", "ions.tpr", "-o", "solv_ions.gro", "-p", "topol.top", "-pname", "NA",
                 "-nname", "CL", "-neutral", "-conc", f'{p["conc"]:g}'], stdin="SOL\n", label="genion (iones)")
            ndx_in = '"Protein" | "LIG"\nq\n' if has_ligand else "q\n"
            run(["make_ndx", "-f", "solv_ions.gro", "-o", "index.ndx"], stdin=ndx_in, label="make_ndx (índices)")
            run(["grompp", "-f", "mdp/em.mdp", "-c", "solv_ions.gro", "-p", "topol.top", "-o", "em.tpr"],
                label="grompp (minimización)")
            for ph in ("nvt", "npt", "md"):    # validación: construye los .tpr con las coordenadas solvatadas
                run(["grompp", "-f", f"mdp/{ph}.mdp", "-c", "solv_ions.gro", "-r", "solv_ions.gro",
                     "-p", "topol.top", "-n", "index.ndx", "-o", f"_check_{ph}.tpr"], label=f"grompp ({ph}, validación)")
        except Exception as e:
            return False, {"error": str(e)}, {}, "\n\n".join(log)

        top_final = open(os.path.join(td, "topol.top")).read()
        mols = re.findall(r"^\s*([A-Za-z0-9_+\-]+)\s+(\d+)\s*$", top_final.split("[ molecules ]")[-1], flags=re.M)
        gro = open(os.path.join(td, "solv_ions.gro")).read().splitlines()
        summary = {"atoms": int(gro[1]), "box": gro[-1].split()[:3], "molecules": {n: int(c) for n, c in mols}}
        built = {}
        for rel in ("topol.top", "solv_ions.gro", "index.ndx", "em.tpr", "posre.itp", "protein.gro",
                    "complex.gro", "newbox.gro", "solv.gro", "ions.tpr", "LIG_GMX.itp", "posre_LIG.itp", "LIG_GMX.gro"):
            path = os.path.join(td, rel)
            if os.path.exists(path):
                built[rel] = open(path, "rb").read()
        for f_ in os.listdir(td):      # topologías que genera pdb2gmx para cadenas
            if f_.startswith("topol_") and f_.endswith(".itp"):
                built[f_] = open(os.path.join(td, f_), "rb").read()
    return True, summary, built, "\n\n".join(log)


def zip_bytes(files):
    """{ruta: str|bytes} -> ZIP en memoria."""
    buf = io.BytesIO()
    with zipfile.ZipFile(buf, "w", zipfile.ZIP_DEFLATED) as z:
        for fname, content in files.items():
            if fname.endswith(".sh"):
                zi = zipfile.ZipInfo(fname)
                zi.external_attr = 0o755 << 16
                z.writestr(zi, content)
            else:
                z.writestr(fname, content)
    return buf.getvalue()


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
            with st.spinner("Descargando estructura..."):
                data = fetch_pdb_from_web(pdb_code)
                if data:
                    st.session_state["raw_pdb"] = data
                    st.session_state["pdb_name"] = pdb_code.upper()
                    # Extraer el nombre de la molécula directamente del texto PDB
                    protein_name = extract_molecule_name(data, default_code=pdb_code)
                    st.success(f"¡Estructura de '{protein_name}' cargada con éxito!")
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
    tab1, tab2, tab3, tab4, tab5, tab6 = st.tabs(["👁️ Visor 3D & Preparación", "📊 Propiedades Fisicoquímicas & Secuencia",
                                                  "🧪 Preparación Docking (Vina)", "🚀 Docking (Vina)",
                                                  "📈 Resultados de Docking", "🧬 GROMACS"])

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
                
                # --- DETECCIÓN AUTOMÁTICA DE HETEROGÉNEOS ---
                all_hets = get_structure_heterogens(raw_pdb)
                
                # Definir qué consideramos cofactores o iones comunes para marcar por defecto
                common_cofts = {"NAP", "NAD", "NDP", "FAD", "HEM", "CA", "ZN", "MG", "MN", "FE", "K", "NA", "CL", "PO4", "SO4"}
                default_retained = [h for h in all_hets if h in common_cofts]
                
                st.info("💡 **Consejo:** Marca abajo únicamente los cofactores (ej. NAP) e iones (ej. CA) que quieras mantener fijos en la proteína. **Deja desmarcado el ligando que vayas a extraer para hacer docking**.")
                
                selected_heterogens = st.multiselect(
                    "Heterogéneos encontrados en el PDB (selecciona cuáles conservar en el receptor):",
                    options=all_hets,
                    default=default_retained,
                    key=f"auto_hets_{pdb_name}",
                )
                
                # Campo libre por si hay cofactores que no figuran como HETATM estándar
                manual_hets_input = st.text_input(
                    "Otros cofactores / iones adicionales a conservar (códigos PDB, separados por comas):",
                    value="",
                    key=f"manual_hets_{pdb_name}",
                    help="Úsalo para introducir manualmente códigos de residuos adicionales."
                )
                
                # Unir selecciones automáticas y manuales
                manual_list = [c.strip().upper() for c in manual_hets_input.split(",") if c.strip()]
                retained_tuple = tuple(sorted(list(set(selected_heterogens + manual_list))))

            if not selected_chains:
                st.warning("Selecciona al menos una cadena.")
                st.stop()

            # Pasar la tupla de elementos retenidos a la función de procesado
            result = process_pdb_structure(raw_pdb, tuple(selected_chains), remove_water, remove_ions, retained_tuple)
            
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
                cofactors=retained_tuple
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
                    pm_custom1 = pymol_ui("t1")
                    pse_key = f"pse_{pdb_name}_{sel_lig}"
                    if st.button("🧩 Generar sesión de PyMOL (.pse)", key=f"btn_{pse_key}"):
                        with st.spinner("Generando la sesión con PLIP + PyMOL..."):
                            st.session_state[pse_key] = export_pymol_session(
                                cpx_dl, f"{i_res}:{i_chain}:{int(i_num)}", i_items,
                                f"resn {i_res} and chain {i_chain} and resi {int(i_num)}", pymol_exe=pm_custom1)
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
                    
                    # PINTAR COFACTORES EN STICKS DENTRO DE LA CAJA
                    if retained_tuple:
                        bv.addStyle({'model': 0, 'resn': list(retained_tuple)},
                                    {'stick': {'colorscheme': 'yellowCarbon', 'radius': 0.20}})
                    
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
            st.error("No hay motor de Vina disponible aquí: el paquete `vina` de pip solo existe para Linux y "
                     "Python ≤ 3.12, y no se encontró el ejecutable oficial.")
            if vina_asset_name():
                if st.button(f"⬇️ Descargar Vina {VINA_VERSION} ({vina_asset_name()}, ≈ 4 MB) desde GitHub",
                             type="primary", key="vina_dl"):
                    try:
                        with st.spinner("Descargando el ejecutable oficial de Vina..."):
                            download_vina_exe()
                        st.rerun()
                    except Exception as e:
                        st.error(f"No se pudo descargar Vina: {e}")
                st.caption("Descarga el binario oficial de la release de AutoDock-Vina (GitHub) a una carpeta "
                           "temporal del servidor/equipo. También puedes dejarlo junto a este script o indicar su ruta arriba.")
            else:
                st.caption("Sistema no reconocido: descarga Vina de https://github.com/ccsb-scripps/AutoDock-Vina/releases "
                           "e indica su ruta arriba.")
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
                search_mode = st.radio(
                    "Modo de búsqueda:",
                    ["Estándar de Vina (más rápido)", "Runs independientes + clústeres (estilo AutoDock4)"],
                    horizontal=True, key="dk_mode",
                    help="Estándar: una sola ejecución de Vina con N búsquedas internas; no informa de clústeres. "
                         "Runs independientes: N ejecuciones con semillas distintas que se agrupan por RMSD, "
                         "lo que da el nº de conformaciones por clúster y el nº total de runs (como AutoDock4).")
                runs_mode = search_mode.startswith("Runs")
                exh, npose, n_runs, tol_c, erange = 8, 9, 20, 2.0, 3.0
                d1, d2, d3 = st.columns(3)
                with d1:
                    scoring = st.selectbox("Función de puntuación", ["vina", "vinardo"], key="dk_sf")
                with d2:
                    if runs_mode:
                        n_runs = st.slider("Nº de runs independientes", 5, 100, 20, key="dk_runs",
                                           help="Cada run es una búsqueda completa (exhaustividad 1) con otra semilla.")
                    else:
                        exh = st.slider("Exhaustividad", 1, 64, 8, key="dk_exh",
                                        help="Más alta = búsqueda más profunda y más lenta.")
                with d3:
                    if runs_mode:
                        tol_c = st.slider("Tolerancia del clustering (RMSD, Å)", 0.5, 5.0, 2.0, 0.5, key="dk_tol",
                                          help="Dos runs caen en el mismo clúster si su RMSD átomo a átomo ≤ esta tolerancia "
                                               "(AutoDock4 usa 2,0 Å).")
                    else:
                        npose = st.slider("Nº de poses", 1, 20, 9, key="dk_np")
                d4, d5, d6 = st.columns(3)
                with d4:
                    if runs_mode:
                        st.caption("Modo runs: cada pose mostrada es el representante (menor energía) de un clúster.")
                    else:
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
                if runs_mode:
                    st.caption("Las runs se reparten entre los núcleos disponibles. Referencia (1 núcleo): ≈ 4 s por run "
                               "para un ligando de 6 torsiones en una proteína de ~5000 átomos.")
                else:
                    st.caption("Referencia: con exhaustividad 8, el redocking de un ligando de 6 torsiones en una "
                               "proteína de ~5000 átomos tardó ≈ 30 s con 1 núcleo.")

                if st.button("🚀 Ejecutar docking", type="primary", use_container_width=True, key="dk_run",
                             disabled=not chosen):
                    out_items = {}
                    bar = st.progress(0.0, text="Preparando...")
                    for k, label in enumerate(chosen):
                        d_ = avail[label]
                        bar.progress(k / len(chosen), text=f"Acoplando {d_['name']} ({k + 1}/{len(chosen)})...")
                        t0 = time.time()
                        try:
                            if runs_mode:
                                def _prog(done, total, _k=k, _n=len(chosen), _nm=d_["name"]):
                                    bar.progress((_k + done / total) / _n, text=f"Acoplando {_nm}: run {done}/{total}")
                                pq, en, clus, run_e = run_vina_independent(
                                    rec_stored[1]["pdbqt"], d_["pdbqt"], box_p[0], box_p[1], n_runs, seed, cpu,
                                    scoring, exe=vina_exe, tol=tol_c, progress=_prog)
                                out_items[d_["name"]] = {"ok": True, "pdbqt": pq, "energies": en,
                                                         "seconds": time.time() - t0, "ref_label": d_["ref_label"],
                                                         "mode": "runs", "runs": n_runs, "clusters": clus,
                                                         "run_energies": run_e}
                            else:
                                pq, en = run_vina(rec_stored[1]["pdbqt"], d_["pdbqt"], box_p[0], box_p[1],
                                                  exh, npose, erange, seed, cpu, scoring, exe=vina_exe)
                                out_items[d_["name"]] = {"ok": True, "pdbqt": pq, "energies": en,
                                                         "seconds": time.time() - t0, "ref_label": d_["ref_label"],
                                                         "mode": "standard"}
                        except Exception as e:
                            out_items[d_["name"]] = {"ok": False, "error": str(e)[:300]}
                    bar.progress(1.0, text="Docking completado")
                    params_txt = (f"{scoring} · {n_runs} runs independientes · clústeres a {tol_c} Å · semilla {seed}"
                                  if runs_mode else
                                  f"{scoring} · exhaustividad {exh} · {npose} poses · semilla {seed}")
                    st.session_state["dock_results"] = {"items": out_items, "pdb_name": pdb_name,
                                                        "params": params_txt}

            dock_store = st.session_state.get("dock_results")
            if dock_store:
                st.markdown("---")
                st.markdown("#### Resultados")
                st.caption(f"Parámetros: {dock_store['params']}")
                rows, rmsd_by_lig, poses_by_lig = [], {}, {}
                for nm, it_ in dock_store["items"].items():
                    if not it_.get("ok"):
                        st.error(f"❌ {nm}: {it_.get('error', 'error desconocido')}")
                        rows.append({"Ligando": nm, "Modo": None, "Mejor afinidad (kcal/mol)": None,
                                     "Ki estimada": None, "Poses": 0, "Tiempo (s)": None,
                                     "RMSD mejor pose (Å)": None, "Estado": "❌ error"})
                        continue
                    poses_d, _, _ = load_poses("out.pdbqt", it_["pdbqt"], "")
                    poses_by_lig[nm] = poses_d
                    best_rmsd = best_any = None
                    if it_["ref_label"] and HAS_RDKIT:
                        try:
                            rr_, rrest_ = it_["ref_label"].split(" ", 1)
                            rn2, rc2 = rrest_.split(":")
                            ref_m = ligand_mol_for_map(extract_single_ligand(raw_pdb, rr_, int(rn2), rc2),
                                                       ensure_smiles(pdb_name, it_["ref_label"], rr_))
                            vals = [ligand_rmsd(pose_mol(p_["molblock"]), ref_m)[0] for p_ in poses_d]
                            rmsd_by_lig[nm] = vals
                            best_rmsd = vals[0] if vals and vals[0] is not None else None
                            best_any = min([v for v in vals if v is not None], default=None)
                        except Exception:
                            pass
                    best_sc = poses_d[0]["score"] if poses_d else (it_["energies"][0][0] if it_["energies"] else None)
                    rows.append({"Ligando": nm, "Modo": "runs + clústeres" if it_.get("mode") == "runs" else "estándar",
                                 "Mejor afinidad (kcal/mol)": best_sc, "Ki estimada": format_ki(best_sc),
                                 "Poses": len(poses_d), "Tiempo (s)": round(it_["seconds"], 1),
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
                    st.dataframe(pose_rows(poses_by_lig[pick], rmsd_by_lig.get(pick)),
                                 hide_index=True, use_container_width=True)
                    st.caption(
                        "**Ki estimada** = exp(ΔG/RT) con ΔG = afinidad de Vina y T = 298,15 K. Es una estimación "
                        "orientativa: la puntuación de Vina no es un ΔG calibrado, así que úsala para comparar ligandos, "
                        "no como valor absoluto. · **RMSD l.b. / u.b. (vs pose 1)**: distancia en Å (átomos pesados) entre "
                        "cada pose y la MEJOR pose de esta misma ejecución (la pose 1 vale 0), no frente al cristal. "
                        "u.b. (cota superior) empareja cada átomo con su homólogo exacto; l.b. (cota inferior) lo empareja "
                        "con el átomo más cercano del mismo elemento, por eso l.b. ≤ u.b.; la diferencia revela simetría "
                        "o grupos intercambiables. Si ambos son pequeños, son la misma pose; si l.b. es bajo y u.b. alto, "
                        "es casi la misma pose con átomos equivalentes permutados. · **RMSD vs cristal**: frente al "
                        "ligando experimental (redocking).")
                    if it_.get("mode") == "runs":
                        with st.expander("📊 Histograma de clústeres (estilo AutoDock4)", expanded=True):
                            st.dataframe([{"Rango": c_["rank"], "Menor energía (kcal/mol)": round(c_["best"], 3),
                                           "Energía media del clúster": round(c_["mean"], 3),
                                           "Nº conformaciones en el clúster": c_["n"],
                                           "Runs totales": it_["runs"]} for c_ in it_["clusters"]],
                                         hide_index=True, use_container_width=True)
                            st.caption("Cada run es una búsqueda independiente; se ordenan por energía y cada una se "
                                       "asigna al primer clúster cuyo representante esté a ≤ la tolerancia (RMSD "
                                       "átomo a átomo). Un clúster muy poblado indica una pose reproducible; uno con "
                                       "pocas conformaciones puede ser un mínimo poco fiable.")
                    elif not any(p_.get("cluster_n") for p_ in poses_by_lig[pick]):
                        st.caption("Modo estándar: Vina no informa de clústeres ni de runs. Usa «Runs independientes + "
                                   "clústeres» para obtener esas columnas.")
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
                    rows5 = pose_rows(poses, rmsds if ref_choice != "Ninguno" else None)
                    st.dataframe(rows5, hide_index=True, use_container_width=True)
                    st.caption("**RMSD l.b. / u.b. (vs pose 1)**: distancia entre cada pose y la mejor pose del archivo "
                               "(no frente al cristal): u.b. empareja cada átomo con su homólogo; l.b. con el átomo más "
                               "cercano del mismo elemento (l.b. ≤ u.b.). **Ki estimada** = exp(ΔG/RT), T = 298,15 K "
                               "(orientativa). **Nº conf. en clúster / Runs totales**: solo en resultados generados "
                               "con el modo «Runs independientes».")
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

                    pv = build_pose_view(protein_pdb, pmol, p_filtered, res_labels, ref_lig_pdb, cofactors=retained_tuple)
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

                    pm_custom5 = pymol_ui("t5")
                    pse5 = f"res_pse_{pose_idx}_{hashlib.md5(cpx.encode()).hexdigest()[:8]}"
                    if st.button("🧩 Generar sesión de PyMOL (.pse) de esta pose", key=f"btn_{pse5}"):
                        with st.spinner("Generando la sesión con PLIP + PyMOL..."):
                            st.session_state[pse5] = export_pymol_session(
                                cpx, "LIG:Z:1", p_items, "resn LIG and chain Z and resi 1", pymol_exe=pm_custom5)
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
                                         "Ki estimada": format_ki(poses[k]["score"]),
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

    with tab6:
        st.markdown("### 🧬 Entrada para GROMACS (dinámica molecular)")
        st.caption("Genera un paquete listo para simular: proteína con los H a pH (PDB2PQR/PROPKA), ligando "
                   "parametrizado (GAFF2 + cargas AM1-BCC con acpype), archivos .mdp y un script reanudable "
                   "con todo el flujo (pdb2gmx → caja → agua → iones → minimización → NVT → NPT → producción).")
        tools = gmx_tools_status()
        lig_tools_ok = bool(tools["acpype"] and tools["antechamber"] and tools["sqm"])
        g1, g2, g3 = st.columns(3)
        g1.write(("✅" if tools["gmx"] else "⚪") + " **GROMACS** " +
                 (tools["gmx_version"] if tools["gmx"] else "no encontrado (opcional: permite construir y validar aquí)"))
        g2.write(("✅" if lig_tools_ok else "⚪") + " **acpype + AmberTools** " +
                 ("disponibles" if lig_tools_ok else "no encontrados (el script del paquete parametrizará el ligando donde se ejecute)"))
        g3.write(("✅" if HAS_PDB2PQR else "⚪") + " **PDB2PQR/PROPKA** " +
                 ("disponible" if HAS_PDB2PQR else "no instalado"))

        if not (HAS_OB and HAS_RDKIT):
            st.error("Esta pestaña necesita Open Babel y RDKit: `pip install openbabel-wheel rdkit`.")
        else:
            dock6 = st.session_state.get("dock_results")
            good6 = {n: v for n, v in (dock6["items"] if dock6 else {}).items() if v.get("ok")}
            cry6 = st.session_state.get("lig_prep")
            kinds = []
            if good6:
                kinds.append("Complejo con una pose de docking")
            if cry6:
                kinds.append("Complejo con el ligando cristalográfico preparado")
            kinds.append("Solo proteína (apo)")
            # clave ligada a las opciones: si aparecen nuevas (tras preparar/acoplar), el selector se reinicia
            kind = st.radio("Sistema a simular:", kinds,
                            key="gm_kind_" + hashlib.md5("|".join(kinds).encode()).hexdigest()[:6])

            mol6, lig_desc = None, ""
            if kind.startswith("Complejo con una pose"):
                nm6 = st.selectbox("Ligando acoplado:", list(good6), key="gm_dock_lig")
                poses6, _, _ = load_poses("out.pdbqt", good6[nm6]["pdbqt"], "")
                k6 = st.selectbox("Pose:", list(range(len(poses6))), key="gm_pose",
                                  format_func=lambda k: f"Pose {k + 1}" + (
                                      f" · {poses6[k]['score']:.2f} kcal/mol" if poses6[k]["score"] is not None else ""))
                mol6, lig_desc = pose_mol(poses6[k6]["molblock"]), f"{nm6} · pose {k6 + 1}"
            elif kind.startswith("Complejo con el ligando cristalográfico"):
                mol6 = Chem.MolFromMolBlock(cry6[1]["sdf"], removeHs=False)
                lig_desc = f"{cry6[1]['name']} (cristalográfico)"
            has_lig6 = mol6 is not None
            if has_lig6 and not any(a.GetAtomicNum() == 1 for a in mol6.GetAtoms()):
                mol6 = Chem.AddHs(mol6, addCoords=True)
                st.warning("El ligando no traía hidrógenos: se añadieron sin considerar el pH. Revisa su estado de "
                           "protonación (el PDBQT de Vina no incluye los H apolares).")

            prepared_rec6 = st.session_state.get("rec_prep")
            if prepared_rec6:
                prot_choice = st.radio("Proteína:", ["Receptor protonado preparado (PDB2PQR/PROPKA)",
                                                     "Proteína limpia (pdb2gmx añade los H)"],
                                       horizontal=True, key="gm_prot")
                protein_text6 = (prepared_rec6[1]["pdb"] if prot_choice.startswith("Receptor protonado")
                                 else result["cleaned"])
            else:
                protein_text6 = result["cleaned"]
                st.info("Prepara el receptor en la pestaña 🧪 (mejor con PDB2PQR + PROPKA) para conservar los estados "
                        "de protonación de His, Asp y Glu. Sin él, `pdb2gmx` añadirá los H con sus valores por defecto.")

            st.markdown("#### Parámetros")
            pa, pb, pc = st.columns(3)
            with pa:
                ff6 = st.selectbox("Campo de fuerzas (proteína)", GMX_AMBER_FFS, key="gm_ff",
                                   help="Familia AMBER: es la compatible con GAFF/GAFF2 para el ligando.")
            with pb:
                water6 = st.selectbox("Modelo de agua", list(GMX_WATERS), key="gm_water")
            with pc:
                box6 = st.selectbox("Caja", ["dodecahedron", "cubic"], key="gm_box",
                                    format_func=lambda b: "Dodecaedro rómbico (menos agua)" if b == "dodecahedron" else "Cúbica")
            pd_, pe, pf = st.columns(3)
            with pd_:
                dist6 = st.number_input("Distancia soluto–borde (nm)", 0.8, 3.0, 1.0, 0.1, key="gm_dist")
            with pe:
                conc6 = st.number_input("NaCl (M)", 0.0, 1.0, 0.15, 0.05, key="gm_conc")
            with pf:
                temp6 = st.number_input("Temperatura (K)", 250.0, 400.0, 300.0, 1.0, key="gm_temp")
            pg_, ph_, pi_ = st.columns(3)
            with pg_:
                pres6 = st.number_input("Presión (bar)", 0.5, 10.0, 1.0, 0.5, key="gm_pres")
            with ph_:
                dt6 = st.selectbox("Paso de integración (fs)", [2, 1, 4], key="gm_dt",
                                   help="2 fs con restricciones h-bonds es lo estándar. 4 fs exige masas de H aumentadas.")
            with pi_:
                em6 = st.number_input("Pasos de minimización", 10, 100000, 5000, 500, key="gm_em")
            pj, pk, pl = st.columns(3)
            with pj:
                nvt6 = st.number_input("NVT (ps)", 1.0, 5000.0, 100.0, 10.0, key="gm_nvt")
            with pk:
                npt6 = st.number_input("NPT (ps)", 1.0, 5000.0, 100.0, 10.0, key="gm_npt")
            with pl:
                prod6 = st.number_input("Producción (ns)", 0.001, 10000.0, 10.0, 1.0, key="gm_prod")
            pm_, pn_ = st.columns(2)
            with pm_:
                out6 = st.number_input("Guardar trayectoria cada (ps)", 1.0, 1000.0, 10.0, 1.0, key="gm_out")
            gaff6, charges6 = "gaff2", "bcc"
            if has_lig6:
                with pn_:
                    gaff6 = st.selectbox("Campo de fuerzas del ligando", ["gaff2", "gaff"], key="gm_gaff")
                charges6 = st.radio("Cargas del ligando", ["bcc", "gas"], horizontal=True, key="gm_chg",
                                    format_func=lambda c: "AM1-BCC (recomendado, más lento)" if c == "bcc"
                                    else "Gasteiger (rápido, menos preciso)")

            build_now6 = False
            if tools["gmx"]:
                build_now6 = st.checkbox("Construir y validar el sistema ahora con GROMACS "
                                         "(pdb2gmx → solvate → genion → grompp)", value=True, key="gm_build")

            p6 = dict(ff=ff6, water=water6, box=box6, dist=float(dist6), conc=float(conc6), temp=float(temp6),
                      pres=float(pres6), dt=dt6 / 1000.0, em_steps=int(em6), nvt_ps=float(nvt6),
                      npt_ps=float(npt6), prod_ns=float(prod6), out_ps=float(out6), gaff=gaff6, charges=charges6)
            sig6 = (kind, lig_desc, tuple(sorted(p6.items())), hashlib.md5(protein_text6.encode()).hexdigest()[:8])

            if st.button("📦 Generar paquete GROMACS", type="primary", use_container_width=True, key="gm_go"):
                try:
                    ligand6, lig_files6, net6 = None, None, 0
                    if has_lig6:
                        mol2_6, sdf6, net6, _nat6 = mol_to_mol2_sdf(mol6)
                        ligand6 = (mol2_6, sdf6, net6)
                        if lig_tools_ok:
                            with st.spinner("Parametrizando el ligando con acpype "
                                            "(AM1-BCC puede tardar de 1 a 3 minutos)..."):
                                lig_files6 = ligand_gromacs_params(mol2_6, net6, gaff6, charges6)
                    files6, info6 = build_gmx_package(p6, protein_text6, ligand6, lig_files6)
                    built_msg, build_log, summ6 = None, "", None
                    build_err = None
                    if build_now6 and tools["gmx"]:
                        if has_lig6 and not lig_files6:
                            build_err = ("No se puede construir aquí sin la topología del ligando "
                                         "(faltan acpype/AmberTools); el script del paquete lo hará.")
                        else:
                            bar6 = st.progress(0.0, text="Construyendo el sistema...")
                            counter = {"n": 0}

                            def _prog6(label):
                                counter["n"] += 1
                                bar6.progress(min(counter["n"] / 11, 0.99), text=label)
                            ok6, summ6, built6, build_log = run_gmx_build(files6, p6, tools["gmx"], has_lig6,
                                                                          progress=_prog6)
                            bar6.empty()
                            if ok6:
                                files6.update(built6)
                                mols6 = summ6["molecules"]
                                built_msg = (f"{summ6['atoms']:,} átomos; " +
                                             ", ".join(f"{n} × {c}" for n, c in mols6.items()) + ".")
                            else:
                                build_err = summ6.get("error", "error desconocido")
                    files6["README.md"] = build_gmx_readme(p6, has_lig6, net6, info6["changes"],
                                                           info6["n_prot_atoms"], built=built_msg)
                    st.session_state["gmx_pkg"] = {
                        "sig": sig6, "zip": zip_bytes(files6), "files": sorted(files6),
                        "script": files6["run_gromacs.sh"], "mdp": {k: v for k, v in files6.items() if k.startswith("mdp/")},
                        "summary": summ6 if built_msg else None, "built_msg": built_msg, "build_err": build_err,
                        "log": build_log, "info": info6, "net": net6, "lig_params": bool(lig_files6)}
                except Exception as e:
                    st.session_state.pop("gmx_pkg", None)
                    st.error(f"No se pudo generar el paquete: {e}")

            pkg = st.session_state.get("gmx_pkg")
            if pkg and pkg["sig"] == sig6:
                name6 = f"{pdb_name}_{safe_name(lig_desc) if has_lig6 else 'apo'}_gromacs"
                if pkg["build_err"]:
                    st.warning("El paquete se generó, pero la construcción con GROMACS falló: " + pkg["build_err"])
                elif pkg["built_msg"]:
                    st.success("✅ Sistema construido y validado con `grompp` (minimización, NVT, NPT y producción): "
                               + pkg["built_msg"])
                    s_ = pkg["summary"]
                    m1, m2, m3, m4 = st.columns(4)
                    m1.metric("Átomos", f"{s_['atoms']:,}")
                    m2.metric("Aguas", f"{s_['molecules'].get('SOL', 0):,}")
                    m3.metric("Na⁺ / Cl⁻", f"{s_['molecules'].get('NA', 0)} / {s_['molecules'].get('CL', 0)}")
                    m4.metric("Caja (nm)", " × ".join(f"{float(x):.1f}" for x in s_["box"]))
                else:
                    st.success("✅ Paquete generado.")
                if has_lig6:
                    st.caption(f"Ligando: carga neta {pkg['net']:+d} · " +
                               ("topología GAFF generada aquí con acpype." if pkg["lig_params"] else
                                "sin topología: el script ejecutará acpype donde corras GROMACS."))
                if pkg["info"]["changes"]:
                    st.caption("Estados de protonación respetados: " +
                               ", ".join(f"{v} × {k}" for k, v in sorted(pkg["info"]["changes"].items())))
                st.download_button("📥 Descargar paquete GROMACS (.zip)", pkg["zip"], file_name=f"{name6}.zip",
                                   mime="application/zip", type="primary", use_container_width=True, key="gm_dl")
                with st.expander("📄 Archivos del paquete"):
                    st.code("\n".join(pkg["files"]), language="text")
                with st.expander("▶️ run_gromacs.sh (flujo completo)"):
                    st.code(pkg["script"], language="bash")
                with st.expander("⚙️ Parámetros .mdp"):
                    mdp_pick = st.selectbox("Fase:", sorted(pkg["mdp"]), key="gm_mdp_pick")
                    st.code(pkg["mdp"][mdp_pick], language="ini")
                if pkg["log"]:
                    with st.expander("🧾 Registro de GROMACS"):
                        st.code(pkg["log"][-6000:], language="text")
            elif pkg:
                st.caption("Los ajustes han cambiado: pulsa de nuevo «Generar paquete GROMACS».")