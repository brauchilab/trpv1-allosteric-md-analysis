08k_v3.3_path_clusters_to_pdb.py
Author: DenyCB
=================================
Construye archivos PDB con cadenas de glicina (GLY) que representan
paths alostéricos como ribbons/flechas en VMD o PyMOL.

DISEÑO
------
  - Cada cluster CSV → un archivo PDB (con un archivo TCL y PML companion)
  - Cada fila del CSV (un cluster) → una cadena GLY independiente
  - Los CA se colocan en las posiciones CB (o CA para GLY) de cada
    residuo del path, interpolando GLYs intermedios cuando la distancia
    entre residuos consecutivos supera INTERP_MAX_STEP
  - El backbone N/CA/C/O se construye con geometría de hoja-β (ángulo
    N-CA-C = 111°, usando el frame local del path y parallel transport
    del vector perpendicular para garantizar ribbon continuo en VMD)
  - Registros SHEET cubren la cadena completa → VMD los lee al cargar
    y asigna estructura E sin recalcular DSSP
  - El TCL companion refuerza la asignación E y crea representaciones
    NewCartoon individuales por cadena (sin llamar mol ssrecalc,
    que resetearía E a coil por ausencia de H-bonds inter-hebra)

ESTRUCTURA DE DIRECTORIOS
--------------------------
  INPUT:  BASE_INPUT_DIR / {sistema} / auto_cluster / clusters_*.csv
  PDB:    PDB_DIR / {7LP9|7LPB|7LPC}.pdb
  OUTPUT: BASE_OUTPUT_DIR / {sistema} / {fwd|rev} / {sistema}_{target}_{dir}.pdb

DEPENDENCIAS
------------
  numpy  (única dependencia externa)
  Python >= 3.7

"""

import os
import re
import csv
import sys
import logging
import numpy as np
from pathlib import Path
from datetime import datetime
from concurrent.futures import ProcessPoolExecutor, as_completed

# ==============================================================
# VERSIÓN
# ==============================================================

VERSION = "08k_v3.3"

# ==============================================================
# ⚙️  CONFIGURACIÓN
#     Modifica esta sección antes de ejecutar.
# ==============================================================

# --- Directorios de entrada ---
BASE_INPUT_DIR = Path(
    r"C:\DinamicasMoleculares\TRPV1_pipeline\output"
    r"\08_Max_Flow\08h_v1.4_paths_MD\method1_TshortestPaths"
)
PDB_DIR = Path(
    r"C:\DinamicasMoleculares\TRPV1_pipeline\input\pdb"
)

# --- Directorio de salida ---
BASE_OUTPUT_DIR = Path(
    r"C:\DinamicasMoleculares\TRPV1_pipeline\output"
    r"\08_Max_Flow\08k_v3.3_paths_pdb"
)

# --- Columna de path a usar ---
# "rep_chainresid"       → path representante del cluster
# "consensus_chainresid" → path consenso del cluster
PATH_COLUMN = "rep_chainresid"

# --- Interpolación ---
# Distancia CB-CB máxima (Å) antes de insertar GLYs intermedios.
# Distancia CA-CA típica en conformación extendida: ~3.8 Å.
# Valores > 5 Å indican un salto topológico entre cadenas o dominios.
# Subir: menos relleno, paths más cortos.
# Bajar: más relleno, paths más suaves (puede distorsionar la trayectoria).
INTERP_MAX_STEP = 5.0

# --- Residuo dummy terminal ---
# VMD no renderiza el último residuo de una cadena: lo usa solo como
# nodo de dirección para el arrowhead. Añadir un GLY fantasma tras el
# último residuo real garantiza que el último nodo del path aparezca
# completamente dibujado con su flecha.
# Desactivar solo si se observan artefactos en la punta del ribbon.
ADD_TERMINAL_DUMMY = False

# --- Rotación de cadenas para homotetrámero ---
# Condición: si el path transita de cadena A a cadena B con el resid
# del nodo A en el rango (CHAIN_ROT_RESID_MIN, CHAIN_ROT_RESID_MAX),
# se aplica la rotación A→D, B→A, C→B, D→C.
# Útil cuando el path cruza subunidades y la cadena "dominante" es B.
# Poner USE_CHAIN_ROTATION = False para desactivar completamente.
USE_CHAIN_ROTATION  = True
CHAIN_ROT_RESID_MIN = 400
CHAIN_ROT_RESID_MAX = 710


# --- Geometría del backbone GLY ---
BOND_N_CA    = 1.458   # Å, N-Cα (estadística PDB)
BOND_CA_C    = 1.523   # Å, Cα-C
BOND_C_O     = 1.231   # Å, C=O (enlace doble, más corto que enlace simple)
ANGLE_N_CA_C = 111.0   # grados, ángulo de enlace Cα en GLY

# --- Procesamiento paralelo ---
# USE_PARALLEL=True procesa los sistemas en paralelo (Windows-compatible).
# N_WORKERS=1 desactiva el paralelismo (útil para depurar).
USE_PARALLEL = True
N_WORKERS    = 4

# --- Archivos de salida opcionales ---
WRITE_TCL   = True   # Script VMD TCL para cada PDB
WRITE_PYMOL = True   # Script PyMOL PML para cada PDB (backup robusto)

# ==============================================================
# 🕒  LOGGING
# ==============================================================

def setup_logging(log_dir: Path) -> logging.Logger:
    """
    Configura logging a consola y a un archivo con timestamp.
    Devuelve el logger raíz del script.
    """
    log_dir.mkdir(parents=True, exist_ok=True)
    stamp    = datetime.now().strftime("%Y%m%d_%H%M%S")
    log_file = log_dir / f"{VERSION}_log_{stamp}.txt"

    logger = logging.getLogger("p2pdb")
    logger.setLevel(logging.DEBUG)

    fmt = logging.Formatter(
        "[%(asctime)s] %(levelname)-8s %(message)s",
        datefmt="%Y-%m-%d %H:%M:%S"
    )

    ch = logging.StreamHandler()
    ch.setLevel(logging.INFO)
    ch.setFormatter(fmt)

    fh = logging.FileHandler(log_file, encoding="utf-8")
    fh.setLevel(logging.DEBUG)
    fh.setFormatter(fmt)

    logger.addHandler(ch)
    logger.addHandler(fh)
    return logger


def get_logger(name: str = "p2pdb") -> logging.Logger:
    """Obtiene (o crea) el logger con nombre dado."""
    return logging.getLogger(name)

# ==============================================================
# 🔧  LECTURA DEL PDB FUENTE
# ==============================================================

def load_cb_map(pdb_file: Path) -> dict:
    """
    Lee el PDB y construye un diccionario de coordenadas CB.

    Usa CB para todos los residuos excepto GLY, para el cual usa CA
    (GLY no tiene cadena lateral, su Cα equivale al centroide).

    Solo conserva la primera conformación alternada (alt loc A o blank).

    Returns:
        dict { "ChainResno" : np.array([x, y, z]) }
        Ejemplo: {"A643": array([10.1, 5.2, 3.4])}
    """
    cb_map = {}
    with open(pdb_file, "r", errors="replace") as f:
        for line in f:
            rec = line[:6].rstrip()
            if rec not in ("ATOM", "HETATM"):
                continue

            atom_name = line[12:16].strip()
            alt_loc   = line[16:17].strip()   # conformación alternada
            res_name  = line[17:20].strip()
            chain     = line[21:22].strip()

            # Saltar conformaciones alternadas que no sean la primera
            if alt_loc not in ("", "A", " "):
                continue

            try:
                res_seq = int(line[22:26])
                x = float(line[30:38])
                y = float(line[38:46])
                z = float(line[46:54])
            except ValueError:
                continue

            is_cb     = (atom_name == "CB")
            is_gly_ca = (res_name == "GLY" and atom_name == "CA")

            if is_cb or is_gly_ca:
                key = f"{chain}{res_seq}"
                if key not in cb_map:  # primera conformación gana
                    cb_map[key] = np.array([x, y, z], dtype=np.float64)

    return cb_map

# ==============================================================
# 🔧  PARSING DE PATHS
# ==============================================================

def parse_path(path_str: str) -> list:
    """
    Parsea un string de path en formato "A511-A551-D590-D637" en
    una lista de tuplas (chain: str, resno: int).

    Acepta:
        "A511-A551-D590"  → [("A", 511), ("A", 551), ("D", 590)]
        "511-551-590"     → [("A", 511), ("A", 551), ("A", 590)]  (sin cadena → A)
    """
    nodes = []
    for token in re.split(r"[-\s]+", path_str.strip()):
        token = token.strip()
        if not token:
            continue
        m = re.match(r"^([A-Za-z])(\d+)$", token)
        if m:
            nodes.append((m.group(1), int(m.group(2))))
            continue
        m2 = re.match(r"^(\d+)$", token)
        if m2:
            nodes.append(("A", int(m2.group(1))))
    return nodes


# ==============================================================
# 🔧  ROTACIÓN DE CADENAS (CORRECCIÓN HOMOTETRÁMÉRO)
# ==============================================================

# Mapa de rotación: A→D, B→A, C→B, D→C
# Equivale a un desplazamiento de -1 en el orden ABCD (módulo 4).
CHAIN_ROTATION_MAP = {'A': 'B', 'B': 'C', 'C': 'D', 'D': 'A'}


def needs_chain_rotation(nodes: list,
                          resid_min: int = CHAIN_ROT_RESID_MIN,
                          resid_max: int = CHAIN_ROT_RESID_MAX) -> bool:
    """
    Detecta si un path necesita rotación de cadenas.

    Condición: existe al menos una transición de cadena A a cadena D
    donde el resid del nodo en A está estrictamente entre resid_min
    y resid_max. Esto identifica paths que "viven" mayoritariamente
    en la subunidad B pero tienen su origen nominal en A.

    Args:
        nodes    : lista de (chain, resid) del path
        resid_min: límite inferior del rango de resid (exclusivo)
        resid_max: límite superior del rango de resid (exclusivo)

    Returns:
        True si se debe aplicar la rotación.
    """
    for i in range(len(nodes) - 1):
        chain_curr, resid_curr = nodes[i]
        chain_next, _          = nodes[i + 1]
        if (chain_curr == 'A'
                and chain_next == 'D'
                and (resid_curr < resid_min or resid_curr > resid_max)):
            return True
    return False


def apply_chain_rotation(nodes: list) -> list:
    """
    Aplica la rotación de cadenas A→D, B→A, C→B, D→C a todos los
    nodos del path. Cadenas no incluidas en CHAIN_ROTATION_MAP
    (E, F, ...) se mantienen sin cambios.

    Args:
        nodes: lista de (chain, resid)

    Returns:
        Nueva lista con las cadenas rotadas.
    """
    return [
        (CHAIN_ROTATION_MAP.get(chain, chain), resid)
        for chain, resid in nodes
    ]


# ==============================================================
# 🔧  INTERPOLACIÓN
# ==============================================================

def interpolate_guides(cb_coords: list, max_step: float = INTERP_MAX_STEP) -> list:
    """
    Expande la lista de posiciones guía insertando puntos intermedios
    equidistantes (~3.8 Å) entre pares que superan max_step.

    Los puntos interpolados permiten que VMD conecte los residuos
    del path sin que aparezcan cortes en el ribbon.

    Args:
        cb_coords : lista de np.array([x, y, z])
        max_step  : distancia máxima permitida antes de interpolar

    Returns:
        Lista expandida de posiciones (incluyendo originales)
    """
    if len(cb_coords) < 2:
        return list(cb_coords)

    result = []
    for i in range(len(cb_coords) - 1):
        p1 = cb_coords[i]
        p2 = cb_coords[i + 1]
        d  = float(np.linalg.norm(p2 - p1))

        result.append(p1.copy())

        if d > max_step:
            # Número de pasos de ~3.8 Å para cubrir la distancia
            n_steps = max(2, int(np.ceil(d / 3.8)))
            for k in range(1, n_steps):
                frac = k / n_steps
                result.append(p1 + frac * (p2 - p1))

    result.append(cb_coords[-1].copy())
    return result

# ==============================================================
# 🔧  CONSTRUCCIÓN DEL BACKBONE GLY CON GEOMETRÍA β
# ==============================================================

def build_beta_backbone(guide_positions: list) -> list:
    """
    Construye N, CA, C, O para cada posición guía.

    CA se coloca exactamente en la posición guía (CB del residuo real).

    N y C se colocan a lo largo de los vectores inter-CA reales,
    escalados para garantizar C(i)–N(i+1) = BOND_C_N = 1.329 Å siempre,
    independientemente de la distancia CA-CA.

    Demostración geométrica:
      C(i)   = CA(i)   + c_frac * t_hat
      N(i+1) = CA(i+1) - n_frac * t_hat
      Con c_frac + n_frac = d - BOND_C_N:
      |C(i) - N(i+1)| = |−d*t + (c_frac+n_frac)*t| = BOND_C_N ✓

    O se coloca usando parallel transport del vector perpendicular
    para mantener la orientación del ribbon coherente a lo largo
    del path (evita flips que cortarían el ribbon en VMD).
    """
    BOND_C_N = 1.329   # Å, enlace peptídico

    n = len(guide_positions)
    if n < 2:
        return []

    residues    = []
    ang_cao     = np.radians(120.8)
    prev_normal = None

    # --- Precomputar vectores y fracciones por segmento inter-CA ---
    # segs[i] = (t_hat, d, c_frac, n_frac) para segmento CA(i)→CA(i+1)
    segs = []
    for i in range(n - 1):
        vec = guide_positions[i + 1] - guide_positions[i]
        d   = float(np.linalg.norm(vec))
        t   = vec / d if d > 1e-8 else np.array([1., 0., 0.])
        # Escalar para que C-N(next) = BOND_C_N exactamente.
        # Si d < BOND_C_N + 0.5, forzar un mínimo para evitar fracciones negativas.
        eff = max(d, BOND_C_N + 0.5)
        total  = BOND_CA_C + BOND_N_CA          # 2.981 Å
        c_frac = (eff - BOND_C_N) * BOND_CA_C / total
        n_frac = (eff - BOND_C_N) * BOND_N_CA / total
        segs.append((t, d, c_frac, n_frac))

    # --- Inicializar vector perpendicular para el átomo O ---
    t0   = segs[0][0]
    ref  = np.array([1., 0., 0.]) if abs(t0[0]) < 0.9 else np.array([0., 1., 0.])
    perp = np.cross(t0, ref)
    perp = perp / np.linalg.norm(perp)

    # --- Construir un residuo por CA ---
    for i in range(n):
        ca = guide_positions[i].copy()

        # ---- Posición de N ----
        # N(i) viene del segmento (i-1 → i): está a n_frac[i-1] Å atrás del CA
        if i > 0:
            t_back, _, _, n_frac = segs[i - 1]
            n_pos = ca - n_frac * t_back
        else:
            # Primer residuo: N colocado en dirección opuesta al primer segmento
            n_pos = ca - BOND_N_CA * segs[0][0]

        # ---- Posición de C ----
        # C(i) va hacia el segmento (i → i+1): está a c_frac[i] Å adelante del CA
        if i < n - 1:
            t_fwd, _, c_frac, _ = segs[i]
            c_pos = ca + c_frac * t_fwd
        else:
            # Último residuo: C colocado en dirección del último segmento
            c_pos = ca + BOND_CA_C * segs[-1][0]

        # ---- Parallel transport del vector perp (para orientación del O) ----
        # Usar el tangente del segmento más cercano
        if 0 < i < n - 1:
            t_bi     = segs[i][0] + segs[i - 1][0]
            t_bi_len = float(np.linalg.norm(t_bi))
            t_curr   = t_bi / t_bi_len if t_bi_len > 1e-8 else segs[i][0]
        else:
            t_curr = segs[min(i, n - 2)][0]

        perp = perp - t_curr * float(np.dot(perp, t_curr))
        p_len = float(np.linalg.norm(perp))
        if p_len < 1e-8:
            ref  = np.array([1., 0., 0.]) if abs(t_curr[0]) < 0.9 else np.array([0., 1., 0.])
            perp = np.cross(t_curr, ref)
            p_len = float(np.linalg.norm(perp))
        perp = perp / p_len

        normal = np.cross(t_curr, perp)

        # Anti-flip de normal: evita q el vector N→O alterne de lado,
        # lo que haría aparecer flechas falsas a mitad del ribbon en VMD.
        if prev_normal is not None and float(np.dot(normal, prev_normal)) < 0.0:
            perp   = -perp
            normal = -normal
        prev_normal = normal.copy()

        # ---- Posición de O ----
        # Calcular normal desde el plano N-CA-C real de este residuo.
        # Más estable que el transporte paralelo de perp en curvas bruscas,
        # porque usa la geometría local exacta en vez de una dirección acumulada.
        nca_vec = n_pos - ca
        cac_vec = c_pos - ca
        nca_len = float(np.linalg.norm(nca_vec))
        cac_len = float(np.linalg.norm(cac_vec))
        if nca_len > 1e-8 and cac_len > 1e-8:
            plane_normal = np.cross(nca_vec / nca_len, cac_vec / cac_len)
            pn_len = float(np.linalg.norm(plane_normal))
            if pn_len > 1e-8:
                plane_normal = plane_normal / pn_len
                # Garantizar consistencia: no permitir flips respecto
                # al residuo anterior (causarían arrowheads falsos en VMD)
                if prev_normal is not None and float(np.dot(plane_normal, prev_normal)) < 0.0:
                    plane_normal = -plane_normal
                prev_normal = plane_normal.copy()
            else:
                plane_normal = prev_normal if prev_normal is not None else normal
        else:
            plane_normal = prev_normal if prev_normal is not None else normal

        c_ca_vec = c_pos - ca
        c_len    = float(np.linalg.norm(c_ca_vec))
        if c_len > 1e-8:
            c_hat = c_ca_vec / c_len
            cos_o = np.cos(np.pi - ang_cao)
            sin_o = np.sin(np.pi - ang_cao)
            o_dir = -cos_o * c_hat + sin_o * plane_normal
            o_pos = c_pos + BOND_C_O * o_dir
        else:
            o_pos = c_pos + BOND_C_O * plane_normal

        residues.append({"N": n_pos, "CA": ca, "C": c_pos, "O": o_pos})

    return residues

# ==============================================================
# 🔧  FORMATO PDB
# ==============================================================

def format_atom(serial: int, atom_name: str, chain: str, res_seq: int,
                x: float, y: float, z: float,
                b_factor: float = 0.0, seg_id: str = "") -> str:
    """
    Formato exacto de registro ATOM del PDB (80 columnas).

    Referencia de columnas (1-indexado):
      1-6   : "ATOM  "
      7-11  : serial (entero, derecha)
      12    : espacio
      13-16 : nombre del átomo (ver atom_fmt abajo)
      17    : indicador de localización alternada (espacio)
      18-20 : nombre del residuo ("GLY")
      21    : espacio
      22    : ID de cadena
      23-26 : número de secuencia del residuo (entero, derecha)
      27    : código de inserción (espacio)
      28-30 : espacios
      31-38 : coordenada X (8.3f)
      39-46 : coordenada Y (8.3f)
      47-54 : coordenada Z (8.3f)
      55-60 : ocupancia (6.2f)
      61-66 : factor de temperatura / B-factor (6.2f)
      67-72 : espacios
      73-76 : identificador de segmento (opcional)
    """
    atom_fmt = {
        "N":  " N  ",
        "CA": " CA ",
        "C":  " C  ",
        "O":  " O  ",
    }
    name_str = atom_fmt.get(atom_name, f" {atom_name:<3}")
    seg_str  = f"{seg_id:<4s}" if seg_id else "    "

    return (
        f"ATOM  "
        f"{serial:5d}"
        f" "
        f"{name_str}"
        f" "
        f"GLY"
        f" "
        f"{chain:1s}"
        f"{res_seq:4d}"
        f" "
        f"   "
        f"{x:8.3f}"
        f"{y:8.3f}"
        f"{z:8.3f}"
        f"  1.00"
        f"{b_factor:6.2f}"
        f"      "
        f"{seg_str}"
    )


def format_ter(serial: int, chain: str, res_seq: int) -> str:
    """Formato de registro TER (fin de cadena)."""
    return f"TER   {serial:5d}      GLY {chain:1s}{res_seq:4d}"


def format_sheet(strand_id: int, sheet_id: str, chain: str,
                 res_start: int, res_end: int) -> str:
    """
    Formato exacto de registro SHEET del PDB.

    VMD lee estos registros al cargar el archivo y asigna estructura E
    (hoja β) a los residuos indicados, SIN necesidad de recalcular DSSP.
    Esto es lo que permite dibujar el ribbon como flecha en NewCartoon.

    Referencia de columnas (1-indexado):
      1-6   : "SHEET "
      7     : espacio
      8-10  : número de hebra (derecha)
      11    : espacio
      12-14 : ID de la hoja (izquierda)
      15-16 : número total de hebras en la hoja
      17    : espacio
      18-20 : nombre residuo inicial ("GLY")
      21    : espacio
      22    : cadena inicial
      23-26 : número de residuo inicial
      27    : código de inserción (espacio)
      28    : espacio
      29-31 : nombre residuo final ("GLY")
      32    : espacio
      33    : cadena final
      34-37 : número de residuo final
      38    : código de inserción (espacio)
      39-40 : sentido (0 = primera hebra de la hoja)
    """
    # Cada cadena/path es su propia hoja de 1 sola hebra
    return (
        f"SHEET "
        f" "
        f"{strand_id:3d}"
        f" "
        f"{sheet_id:<3s}"
        f" 1"
        f" "
        f"GLY"
        f" "
        f"{chain:1s}"
        f"{res_start:4d}"
        f" "
        f" "
        f"GLY"
        f" "
        f"{chain:1s}"
        f"{res_end:4d}"
        f" "
        f" 0"
    )

# ==============================================================
# 🔧  ESCRITURA DEL PDB
# ==============================================================

def write_pdb(output_path: Path, chain_data: list,
              system_name: str, source_pdb: str,
              direction: str, target: str) -> bool:
    """
    Escribe el PDB completo con todos los paths de un archivo CSV.

    Estructura del output:
        REMARK ...
        SHEET records (uno por cadena, cubre residuos 1 a N)
        ATOM records (N, CA, C, O por residuo GLY)
        TER records (separan cadenas)
        END

    El orden SHEET-antes-ATOM es importante para que VMD lea
    la asignación de estructura secundaria antes de cargar los átomos.

    Args:
        chain_data : lista de dicts por cadena/path
        system_name: ej. "7LP9_W426A"
        source_pdb : ej. "7LP9"
        direction  : "fwd" o "rev"
        target     : "filter" o "gate"

    Returns:
        True si se escribió correctamente.
    """
    lines = []

    # --- REMARK header ---
    now = datetime.now().strftime("%Y-%m-%d %H:%M:%S")
    lines += [
        f"REMARK  {VERSION} — {now}",
        f"REMARK  Sistema: {system_name}  PDB fuente: {source_pdb}.pdb",
        f"REMARK  Dirección: {direction}  Target: {target}",
        f"REMARK  Columna de path: {PATH_COLUMN}",
        f"REMARK  Geometría: GLY backbone, ángulo N-CA-C={ANGLE_N_CA_C}°",
        f"REMARK  B-factor = rank del cluster (1 = mayor score)",
        f"REMARK  Cada cadena = un cluster path",
        f"REMARK  -----------------------------------------------",
    ]
    for cd in chain_data:
        if cd["residues"]:
            lines.append(
                f"REMARK  Cadena {cd['chain_id']}: "
                f"cluster {cd['cluster_id']} | "
                f"score {cd['score']:.1f} | "
                f"{cd['n_original']} residuos guía + "
                f"{cd['n_total'] - cd['n_original']} interp = "
                f"{cd['n_total']} total | "
                f"path: {cd['path_str']}"
            )
    lines.append("REMARK  -----------------------------------------------")

    # --- SHEET records ---
    # Deben ir ANTES de los ATOM records para que VMD los lea primero
    for cd in chain_data:
        n_res = len(cd["residues"])
        if n_res == 0:
            continue
        n_real = cd.get("n_real", n_res)   # excluye el dummy del SHEET
        lines.append(format_sheet(
            strand_id = cd["rank"],
            sheet_id  = f"S{cd['chain_id']}",
            chain     = cd["chain_id"],
            res_start = 1,
            res_end   = n_real
        ))

    # --- ATOM records ---
    serial = 1
    for cd in chain_data:
        chain_id = cd["chain_id"]
        residues = cd["residues"]
        b_factor = float(cd["mean_imp_score"])
        seg_id   = f"P{cd['rank']:02d}"

        if not residues:
            continue

        for res_idx, res_atoms in enumerate(residues):
            res_seq = res_idx + 1
            for atom_name in ("N", "CA", "C", "O"):
                xyz = res_atoms[atom_name]
                lines.append(format_atom(
                    serial    = serial,
                    atom_name = atom_name,
                    chain     = chain_id,
                    res_seq   = res_seq,
                    x         = float(xyz[0]),
                    y         = float(xyz[1]),
                    z         = float(xyz[2]),
                    b_factor  = b_factor,
                    seg_id    = seg_id
                ))
                serial += 1

        lines.append(format_ter(serial, chain_id, len(residues)))
        serial += 1
# --- CONECT records: enlaces C(i)–N(i+1) entre residuos consecutivos ---
    #
    # La distancia C(i)–N(i+1) depende de la separación CA–CA del path real.
    # Cuando dos residuos del path están > ~4.1 Å de distancia, el enlace
    # peptídico supera el umbral de detección automática de VMD (~1.8 Å) y
    # el backbone aparece cortado en ribbons, newcartoon y tube.
    #
    # Los registros CONECT especifican conectividad explícita independiente
    # de la distancia, solucionando el problema para todos los representaciones.
    #
    # Fórmula de serial por residuo i (0-indexed) en cadena con inicio S:
    #   N(i)  = S + 4*i
    #   CA(i) = S + 4*i + 1
    #   C(i)  = S + 4*i + 2
    #   O(i)  = S + 4*i + 3
    # Enlace interpeptídico: C(i) → N(i+1)
    #   C(i) serial = S + 4*i + 2
    #   N(i+1) serial = S + 4*i + 4

    conect_lines = []
    running_serial = 1
    for cd in chain_data:
        n_res = len(cd["residues"])
        if n_res >= 1:
            cs = running_serial
            for i in range(n_res):
                n_ser  = cs + 4 * i        # N(i)
                ca_ser = cs + 4 * i + 1    # CA(i)
                c_ser  = cs + 4 * i + 2    # C(i)
                o_ser  = cs + 4 * i + 3    # O(i)
                # Todos los enlaces intra-residuo explícitos
                conect_lines.append(f"CONECT{n_ser:5d}{ca_ser:5d}")
                conect_lines.append(f"CONECT{ca_ser:5d}{c_ser:5d}")
                conect_lines.append(f"CONECT{c_ser:5d}{o_ser:5d}")
                # Enlace peptídico inter-residuo
                if i < n_res - 1:
                    n_next = cs + 4 * (i + 1)
                    conect_lines.append(f"CONECT{c_ser:5d}{n_next:5d}")
        running_serial += 4 * n_res + 1
    lines += conect_lines
    # lines.append("END")  ← esta línea ya está a continuación


    lines.append("END")

    try:
        output_path.parent.mkdir(parents=True, exist_ok=True)
        output_path.write_text("\n".join(lines) + "\n", encoding="utf-8")
        return True
    except Exception as e:
        get_logger().error(f"    Error escribiendo PDB {output_path}: {e}")
        return False

# ==============================================================
# 🔧  TCL PARA VMD
# ==============================================================

def write_tcl(tcl_path: Path, pdb_path: Path, chain_data: list,
              system_name: str, direction: str, target: str) -> None:
    """
    Escribe el script TCL companion para VMD.

    Estrategia VMD:
      1. El PDB incluye registros SHEET → VMD los lee al cargar y
         asigna estructura E automáticamente.
      2. El TCL refuerza la asignación E por si VMD recalcula SS.
      3. CRÍTICO: NO se llama 'mol ssrecalc'. Esa función ejecuta
         DSSP, que asigna coil (C) a hebras aisladas sin H-bonds
         inter-hebra, borrando la asignación E.
      4. Se crea una representación NewCartoon por cadena para
         permitir mostrar/ocultar paths individuales.

    Los colores son ColorIDs fijos de VMD (0-32, sin blanco=8).
    """
    chains_active = [cd for cd in chain_data if cd["residues"]]
#    color_cycle   = [0,1,2,3,4,5,6,7,9,10,11,12,13,14,15,17,18,19,20,21,22,23]
    pdb_stem      = pdb_path.stem
    pdb_name      = pdb_path.name

    with open(tcl_path, "w", encoding="utf-8") as con:
        def w(line=""):
            con.write(line + "\n")

        w("# " + "=" * 57)
        w(f"# VMD TCL — {system_name} | {direction} | {target}")
        w(f"# Generado por {VERSION}")
        w("# " + "=" * 57)
        w("# INSTRUCCIONES DE USO:")
        w("# 1. Abrir VMD")
        w("# 2. File > New Molecule > cargar PDB original (7LP9.pdb, etc.)")
        w("# 3. File > New Molecule > cargar el ribbon PDB como")
        w(f"#    molecula SEPARADA: {pdb_name}")
        w("# 4. Extensions > Tk Console")
        w("# 5. En la consola, escribir:")
        tcl_path_fwd = str(tcl_path).replace("\\", "/")
        w(f"#       source {{{tcl_path_fwd}}}")
        w("# " + "=" * 57)
        w()

        # Detectar molécula ribbon por nombre de archivo
        w("# --- Detectar la molécula ribbon por nombre de archivo ---")
        w("set path_mol -1")
        w("foreach mid [molinfo list] {")
        w("  set fname [lindex [molinfo $mid get filename] 0]")
        w(f"  if {{[string match *{pdb_stem}* $fname]}} {{")
        w("    set path_mol $mid")
        w("  }")
        w("}")
        w("if {$path_mol == -1} {")
        w("  puts {AVISO: no se encontro la molecula ribbon por nombre.}")
        w("  puts {       Usando la molecula top. Ajusta path_mol si es incorrecto.}")
        w("  set path_mol [molinfo top]")
        w("}")
        w("puts [format {Molecula ribbon detectada: mol %d} $path_mol]")
        w()

        # Reforzar estructura E SIN ssrecalc
        w("# --- Reforzar estructura E (hoja-beta) en todos los residuos ---")
        w("# Los registros SHEET del PDB ya asignan E al cargar.")
        w("# Esta linea lo refuerza por si acaso.")
        w("# CRITICO: NO llamar 'mol ssrecalc' despues de esto.")
        w("# mol ssrecalc ejecuta DSSP y reemplaza E con coil (C)")
        w("# en hebras aisladas sin pares H-bond inter-hebra.")
      # --- Eliminar representacion Lines por defecto ---
        w("mol delrep 0 $path_mol")
        w()

        # --- Representacion NewCartoon por cadena/path ---
        w("# --- Una representacion NewCartoon para todos los paths ---")
        w("# Coloreado por Chain: cada cadena (= cada cluster) recibe")
        w("# un color VMD distinto automaticamente.")
        w("mol representation NewCartoon 0.3 10 4.1 0")
        w("mol color Chain")
#        w("mol selection {all}")

        # Generar selección que excluye el residuo dummy terminal de cada cadena.
        # El dummy PERMANECE en la cadena del PDB (para dar dirección al spline),
        # pero al excluirlo de la selección VMD trata el último residuo real (643)
        # como el fin del segmento E → dibuja flecha → en lugar de nail-head.
        # Si se asignara como C (coil), el segmento E terminaría en E→C → nail-head.
        # Si se incluyera en la selección como E, VMD dibujaría la flecha
        # en el dummy (más allá del target).
        if ADD_TERMINAL_DUMMY and any(
            cd["n_total"] > cd.get("n_real", cd["n_total"])
            for cd in chains_active
        ):
            parts = [
                f"not (chain {cd['chain_id']} and resid {cd.get('n_real', cd['n_total']) + 1})"
                for cd in chains_active
                if cd["n_total"] > cd.get("n_real", cd["n_total"])
            ]
            sel_str = " and ".join(parts)
        else:
            sel_str = "all"
        w(f"mol selection {{{sel_str}}}")
 
        w("mol material Opaque")
        w("mol addrep $path_mol")
        w()
        w("# --- Asignar estructura E DESPUES de mol addrep ---")
        w("# mol addrep puede disparar STRIDE internamente en algunas")
        w("# versiones de VMD, sobrescribiendo la asignacion E previa.")
        w("# CRITICO: NO llamar mol ssrecalc ni mol modrep despues de esto.")
        w("# mol modrep sobreescribiria la representacion con el estado")
        w("# global actual de VMD (ultimo chain/color usado), borrando")
        w("# la seleccion 'all' y el color Chain.")
        w("set all_sel [atomselect $path_mol all]")
        w("$all_sel set structure E")
        w("$all_sel delete")
        w()

        w("# display update fuerza un redibujado sin tocar ningun parametro")
        w("# de representacion.")
        w("display update")
        w()
        w("puts {Estructura E asignada. Ribbon coloreado por chain.}")
        w()

        # Comandos útiles
        w("# " + "=" * 57)
        w("# REFERENCIA DE CHAINS:")
        for cd in chains_active:
            w(f"# chain {cd['chain_id']} → cluster {cd['cluster_id']} | {cd['path_str']}")
        w("# " + "=" * 57)
        w("#")
        w("# COMANDOS UTILES EN LA CONSOLA:")
        w("# Ocultar/mostrar un path (N = indice de representacion, desde 0):")
        w("#   mol showrep $path_mol 0 0   <- ocultar rep 0 (chain A)")
        w("#   mol showrep $path_mol 0 1   <- mostrar rep 0")
        w("#")
        w("# Cambiar grosor del ribbon (para todas las reps):")
        w("#   for {set i 0} {$i < [molinfo $path_mol get numreps]} {incr i} {")
        w("#     mol representation NewCartoon 0.5 10 4.1 0")
        w("#     mol modrep $i $path_mol")
        w("#   }")
        w("#")
        w("# Colorear por rank del cluster (B-factor, 1=mejor):")
        w("#   for {set i 0} {$i < [molinfo $path_mol get numreps]} {incr i} {")
        w("#     mol color Beta")
        w("#     mol modrep $i $path_mol")
        w("#   }")
        w("# " + "=" * 57)
        w()
        w(f"puts {{Listo: {len(chains_active)} paths como NewCartoon. Ver comentarios arriba.}}")


# ==============================================================
# 🔧  PML PARA PYMOL (BACKUP)
# ==============================================================

def write_pymol(pml_path: Path, pdb_path: Path, chain_data: list,
                system_name: str, direction: str, target: str) -> None:
    """
    Escribe un script PyMOL como alternativa robusta a VMD.

    PyMOL con 'cartoon_trace_atoms=1' dibuja el cartoon directamente
    desde los CA, sin depender de la asignación de estructura secundaria.
    Esto garantiza que todos los paths aparezcan correctamente.
    """
    chains_active = [cd for cd in chain_data if cd["residues"]]
    obj           = pdb_path.stem
    pdb_fwd       = str(pdb_path).replace("\\", "/")

    with open(pml_path, "w", encoding="utf-8") as f:
        def w(line=""):
            f.write(line + "\n")

        w(f"# PyMOL script — {system_name} | {direction} | {target}")
        w(f"# Generado por {VERSION}")
        w()
        w(f"load {pdb_fwd}, {obj}")
        w(f"hide everything, {obj}")
        w(f"set cartoon_trace_atoms, 1, {obj}")
        w(f"set cartoon_flat_sheets, 1, {obj}")
        w(f"set cartoon_tube_radius, 0.3, {obj}")
        w(f"show cartoon, {obj}")
        w(f"spectrum chain, rainbow, {obj}")
        w()
        w("# Referencia de chains:")
        for cd in chains_active:
            w(f"# chain {cd['chain_id']}: cluster {cd['cluster_id']} | {cd['path_str']}")
        w()
        w("# Para mostrar/ocultar un path individual:")
        for cd in chains_active[:3]:
            ch = cd["chain_id"]
            w(f"# hide cartoon, {obj} and chain {ch}")
            w(f"# show cartoon, {obj} and chain {ch}")

# ==============================================================
# 🔧  PROCESAMIENTO DE UN ARCHIVO CSV
# ==============================================================

CHAIN_LETTERS = list(
    "ABCDEFGHIJKLMNOPQRSTUVWXYZabcdefghijklmnopqrstuvwxyz0123456789"
)


def process_csv_file(csv_path: Path, cb_map: dict, system_name: str,
                     sys_output: Path, logger: logging.Logger) -> dict:
    """
    Procesa un archivo CSV de clusters → PDB + TCL + PML.

    Naming de archivos de entrada:
        clusters_7LP9_W426A_fwd_filter.csv

    Naming de archivos de salida:
        {sys_output}/fwd/7LP9_W426A_filter_fwd.pdb
        {sys_output}/fwd/7LP9_W426A_filter_fwd.tcl
        {sys_output}/fwd/7LP9_W426A_filter_fwd.pml

    Returns:
        dict con estadísticas del procesamiento
    """
    fname     = csv_path.stem  # ej. "clusters_7LP9_W426A_fwd_filter"
    direction = ("fwd" if "_fwd_" in fname
                 else "rev" if "_rev_" in fname else "unk")
    target    = ("filter" if "_filter" in fname
                 else "gate" if "_gate" in fname else "unk")

    out_dir  = sys_output / direction
    out_stem = f"{system_name}_{target}_{direction}"
    out_pdb  = out_dir / f"{out_stem}.pdb"
    out_tcl  = out_dir / f"{out_stem}.tcl"
    out_pml  = out_dir / f"{out_stem}.pml"
    out_dir.mkdir(parents=True, exist_ok=True)

    logger.info(f"  CSV: {csv_path.name}  →  {out_pdb.name}")

    # --- Leer CSV ---
    try:
        with open(csv_path, newline="", encoding="utf-8", errors="replace") as f:
            rows = list(csv.DictReader(f))
    except Exception as e:
        logger.error(f"  No se pudo leer {csv_path.name}: {e}")
        return {}

    if not rows:
        logger.warning(f"  CSV vacío: {csv_path.name}")
        return {}

    if PATH_COLUMN not in rows[0]:
        available = list(rows[0].keys())
        logger.error(f"  Columna '{PATH_COLUMN}' no encontrada en {csv_path.name}")
        logger.error(f"  Columnas disponibles: {available}")
        return {}

    # --- Procesar cada cluster ---
    chain_data = []
    n_ok = n_skip = 0

    for rank, row in enumerate(rows):
        cluster_id = row.get("cluster_id", str(rank + 1))
        path_str   = row.get(PATH_COLUMN, "").strip()

        try:
            score = float(row.get("total_imp_score",
                          row.get("n_paths_total", 0)) or 0)
        except (ValueError, TypeError):
            score = 0.0

        try:
            mean_imp_score = float(row.get("mean_imp_score", 0) or 0)
        except (ValueError, TypeError):
            mean_imp_score = 0.0

        if not path_str:
            logger.warning(f"  Cluster {cluster_id}: path vacío, omitido")
            n_skip += 1
            continue

        # Parsear nodos del path
        nodes = parse_path(path_str)
        # Rotación de cadenas (homotetrámero): detectar y aplicar si procede
        if USE_CHAIN_ROTATION and len(nodes) >= 2:
            if needs_chain_rotation(nodes):
                nodes = apply_chain_rotation(nodes)
                logger.debug(f"  Cluster {cluster_id}: rotación de cadenas A→D aplicada")
                
                
        if len(nodes) < 2:
            logger.warning(f"  Cluster {cluster_id}: < 2 nodos válidos, omitido")
            n_skip += 1
            continue

        # Obtener coordenadas CB
        cb_coords = []
        missing   = []
        for chain, resno in nodes:
            key = f"{chain}{resno}"
            if key in cb_map:
                cb_coords.append(cb_map[key])
            else:
                missing.append(key)

        if missing:
            logger.warning(f"  Cluster {cluster_id}: {len(missing)} nodos sin CB en PDB: {missing}")

        if len(cb_coords) < 2:
            logger.warning(f"  Cluster {cluster_id}: < 2 posiciones CB, omitido")
            n_skip += 1
            continue

        # Interpolar gaps grandes
        guide_pos  = interpolate_guides(cb_coords, INTERP_MAX_STEP)
        n_original = len(cb_coords)
        n_total    = len(guide_pos)
        n_interp   = n_total - n_original

        if n_interp > 0:
            logger.debug(f"  Cluster {cluster_id}: {n_interp} GLYs interpolados")
            
            
# Añadir residuo dummy terminal para que VMD dibuje el último
        # nodo del path con su arrowhead completo.
        # VMD necesita un nodo "más allá" del último residuo para
        # calcular la dirección de salida del ribbon; sin él, el último
        # residuo real no se renderiza visualmente.
        if ADD_TERMINAL_DUMMY and len(guide_pos) >= 2:
            last_step = guide_pos[-1] - guide_pos[-2]
            step_len  = float(np.linalg.norm(last_step))
            step_dir  = (last_step / step_len
                         if step_len > 1e-8
                         else np.array([1.0, 0.0, 0.0]))
            guide_pos = guide_pos + [guide_pos[-1] + step_dir * 3.8]
            n_total += 1
            
            
        # Construir backbone
        residues = build_beta_backbone(guide_pos)
        if not residues:
            logger.warning(f"  Cluster {cluster_id}: backbone vacío, omitido")
            n_skip += 1
            continue

        # Asignar letra de cadena
        if rank >= len(CHAIN_LETTERS):
            logger.warning(f"  Cluster {cluster_id}: supera máximo de cadenas, omitido")
            n_skip += 1
            continue

        chain_data.append({
            "chain_id":  CHAIN_LETTERS[rank],
            "rank":      rank + 1,
            "cluster_id": cluster_id,
            "path_str":  path_str,
            "score":     score,
            "mean_imp_score": mean_imp_score,
            "residues":  residues,
            "n_original": n_original,
            "n_total":   n_total,
            'n_real': (n_total - 1) if ADD_TERMINAL_DUMMY else n_total,   # ← agregar esta línea
        })
        logger.debug(
            f"  Cluster {cluster_id} → chain {CHAIN_LETTERS[rank]}: "
            f"{n_original} guía + {n_interp} interp = {n_total} residuos"
        )
        n_ok += 1

    if not chain_data:
        logger.warning(f"  Sin cadenas válidas para {csv_path.name}")
        return {}

    # --- Escribir PDB ---
    source_pdb = system_name[:4]
    ok = write_pdb(out_pdb, chain_data, system_name, source_pdb, direction, target)
    if ok:
        logger.info(f"  PDB OK: {out_pdb.name}  ({n_ok} cadenas, {n_skip} omitidas)")
    else:
        logger.error(f"  FALLO escritura PDB: {out_pdb}")

    # --- Escribir TCL ---
    if WRITE_TCL:
        write_tcl(out_tcl, out_pdb, chain_data, system_name, direction, target)
        logger.info(f"  TCL OK: {out_tcl.name}")

    # --- Escribir PyMOL ---
    if WRITE_PYMOL:
        write_pymol(out_pml, out_pdb, chain_data, system_name, direction, target)
        logger.info(f"  PML OK: {out_pml.name}")

    return {
        "system":    system_name,
        "direction": direction,
        "target":    target,
        "n_ok":      n_ok,
        "n_skip":    n_skip,
        "pdb":       str(out_pdb),
    }

# ==============================================================
# 🔧  PROCESAMIENTO DE UN SISTEMA COMPLETO
# ==============================================================

def process_system(args: tuple) -> list:
    """
    Procesa todos los CSV de un directorio de sistema.

    Esta función se ejecuta en procesos hijos (parallel). El logger
    se reconstruye aquí porque los objetos logging no son serializables
    para multiprocessing en Windows.

    Args:
        args: (system_dir: Path, base_output: Path)

    Returns:
        Lista de dicts resultado por archivo CSV procesado.
    """
    system_dir, base_output = args

    # Logger del subproceso (solo consola, para no colisionar archivos)
    logger = logging.getLogger(f"p2pdb.{system_dir.name}")
    if not logger.handlers:
        h = logging.StreamHandler()
        h.setFormatter(logging.Formatter(
            "[%(asctime)s] %(levelname)-8s [%(name)s] %(message)s",
            datefmt="%Y-%m-%d %H:%M:%S"
        ))
        logger.addHandler(h)
    logger.setLevel(logging.DEBUG)

    system_name = system_dir.name
    pdb_code    = system_name[:4].upper()
    pdb_file    = Path(PDB_DIR) / f"{pdb_code}.pdb"

    logger.info(f"SISTEMA: {system_name}  |  PDB: {pdb_file.name}")

    if not pdb_file.exists():
        logger.error(f"  PDB no encontrado: {pdb_file}")
        return []

    try:
        cb_map = load_cb_map(pdb_file)
        logger.info(f"  CB map: {len(cb_map)} átomos cargados")
    except Exception as e:
        logger.error(f"  Error cargando PDB: {e}")
        return []

    auto_cluster_dir = system_dir / "auto_cluster"
    if not auto_cluster_dir.exists():
        logger.warning(f"  Directorio auto_cluster/ no encontrado en {system_name}")
        return []

    csv_files = sorted(auto_cluster_dir.glob("clusters_*.csv"))
    if not csv_files:
        logger.warning(f"  No se encontraron CSV en {auto_cluster_dir}")
        return []

    logger.info(f"  {len(csv_files)} archivos CSV encontrados")

    sys_output = base_output / system_name
    results    = []

    for csv_path in csv_files:
        try:
            result = process_csv_file(csv_path, cb_map, system_name, sys_output, logger)
            if result:
                results.append(result)
        except Exception as e:
            logger.error(f"  ERROR en {csv_path.name}: {e}", exc_info=True)

    logger.info(f"SISTEMA {system_name}: {len(results)}/{len(csv_files)} CSV procesados")
    return results

# ==============================================================
# 🚀  MAIN
# ==============================================================

def main():

    # Logging principal
    log_dir = BASE_OUTPUT_DIR / "logs"
    logger  = setup_logging(log_dir)

    t_start = datetime.now()

    logger.info("=" * 60)
    logger.info(f"INICIO  {VERSION}  —  {t_start:%Y-%m-%d %H:%M:%S}")
    logger.info("=" * 60)
    logger.info(f"Input dir    : {BASE_INPUT_DIR}")
    logger.info(f"PDB dir      : {PDB_DIR}")
    logger.info(f"Output dir   : {BASE_OUTPUT_DIR}")
    logger.info(f"Columna path : {PATH_COLUMN}")
    logger.info(f"Interp. max  : {INTERP_MAX_STEP} Å")
    logger.info(f"Paralelo     : {USE_PARALLEL}  ({N_WORKERS} workers)")
    logger.info(f"Escribe TCL  : {WRITE_TCL}   |  Escribe PML: {WRITE_PYMOL}")
    logger.info("=" * 60)

    # Descubrir directorios de sistema
    if not BASE_INPUT_DIR.exists():
        logger.error(f"Directorio de entrada no encontrado: {BASE_INPUT_DIR}")
        sys.exit(1)

    system_dirs = sorted([
        d for d in BASE_INPUT_DIR.iterdir()
        if d.is_dir() and (d / "auto_cluster").exists()
    ])

    if not system_dirs:
        logger.error("No se encontraron directorios de sistema con auto_cluster/")
        sys.exit(1)

    logger.info(f"Sistemas detectados: {len(system_dirs)}")
    for d in system_dirs:
        logger.info(f"  {d.name}")

    # Preparar argumentos para ejecución paralela
    args_list = [(d, BASE_OUTPUT_DIR) for d in system_dirs]

    # Ejecutar
    all_results = []

    if USE_PARALLEL and N_WORKERS > 1 and len(system_dirs) > 1:
        logger.info(f"Ejecutando en paralelo ({N_WORKERS} workers)...")
        with ProcessPoolExecutor(max_workers=N_WORKERS) as executor:
            futures = {
                executor.submit(process_system, args): args[0].name
                for args in args_list
            }
            for future in as_completed(futures):
                sys_name = futures[future]
                try:
                    results = future.result()
                    all_results.extend(results)
                    logger.info(f"Completado: {sys_name}  ({len(results)} archivos)")
                except Exception as e:
                    logger.error(f"FALLO: {sys_name} — {e}", exc_info=True)
    else:
        logger.info("Ejecutando en modo secuencial...")
        for args in args_list:
            try:
                results = process_system(args)
                all_results.extend(results)
            except Exception as e:
                logger.error(f"FALLO: {args[0].name} — {e}", exc_info=True)

    # Resumen final
    t_end   = datetime.now()
    elapsed = (t_end - t_start).total_seconds()

    logger.info("=" * 60)
    logger.info(f"FIN  {VERSION}  —  {t_end:%Y-%m-%d %H:%M:%S}")
    logger.info(f"Tiempo total: {elapsed:.1f} s")
    logger.info(f"PDB escritos: {len(all_results)}")
    for r in all_results:
        logger.info(
            f"  {r.get('system')} | {r.get('direction')} | {r.get('target')}"
            f" → {r.get('n_ok')} cadenas ({r.get('n_skip')} omitidas)"
        )
    logger.info("=" * 60)


# ==============================================================
# PUNTO DE ENTRADA
# Necesario en Windows para ProcessPoolExecutor
# ==============================================================

if __name__ == "__main__":
    main()
# -*- coding: utf-8 -*-
"""
Created on Sun May 31 18:53:54 2026

@author: denyc
"""
