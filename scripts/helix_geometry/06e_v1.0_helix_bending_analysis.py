"""
================================================================================
helix_bending_analysis.py
================================================================================
Version     : 1.0.0
Author      : DenyCB
Description: Analysis of bending angles of transmembrane helices
             in molecular dynamics simulations of TRPV1 (WT and mutants).
             Calculates local per-residue angles (sliding window of Cα atoms) and
             global per-helix angles (triangle and/or half-vector methods),
             for each frame, chain, replicate, and system. Generates heatmaps (resid × frame)
             and global-angle time-series plots, using a color/Y-axis scale normalized
             to the global maximum for each helix across all systems. Intermediate data
             are saved as CSV files to avoid recalculation when only the visualization
             is changed.

Dependencies: MDAnalysis, numpy, pandas, matplotlib, multiprocessing
              (all compatible with Windows 10 and Spyder)
================================================================================
"""

# ==============================================================================
# SECCIÓN 0 — IMPORTACIONES
# ==============================================================================
import os
import sys
import time
import logging
import warnings
import itertools
import traceback
import multiprocessing
from datetime import datetime
from functools import partial

import numpy as np
import pandas as pd
import matplotlib
matplotlib.use('Agg')          # backend sin GUI, compatible con multiprocessing
import matplotlib.pyplot as plt
import matplotlib.colors as mcolors
from matplotlib.backends.backend_pdf import PdfPages

import MDAnalysis as mda
from MDAnalysis.analysis import align

warnings.filterwarnings('ignore')   # suprimir avisos menores de MDAnalysis


# ==============================================================================
# SECCIÓN 1 — CONFIGURACIÓN PRINCIPAL (EDITAR AQUÍ)
# ==============================================================================

# ------------------------------------------------------------------------------
# 1.1  Rutas de carpetas
# ------------------------------------------------------------------------------
INPUT_DIR  = r"C:\DinamicasMoleculares\TRPV1_pipeline\input\md"
OUTPUT_DIR = r"C:\DinamicasMoleculares\TRPV1_pipeline\output\06e_helix_analysis"

# ------------------------------------------------------------------------------
# 1.2  Definición de sistemas y réplicas
#       Formato: {"PDB_MUT": ("PDB_ID", "MUT_ID")}
#       El script buscará archivos tipo  PDB_MUT_repN.pdb  y  PDB_MUT_repN.dcd
# ------------------------------------------------------------------------------
SYSTEMS = {
    "7LP9_WT"    : ("7LP9", "WT"),
    "7LP9_W426A" : ("7LP9", "W426A"),
    "7LP9_W697A" : ("7LP9", "W697A"),
    "7LP9_Y441A" : ("7LP9", "Y441A"),
    "7LPB_WT"    : ("7LPB", "WT"),
    "7LPB_W426A" : ("7LPB", "W426A"),
    "7LPB_W697A" : ("7LPB", "W697A"),
    "7LPB_Y441A" : ("7LPB", "Y441A"),
    "7LPC_WT"    : ("7LPC", "WT"),
    "7LPC_W426A" : ("7LPC", "W426A"),
    "7LPC_W697A" : ("7LPC", "W697A"),
    "7LPC_Y441A" : ("7LPC", "Y441A"),
}

REPLICAS = [1, 2, 3]           # índices de réplica

# ------------------------------------------------------------------------------
# 1.3  Cadenas (chainID en el PDB)
# ------------------------------------------------------------------------------
CHAINS = ["A", "B", "C", "D"]

# ------------------------------------------------------------------------------
# 1.4  Helix definition
# Format: {"Name": (start_resid, end_resid)}
# Residues are absolute and identical across all chains.
# Modify/add/remove helices here according to your system.
# ------------------------------------------------------------------------------
HELICES = {
    "TM1"  : (432, 450),
    "TM2"  : (460, 480),
    "TM3"  : (505, 525),
    "TM4"  : (535, 555),
    "TM5"  : (572, 592),
    "TM6"  : (638, 658),
    "TRPh" : (683, 700),
}

# ------------------------------------------------------------------------------
# 1.5  Angle calculation parameters
# ------------------------------------------------------------------------------
# Window for local per-residue angle: number of neighbors on each side
# (total window = 2*LOCAL_WINDOW + 1 residues)
LOCAL_WINDOW = 2                # window of 5 (±2)

# Reference atom for the helix axis
BACKBONE_ATOM = "CA"           # alpha carbon

# Global angle methods (switch): at least one must be True
USE_TRIANGLE_METHOD = True     # triangle: end-to-end-maximum deviation
USE_MIDVECTOR_METHOD = True    # first/second half helix vectors

# ------------------------------------------------------------------------------
# 1.6  Parallelization parameters
# ------------------------------------------------------------------------------
N_WORKERS = 6                  # number of parallel processes (win10 compatible)

# ------------------------------------------------------------------------------
# 1.7  Simulation parameters
# ------------------------------------------------------------------------------
N_FRAMES = 500                 # number of frames per trajectory
DT_NS    = 0.1                 # time step between saved frames (ns)

# ------------------------------------------------------------------------------
# 1.8  Processing switches
# ------------------------------------------------------------------------------
RECALCULATE_ANGLES = True      # True: recalculate even if previous CSVs exist
                               # False: load existing CSVs if available
SAVE_CSV   = True              # save angle data to CSV
SAVE_PNG   = True              # save each figure as an individual PNG
SAVE_PDF   = True              # save all figures in a summary PDF

# ------------------------------------------------------------------------------
# 1.9  Visualization parameters
# ------------------------------------------------------------------------------
# Replica colors in time-series plots
REPLICA_COLORS = {1: "black", 2: "blue", 3: "red"}

# DPI of exported figures
FIGURE_DPI = 150

# Script version (for the log)
SCRIPT_VERSION = "1.0.0"
SCRIPT_NAME    = "helix_bending_analysis.py"

# ==============================================================================
# SECCIÓN 2 — CONFIGURACIÓN DE LOGGING
# ==============================================================================

def setup_logger(output_dir: str) -> logging.Logger:
    """
    Configura el logger para escribir simultáneamente en consola y en un
    archivo .log con timestamp en el nombre.

    Parámetros
    ----------
    output_dir : str
        Directorio donde se guardará el archivo de log.

    Retorna
    -------
    logging.Logger
        Logger configurado.
    """
    os.makedirs(output_dir, exist_ok=True)
    log_dir = os.path.join(output_dir, "logs")
    os.makedirs(log_dir, exist_ok=True)

    timestamp = datetime.now().strftime("%Y%m%d_%H%M%S")
    log_filename = os.path.join(log_dir, f"{SCRIPT_NAME.replace('.py','')}_{timestamp}.log")

    logger = logging.getLogger("helix_analysis")
    logger.setLevel(logging.DEBUG)

    fmt = logging.Formatter("[%(asctime)s] %(levelname)s — %(message)s",
                             datefmt="%Y-%m-%d %H:%M:%S")

    # Handler consola
    ch = logging.StreamHandler(sys.stdout)
    ch.setLevel(logging.INFO)
    ch.setFormatter(fmt)

    # Handler archivo
    fh = logging.FileHandler(log_filename, encoding="utf-8")
    fh.setLevel(logging.DEBUG)
    fh.setFormatter(fmt)

    logger.addHandler(ch)
    logger.addHandler(fh)

    return logger, log_filename


def log_config(logger: logging.Logger) -> None:
    """
    Escribe en el log un resumen de las variables de configuración principales.
    """
    logger.info("=" * 70)
    logger.info(f"  Script  : {SCRIPT_NAME}  v{SCRIPT_VERSION}")
    logger.info(f"  Inicio  : {datetime.now().strftime('%Y-%m-%d %H:%M:%S')}")
    logger.info("=" * 70)
    logger.info(f"  INPUT_DIR       : {INPUT_DIR}")
    logger.info(f"  OUTPUT_DIR      : {OUTPUT_DIR}")
    logger.info(f"  Sistemas        : {list(SYSTEMS.keys())}")
    logger.info(f"  Réplicas        : {REPLICAS}")
    logger.info(f"  Cadenas         : {CHAINS}")
    logger.info(f"  Hélices         : {list(HELICES.keys())}")
    logger.info(f"  LOCAL_WINDOW    : ±{LOCAL_WINDOW}  (ventana {2*LOCAL_WINDOW+1} residuos)")
    logger.info(f"  BACKBONE_ATOM   : {BACKBONE_ATOM}")
    logger.info(f"  USE_TRIANGLE    : {USE_TRIANGLE_METHOD}")
    logger.info(f"  USE_MIDVECTOR   : {USE_MIDVECTOR_METHOD}")
    logger.info(f"  N_WORKERS       : {N_WORKERS}")
    logger.info(f"  N_FRAMES        : {N_FRAMES}")
    logger.info(f"  RECALCULATE     : {RECALCULATE_ANGLES}")
    logger.info(f"  SAVE_CSV        : {SAVE_CSV}")
    logger.info(f"  SAVE_PNG        : {SAVE_PNG}")
    logger.info(f"  SAVE_PDF        : {SAVE_PDF}")
    logger.info("=" * 70)


# ==============================================================================
# SECCIÓN 3 — FUNCIONES DE CÁLCULO DE ÁNGULOS
# ==============================================================================

def angle_between_vectors(v1: np.ndarray, v2: np.ndarray) -> float:
    """
    Calcula el ángulo (en grados) entre dos vectores 3D.

    Parámetros
    ----------
    v1, v2 : np.ndarray de shape (3,)
        Vectores en coordenadas cartesianas.

    Retorna
    -------
    float
        Ángulo en grados [0°, 180°].
    """
    n1 = np.linalg.norm(v1)
    n2 = np.linalg.norm(v2)
    if n1 == 0 or n2 == 0:
        return np.nan
    cos_a = np.clip(np.dot(v1, v2) / (n1 * n2), -1.0, 1.0)
    return np.degrees(np.arccos(cos_a))


def local_bending_angles(ca_coords: np.ndarray, window: int = 2) -> np.ndarray:
    """
    Calcula el ángulo local de flexión (bending) para cada residuo de una
    hélice usando una ventana deslizante de Cα.

    Método: para el residuo central i, se forman dos vectores:
        v1 = CA[i]   - CA[i - window]
        v2 = CA[i + window] - CA[i]
    El ángulo de bending es el suplemento del ángulo entre v1 y v2
    (0° = hélice perfectamente recta en ese punto; 180° sería doblada 180°).
    Los residuos en los extremos (sin suficientes vecinos) reciben NaN.

    Parámetros
    ----------
    ca_coords : np.ndarray de shape (N_resid, 3)
        Coordenadas de los Cα de la hélice, en orden secuencial de resid.
    window : int
        Número de residuos a cada lado del residuo central.

    Retorna
    -------
    np.ndarray de shape (N_resid,)
        Ángulo local de bending por residuo (grados). NaN en los extremos.
    """
    n = len(ca_coords)
    angles = np.full(n, np.nan)
    for i in range(window, n - window):
        v1 = ca_coords[i]     - ca_coords[i - window]
        v2 = ca_coords[i + window] - ca_coords[i]
        # Ángulo de bending = desviación de la linealidad (180° - ángulo_v1_v2)
        angles[i] = 180.0 - angle_between_vectors(v1, v2)
    return angles


def global_angle_triangle(ca_coords: np.ndarray) -> float:
    """
    Calcula el ángulo global de bending de una hélice mediante el método
    del triángulo: usa los dos Cα extremos como base y el Cα más alejado
    del eje que los une como vértice.

    El eje de referencia es el vector que une el primer y último Cα de la
    hélice. Se proyecta cada Cα sobre ese eje y se calcula la distancia
    perpendicular. El Cα con mayor distancia perpendicular es el vértice.
    El ángulo reportado es el ángulo en ese vértice (∠ extremo1-vértice-extremo2),
    que representa la magnitud del doblez.

    Parámetros
    ----------
    ca_coords : np.ndarray de shape (N_resid, 3)
        Coordenadas de los Cα de la hélice.

    Retorna
    -------
    float
        Ángulo de bending global en grados. 0° = hélice perfectamente recta.
        NaN si hay menos de 3 residuos.
    """
    if len(ca_coords) < 3:
        return np.nan

    p0 = ca_coords[0]
    p1 = ca_coords[-1]
    axis = p1 - p0
    axis_len = np.linalg.norm(axis)

    if axis_len < 1e-8:
        return np.nan

    axis_unit = axis / axis_len

    # Distancia perpendicular de cada Cα al eje p0→p1
    perp_dists = []
    for p in ca_coords:
        proj = np.dot(p - p0, axis_unit) * axis_unit
        perp = (p - p0) - proj
        perp_dists.append(np.linalg.norm(perp))

    vertex_idx = int(np.argmax(perp_dists))
    vertex = ca_coords[vertex_idx]

    # Ángulo en el vértice: vectores vértice→extremo1 y vértice→extremo2
    v1 = p0 - vertex
    v2 = p1 - vertex
    angle_vertex = angle_between_vectors(v1, v2)

    # Convertir: 180° = recto, 0° = doblado al máximo
    # Retornamos desviación de la linealidad = 180° - ángulo_en_vértice
    return 180.0 - angle_vertex


def global_angle_midvector(ca_coords: np.ndarray) -> float:
    """
    Calcula el ángulo global de bending de una hélice mediante el método
    de los vectores mitad (half-helix vectors).

    La hélice se divide en dos mitades iguales. Se calcula el vector
    promedio de la primera mitad (N→C dentro de esa mitad) y el de la
    segunda mitad. El ángulo entre estos dos vectores es el ángulo de
    bending global: 0° = perfectamente recto.

    Parámetros
    ----------
    ca_coords : np.ndarray de shape (N_resid, 3)
        Coordenadas de los Cα de la hélice.

    Retorna
    -------
    float
        Ángulo de bending global en grados.
    """
    n = len(ca_coords)
    if n < 4:
        return np.nan

    mid = n // 2
    # Vector de la primera mitad: del primer al último Cα de esa mitad
    v1 = ca_coords[mid - 1] - ca_coords[0]
    # Vector de la segunda mitad: del primer al último Cα de esa mitad
    v2 = ca_coords[-1] - ca_coords[mid]

    return 180.0 - angle_between_vectors(v1, v2)


# ==============================================================================
# SECCIÓN 4 — FUNCIÓN DE ANÁLISIS DE UN SISTEMA/RÉPLICA
#             (ejecutada en paralelo por cada worker)
# ==============================================================================

def analyze_replica(task: dict) -> dict:
    """
    Carga la trayectoria de un sistema/réplica y calcula los ángulos locales
    y globales para todas las hélices y cadenas definidas.

    Esta función es la unidad de trabajo paralelo. Recibe un diccionario
    con los parámetros del sistema y devuelve un diccionario con los
    resultados crudos en arrays de NumPy.

    Parámetros
    ----------
    task : dict
        Claves requeridas:
        - system_name : str     (ej. "7LP9_WT")
        - replica     : int     (ej. 1)
        - pdb_path    : str     (ruta al PDB)
        - dcd_path    : str     (ruta al DCD)
        - helices     : dict    (nombre → (resid_ini, resid_fin))
        - chains      : list    (ej. ["A","B","C","D"])
        - atom        : str     (ej. "CA")
        - window      : int     (ventana local)
        - use_triangle: bool
        - use_midvec  : bool
        - n_frames    : int

    Retorna
    -------
    dict con claves:
        - system_name, replica
        - local_angles  : dict[helix][chain] = array (n_frames, n_resid)
        - global_triangle : dict[helix][chain] = array (n_frames,)   si use_triangle
        - global_midvec   : dict[helix][chain] = array (n_frames,)   si use_midvec
        - resids          : dict[helix] = lista de resids
        - error           : str o None
    """
    system_name = task["system_name"]
    replica     = task["replica"]
    pdb_path    = task["pdb_path"]
    dcd_path    = task["dcd_path"]
    helices     = task["helices"]
    chains      = task["chains"]
    atom        = task["atom"]
    window      = task["window"]
    use_tri     = task["use_triangle"]
    use_mid     = task["use_midvec"]

    result = {
        "system_name"    : system_name,
        "replica"        : replica,
        "local_angles"   : {},
        "global_triangle": {},
        "global_midvec"  : {},
        "resids"         : {},
        "error"          : None,
    }

    try:
        # ------------------------------------------------------------------
        # 4.1  Cargar universo MDAnalysis
        # ------------------------------------------------------------------
        u = mda.Universe(pdb_path, dcd_path)
        n_frames_traj = len(u.trajectory)

        for helix_name, (res_start, res_end) in helices.items():
            result["local_angles"][helix_name]    = {}
            result["global_triangle"][helix_name] = {}
            result["global_midvec"][helix_name]   = {}

            for chain in chains:
                # Selección: Cα de la cadena en el rango de resid de la hélice
                sel_str = (f"chainID {chain} and resid {res_start}:{res_end} "
                           f"and name {atom}")
                ag = u.select_atoms(sel_str)

                if len(ag) == 0:
                    # No se encontraron átomos: marcar como NaN
                    result["local_angles"][helix_name][chain]    = None
                    result["global_triangle"][helix_name][chain] = None
                    result["global_midvec"][helix_name][chain]   = None
                    continue

                # Resids reales presentes en la selección (para el eje X del heatmap)
                resids = ag.resids
                if helix_name not in result["resids"]:
                    result["resids"][helix_name] = resids.tolist()

                n_res    = len(ag)
                loc_mat  = np.full((n_frames_traj, n_res), np.nan)
                glob_tri = np.full(n_frames_traj, np.nan)
                glob_mid = np.full(n_frames_traj, np.nan)

                # ----------------------------------------------------------
                # 4.2  Iterar sobre frames
                # ----------------------------------------------------------
                for fi, ts in enumerate(u.trajectory):
                    ca_coords = ag.positions.copy()  # shape (n_res, 3)

                    # Ángulos locales (ventana deslizante)
                    loc_mat[fi, :] = local_bending_angles(ca_coords, window)

                    # Ángulo global — método triángulo
                    if use_tri:
                        glob_tri[fi] = global_angle_triangle(ca_coords)

                    # Ángulo global — método vectores mitad
                    if use_mid:
                        glob_mid[fi] = global_angle_midvector(ca_coords)

                result["local_angles"][helix_name][chain]    = loc_mat
                result["global_triangle"][helix_name][chain] = glob_tri
                result["global_midvec"][helix_name][chain]   = glob_mid

    except Exception as e:
        result["error"] = traceback.format_exc()

    return result


# ==============================================================================
# SECCIÓN 5 — GUARDADO Y CARGA DE CSV
# ==============================================================================

def save_angles_to_csv(result: dict, csv_dir: str) -> None:
    """
    Guarda los ángulos calculados de un sistema/réplica en archivos CSV.

    Estructura de archivos generados:
        csv_dir/
            {system}_{rep}_local_{helix}_{chain}.csv
            {system}_{rep}_global_triangle_{helix}_{chain}.csv
            {system}_{rep}_global_midvec_{helix}_{chain}.csv

    Cada CSV tiene filas = frames, columnas = resids (para local)
    o una sola columna "angle" (para global).

    Parámetros
    ----------
    result  : dict   Resultado devuelto por analyze_replica().
    csv_dir : str    Directorio de destino.
    """
    os.makedirs(csv_dir, exist_ok=True)
    sys_name = result["system_name"]
    rep      = result["replica"]

    for helix_name in result["local_angles"]:
        resids = result["resids"].get(helix_name, [])

        for chain in result["local_angles"][helix_name]:
            prefix = f"{sys_name}_rep{rep}"

            # Ángulos locales (matriz frames × resids)
            loc_data = result["local_angles"][helix_name][chain]
            if loc_data is not None:
                df_loc = pd.DataFrame(loc_data, columns=resids)
                df_loc.index.name = "frame"
                path_loc = os.path.join(csv_dir,
                    f"{prefix}_local_{helix_name}_{chain}.csv")
                df_loc.to_csv(path_loc)

            # Ángulo global triángulo
            tri_data = result["global_triangle"][helix_name][chain]
            if tri_data is not None:
                df_tri = pd.DataFrame(tri_data, columns=["angle"])
                df_tri.index.name = "frame"
                path_tri = os.path.join(csv_dir,
                    f"{prefix}_global_triangle_{helix_name}_{chain}.csv")
                df_tri.to_csv(path_tri)

            # Ángulo global vectores mitad
            mid_data = result["global_midvec"][helix_name][chain]
            if mid_data is not None:
                df_mid = pd.DataFrame(mid_data, columns=["angle"])
                df_mid.index.name = "frame"
                path_mid = os.path.join(csv_dir,
                    f"{prefix}_global_midvec_{helix_name}_{chain}.csv")
                df_mid.to_csv(path_mid)


def load_angles_from_csv(system_name: str, replica: int, csv_dir: str,
                         helices: dict, chains: list) -> dict:
    """
    Reconstruye el diccionario de resultados cargando los CSV previamente
    guardados por save_angles_to_csv().

    Parámetros
    ----------
    system_name : str
    replica     : int
    csv_dir     : str
    helices     : dict  (nombre → (resid_ini, resid_fin))
    chains      : list

    Retorna
    -------
    dict  con la misma estructura que analyze_replica()
    """
    result = {
        "system_name"    : system_name,
        "replica"        : replica,
        "local_angles"   : {},
        "global_triangle": {},
        "global_midvec"  : {},
        "resids"         : {},
        "error"          : None,
    }
    prefix = f"{system_name}_rep{replica}"

    for helix_name in helices:
        result["local_angles"][helix_name]    = {}
        result["global_triangle"][helix_name] = {}
        result["global_midvec"][helix_name]   = {}

        for chain in chains:
            # Local
            path_loc = os.path.join(csv_dir,
                f"{prefix}_local_{helix_name}_{chain}.csv")
            if os.path.exists(path_loc):
                df = pd.read_csv(path_loc, index_col=0)
                result["local_angles"][helix_name][chain] = df.values
                if helix_name not in result["resids"]:
                    result["resids"][helix_name] = [int(c) for c in df.columns]
            else:
                result["local_angles"][helix_name][chain] = None

            # Global triángulo
            path_tri = os.path.join(csv_dir,
                f"{prefix}_global_triangle_{helix_name}_{chain}.csv")
            if os.path.exists(path_tri):
                result["global_triangle"][helix_name][chain] = \
                    pd.read_csv(path_tri, index_col=0)["angle"].values
            else:
                result["global_triangle"][helix_name][chain] = None

            # Global vectores mitad
            path_mid = os.path.join(csv_dir,
                f"{prefix}_global_midvec_{helix_name}_{chain}.csv")
            if os.path.exists(path_mid):
                result["global_midvec"][helix_name][chain] = \
                    pd.read_csv(path_mid, index_col=0)["angle"].values
            else:
                result["global_midvec"][helix_name][chain] = None

    return result


# ==============================================================================
# SECCIÓN 6 — CÁLCULO DE RANGOS GLOBALES (para escalas unificadas)
# ==============================================================================

def compute_global_ranges(all_results: list, helices: dict,
                          chains: list, use_tri: bool, use_mid: bool) -> dict:
    """
    Recorre todos los resultados calculados y determina, para cada hélice,
    el valor máximo global de ángulo local y de ángulo global. Estos valores
    se usan para unificar las escalas de color (heatmaps) y de los ejes Y
    (gráficos de series temporales) entre todos los sistemas y réplicas.

    Parámetros
    ----------
    all_results : list   Lista de dicts devueltos por analyze_replica().
    helices     : dict
    chains      : list
    use_tri     : bool
    use_mid     : bool

    Retorna
    -------
    dict con estructura:
        {helix_name: {
            "local_max"   : float,
            "tri_max"     : float,
            "mid_max"     : float,
        }}
    """
    ranges = {h: {"local_max": 0.0, "tri_max": 0.0, "mid_max": 0.0}
              for h in helices}

    for res in all_results:
        if res["error"] is not None:
            continue
        for helix_name in helices:
            for chain in chains:
                # Ángulo local
                loc = res["local_angles"].get(helix_name, {}).get(chain)
                if loc is not None:
                    val = np.nanmax(loc)
                    if val > ranges[helix_name]["local_max"]:
                        ranges[helix_name]["local_max"] = val

                # Global triángulo
                if use_tri:
                    tri = res["global_triangle"].get(helix_name, {}).get(chain)
                    if tri is not None:
                        val = np.nanmax(tri)
                        if val > ranges[helix_name]["tri_max"]:
                            ranges[helix_name]["tri_max"] = val

                # Global vectores mitad
                if use_mid:
                    mid = res["global_midvec"].get(helix_name, {}).get(chain)
                    if mid is not None:
                        val = np.nanmax(mid)
                        if val > ranges[helix_name]["mid_max"]:
                            ranges[helix_name]["mid_max"] = val

    return ranges


# ==============================================================================
# SECCIÓN 7 — GENERACIÓN DE FIGURAS
# ==============================================================================

# Colormap personalizado: blanco → azul → rojo
HELIX_CMAP = mcolors.LinearSegmentedColormap.from_list(
    "helix_bending",
    [(0.0, "white"), (0.5, "royalblue"), (1.0, "red")]
)


def plot_heatmap(result: dict, helix_name: str, global_ranges: dict,
                 chains: list, dt_ns: float,
                 plot_dir: str, pdf_pages, save_png: bool) -> None:
    """
    Genera la figura de heatmap de ángulos locales para un sistema/réplica/hélice.

    Diseño:
        - Una figura con 4 paneles (uno por cadena, A–D), disposición 1×4
        - Eje X = resid (número de residuo dentro de la hélice)
        - Eje Y = tiempo (ns), convertido desde frame × dt_ns
        - Color = ángulo local de bending (grados)
        - Escala de color: 0° (blanco) → máximo global de esa hélice (rojo)
        - Barra de color compartida

    Parámetros
    ----------
    result       : dict   Resultado de analyze_replica()
    helix_name   : str
    global_ranges: dict   Rangos globales devueltos por compute_global_ranges()
    chains       : list
    dt_ns        : float  Paso temporal entre frames en nanosegundos
    plot_dir     : str    Carpeta de destino para PNG
    pdf_pages    : PdfPages o None
    save_png     : bool
    """
    sys_name  = result["system_name"]
    rep       = result["replica"]
    vmax      = global_ranges[helix_name]["local_max"]
    resids    = result["resids"].get(helix_name, [])
    n_frames  = None

    fig, axes = plt.subplots(1, 4, figsize=(16, 5), sharey=True)
    fig.suptitle(
        f"Ángulo local de bending — {helix_name} | {sys_name} | Réplica {rep}",
        fontsize=13, fontweight="bold"
    )

    im = None
    for ax, chain in zip(axes, chains):
        loc_data = result["local_angles"].get(helix_name, {}).get(chain)

        if loc_data is None or len(resids) == 0:
            ax.set_title(f"Cadena {chain}")
            ax.text(0.5, 0.5, "Sin datos", ha="center", va="center",
                    transform=ax.transAxes)
            continue

        # loc_data shape: (n_frames, n_resid)
        # Para imshow: filas = Y (frame/tiempo), columnas = X (resid)
        n_frames = loc_data.shape[0]
        time_axis = np.arange(n_frames) * dt_ns    # eje Y en ns

        im = ax.imshow(
            loc_data,                               # shape (n_frames, n_resid)
            aspect="auto",
            origin="lower",
            cmap=HELIX_CMAP,
            vmin=0.0,
            vmax=vmax if vmax > 0 else 1.0,
            extent=[resids[0] - 0.5, resids[-1] + 0.5,
                    time_axis[0],    time_axis[-1]],
            interpolation="nearest",
        )

        ax.set_title(f"Cadena {chain}", fontsize=11)
        ax.set_xlabel("Residuo (resid)", fontsize=9)

    axes[0].set_ylabel("Tiempo (ns)", fontsize=10)

    # Barra de color compartida
    if im is not None:
        cbar = fig.colorbar(im, ax=axes, orientation="vertical",
                            fraction=0.02, pad=0.02)
        cbar.set_label("Ángulo de bending (°)", fontsize=10)

    plt.tight_layout()

    # Guardar
    fname = f"{sys_name}_rep{rep}_{helix_name}_heatmap"
    if save_png:
        png_path = os.path.join(plot_dir, fname + ".png")
        fig.savefig(png_path, dpi=FIGURE_DPI, bbox_inches="tight")

    if pdf_pages is not None:
        pdf_pages.savefig(fig, bbox_inches="tight")

    plt.close(fig)


def plot_global_timeseries(system_name: str, helix_name: str,
                           system_results: list, global_ranges: dict,
                           chains: list, replicas: list,
                           dt_ns: float, method: str,
                           plot_dir: str, pdf_pages, save_png: bool) -> None:
    """
    Genera el gráfico de series temporales del ángulo global de bending para
    un sistema y hélice, con 12 líneas (4 cadenas × 3 réplicas).

    Diseño:
        - Una sola figura con un panel
        - Eje X = tiempo (ns)
        - Eje Y = ángulo global de bending (grados), escala unificada al
                  máximo global de esa hélice en todos los sistemas
        - Colores de réplica: {1: negro, 2: azul, 3: rojo}
        - Estilo de línea diferenciado por cadena (solid, dashed, dotted, dashdot)
        - Leyenda con réplica y cadena

    Parámetros
    ----------
    system_name    : str
    helix_name     : str
    system_results : list   Lista de resultados del sistema (una entrada por réplica)
    global_ranges  : dict
    chains         : list
    replicas       : list
    dt_ns          : float
    method         : str    "triangle" o "midvec"
    plot_dir       : str
    pdf_pages      : PdfPages o None
    save_png       : bool
    """
    ymax_key = "tri_max" if method == "triangle" else "mid_max"
    ymax     = global_ranges[helix_name][ymax_key]
    method_key = "global_triangle" if method == "triangle" else "global_midvec"
    method_label = "Triángulo" if method == "triangle" else "Vectores mitad"

    linestyles = ["-", "--", ":", "-."]   # uno por cadena

    fig, ax = plt.subplots(figsize=(10, 4))
    ax.set_title(
        f"Ángulo global de bending ({method_label}) — {helix_name} | {system_name}",
        fontsize=12, fontweight="bold"
    )

    for rep_result in system_results:
        rep = rep_result["replica"]
        color = REPLICA_COLORS.get(rep, "gray")
        n_frames = None

        for ci, chain in enumerate(chains):
            data = rep_result[method_key].get(helix_name, {}).get(chain)
            if data is None:
                continue
            if n_frames is None:
                n_frames = len(data)
                time_axis = np.arange(n_frames) * dt_ns

            ax.plot(
                time_axis, data,
                color=color,
                linestyle=linestyles[ci % len(linestyles)],
                linewidth=0.8,
                alpha=0.85,
                label=f"Rep{rep} Chain{chain}",
            )

    ax.set_xlabel("Tiempo (ns)", fontsize=10)
    ax.set_ylabel("Ángulo de bending (°)", fontsize=10)
    ax.set_ylim(0, ymax * 1.05 if ymax > 0 else 1.0)

    # Leyenda resumida (solo muestra una entrada por réplica y una por cadena)
    handles, labels = ax.get_legend_handles_labels()
    ax.legend(handles, labels, fontsize=6, ncol=3, loc="upper right",
              framealpha=0.6)

    plt.tight_layout()

    fname = f"{system_name}_{helix_name}_global_{method}"
    if save_png:
        png_path = os.path.join(plot_dir, fname + ".png")
        fig.savefig(png_path, dpi=FIGURE_DPI, bbox_inches="tight")

    if pdf_pages is not None:
        pdf_pages.savefig(fig, bbox_inches="tight")

    plt.close(fig)


# ==============================================================================
# SECCIÓN 8 — FUNCIÓN PRINCIPAL
# ==============================================================================

def main():
    # --------------------------------------------------------------------------
    # 8.1  Inicializar logger y directorios de salida
    # --------------------------------------------------------------------------
    logger, log_path = setup_logger(OUTPUT_DIR)
    log_config(logger)
    logger.info(f"Log guardado en: {log_path}")

    csv_dir  = os.path.join(OUTPUT_DIR, "csv")
    plot_dir = os.path.join(OUTPUT_DIR, "plots")
    os.makedirs(csv_dir,  exist_ok=True)
    os.makedirs(plot_dir, exist_ok=True)

    t_start = time.time()

    # --------------------------------------------------------------------------
    # 8.2  Construir lista de tareas (una por sistema × réplica)
    # --------------------------------------------------------------------------
    tasks = []
    for system_name in SYSTEMS:
        for rep in REPLICAS:
            pdb_path = os.path.join(INPUT_DIR, f"{system_name}_rep{rep}.pdb")
            dcd_path = os.path.join(INPUT_DIR, f"{system_name}_rep{rep}.dcd")

            if not os.path.exists(pdb_path):
                logger.warning(f"PDB no encontrado: {pdb_path}  — omitiendo")
                continue
            if not os.path.exists(dcd_path):
                logger.warning(f"DCD no encontrado: {dcd_path}  — omitiendo")
                continue

            tasks.append({
                "system_name" : system_name,
                "replica"     : rep,
                "pdb_path"    : pdb_path,
                "dcd_path"    : dcd_path,
                "helices"     : HELICES,
                "chains"      : CHAINS,
                "atom"        : BACKBONE_ATOM,
                "window"      : LOCAL_WINDOW,
                "use_triangle": USE_TRIANGLE_METHOD,
                "use_midvec"  : USE_MIDVECTOR_METHOD,
                "n_frames"    : N_FRAMES,
            })

    logger.info(f"Tareas generadas: {len(tasks)} "
                f"({len(SYSTEMS)} sistemas × {len(REPLICAS)} réplicas)")

    # --------------------------------------------------------------------------
    # 8.3  Ejecutar análisis (paralelo o desde CSV)
    # --------------------------------------------------------------------------
    all_results = []

    for task in tasks:
        sys_name = task["system_name"]
        rep      = task["replica"]

        # Verificar si ya existen CSV y RECALCULATE_ANGLES == False
        sample_csv = os.path.join(
            csv_dir,
            f"{sys_name}_rep{rep}_local_{list(HELICES.keys())[0]}_{CHAINS[0]}.csv"
        )
        if not RECALCULATE_ANGLES and os.path.exists(sample_csv):
            logger.info(f"  Cargando CSV existentes: {sys_name} rep{rep}")
            res = load_angles_from_csv(sys_name, rep, csv_dir, HELICES, CHAINS)
            all_results.append(res)
        else:
            all_results.append(None)   # placeholder para el proceso paralelo

    # Filtrar las que necesitan calcularse
    tasks_to_run = [t for t, r in zip(tasks, all_results) if r is None]
    indices_to_run = [i for i, r in enumerate(all_results) if r is None]

    if tasks_to_run:
        logger.info(f"Calculando ángulos para {len(tasks_to_run)} réplicas "
                    f"usando {N_WORKERS} workers...")

        # multiprocessing.Pool compatible con Windows (requiere if __name__...)
        with multiprocessing.Pool(processes=N_WORKERS) as pool:
            computed = pool.map(analyze_replica, tasks_to_run)

        for idx, res in zip(indices_to_run, computed):
            all_results[idx] = res
            sys_name = res["system_name"]
            rep      = res["replica"]

            if res["error"]:
                logger.error(f"  ERROR en {sys_name} rep{rep}:\n{res['error']}")
            else:
                logger.info(f"  Completado: {sys_name} rep{rep}")
                if SAVE_CSV:
                    save_angles_to_csv(res, csv_dir)
                    logger.info(f"    CSV guardados: {sys_name} rep{rep}")
    else:
        logger.info("Todos los datos cargados desde CSV existentes.")

    t_calc = time.time()
    logger.info(f"Tiempo de cálculo: {(t_calc - t_start):.1f} s")

    # --------------------------------------------------------------------------
    # 8.4  Calcular rangos globales para escalas unificadas
    # --------------------------------------------------------------------------
    logger.info("Calculando rangos globales de ángulos (escalas unificadas)...")
    valid_results = [r for r in all_results if r is not None and r["error"] is None]
    global_ranges = compute_global_ranges(
        valid_results, HELICES, CHAINS, USE_TRIANGLE_METHOD, USE_MIDVECTOR_METHOD
    )

    for h, rng in global_ranges.items():
        logger.info(f"  {h}: local_max={rng['local_max']:.2f}°  "
                    f"tri_max={rng['tri_max']:.2f}°  "
                    f"mid_max={rng['mid_max']:.2f}°")

    # --------------------------------------------------------------------------
    # 8.5  Organizar resultados por sistema
    # --------------------------------------------------------------------------
    # system_map[system_name] = lista de resultados (uno por réplica)
    system_map = {s: [] for s in SYSTEMS}
    for res in valid_results:
        if res["system_name"] in system_map:
            system_map[res["system_name"]].append(res)

    # --------------------------------------------------------------------------
    # 8.6  Generar figuras
    # --------------------------------------------------------------------------
    logger.info("Generando figuras...")

    # Abrir PDF si corresponde
    pdf_path = os.path.join(plot_dir, "helix_bending_summary.pdf")
    pdf_ctx  = PdfPages(pdf_path) if SAVE_PDF else None

    # --- 8.6a  Heatmaps (por sistema × réplica × hélice) ---
    logger.info("  Generando heatmaps de ángulo local...")
    for sys_name, rep_results in system_map.items():
        for res in rep_results:
            for helix_name in HELICES:
                try:
                    plot_heatmap(
                        result       = res,
                        helix_name   = helix_name,
                        global_ranges= global_ranges,
                        chains       = CHAINS,
                        dt_ns        = DT_NS,
                        plot_dir     = plot_dir,
                        pdf_pages    = pdf_ctx,
                        save_png     = SAVE_PNG,
                    )
                except Exception as e:
                    logger.error(f"    Error heatmap {sys_name} rep{res['replica']} "
                                 f"{helix_name}: {e}")

    logger.info("  Heatmaps completados.")

    # --- 8.6b  Series temporales del ángulo global (por sistema × hélice) ---
    logger.info("  Generando gráficos de ángulo global (series temporales)...")
    for sys_name, rep_results in system_map.items():
        for helix_name in HELICES:
            for method in (["triangle"] if USE_TRIANGLE_METHOD else []) + \
                          (["midvec"]   if USE_MIDVECTOR_METHOD  else []):
                try:
                    plot_global_timeseries(
                        system_name    = sys_name,
                        helix_name     = helix_name,
                        system_results = rep_results,
                        global_ranges  = global_ranges,
                        chains         = CHAINS,
                        replicas       = REPLICAS,
                        dt_ns          = DT_NS,
                        method         = method,
                        plot_dir       = plot_dir,
                        pdf_pages      = pdf_ctx,
                        save_png       = SAVE_PNG,
                    )
                except Exception as e:
                    logger.error(f"    Error timeseries {sys_name} {helix_name} "
                                 f"{method}: {e}")

    if pdf_ctx is not None:
        pdf_ctx.close()
        logger.info(f"  PDF resumen guardado: {pdf_path}")

    # --------------------------------------------------------------------------
    # 8.7  Cierre del script
    # --------------------------------------------------------------------------
    t_end = time.time()
    t_total = t_end - t_start
    logger.info("=" * 70)
    logger.info(f"  Script finalizado exitosamente.")
    logger.info(f"  Tiempo total de ejecución: {t_total:.1f} s "
                f"({t_total/60:.1f} min)")
    logger.info(f"  Outputs en: {OUTPUT_DIR}")
    logger.info("=" * 70)


# ==============================================================================
# SECCIÓN 9 — PUNTO DE ENTRADA
#   IMPORTANTE: el bloque if __name__ == "__main__" es OBLIGATORIO en Windows
#   para que multiprocessing.Pool funcione correctamente desde Spyder.
# ==============================================================================
if __name__ == "__main__":
    main()
