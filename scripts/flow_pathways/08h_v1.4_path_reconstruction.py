# ==============================================================
# 08h_v1.4_path_reconstruction.py
# Author: DenyCB
#
# CHANGES v1.3 -> v1.4:
#   1. CLUSTERING FIX (equal-height dendrogram):
#      _determine_n_clusters_gap now detects when all merge heights
#      are equal (np.all(gaps == 0.0)) and returns min(max_k, n_paths)
#      instead of falling back to 2. This fixes the case where N paths
#      all merge at the same dendrogram level and MAX_CLUSTERS < N
#      caused the function to return 1 cluster via silhouette fallback.
#
#   2. NEW FOLDER STRUCTURE FOR MATRICES AND GEPHI:
#      MATRICES_DIR (single flat folder) replaced by:
#        OUTPUT_DIR/matrices/Method1/
#        OUTPUT_DIR/matrices/Method2/
#        OUTPUT_DIR/gephi/Method1/
#        OUTPUT_DIR/gephi/Method2/
#      method_matrices_dirs and method_gephi_dirs dicts map
#      method_label -> directory, passed through task_args.
#
#   3. TARGETS CONSOLIDATED IN MATRICES AND GEPHI:
#      process_one no longer saves matrices/gephi directly.
#      Instead it returns E_acc in the result dict.
#      _run_one_system accumulates E_acc per method_label by summing
#      across all directions (fwd/rev) and targets (filter/gate).
#      One consolidated matrix and one gephi file are saved per
#      system per method after all combinations are processed.
#
# CHANGES v1.2 -> v1.3 (preserved):
#   1. THREE PDBs (one per condition)
#   2. VECTORIZED compute_node_importance
#   3. PRECOMPUTED Cb SIM CACHE
#
# CHANGES v1.1 -> v1.2 (preserved):
#   1. Cb-AWARE NW DISTANCE (sigma=8.0 A)
#   2. EDGE ACCUMULATION with numpy array
#   3. DUAL CLUSTERING MODES (auto_cluster + fixed_cluster)
#   4. SILHOUETTE AS PRIMARY FOR AUTO MODE
#
# ALL V1.1 CHANGES PRESERVED:
#   - deduplicate_paths() consolidation with n_duplicates
#   - FASTA chain format only
#   - Dendrogram height = 60 inches
#   - ProcessPoolExecutor parallelization
#
# PURPOSE:
#   Reconstruct allosteric communication paths from Markov
#   transition (T) and fundamental (N) matrices produced by
#   08g_v2.1. Implements two path-finding methods:
#
#   Method 1 - T-shortest paths (method1_TshortestPaths):
#     Edge weight = -log(T_mean[i,j] + eps)
#     Finds paths maximising cumulative transition probability.
#     Same graph for all seeds; seed-independent network structure.
#
#   Method 2 - N-importance paths (method2_NimportancePaths):
#     node_importance[i] = sum_k( N_mean[k,i] * w_k ) over seeds
#     Edge weight = -log(node_importance[i] * T_mean[i,j] + eps)
#     One graph per system (pre-computed node_importance).
#
# INPUTS:
#   - consolidated/ CSVs from 08g_v2.1
#   - edge_importance/ CSVs from 08g_v2.1
#   - PDB files for Cb distance computation (one per condition)
#
# OUTPUTS:
#   OUTPUT_DIR/matrices/Method1/
#     matrix_{sys}_method1.csv  (consolidated across targets+directions)
#   OUTPUT_DIR/matrices/Method2/
#     matrix_{sys}_method2.csv
#   OUTPUT_DIR/gephi/Method1/
#     {sys}_method1_edges_gephi.csv
#   OUTPUT_DIR/gephi/Method2/
#     {sys}_method2_edges_gephi.csv
#   OUTPUT_DIR/method1_TshortestPaths/{sys}/auto_cluster/
#   OUTPUT_DIR/method1_TshortestPaths/{sys}/fixed_cluster/
#   OUTPUT_DIR/method2_NimportancePaths/{sys}/auto_cluster/
#   OUTPUT_DIR/method2_NimportancePaths/{sys}/fixed_cluster/
#     clusters_{sys}_{dir}_{target}.csv
#     alignment_{sys}_{dir}_{target}.fasta
#     dendrogram_{sys}_{dir}_{target}.png
#   OUTPUT_DIR/logs/
# ==============================================================

import os
import re
import glob
import itertools
import logging
import warnings
from collections import Counter, defaultdict
from concurrent.futures import ProcessPoolExecutor, as_completed
from datetime import datetime

import numpy as np
import pandas as pd
import networkx as nx
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt

from scipy.cluster.hierarchy import linkage, fcluster, dendrogram
from scipy.spatial.distance import squareform
from sklearn.metrics import silhouette_score

warnings.filterwarnings('ignore')

# ==============================================================
# SECTION 1: ANALYSIS SWITCHES
# ==============================================================

IS_NMA            = False    # True = NMA; False = MD (consolidated mean)

RUN_METHOD1       = True    # T-shortest paths
RUN_METHOD2       = True    # N-importance paths
RUN_FWD           = True    # forward direction (seed -> target)
RUN_REV           = True    # reverse direction (target -> seed)
RUN_WT_COMPARISON = False   # optional WT vs mutant comparison block

# Clustering modes — both can be True; linkage computed once
RUN_CLUSTER_AUTO  = True    # silhouette-primary auto k selection
RUN_CLUSTER_FIXED = True    # fixed k = FIXED_K_CLUSTERS

# Parallel workers at system level (1 = sequential)
N_CORES = 1

# ==============================================================
# SECTION 2: INPUT / OUTPUT CONFIGURATION
# ==============================================================

# --- NMA (uncomment to use) ---#INPUT_08G_DIR = (
#    "C:/DinamicasMoleculares/TRPV1_pipeline/output/"
#    "08_Max_Flow/08g_v2.1_markov_NMA"
#)
#MAP_FILE = (
#    "C:/DinamicasMoleculares/TRPV1_pipeline/config_files/"
#    "trimmed_residue_map.csv"
#)
#OUTPUT_DIR = (
#    "C:/DinamicasMoleculares/TRPV1_pipeline/output/"
#    "08_Max_Flow/08h_v1.4_paths_NMA"
#)

# --- MD (comment out NMA block above and uncomment below) ---
INPUT_08G_DIR = (
     "C:/DinamicasMoleculares/TRPV1_pipeline/output/"
     "08_Max_Flow/08g_v2.1_markov_MD"
)
MAP_FILE = (
     "C:/DinamicasMoleculares/TRPV1_pipeline/output/residue_map_MD.csv"
)
OUTPUT_DIR = (
     "C:/DinamicasMoleculares/TRPV1_pipeline/output/"
     "08_Max_Flow/08h_v1.4_paths_MD"
)

# PDB files for Cb distance computation — one per condition.
# Key = PDB prefix used in system names (e.g. '7LP9' from '7LP9_WT').
# Value = filename without .pdb extension.
# Mutant systems (e.g. 7LP9_W426A) use the same PDB as their WT.
PDB_DIR   = "C:/DinamicasMoleculares/TRPV1_pipeline/input/pdb"
PDB_FILES = {
    '7LP9': '7LP9',
    '7LPB': '7LPB',
    '7LPC': '7LPC',
}

CONSOLIDATED_DIR      = os.path.join(INPUT_08G_DIR, "consolidated")
EDGE_IMP_DIR          = os.path.join(INPUT_08G_DIR, "edge_importance")
METHOD1_DIR           = os.path.join(OUTPUT_DIR, "method1_TshortestPaths")
METHOD2_DIR           = os.path.join(OUTPUT_DIR, "method2_NimportancePaths")
MATRICES_METHOD1_DIR  = os.path.join(OUTPUT_DIR, "matrices", "Method1")
MATRICES_METHOD2_DIR  = os.path.join(OUTPUT_DIR, "matrices", "Method2")
GEPHI_METHOD1_DIR     = os.path.join(OUTPUT_DIR, "gephi", "Method1")
GEPHI_METHOD2_DIR     = os.path.join(OUTPUT_DIR, "gephi", "Method2")
LOGS_DIR              = os.path.join(OUTPUT_DIR, "logs")

for _d in [OUTPUT_DIR, METHOD1_DIR, METHOD2_DIR,
           MATRICES_METHOD1_DIR, MATRICES_METHOD2_DIR,
           GEPHI_METHOD1_DIR, GEPHI_METHOD2_DIR, LOGS_DIR]:
    os.makedirs(_d, exist_ok=True)

# ==============================================================
# SECTION 3: ANALYSIS PARAMETERS
# ==============================================================

TARGET_FILTER = 643
TARGET_GATE   = 679

K_PATHS    = 5       # k shortest paths per seed per target node set
EPS        = 1e-10   # avoid log(0)

# Needleman-Wunsch gap penalty
GAP_PENALTY = 0.5

# Cb-aware NW substitution parameters
USE_CB_DISTANCES  = True   # False = fall back to binary 0/1
CB_SIGMA_ANGSTROM = 8.0    # sigma in Angstroms: exp(-d/sigma)
                            # d=0: sim=1.0; d=4A: ~0.61; d=8A: ~0.37

# Clustering
MAX_CLUSTERS      = 10
MIN_CLUSTERS      = 2
FIXED_K_CLUSTERS  = 10    # k used in fixed clustering mode
SIL_MIN_THRESHOLD = 0.10  # below this, silhouette is inconclusive
                           # -> fall back to max-gap in auto mode

# Chain constants
TARGET_CHAINS       = ['A', 'B', 'C', 'D']
TARGET_FILTER_NODES = [f"{c}{TARGET_FILTER}" for c in TARGET_CHAINS]
TARGET_GATE_NODES   = [f"{c}{TARGET_GATE}"   for c in TARGET_CHAINS]
CHAIN_MAP           = {'A': 0, 'B': 1, 'C': 2, 'D': 3}
REVERSE_MAP         = {0: 'A', 1: 'B', 2: 'C', 3: 'D'}

# ==============================================================
# SECTION 4: SEEDS
# ==============================================================

SEEDS = {
    "C289": 1.0, "B289": 1.0,
    "C296": 1.0, "B296": 1.0,
    "A297": 1.0, "C297": 1.0, "B297": 1.0, "D297": 1.0,
    "C298": 1.0, "A301": 1.0,
    "B337": 1.0, "D337": 1.0, "A337": 1.0,
    "A338": 1.0, "D338": 1.0, "B338": 1.0,
    "B341": 1.0, "C341": 1.0, "D341": 1.0,
    "D344": 1.0, "B344": 1.0,
    "C345": 1.0, "B345": 1.0, "A345": 1.0,
    "D346": 1.0, "B346": 1.0, "C346": 1.0, "A346": 1.0,
    "B347": 1.0, "A347": 1.0, "C347": 1.0, "D347": 1.0,
    "C349": 1.0, "A349": 1.0, "B349": 1.0, "D349": 1.0,
    "C350": 1.0, "A350": 1.0, "D350": 1.0, "B350": 1.0,
    "A351": 1.0, "A352": 1.0, "C352": 1.0, "D354": 1.0,
    "D378": 1.0, "D382": 1.0,
    "A387": 1.0, "C387": 1.0, "A388": 1.0,
    "D395": 1.0, "C395": 1.0, "C396": 1.0,
    "D399": 1.0, "C399": 1.0, "A399": 1.0, "B399": 1.0,
    "A400": 1.0, "B400": 1.0, "D400": 1.0, "C400": 1.0,
    "A408": 1.0, "C408": 1.0, "B408": 1.0, "D408": 1.0,
    "C409": 1.0, "B409": 1.0, "D409": 1.0, "A409": 1.0,
    "B410": 1.0, "C410": 1.0, "A410": 1.0, "D410": 1.0,
    "C411": 1.0, "A411": 1.0, "D411": 1.0, "B411": 1.0,
    "D412": 1.0, "B412": 1.0, "A412": 1.0, "C412": 1.0,
    "C413": 1.0, "A413": 1.0, "B413": 1.0, "D413": 1.0,
    "A414": 1.0, "D414": 1.0, "C414": 1.0, "B414": 1.0,
    "D418": 1.0, "A418": 1.0, "B418": 1.0, "C418": 1.0,
    "D419": 1.0, "A419": 1.0, "B419": 1.0, "C419": 1.0,
    "B421": 1.0, "C426": 1.0, "B430": 1.0, "A430": 1.0,
    "D745": 1.0, "A750": 1.0, "D746": 1.0, "A749": 1.0,
    "D747": 1.0, "B746": 1.0, "A746": 1.0, "C745": 1.0,
    "B731": 1.0, "D726": 1.0, "A747": 1.0, "C726": 1.0,
    "B729": 1.0, "C746": 1.0, "B730": 1.0, "D744": 1.0,
    "B747": 1.0, "A730": 1.0, "A731": 1.0, "D728": 1.0,
    "A551": 1.0, "B551": 1.0, "C551": 1.0, "D551": 1.0,
    "A550": 1.0, "B550": 1.0, "C550": 1.0, "D550": 1.0,
    "A570": 1.0, "B570": 1.0, "C570": 1.0, "D570": 1.0,
    "A571": 1.0, "B571": 1.0, "C571": 1.0, "D571": 1.0,
    "A510": 1.0, "B510": 1.0, "C510": 1.0, "D510": 1.0,
    "A511": 1.0, "B511": 1.0, "C511": 1.0, "D511": 1.0,
    "A512": 1.0, "B512": 1.0, "C512": 1.0, "D512": 1.0,
    "A516": 1.0, "B516": 1.0, "C516": 1.0, "D516": 1.0,
    "A462": 1.0, "C461": 1.0, "B461": 1.0, "D461": 1.0,
    "A457": 1.0, "B457": 1.0, "D462": 1.0, "A455": 1.0,
    "A461": 1.0, "B455": 1.0, "B460": 1.0, "C460": 1.0,
    "D460": 1.0, "A456": 1.0, "A458": 1.0, "C457": 1.0,
    "B536": 1.0, "A536": 1.0,
}

# ==============================================================
# SECTION 5: LOGGING
# ==============================================================

_timestamp  = datetime.now().strftime("%Y%m%d_%H%M%S")
_master_log = os.path.join(LOGS_DIR, f"08h_v1.4_master_log_{_timestamp}.txt")

logging.basicConfig(
    level=logging.INFO,
    format="[%(asctime)s] 08h_v1.4 | %(message)s",
    datefmt="%Y-%m-%d %H:%M:%S",
    handlers=[
        logging.FileHandler(_master_log),
        logging.StreamHandler(),
    ],
)
log = logging.getLogger()


def log_msg(msg, sys_log=None):
    log.info(msg)
    if sys_log is not None:
        try:
            with open(sys_log, 'a') as _fh:
                _fh.write(
                    f"[{datetime.now().strftime('%Y-%m-%d %H:%M:%S')}]"
                    f" 08h_v1.4 | {msg}\n"
                )
        except Exception:
            pass

# ==============================================================
# SECTION 6: CB COORDINATE LOADING
# ==============================================================

def load_cb_coords_from_pdb(pdb_path, chains=None):
    """Parse PDB ATOM records and extract Cb coordinates (Ca for Gly).
    Pure Python parser — no external molecular library required.
    Handles alternate location codes (keeps first occurrence).

    Parameters
    ----------
    pdb_path : str
        Path to PDB file.
    chains : list of str or None
        Chain IDs to include. None = ['A','B','C','D'].

    Returns
    -------
    dict {chainresid_str: np.array([x,y,z], float32)}
        e.g. {'A300': array([10.1, 5.2, 3.3]), ...}
    """
    if chains is None:
        chains = ['A', 'B', 'C', 'D']
    coords = {}
    try:
        with open(pdb_path, 'r') as fh:
            for line in fh:
                if not line.startswith('ATOM'):
                    continue
                if len(line) < 54:
                    continue
                atom_name = line[12:16].strip()
                alt_loc   = line[16]
                if alt_loc not in (' ', 'A'):
                    continue
                res_name = line[17:20].strip()
                chain_id = line[21]
                if chain_id not in chains:
                    continue
                resnum_str = line[22:26].strip()
                if not resnum_str.lstrip('-').isdigit():
                    continue
                resnum = int(resnum_str)
                if resnum <= 0:
                    continue
                # Use CB for non-Gly, CA for Gly
                want_atom = ('CB' if res_name != 'GLY' else 'CA')
                if atom_name != want_atom:
                    continue
                label = f"{chain_id}{resnum}"
                if label in coords:
                    continue   # keep first alt loc
                try:
                    x = float(line[30:38])
                    y = float(line[38:46])
                    z = float(line[46:54])
                    coords[label] = np.array([x, y, z], dtype=np.float32)
                except ValueError:
                    pass
    except Exception as e:
        log_msg(f"WARNING: could not load PDB {pdb_path}: {e}")
    return coords


def precompute_cb_sim_cache(cb_coords, sigma):
    """Precompute all pairwise Cb similarity values from coordinate dict.
    Returns dict {(a, b): sim_value} for all label pairs a != b.
    Both (a,b) and (b,a) stored for O(1) symmetric lookup.
    Computed once per PDB at startup; avoids repeated np.linalg.norm
    calls inside the NW distance matrix computation."""
    labels     = list(cb_coords.keys())
    n          = len(labels)
    coords_arr = np.array([cb_coords[l] for l in labels], dtype=np.float32)
    cache      = {}
    for i in range(n):
        for j in range(i + 1, n):
            d   = float(np.linalg.norm(coords_arr[i] - coords_arr[j]))
            sim = float(np.exp(-d / sigma))
            cache[(labels[i], labels[j])] = sim
            cache[(labels[j], labels[i])] = sim
    return cache

# ==============================================================
# SECTION 7: HELPER FUNCTIONS — I/O
# ==============================================================

def load_matrix_csv(path):
    """Load CSV written by R write.csv(row.names=TRUE)."""
    df = pd.read_csv(path, index_col=0)
    df.index   = df.index.astype(str)
    df.columns = df.columns.astype(str)
    return df


def discover_systems(consolidated_dir):
    """Discover system_base names from N_fwd_filter_*_mean.csv."""
    pattern = os.path.join(consolidated_dir, "N_fwd_filter_*_mean.csv")
    systems = set()
    for f in glob.glob(pattern):
        m = re.match(r"N_fwd_filter_(.+)_mean\.csv$", os.path.basename(f))
        if m:
            systems.add(m.group(1))
    return sorted(systems)


def load_system_matrices(consolidated_dir, sys_base, direction, target):
    """Load T_mean and N_mean DataFrames. Returns (T_df, N_df) or (None, None)."""
    t_path = os.path.join(
        consolidated_dir,
        f"T_{direction}_{target}_{sys_base}_mean.csv"
    )
    n_path = os.path.join(
        consolidated_dir,
        f"N_{direction}_{target}_{sys_base}_mean.csv"
    )
    if not os.path.exists(t_path) or not os.path.exists(n_path):
        return None, None
    return load_matrix_csv(t_path), load_matrix_csv(n_path)


def load_edge_importance(edge_imp_dir, sys_base, direction, target):
    """Load edge importance CSV from 08g edge_importance/ folder."""
    path = os.path.join(
        edge_imp_dir,
        f"matrix_{sys_base}_{direction}_{target}.csv"
    )
    if not os.path.exists(path):
        return None
    return load_matrix_csv(path)

# ==============================================================
# SECTION 8: HELPER FUNCTIONS — GRAPH CONSTRUCTION
# ==============================================================

def compute_node_importance(N_df, seeds_dict):
    """Seed-weighted node importance: sum_k(N[k,j]*w_k) over all seeds k.
    Vectorized: uses numpy matrix multiply instead of Python double loop."""
    n_labels  = set(N_df.index.tolist())
    all_cols  = N_df.columns.tolist()
    valid     = [(s, w) for s, w in seeds_dict.items() if s in n_labels]
    node_imp  = pd.Series(0.0, index=all_cols)
    if not valid:
        return node_imp
    total_w   = sum(w for _, w in valid)
    seed_lbls = [s for s, w in valid]
    weights   = np.array([w / total_w for s, w in valid], dtype=np.float64)
    N_sub     = N_df.loc[seed_lbls].values.astype(np.float64)
    node_imp_vals = (weights[:, None] * N_sub).sum(axis=0)
    return pd.Series(node_imp_vals, index=all_cols)


def build_graph(T_df, node_importance=None, eps=EPS):
    """Build weighted directed DiGraph from T_mean.
    Method 1 (node_importance=None): weight = -log(T[i,j] + eps)
    Method 2: weight = -log(node_importance[i] * T[i,j] + eps)"""
    G      = nx.DiGraph()
    labels = T_df.index.tolist()
    G.add_nodes_from(labels)
    T_arr  = T_df.values
    idx    = {lbl: i for i, lbl in enumerate(labels)}

    for i_lbl in labels:
        i    = idx[i_lbl]
        ni_i = (1.0 if node_importance is None
                else float(node_importance.get(i_lbl, 0.0)))
        for j_lbl in labels:
            j    = idx[j_lbl]
            t_ij = float(T_arr[i, j])
            if t_ij <= 0.0:
                continue
            w = (-np.log(t_ij + eps) if node_importance is None
                 else -np.log(ni_i * t_ij + eps))
            G.add_edge(i_lbl, j_lbl, weight=float(w))
    return G

# ==============================================================
# SECTION 9: HELPER FUNCTIONS — PATH FINDING
# ==============================================================

def find_k_paths(G, source, target_nodes, k, syslog=None):
    """Find top-k shortest simple paths (Yen's) from source to any
    node in target_nodes. Returns list of (total_weight, path)."""
    if source not in G:
        return []
    candidates = []
    for tgt in target_nodes:
        if tgt not in G or tgt == source:
            continue
        try:
            gen = nx.shortest_simple_paths(G, source, tgt, weight='weight')
            for path in itertools.islice(gen, k):
                w = sum(G[u][v]['weight'] for u, v in zip(path[:-1], path[1:]))
                candidates.append((w, path))
        except (nx.NetworkXNoPath, nx.NodeNotFound):
            pass
        except Exception as e:
            if syslog:
                log_msg(f"    WARNING {source}->{tgt}: {e}", syslog)
    candidates.sort(key=lambda x: x[0])
    return candidates[:k]

# ==============================================================
# SECTION 10: HELPER FUNCTIONS — CHAIN RELATIVISATION
# ==============================================================

def relativize_path(path):
    """Shift chains so first node is always chain A.
    Returns (relativized_path, chain_origin_letter)."""
    if not path:
        return [], 'A'
    first_chain = path[0][0]
    shift       = CHAIN_MAP.get(first_chain, 0)
    rel = [
        REVERSE_MAP[(CHAIN_MAP[node[0]] - shift) % 4] + node[1:]
        for node in path
    ]
    return rel, first_chain

# ==============================================================
# SECTION 11: HELPER FUNCTIONS — PATH SCORING
# ==============================================================

def path_log_prob(path, T_df):
    """Sum of log(T[i,j]) along path edges."""
    total = 0.0
    for u, v in zip(path[:-1], path[1:]):
        if u in T_df.index and v in T_df.columns:
            t = float(T_df.loc[u, v])
            total += np.log(max(t, EPS))
    return float(total)


def path_importance_score(path, E_df):
    """Sum of edge_importance[i,j] along path. Returns 0 if E_df None."""
    if E_df is None:
        return 0.0
    total = 0.0
    for u, v in zip(path[:-1], path[1:]):
        if u in E_df.index and v in E_df.columns:
            total += float(E_df.loc[u, v])
    return float(total)

# ==============================================================
# SECTION 12: HELPER FUNCTIONS — DEDUPLICATION
# ==============================================================

def deduplicate_paths(raw_paths):
    """Consolidate identical relativized paths into one entry.
    Paths that were identical before relativization OR that become
    identical after chain-shift are collapsed. Each unique path retains:
      n_duplicates    : count of original paths with this sequence
      imp_score       : sum of imp_scores across all duplicates
      log_prob        : mean log_prob across all duplicates
      chain_origins   : set of original chain letters
      n_chain_origins : cardinality of chain_origins
      seeds_in_path   : set of seed labels that generated this path
    This is CONSOLIDATION (sum), not elimination: n_duplicates
    records how many original paths collapsed here."""
    seen = defaultdict(list)
    for p in raw_paths:
        key = tuple(p['path_rel'])
        seen[key].append(p)

    deduped = []
    for key, group in seen.items():
        rep                    = {k: v for k, v in group[0].items()}
        rep['n_duplicates']    = len(group)
        rep['imp_score']       = sum(g['imp_score'] for g in group)
        rep['log_prob']        = float(np.mean([g['log_prob'] for g in group]))
        rep['chain_origins']   = set(g['chain_origin'] for g in group)
        rep['n_chain_origins'] = len(rep['chain_origins'])
        rep['seeds_in_path']   = set(g['seed'] for g in group)
        deduped.append(rep)
    return deduped

# ==============================================================
# SECTION 13: HELPER FUNCTIONS — NW DISTANCE + CLUSTERING
# ==============================================================

def _cb_sim(a, b, cb_coords, sigma, cb_cache=None):
    """Substitution similarity for NW using Cb distances.
    Returns 1.0 for identical nodes, exp(-d/sigma) if both in
    cb_coords, 0.0 otherwise. Uses cb_cache for O(1) lookup
    when provided (avoids np.linalg.norm per pair)."""
    if a == b:
        return 1.0
    if cb_cache is not None:
        return cb_cache.get((a, b), 0.0)
    if cb_coords and a in cb_coords and b in cb_coords:
        d = float(np.linalg.norm(cb_coords[a] - cb_coords[b]))
        return float(np.exp(-d / sigma))
    return 0.0


def nw_distance(seq1, seq2, gap_penalty=GAP_PENALTY,
                cb_coords=None, cb_sigma=CB_SIGMA_ANGSTROM,
                cb_cache=None):
    """Global NW alignment distance in [0, 1].
    Substitution score: _cb_sim(a,b) if cb_coords or cb_cache provided,
    else binary. distance = 1 - max(0, NW_score / max(len1, len2))."""
    n, m = len(seq1), len(seq2)
    if n == 0 and m == 0:
        return 0.0
    if n == 0 or m == 0:
        return 1.0

    dp = np.zeros((n + 1, m + 1))
    for i in range(1, n + 1):
        dp[i][0] = -i * gap_penalty
    for j in range(1, m + 1):
        dp[0][j] = -j * gap_penalty

    for i in range(1, n + 1):
        for j in range(1, m + 1):
            if cb_coords is not None or cb_cache is not None:
                match = _cb_sim(seq1[i - 1], seq2[j - 1],
                                cb_coords, cb_sigma, cb_cache)
            else:
                match = 1.0 if seq1[i - 1] == seq2[j - 1] else 0.0
            dp[i][j] = max(
                dp[i - 1][j - 1] + match,
                dp[i - 1][j]     - gap_penalty,
                dp[i][j - 1]     - gap_penalty,
            )

    max_score = float(max(n, m))
    return float(np.clip(1.0 - max(0.0, dp[n][m] / max_score), 0.0, 1.0))


def build_distance_matrix(paths_list, gap_penalty=GAP_PENALTY,
                          cb_coords=None, cb_sigma=CB_SIGMA_ANGSTROM,
                          cb_cache=None):
    """Symmetric NW distance matrix for a list of paths."""
    n = len(paths_list)
    D = np.zeros((n, n))
    for i in range(n):
        for j in range(i + 1, n):
            d       = nw_distance(paths_list[i], paths_list[j],
                                  gap_penalty, cb_coords, cb_sigma,
                                  cb_cache)
            D[i, j] = d
            D[j, i] = d
    return D


def _determine_n_clusters_gap(Z, n_paths, max_k, min_k):
    """Optimal k from max gap in linkage merge heights.
    When all merge heights are equal (all gaps == 0), returns
    min(max_k, n_paths) to give as many clusters as possible."""
    if n_paths <= min_k:
        return min(n_paths, min_k)
    heights = Z[:, 2]
    if len(heights) < 2:
        return min_k
    gaps     = np.diff(heights)
    # All heights equal: each path is equidistant from all others.
    # Return as many clusters as allowed up to n_paths.
    if np.all(gaps == 0.0):
        return min(max_k, n_paths)
    rev_gaps = gaps[::-1]
    gap_idx  = int(np.argmax(rev_gaps))
    return int(max(min_k, min(max_k, gap_idx + 2)))


def _determine_n_clusters_auto(Z, D, n_paths, max_k, min_k,
                                sil_threshold=SIL_MIN_THRESHOLD):
    """Determine optimal k with silhouette as primary criterion.
    Falls back to max-gap if max silhouette < sil_threshold.
    Returns (n_clust, sil_scores_dict, method_used_str)."""
    if n_paths <= min_k:
        return min(n_paths, min_k), {}, 'trivial'

    sil_scores = {}
    for k in range(min_k, min(max_k + 1, n_paths)):
        labels_k = fcluster(Z, k, criterion='maxclust')
        if len(set(labels_k)) < 2:
            continue
        sil = float(silhouette_score(D, labels_k, metric='precomputed'))
        sil_scores[k] = round(sil, 4)

    if sil_scores and max(sil_scores.values()) >= sil_threshold:
        best_k = max(sil_scores, key=sil_scores.get)
        return best_k, sil_scores, 'silhouette'

    # Fallback: max-gap
    k_gap = _determine_n_clusters_gap(Z, n_paths, max_k, min_k)
    return k_gap, sil_scores, 'max_gap_fallback'


def _compute_linkage(D):
    """Compute complete linkage from NW distance matrix.
    Returns (Z, valid) where valid=False if matrix is degenerate."""
    condensed = squareform(D, checks=False)
    if condensed.max() == 0.0:
        return None, False
    Z = linkage(condensed, method='complete')
    return Z, True


def find_medoid(global_indices, D):
    """Return index of medoid within a cluster."""
    if len(global_indices) == 1:
        return global_indices[0]
    sub    = D[np.ix_(global_indices, global_indices)]
    mean_d = sub.mean(axis=1)
    return global_indices[int(np.argmin(mean_d))]


def get_cut_height(Z, n_paths, n_clusters):
    """Return dendrogram cut height for given n_clusters."""
    if n_clusters >= n_paths or len(Z) == 0:
        return Z[0, 2] / 2.0 if len(Z) > 0 else 0.0
    n_merges = min(n_paths - n_clusters, len(Z))
    if n_merges <= 0:
        return Z[0, 2] / 2.0
    h_below = Z[n_merges - 1, 2]
    h_above = Z[n_merges, 2] if n_merges < len(Z) else h_below + 1e-6
    return float((h_below + h_above) / 2.0)

# ==============================================================
# SECTION 14: HELPER FUNCTIONS — CONSENSUS SEQUENCE
# ==============================================================

def _align_pair_traceback(seq1, seq2, gap_penalty=GAP_PENALTY,
                          cb_coords=None, cb_sigma=CB_SIGMA_ANGSTROM,
                          cb_cache=None):
    """NW with traceback using same substitution scores as nw_distance.
    Returns (aligned_seq1, aligned_seq2) with None for gap positions."""
    n, m = len(seq1), len(seq2)
    dp   = np.zeros((n + 1, m + 1))
    for i in range(1, n + 1):
        dp[i][0] = -i * gap_penalty
    for j in range(1, m + 1):
        dp[0][j] = -j * gap_penalty

    for i in range(1, n + 1):
        for j in range(1, m + 1):
            if cb_coords is not None or cb_cache is not None:
                match = _cb_sim(seq1[i-1], seq2[j-1],
                                cb_coords, cb_sigma, cb_cache)
            else:
                match = 1.0 if seq1[i-1] == seq2[j-1] else 0.0
            dp[i][j] = max(
                dp[i-1][j-1] + match,
                dp[i-1][j]   - gap_penalty,
                dp[i][j-1]   - gap_penalty,
            )

    al1, al2 = [], []
    i, j = n, m
    while i > 0 or j > 0:
        if i > 0 and j > 0:
            if cb_coords is not None or cb_cache is not None:
                match = _cb_sim(seq1[i-1], seq2[j-1],
                                cb_coords, cb_sigma, cb_cache)
            else:
                match = 1.0 if seq1[i-1] == seq2[j-1] else 0.0
            if np.isclose(dp[i][j], dp[i-1][j-1] + match):
                al1.append(seq1[i - 1])
                al2.append(seq2[j - 1])
                i -= 1
                j -= 1
                continue
        if i > 0 and (j == 0 or
                      np.isclose(dp[i][j], dp[i-1][j] - gap_penalty)):
            al1.append(seq1[i - 1])
            al2.append(None)
            i -= 1
        elif j > 0:
            al1.append(None)
            al2.append(seq2[j - 1])
            j -= 1
        else:
            break
    return al1[::-1], al2[::-1]


def consensus_from_cluster(member_rel_paths, medoid_sub_idx,
                           gap_penalty=GAP_PENALTY,
                           cb_coords=None, cb_sigma=CB_SIGMA_ANGSTROM,
                           cb_cache=None):
    """Progressive alignment of all members against medoid.
    Returns consensus as list of residue labels."""
    ref = member_rel_paths[medoid_sub_idx]
    if len(member_rel_paths) == 1:
        return list(ref)
    projections = []
    for idx, path in enumerate(member_rel_paths):
        if idx == medoid_sub_idx:
            projections.append(list(ref))
        else:
            al_ref, al_path = _align_pair_traceback(
                ref, path, gap_penalty, cb_coords, cb_sigma, cb_cache
            )
            projected = [p for r, p in zip(al_ref, al_path) if r is not None]
            projections.append(projected)
    consensus = []
    for pos in range(len(ref)):
        votes = [proj[pos] for proj in projections
                 if pos < len(proj) and proj[pos] is not None]
        if votes:
            consensus.append(Counter(votes).most_common(1)[0][0])
    return consensus

# ==============================================================
# SECTION 15: HELPER FUNCTIONS — FORMATTING & OUTPUT
# ==============================================================

def fmt_chainresid(path):
    """Format path as 'A300-A350-D420'."""
    return '-'.join(path)


def fmt_resid_only(path):
    """Format path residue numbers as '300 350 420'."""
    return ' '.join(node[1:] for node in path)


def build_edge_accumulation(deduped_paths, all_labels, E_df=None, syslog=None):
    n            = len(all_labels)
    label_to_idx = {lbl: i for i, lbl in enumerate(all_labels)}
    E            = np.zeros((n, n), dtype=np.float64)

    # Pre-extraer E_df como array numpy para lookup rápido
    edf_arr = None
    edf_idx = None
    if E_df is not None:
        edf_labels = E_df.index.tolist()
        edf_idx    = {lbl: i for i, lbl in enumerate(edf_labels)}
        edf_arr    = E_df.values.astype(np.float64)

    n_edges_acc = 0
    for p in deduped_paths:
        path  = p['path_rel']
        n_dup = float(p.get('n_duplicates', 1))

        for u, v in zip(path[:-1], path[1:]):
            i = label_to_idx.get(u)
            j = label_to_idx.get(v)
            if i is None or j is None:
                continue

            edge_val = 0.0
            if edf_arr is not None:
                ei = edf_idx.get(u)
                ej = edf_idx.get(v)
                if ei is not None and ej is not None:
                    v_tmp = float(edf_arr[ei, ej])
                    if not np.isnan(v_tmp):
                        edge_val = v_tmp

            weight    = (edge_val * n_dup) if edge_val > 0.0 else n_dup
            E[i, j]  += weight
            n_edges_acc += 1

    if syslog:
        log_msg(
            f"      Edge accumulation: {n_edges_acc} edge traversals, "
            f"{int((E > 0).sum())} non-zero cells", syslog
        )

    E_sym = (E + E.T) / 2.0
    return pd.DataFrame(E_sym, index=all_labels, columns=all_labels)


def save_distogram_csv(E_sym, out_path):
    """Save matrix as 08d-compatible CSV (row.names=TRUE equivalent)."""
    E_sym.to_csv(out_path)


def save_gephi_csv(E_sym, out_path, type_label):
    """Save Gephi edge list: Source, Target, type, weight, weight_abs."""
    max_val = float(E_sym.values.max())
    if max_val <= 0.0:
        max_val = 1.0
    rows   = []
    labels = E_sym.index.tolist()
    for i, src in enumerate(labels):
        for j, tgt in enumerate(labels):
            if i >= j:
                continue
            val = float(E_sym.loc[src, tgt])
            if val <= 0.0:
                continue
            rows.append({
                'Source':     src,
                'Target':     tgt,
                'type':       type_label,
                'weight':     val / max_val,
                'weight_abs': val,
            })
    pd.DataFrame(
        rows,
        columns=['Source', 'Target', 'type', 'weight', 'weight_abs']
    ).to_csv(out_path, index=False)


def save_dendrogram_png(Z, n_unique, n_clust, out_path, title,
                        gap_penalty, mode_label):
    """Save dendrogram PNG with cut line. Height = 60 inches."""
    fig_w = max(10, n_unique * 0.20)
    fig_h = 60
    fig, ax = plt.subplots(figsize=(fig_w, fig_h))
    dendrogram(Z, ax=ax, labels=None,
               leaf_rotation=90, leaf_font_size=6,
               color_threshold=0)
    cut_h = get_cut_height(Z, n_unique, n_clust)
    ax.axhline(y=cut_h, color='red', linestyle='--', linewidth=2,
               label=f'Cut  k={n_clust}  [{mode_label}]')
    ax.legend(fontsize=9)
    ax.set_title(title, fontsize=9)
    ax.set_xlabel('Unique paths after deduplication', fontsize=8)
    ax.set_ylabel(f'NW distance  (gap={gap_penalty}  Cb-aware={USE_CB_DISTANCES})',
                  fontsize=8)
    plt.tight_layout()
    fig.savefig(out_path, dpi=150)
    plt.close(fig)

# ==============================================================
# SECTION 16: CLUSTERING OUTPUT HELPER
# ==============================================================

def _save_cluster_outputs(deduped, rel_list, D, labels, n_clust, Z,
                          out_dir, tag, mode_label,
                          gap_penalty, cb_coords, cb_sigma, syslog,
                          cb_cache=None):
    """Save cluster CSV, FASTA, and dendrogram for one clustering mode.
    Called once for 'auto_cluster' and once for 'fixed_cluster'."""
    os.makedirs(out_dir, exist_ok=True)

    cluster_rows = []
    fasta_lines  = []

    for cid in range(1, n_clust + 1):
        mem_idx = [i for i, l in enumerate(labels) if l == cid]
        members = [deduped[i] for i in mem_idx]
        if not members:
            continue

        imp_scores      = [m['imp_score']       for m in members]
        log_probs       = [m['log_prob']        for m in members]
        lengths         = [m['path_length']     for m in members]
        n_dups          = [m['n_duplicates']    for m in members]
        n_chain_origs   = [m['n_chain_origins'] for m in members]
        all_chain_origs = set.union(*[m['chain_origins'] for m in members])
        all_seeds       = set.union(*[m['seeds_in_path'] for m in members])

        n_members           = len(members)
        n_paths_total       = int(sum(n_dups))
        n_chains            = len(all_chain_origs)
        chain_freq          = round(n_chains / 4.0, 3)
        mean_chains_unique  = round(float(np.mean(n_chain_origs)), 2)
        mean_imp            = float(np.mean(imp_scores))
        total_imp           = float(np.sum(imp_scores))
        mean_logprob        = float(np.mean(log_probs))
        mean_length         = float(np.mean(lengths))

        # Medoid
        sub_arr = np.array(mem_idx)
        med_sub = find_medoid(
            list(range(len(mem_idx))),
            D[np.ix_(sub_arr, sub_arr)]
        )
        med_gi   = mem_idx[med_sub]
        rep_path = deduped[med_gi]['path_rel']

        # Consensus
        mem_rel   = [rel_list[i] for i in mem_idx]
        consensus = consensus_from_cluster(
            mem_rel, med_sub, gap_penalty, cb_coords, cb_sigma, cb_cache
        )

        rep_cr = fmt_chainresid(rep_path)
        rep_ro = fmt_resid_only(rep_path)
        con_cr = fmt_chainresid(consensus)
        con_ro = fmt_resid_only(consensus)

        cluster_rows.append({
            'cluster_mode':           mode_label,
            'cluster_id':             cid,
            'n_members':              n_members,
            'n_paths_total':          n_paths_total,
            'n_chains':               n_chains,
            'chain_freq':             chain_freq,
            'mean_chains_per_unique': mean_chains_unique,
            'n_seeds':                len(all_seeds),
            'mean_imp_score':         round(mean_imp, 6),
            'total_imp_score':        round(total_imp, 6),
            'mean_log_prob':          round(mean_logprob, 6),
            'mean_path_length':       round(mean_length, 2),
            'rep_chainresid':         rep_cr,
            'rep_resid_only':         rep_ro,
            'consensus_chainresid':   con_cr,
            'consensus_resid_only':   con_ro,
        })

        # FASTA block — chain format only
        fasta_lines += [
            f">cluster_{cid}_{mode_label}"
            f"  n_members={n_members}"
            f"  n_paths_total={n_paths_total}"
            f"  mean_imp={mean_imp:.4f}"
            f"  n_chains={n_chains}",
            f"#consensus:      {con_cr}",
            f"#representative: {rep_cr}",
        ]
        for mi, m in enumerate(members):
            fasta_lines += [
                f">cluster_{cid}_member_{mi+1}"
                f"  seed={','.join(sorted(m['seeds_in_path']))}"
                f"  chain_origin={','.join(sorted(m['chain_origins']))}"
                f"  n_dup={m['n_duplicates']}"
                f"  imp={m['imp_score']:.4f}"
                f"  logp={m['log_prob']:.4f}"
                f"  len={m['path_length']}",
                fmt_chainresid(m['path_rel']),
            ]
        fasta_lines.append("")

    # ---- Save cluster CSV ----
    csv_path = os.path.join(out_dir, f"clusters_{tag}.csv")
    pd.DataFrame(cluster_rows).to_csv(csv_path, index=False)
    log_msg(f"      [{mode_label}] clusters saved: {os.path.basename(csv_path)}",
            syslog)

    # ---- Save FASTA ----
    fasta_path = os.path.join(out_dir, f"alignment_{tag}.fasta")
    with open(fasta_path, 'w') as fh:
        fh.write('\n'.join(fasta_lines) + '\n')
    log_msg(f"      [{mode_label}] alignment saved: {os.path.basename(fasta_path)}",
            syslog)

    # ---- Save dendrogram ----
    if Z is not None and len(deduped) > 1:
        png_path = os.path.join(out_dir, f"dendrogram_{tag}.png")
        save_dendrogram_png(
            Z, len(deduped), n_clust, png_path,
            title=f"Path clustering: {tag}  [{mode_label}]  k={n_clust}",
            gap_penalty=gap_penalty,
            mode_label=mode_label,
        )
        log_msg(f"      [{mode_label}] dendrogram saved: {os.path.basename(png_path)}",
                syslog)

    return cluster_rows

# ==============================================================
# SECTION 17: MAIN PROCESSING FUNCTION
# ==============================================================

def process_one(sys_base, direction, target_str,
                T_df, N_df, E_df,
                method_label, method_dir,
                seeds_dict, gap_penalty, k_paths,
                max_clusters, min_clusters, fixed_k,
                run_auto, run_fixed,
                cb_coords, cb_sigma,
                syslog,
                cb_cache=None):
    """Full pipeline for one system x direction x target x method.
    Returns result dict including E_acc (edge accumulation DataFrame)
    for consolidation across targets/directions in _run_one_system."""
    sys_out = os.path.join(method_dir, sys_base)
    os.makedirs(sys_out, exist_ok=True)
    tag = f"{sys_base}_{direction}_{target_str}"
    log_msg(f"    [{method_label}] {tag}", syslog)

    target_nodes = (TARGET_FILTER_NODES if target_str == 'filter'
                    else TARGET_GATE_NODES)

    # ---- node_importance (Method 2 only) ----
    node_importance = None
    if method_label == 'method2':
        node_importance = compute_node_importance(N_df, seeds_dict)
        log_msg(
            f"      node_importance: "
            f"{int((node_importance > 0).sum())} non-zero nodes", syslog
        )

    # ---- Build graph ----
    G = build_graph(T_df, node_importance)
    log_msg(
        f"      Graph: {G.number_of_nodes()} nodes, "
        f"{G.number_of_edges()} edges", syslog
    )

    # ---- Generate paths ----
    valid_seeds = [s for s in seeds_dict if s in G.nodes]
    log_msg(
        f"      Seeds in graph: {len(valid_seeds)} / {len(seeds_dict)}",
        syslog
    )

    raw_paths = []
    for seed in valid_seeds:
        candidates = find_k_paths(G, seed, target_nodes, k_paths, syslog)
        for rank, (weight, path) in enumerate(candidates):
            rel_path, chain_orig = relativize_path(path)
            raw_paths.append({
                'seed':          seed,
                'rank_in_seed':  rank,
                'chain_origin':  chain_orig,
                'path_original': path,
                'path_rel':      rel_path,
                'path_length':   len(path),
                'log_prob':      path_log_prob(path, T_df),
                'imp_score':     path_importance_score(path, E_df),
                'weight':        weight,
                'target_hit':    path[-1] if path else '',
            })

    log_msg(f"      Raw paths: {len(raw_paths)}", syslog)
    if not raw_paths:
        log_msg(f"      No paths — skipping {tag}", syslog)
        return None

    # ---- Deduplicate / consolidate identical relativized paths ----
    deduped  = deduplicate_paths(raw_paths)
    n_unique = len(deduped)
    log_msg(
        f"      After deduplication: {n_unique} unique paths "
        f"(from {sum(p['n_duplicates'] for p in deduped)} original paths)",
        syslog
    )

    # ---- Edge accumulation (returned for consolidation in _run_one_system) ----
    # Weight = E_df[u_rel, v_rel] * n_duplicates  (or n_dup fallback)
    all_labels_T = T_df.index.tolist()
    E_acc = build_edge_accumulation(deduped, all_labels_T, E_df, syslog)

    # ---- NW distance matrix (computed once, reused by both modes) ----
    rel_list = [p['path_rel'] for p in deduped]
    log_msg(
        f"      NW distance matrix {n_unique}x{n_unique} "
        f"(Cb-aware={cb_cache is not None or (cb_coords is not None and USE_CB_DISTANCES)})...",
        syslog
    )
    D = build_distance_matrix(rel_list, gap_penalty, cb_coords, cb_sigma,
                              cb_cache)

    # ---- Linkage (computed once, reused by both modes) ----
    Z, valid_linkage = _compute_linkage(D)
    if not valid_linkage:
        log_msg(
            f"      WARNING: all-zero distance matrix (all paths identical"
            f" after dedup). 1 cluster for all modes.", syslog
        )
        labels_degenerate = np.ones(n_unique, dtype=int)
        Z_use             = None
        n_auto            = 1
        n_fixed           = 1
        sil_scores        = {}
        sil_method        = 'degenerate'
    else:
        n_auto, sil_scores, sil_method = _determine_n_clusters_auto(
            Z, D, n_unique, max_clusters, min_clusters
        )
        n_fixed = min(fixed_k, n_unique) if n_unique >= 2 else 1
        Z_use   = Z
        log_msg(
            f"      auto k={n_auto} [{sil_method}] | "
            f"fixed k={n_fixed} | sil={sil_scores}",
            syslog
        )

    # ---- AUTO clustering ----
    all_results = []
    if run_auto:
        if valid_linkage:
            labels_auto = fcluster(Z_use, n_auto, criterion='maxclust')
        else:
            labels_auto = labels_degenerate
        auto_dir = os.path.join(sys_out, "auto_cluster")
        _save_cluster_outputs(
            deduped, rel_list, D, labels_auto, n_auto, Z_use,
            auto_dir, tag, 'auto_cluster',
            gap_penalty, cb_coords, cb_sigma, syslog, cb_cache
        )
        all_results.append(('auto', n_auto, sil_method))

    # ---- FIXED clustering ----
    if run_fixed:
        if valid_linkage and n_unique >= 2:
            labels_fixed = fcluster(Z_use, n_fixed, criterion='maxclust')
        elif valid_linkage:
            labels_fixed = np.ones(n_unique, dtype=int)
        else:
            labels_fixed = labels_degenerate
        fixed_dir = os.path.join(sys_out, "fixed_cluster")
        _save_cluster_outputs(
            deduped, rel_list, D, labels_fixed, n_fixed, Z_use,
            fixed_dir, tag, 'fixed_cluster',
            gap_penalty, cb_coords, cb_sigma, syslog, cb_cache
        )
        all_results.append(('fixed', n_fixed, 'fixed'))

    return {
        'system':    sys_base,
        'direction': direction,
        'target':    target_str,
        'method':    method_label,
        'n_raw':     len(raw_paths),
        'n_unique':  n_unique,
        'n_auto':    n_auto if valid_linkage else 1,
        'n_fixed':   n_fixed if valid_linkage else 1,
        'sil_method': sil_method if valid_linkage else 'degenerate',
        'sil_scores': sil_scores,
        'E_acc':      E_acc,
    }

# ==============================================================
# SECTION 18: WORKER WRAPPER
# ==============================================================

def _run_one_system(args):
    """Top-level worker: all direction/target/method combos for one system.
    Accumulates E_acc per method_label across all directions and targets,
    then saves one consolidated matrix and one Gephi file per method."""
    (sys_base, directions, targets, methods_cfg,
     consolidated_dir, edge_imp_dir,
     method_matrices_dirs, method_gephi_dirs,
     seeds_dict, gap_penalty, k_paths,
     max_clusters, min_clusters, fixed_k,
     run_auto, run_fixed,
     all_cb_caches, cb_sigma,
     logs_dir, timestamp) = args

    sys_log_path = os.path.join(
        logs_dir, f"08h_{sys_base}_{timestamp}.txt"
    )

    # Select Cb sim cache for this system's condition PDB
    pdb_key      = sys_base.split('_')[0]
    cb_cache_sys = all_cb_caches.get(pdb_key) if all_cb_caches else None

    results = []
    # Accumulate E_acc per method_label (sum across directions and targets)
    e_acc_by_method = {}

    for direction in directions:
        for target_str in targets:
            T_df, N_df = load_system_matrices(
                consolidated_dir, sys_base, direction, target_str
            )
            if T_df is None:
                continue
            E_df = load_edge_importance(
                edge_imp_dir, sys_base, direction, target_str
            )

            for method_label, method_dir in methods_cfg:
                try:
                    result = process_one(
                        sys_base      = sys_base,
                        direction     = direction,
                        target_str    = target_str,
                        T_df          = T_df,
                        N_df          = N_df,
                        E_df          = E_df,
                        method_label  = method_label,
                        method_dir    = method_dir,
                        seeds_dict    = seeds_dict,
                        gap_penalty   = gap_penalty,
                        k_paths       = k_paths,
                        max_clusters  = max_clusters,
                        min_clusters  = min_clusters,
                        fixed_k       = fixed_k,
                        run_auto      = run_auto,
                        run_fixed     = run_fixed,
                        cb_coords     = None,
                        cb_sigma      = cb_sigma,
                        syslog        = sys_log_path,
                        cb_cache      = cb_cache_sys,
                    )
                    if result is not None:
                        # Extract and accumulate E_acc; remove from result
                        # before appending to results list
                        e_acc_new = result.pop('E_acc', None)
                        results.append(result)
                        if e_acc_new is not None:
                            if method_label not in e_acc_by_method:
                                e_acc_by_method[method_label] = e_acc_new.copy()
                            else:
                                e_acc_by_method[method_label] = (
                                    e_acc_by_method[method_label]
                                    .add(e_acc_new, fill_value=0)
                                )
                except Exception as exc:
                    try:
                        with open(sys_log_path, 'a') as fh:
                            fh.write(
                                f"ERROR {sys_base}/{direction}/"
                                f"{target_str}/{method_label}: {exc}\n"
                            )
                    except Exception:
                        pass

    # ---- Save consolidated matrices and Gephi per method ----
    for method_label, E_acc_total in e_acc_by_method.items():
        mat_dir   = method_matrices_dirs.get(method_label, '')
        gephi_dir = method_gephi_dirs.get(method_label, '')
        os.makedirs(mat_dir,   exist_ok=True)
        os.makedirs(gephi_dir, exist_ok=True)

        mat_fname   = f"matrix_{sys_base}_{method_label}.csv"
        gephi_fname = f"{sys_base}_{method_label}_edges_gephi.csv"

        save_distogram_csv(
            E_acc_total,
            os.path.join(mat_dir, mat_fname)
        )
        save_gephi_csv(
            E_acc_total,
            os.path.join(gephi_dir, gephi_fname),
            type_label=method_label,
        )
        try:
            with open(sys_log_path, 'a') as fh:
                fh.write(
                    f"[{datetime.now().strftime('%Y-%m-%d %H:%M:%S')}]"
                    f" 08h_v1.4 | Consolidated matrix saved:"
                    f" {mat_fname}\n"
                )
                fh.write(
                    f"[{datetime.now().strftime('%Y-%m-%d %H:%M:%S')}]"
                    f" 08h_v1.4 | Consolidated gephi saved:"
                    f" {gephi_fname}\n"
                )
        except Exception:
            pass

    return results

# ==============================================================
# SECTION 19: MAIN BLOCK (Windows-safe multiprocessing guard)
# ==============================================================

if __name__ == '__main__':

    log_msg("=== STARTING 08h_v1.4 PATH RECONSTRUCTION ===")
    log_msg(f"IS_NMA: {IS_NMA} | METHOD1: {RUN_METHOD1} | METHOD2: {RUN_METHOD2}")
    log_msg(f"RUN_FWD: {RUN_FWD} | RUN_REV: {RUN_REV}")
    log_msg(f"K_PATHS: {K_PATHS} | GAP_PENALTY: {GAP_PENALTY}")
    log_msg(f"USE_CB_DISTANCES: {USE_CB_DISTANCES} | CB_SIGMA: {CB_SIGMA_ANGSTROM} A")
    log_msg(f"RUN_CLUSTER_AUTO: {RUN_CLUSTER_AUTO} | RUN_CLUSTER_FIXED: {RUN_CLUSTER_FIXED}")
    log_msg(f"FIXED_K: {FIXED_K_CLUSTERS} | MAX_CLUSTERS: {MAX_CLUSTERS}")
    log_msg(f"N_CORES: {N_CORES} | Seeds defined: {len(SEEDS)}")
    log_msg(f"INPUT_08G_DIR: {INPUT_08G_DIR}")
    log_msg(f"OUTPUT_DIR:    {OUTPUT_DIR}")

    # ---- Load Cb coordinates and precompute sim caches (one per PDB) ----
    all_cb_caches = {}
    if USE_CB_DISTANCES:
        for pdb_key, pdb_name in PDB_FILES.items():
            pdb_path = os.path.join(PDB_DIR, f"{pdb_name}.pdb")
            log_msg(f"Loading Cb coordinates from: {pdb_path}")
            cb_coords_pdb = load_cb_coords_from_pdb(pdb_path)
            if cb_coords_pdb:
                log_msg(
                    f"  {pdb_key}: {len(cb_coords_pdb)} residues loaded. "
                    f"Precomputing sim cache..."
                )
                all_cb_caches[pdb_key] = precompute_cb_sim_cache(
                    cb_coords_pdb, CB_SIGMA_ANGSTROM
                )
                log_msg(
                    f"  {pdb_key}: sim cache ready "
                    f"({len(all_cb_caches[pdb_key])} pairs)"
                )
            else:
                log_msg(
                    f"  WARNING: no Cb coordinates loaded for {pdb_key} "
                    f"— will fall back to binary NW for this condition"
                )
        if not all_cb_caches:
            log_msg(
                "WARNING: no Cb caches built — "
                "falling back to binary NW substitution for all systems"
            )

    # ---- Discover systems ----
    log_msg("Discovering systems from consolidated/...")
    all_systems = discover_systems(CONSOLIDATED_DIR)
    log_msg(f"Systems found: {len(all_systems)} -> {all_systems}")

    # ---- Assemble task list ----
    directions = []
    if RUN_FWD:
        directions.append('fwd')
    if RUN_REV:
        directions.append('rev')

    targets = ['filter', 'gate']

    methods_cfg = []
    if RUN_METHOD1:
        methods_cfg.append(('method1', METHOD1_DIR))
    if RUN_METHOD2:
        methods_cfg.append(('method2', METHOD2_DIR))

    # Method-specific output directories for consolidated matrices and gephi
    method_matrices_dirs = {
        'method1': MATRICES_METHOD1_DIR,
        'method2': MATRICES_METHOD2_DIR,
    }
    method_gephi_dirs = {
        'method1': GEPHI_METHOD1_DIR,
        'method2': GEPHI_METHOD2_DIR,
    }

    task_args = [
        (sys_base, directions, targets, methods_cfg,
         CONSOLIDATED_DIR, EDGE_IMP_DIR,
         method_matrices_dirs, method_gephi_dirs,
         SEEDS, GAP_PENALTY, K_PATHS,
         MAX_CLUSTERS, MIN_CLUSTERS, FIXED_K_CLUSTERS,
         RUN_CLUSTER_AUTO, RUN_CLUSTER_FIXED,
         all_cb_caches, CB_SIGMA_ANGSTROM,
         LOGS_DIR, _timestamp)
        for sys_base in all_systems
    ]

    # ---- Execute ----
    all_results = []

    if N_CORES > 1 and len(all_systems) > 1:
        log_msg(f"Running parallel: {N_CORES} workers, {len(all_systems)} systems")
        with ProcessPoolExecutor(max_workers=N_CORES) as executor:
            futures = {
                executor.submit(_run_one_system, args): args[0]
                for args in task_args
            }
            for future in as_completed(futures):
                sys_name = futures[future]
                try:
                    sys_results = future.result()
                    all_results.extend(sys_results)
                    log_msg(
                        f"  Completed: {sys_name} "
                        f"({len(sys_results)} combinations)"
                    )
                except Exception as exc:
                    log_msg(f"  ERROR {sys_name}: {exc}")
    else:
        log_msg("Running sequential...")
        for args in task_args:
            sys_results = _run_one_system(args)
            all_results.extend(sys_results)
            log_msg(f"  Completed: {args[0]}")

    # ---- Optional WT comparison block ----
    if RUN_WT_COMPARISON:
        log_msg("=== WT COMPARISON BLOCK (not yet implemented) ===")

    # ---- Global summary ----
    log_msg("\n=== GLOBAL SUMMARY ===")
    log_msg(f"Systems: {len(all_systems)} | Combinations: {len(all_results)}")
    for r in all_results:
        log_msg(
            f"  {r['system']:30s}  {r['direction']:3s}  {r['target']:6s}"
            f"  {r['method']:7s}"
            f"  raw={r['n_raw']:4d}  unique={r['n_unique']:3d}"
            f"  auto_k={r['n_auto']}  fixed_k={r['n_fixed']}"
            f"  [{r['sil_method']}]"
        )

    summary_df   = pd.DataFrame(
        [{k: v for k, v in r.items() if k != 'sil_scores'} for r in all_results]
    )
    summary_path = os.path.join(OUTPUT_DIR, "08h_global_summary.csv")
    if not summary_df.empty:
        summary_df.to_csv(summary_path, index=False)
        log_msg(f"Global summary: {summary_path}")

    log_msg("=== 08h_v1.4 FINISHED ===")
