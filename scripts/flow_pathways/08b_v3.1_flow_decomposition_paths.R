# ==========================================================
# 08b_v3.1_flow_decomposition_paths.R
# Author: DenyCB
#
# PURPOSE:
#   Decompose accumulated edge flows (from 08a max-flow analysis)
#   into allosteric paths, relativize chain labels, and consolidate
#   highly similar paths into representative "river" families using
#   Smith-Waterman alignment with spatial distance penalties and
#   Neighbor-Joining clustering.
#
#   Processing pipeline:
#   STEP 1  — Load flow RDS objects from 08a output
#   STEP 2  — Consolidate replicas (MD only): sum edge flows
#   STEP 3  — Load PDB and build Cα spatial distance matrix
#   STEP 4  — Flow decomposition (Dijkstra-based, iterative)
#   STEP 5  — Save raw paths (pre-relativization)
#   STEP 6  — Relativize chain labels (source → always chain A)
#   STEP 7  — Iterative edge reaggregation + re-decomposition loop
#   STEP 8  — Build final edge list and interaction matrix
#   STEP 9  — Smith-Waterman pairwise alignment of paths
#   STEP 10 — Neighbor-Joining tree and clustering into N clados
#   STEP 11 — Save consolidated river families + all intermediate outputs
#
# ANALYSIS MODES (controlled by switches):
#   - Undirected: uses acc_abs (total activity per edge, symmetric)
#   - Directed:   uses acc_fwd / acc_rev (net flow per direction)
#     Both can be run independently via RUN_UNDIRECTED / RUN_DIRECTED
#
# INPUTS:
#   - RDS files from 08a (intracellular, LBD, combined per system)
#   - Residue mapping CSV (node_id -> chain, residue, system)
#   - PDB files (for Cα spatial distance matrix used in alignment)
#
# OUTPUTS (per system, per mode):
#   - paths_raw_*        : paths before relativization (CSV + RDS)
#   - paths_individual_* : all individual paths with path_chainres_rel,
#                          before summarise collapse (CSV + RDS)
#   - paths_rel_*        : paths after relativization (CSV + RDS)
#   - edges_*            : Gephi-compatible edge CSV
#   - matrix_*           : residue x residue flow matrix (CSV)
#   - dist_matrix_*      : SW pairwise distance matrix between paths (CSV)
#   - sw_alignments_*    : pairwise SW alignments with gaps, long format (RDS)
#   - nj_tree_*          : Newick-format NJ tree (TXT)
#   - dendro_*           : dendrogram plot for manual threshold selection (PNG)
#   - rivers_*           : consolidated river families (CSV + RDS)
#
# NOTES:
#   - Compatible with Windows 10 via PSOCK parallel cluster
#   - NMA (no replicas) and MD (with replicas) handled via IS_NMA switch
#   - Smith-Waterman uses discrete distance bins based on Cα-Cα distance:
#     match: <5Å (score +2), near: 5-10Å (score +1), far: 10-15Å (score -1),
#     mismatch: >15Å (score -2), gap open: -3, gap extend: -1
#   - NJ tree cut at N_RIVERS clados; if not achievable, outputs full tree
#   - SW pairwise alignments saved as RDS (more efficient than CSV for
#     long-format data with ~11,175 pairs × ~20 positions per pair)
#
# ==========================================================

suppressPackageStartupMessages({
  library(dplyr)
  library(readr)
  library(stringr)
  library(purrr)
  library(parallel)
  library(igraph)
  library(bio3d)
  library(ape)        # for NJ tree and clustering
  library(ggplot2)
  library(ggdendro)
})

# ==========================================================
# SECTION 1: ANALYSIS SWITCHES
# ==========================================================

# Which flow accumulations to decompose
RUN_TOTAL          <- FALSE   # decompose acc_total (filter + gate combined)
RUN_FILTER_SEP     <- TRUE   # decompose acc_filter separately
RUN_GATE_SEP       <- TRUE   # decompose acc_gate separately

# Analysis modes
RUN_UNDIRECTED     <- TRUE   # use acc_abs (symmetric, direction-agnostic)
RUN_DIRECTED       <- FALSE  # use acc_fwd / acc_rev (requires 08a v1.3 RDS)

# Iterative consolidation loop
RUN_ITERATIVE_CONSOLIDATION <- TRUE
MAX_CONSOLIDATION_ITER      <- 2

# NJ clustering
RUN_NJ_CLUSTERING  <- TRUE   # run SW alignment + NJ tree + river clustering

# NMA vs MD mode
IS_NMA             <- FALSE   # TRUE = NMA (no replicas), FALSE = MD


MAX_PATHS_FOR_SW <- 150

# ==========================================================
# SECTION 2: INPUT/OUTPUT CONFIGURATION
# ==========================================================

# --- NMA (uncomment to use) ---
#INPUT_RDS_DIR <- "C:/DinamicasMoleculares/TRPV1_pipeline/output/08_Max_Flow/08a_v1.3_allosteric_edges_NMA/rds"
#MAP_FILE      <- "C:/DinamicasMoleculares/TRPV1_pipeline/config_files/trimmed_residue_map.csv"
#PDB_DIR       <- "C:/DinamicasMoleculares/TRPV1_pipeline/input/pdb/"
#OUTPUT_DIR    <- "C:/DinamicasMoleculares/TRPV1_pipeline/output/08_Max_Flow/08b_v3.1_flow_paths_NMA_score_not_sqrt_min9"

# --- MD (comment out NMA block above and uncomment below) ---
 INPUT_RDS_DIR <- "C:/DinamicasMoleculares/TRPV1_pipeline/output/08_Max_Flow/08a_v1.3_allosteric_edges_MD/rds"
 MAP_FILE      <- "C:/DinamicasMoleculares/TRPV1_pipeline/output/residue_map_MD.csv"
 PDB_DIR       <- "C:/DinamicasMoleculares/TRPV1_pipeline/input/pdb/"
 OUTPUT_DIR    <- "C:/DinamicasMoleculares/TRPV1_pipeline/output/08_Max_Flow/08b_v3.1_flow_paths_MD_min8"

dir.create(OUTPUT_DIR,                            recursive = TRUE,  showWarnings = FALSE)
dir.create(file.path(OUTPUT_DIR, "raw"),          showWarnings = FALSE)
dir.create(file.path(OUTPUT_DIR, "relativized"),  showWarnings = FALSE)
dir.create(file.path(OUTPUT_DIR, "edges"),        showWarnings = FALSE)
dir.create(file.path(OUTPUT_DIR, "matrices"),     showWarnings = FALSE)
dir.create(file.path(OUTPUT_DIR, "final"),        showWarnings = FALSE)
dir.create(file.path(OUTPUT_DIR, "rivers"),       showWarnings = FALSE)
dir.create(file.path(OUTPUT_DIR, "trees"),        showWarnings = FALSE)
dir.create(file.path(OUTPUT_DIR, "logs"),         showWarnings = FALSE)

# ==========================================================
# SECTION 3: ANALYSIS PARAMETERS
# ==========================================================

# Target residues
TARGET_FILTER <- 643
TARGET_GATE   <- 679

# Minimum path length in number of residues
MIN_PATH_LENGTH <- 8

# Minimum flow threshold for edge inclusion
MIN_FLOW_THRESHOLD <- 1e-6

# Maximum number of paths to extract per decomposition
MAX_PATHS <- 500

# Similarity threshold (kept for reference; NJ used instead)
SIMILARITY_THRESHOLD <- 0.75

# Number of river families (clados) to cut the NJ tree into
# Set to NULL to skip automatic cut and output full tree only
N_RIVERS <- 10

# Smith-Waterman alignment scoring (discrete Cα distance bins)
# Scores for substitution based on spatial distance between residues:
SW_MATCH_SCORE   <-  2   # Cα-Cα distance < 5 Å  (structurally equivalent)
SW_NEAR_SCORE    <-  1   # Cα-Cα distance 5-10 Å (structurally close)
SW_FAR_SCORE     <- -1   # Cα-Cα distance 10-15 Å (structurally distant)
SW_MISMATCH_SCORE <- -2  # Cα-Cα distance > 15 Å (structurally unrelated)
SW_GAP_OPEN      <- -3   # affine gap opening penalty
SW_GAP_EXTEND    <- -1   # affine gap extension penalty

# Distance bins for spatial substitution scoring (Å)
SW_DIST_NEAR  <-  5.0   # below this → match
SW_DIST_FAR   <- 10.0   # below this → near
SW_DIST_MAX   <- 15.0   # below this → far; above → mismatch

# Parallel cores (PSOCK cluster, Windows 10 compatible)
N_CORES <- 6

# Maximum candidate source nodes for Dijkstra (top-K by node_flow).
# Limits shortest_paths() call from O(all nodes) to O(K) destinations,
# drastically reducing per-iteration compute time.
# K_SOURCES=600 ensures the full intracelullar domain (res. 277-548) is
# reachable as candidate seeds; increase if deep intracelullar paths are
# still missing, decrease if runtime is prohibitive.
K_SOURCES <- 600

# ==========================================================
# SECTION 4: LOGGING SETUP
# ==========================================================

master_log <- file.path(
  OUTPUT_DIR, "logs",
  paste0("08b_v3.1_master_log_",
         format(Sys.time(), "%Y%m%d_%H%M%S"), ".txt")
)

log_msg <- function(msg, log_file = NULL) {
  txt <- sprintf("[%s] 08b_v3.1 | %s",
                 format(Sys.time(), "%Y-%m-%d %H:%M:%OS3"), msg)
  cat(txt, "\n")
  if (!is.null(log_file)) cat(txt, "\n", file = log_file, append = TRUE)
}

log_msg("=== STARTING 08b_v3.1 FLOW DECOMPOSITION + NJ CLUSTERING ===",
        master_log)
log_msg(paste("INPUT_RDS_DIR:", INPUT_RDS_DIR), master_log)
log_msg(paste("PDB_DIR:", PDB_DIR), master_log)
log_msg(paste("IS_NMA:", IS_NMA,
              "| RUN_UNDIRECTED:", RUN_UNDIRECTED,
              "| RUN_DIRECTED:", RUN_DIRECTED), master_log)
log_msg(paste("RUN_NJ_CLUSTERING:", RUN_NJ_CLUSTERING,
              "| N_RIVERS:", N_RIVERS), master_log)
log_msg(paste("SW_MATCH:", SW_MATCH_SCORE,
              "| SW_NEAR:", SW_NEAR_SCORE,
              "| SW_FAR:", SW_FAR_SCORE,
              "| SW_MISMATCH:", SW_MISMATCH_SCORE,
              "| GAP_OPEN:", SW_GAP_OPEN,
              "| GAP_EXTEND:", SW_GAP_EXTEND), master_log)

# ==========================================================
# SECTION 5: LOAD RESIDUE MAPPING
# ==========================================================

log_msg("Loading residue mapping...", master_log)
mapping <- read_csv(MAP_FILE, show_col_types = FALSE)
log_msg(paste("Mapping loaded:", nrow(mapping), "rows,",
              length(unique(mapping$system)), "systems"), master_log)

# ==========================================================
# SECTION 6: HELPER FUNCTIONS
# ==========================================================

# --- Extract system name from RDS filename ---
get_system_from_rds <- function(filename) {
  nm <- basename(filename) %>% str_remove("\\.rds$")
  nm %>% str_remove("_intracellular$|_LBD$|_combined$")
}

get_structure_from_system <- function(system_name) {
  str_split(system_name, "_")[[1]][1]
}

# --- Load PDB and build Cα spatial distance matrix ---
# Returns a named matrix where rownames/colnames are "chainResno" labels
# (e.g. "A300", "B426") matching the map_df format.
# Only residues present in map_df are included.
build_ca_distance_matrix <- function(structure_name, pdb_dir, map_df,
                                     log_file) {
  
  pdb_file <- file.path(pdb_dir,
                        paste0(tolower(structure_name), ".pdb"))
  
  if (!file.exists(pdb_file)) {
    log_msg(paste("  WARNING: PDB not found:", pdb_file,
                  "— spatial scoring unavailable"), log_file)
    return(NULL)
  }
  
  pdb <- tryCatch(read.pdb(pdb_file), error = function(e) {
    log_msg(paste("  ERROR reading PDB:", e$message), log_file)
    NULL
  })
  if (is.null(pdb)) return(NULL)
  
  ca_atoms <- pdb$atom[pdb$atom$elety == "CA", ]
  
  # Build label map: "chainResno" -> row index in ca_atoms
  ca_labels <- paste0(ca_atoms$chain, ca_atoms$resno)
  
  # Keep only residues present in map_df
  map_labels <- paste0(map_df$chain, map_df$residue)
  keep       <- ca_labels %in% map_labels
  
  ca_sub    <- ca_atoms[keep, ]
  ca_labels <- ca_labels[keep]
  
  if (nrow(ca_sub) == 0) {
    log_msg("  WARNING: no CA atoms matched map_df labels", log_file)
    return(NULL)
  }
  
  coords <- as.matrix(ca_sub[, c("x", "y", "z")])
  dmat   <- as.matrix(dist(coords))
  rownames(dmat) <- ca_labels
  colnames(dmat) <- ca_labels
  
  log_msg(paste("  Cα distance matrix built:", nrow(dmat), "x", ncol(dmat)),
          log_file)
  dmat
}

# --- Smith-Waterman alignment score between two paths ---
# path1, path2: character vectors of chain+residue labels (e.g. "A300")
# dmat: Cα distance matrix with matching rownames/colnames
# Returns: alignment score (higher = more similar)
sw_align_score <- function(path1, path2, dmat,
                           match_score, near_score, far_score,
                           mismatch_score, gap_open, gap_extend,
                           dist_near, dist_far, dist_max) {
  
  n <- length(path1)
  m <- length(path2)
  
  if (n == 0 || m == 0) return(0)
  
  # Substitution score between two residue labels
  sub_score <- function(r1, r2) {
    if (r1 == r2) return(match_score)
    if (is.null(dmat) || !r1 %in% rownames(dmat) ||
        !r2 %in% colnames(dmat)) return(mismatch_score)
    d <- dmat[r1, r2]
    if (d < dist_near) return(match_score)
    if (d < dist_far)  return(near_score)
    if (d < dist_max)  return(far_score)
    return(mismatch_score)
  }
  
  # Affine gap Smith-Waterman DP matrices
  # H: best score ending here
  # E: best score ending with gap in path2 (horizontal)
  # F: best score ending with gap in path1 (vertical)
  H <- matrix(0, n + 1, m + 1)
  E <- matrix(-Inf, n + 1, m + 1)
  F <- matrix(-Inf, n + 1, m + 1)
  
  max_score <- 0
  
  for (i in 1:n) {
    for (j in 1:m) {
      s <- sub_score(path1[i], path2[j])
      
      E[i+1, j+1] <- max(H[i+1, j] + gap_open + gap_extend,
                         E[i+1, j] + gap_extend)
      F[i+1, j+1] <- max(H[i, j+1] + gap_open + gap_extend,
                         F[i, j+1] + gap_extend)
      H[i+1, j+1] <- max(0,
                         H[i, j] + s,
                         E[i+1, j+1],
                         F[i+1, j+1])
      
      if (H[i+1, j+1] > max_score) max_score <- H[i+1, j+1]
    }
  }
  
  max_score
}

# --- Smith-Waterman traceback: return explicit pairwise alignment with gaps ---
# path1, path2: character vectors of chain+residue labels (e.g. "A300")
# dmat: Cα distance matrix with matching rownames/colnames
# Returns: data.frame with columns pos, residue_i, residue_j
#   where "-" indicates a gap. pos is the alignment column index.
sw_align_traceback <- function(path1, path2, dmat,
                               match_score, near_score, far_score,
                               mismatch_score, gap_open, gap_extend,
                               dist_near, dist_far, dist_max) {
  
  n <- length(path1)
  m <- length(path2)
  
  if (n == 0 || m == 0) {
    return(data.frame(pos = integer(0),
                      residue_i = character(0),
                      residue_j = character(0),
                      stringsAsFactors = FALSE))
  }
  
  # Substitution score between two residue labels
  sub_score <- function(r1, r2) {
    if (r1 == r2) return(match_score)
    if (is.null(dmat) || !r1 %in% rownames(dmat) ||
        !r2 %in% colnames(dmat)) return(mismatch_score)
    d <- dmat[r1, r2]
    if (d < dist_near) return(match_score)
    if (d < dist_far)  return(near_score)
    if (d < dist_max)  return(far_score)
    return(mismatch_score)
  }
  
  # Fill DP matrices (identical to sw_align_score)
  H <- matrix(0, n + 1, m + 1)
  E <- matrix(-Inf, n + 1, m + 1)
  F <- matrix(-Inf, n + 1, m + 1)
  # Traceback pointer matrix: 0=stop, 1=diag, 2=up(gap in j), 3=left(gap in i)
  ptr <- matrix(0L, n + 1, m + 1)
  
  max_score <- 0
  max_i <- 0L
  max_j <- 0L
  
  for (i in 1:n) {
    for (j in 1:m) {
      s <- sub_score(path1[i], path2[j])
      
      E[i+1, j+1] <- max(H[i+1, j] + gap_open + gap_extend,
                         E[i+1, j] + gap_extend)
      F[i+1, j+1] <- max(H[i, j+1] + gap_open + gap_extend,
                         F[i, j+1] + gap_extend)
      
      diag_val <- H[i, j] + s
      up_val   <- F[i+1, j+1]
      left_val <- E[i+1, j+1]
      
      best <- max(0, diag_val, up_val, left_val)
      H[i+1, j+1] <- best
      
      if (best == 0) {
        ptr[i+1, j+1] <- 0L
      } else if (best == diag_val) {
        ptr[i+1, j+1] <- 1L
      } else if (best == up_val) {
        ptr[i+1, j+1] <- 2L
      } else {
        ptr[i+1, j+1] <- 3L
      }
      
      if (best > max_score) {
        max_score <- best
        max_i <- i + 1L
        max_j <- j + 1L
      }
    }
  }
  
  # Traceback from max cell
  aln_i <- character(0)
  aln_j <- character(0)
  ci <- max_i
  cj <- max_j
  
  while (ci > 1 && cj > 1 && ptr[ci, cj] != 0L) {
    p <- ptr[ci, cj]
    if (p == 1L) {
      aln_i <- c(path1[ci - 1], aln_i)
      aln_j <- c(path2[cj - 1], aln_j)
      ci <- ci - 1L
      cj <- cj - 1L
    } else if (p == 2L) {
      aln_i <- c(path1[ci - 1], aln_i)
      aln_j <- c("-", aln_j)
      ci <- ci - 1L
    } else {
      aln_i <- c("-", aln_i)
      aln_j <- c(path2[cj - 1], aln_j)
      cj <- cj - 1L
    }
  }
  
  if (length(aln_i) == 0) {
    return(data.frame(pos = integer(0),
                      residue_i = character(0),
                      residue_j = character(0),
                      stringsAsFactors = FALSE))
  }
  
  data.frame(pos       = seq_along(aln_i),
             residue_i = aln_i,
             residue_j = aln_j,
             stringsAsFactors = FALSE)
}

# --- Compute pairwise SW distance matrix for a set of paths ---
# paths_vec: character vector of path strings (e.g. "A300-A301-B500")
# dmat: Cα distance matrix
# Returns: list with two elements:
#   $dist_mat   — symmetric distance matrix (1 - normalized_similarity)
#   $alignments — list of pairwise alignment data.frames, each with columns
#                 path_i, path_j, pos, residue_i, residue_j.
#                 path_i and path_j are 1-based indices into paths_vec.
#                 Stored as RDS (not CSV) for efficiency.
compute_sw_distance_matrix <- function(paths_vec, dmat,
                                       match_score, near_score,
                                       far_score, mismatch_score,
                                       gap_open, gap_extend,
                                       dist_near, dist_far, dist_max,
                                       log_file) {
  
  n      <- length(paths_vec)
  parsed <- lapply(paths_vec, function(p) str_split(p, "-")[[1]])
  
  # Self-alignment scores for normalization
  self_scores <- sapply(seq_len(n), function(i) {
    sw_align_score(parsed[[i]], parsed[[i]], dmat,
                   match_score, near_score, far_score, mismatch_score,
                   gap_open, gap_extend, dist_near, dist_far, dist_max)
  })
  
  dist_mat   <- matrix(0, n, n)
  alignments <- vector("list", n * (n - 1) / 2)
  aln_idx    <- 0L
  
  n_pairs <- n * (n - 1) / 2
  done    <- 0L
  
  for (i in 1:(n - 1)) {
    for (j in (i + 1):n) {
      
      score_ij <- sw_align_score(parsed[[i]], parsed[[j]], dmat,
                                 match_score, near_score, far_score,
                                 mismatch_score,
                                 gap_open, gap_extend,
                                 dist_near, dist_far, dist_max)
      
      # Normalize by geometric mean of self-scores
      denom <- sqrt(self_scores[i] * self_scores[j])
      sim   <- if (denom > 0) score_ij / denom else 0
      sim   <- min(sim, 1)
      
      dist_mat[i, j] <- 1 - sim
      dist_mat[j, i] <- 1 - sim
      
      # Compute and store traceback alignment for this pair
      aln_df <- sw_align_traceback(parsed[[i]], parsed[[j]], dmat,
                                   match_score, near_score, far_score,
                                   mismatch_score, gap_open, gap_extend,
                                   dist_near, dist_far, dist_max)
      aln_idx <- aln_idx + 1L
      if (nrow(aln_df) > 0) {
        aln_df$path_i <- i
        aln_df$path_j <- j
        alignments[[aln_idx]] <- aln_df[, c("path_i", "path_j",
                                            "pos", "residue_i", "residue_j")]
      }
      
      done <- done + 1L
      if (done %% 1000 == 0) {
        log_msg(paste("    SW alignment:", done, "/", n_pairs, "pairs"),
                log_file)
      }
    }
  }
  
  rownames(dist_mat) <- seq_len(n)
  colnames(dist_mat) <- seq_len(n)
  
  # Remove NULL entries from alignments list (pairs with empty traceback)
  alignments <- Filter(Negate(is.null), alignments)
  
  list(dist_mat = dist_mat, alignments = alignments)
}

# --- Build NJ tree and cut into N clados ---
# dist_mat: symmetric distance matrix
# n_rivers: number of groups to cut the tree into (NULL = no cut)
# Returns: list(tree, groups) where groups is integer vector of cluster IDs
build_nj_tree_and_cut <- function(dist_mat, n_rivers, log_file) {
  
  nj_tree <- tryCatch(
    ape::nj(as.dist(dist_mat)),
    error = function(e) {
      log_msg(paste("  NJ tree error:", e$message), log_file)
      NULL
    }
  )
  
  if (is.null(nj_tree)) return(list(tree = NULL, groups = NULL))
  
  groups <- NULL
  
  if (!is.null(n_rivers) && n_rivers >= 2) {
    # Cut tree into n_rivers groups using hclust on the distance matrix
    # (NJ trees don't support direct cutree; use hclust as proxy)
    hc     <- hclust(as.dist(dist_mat), method = "average")
    groups <- cutree(hc, k = min(n_rivers, nrow(dist_mat)))
    log_msg(paste("  NJ cut into", length(unique(groups)), "groups"),
            log_file)
  }
  
  list(tree = nj_tree, groups = groups, hclust = hc)
}

# --- Dijkstra-based flow decomposition ---
# --- Dijkstra-based flow decomposition ---
# Identical to 08b_v1.2 — see that script for full documentation.
decompose_flow_to_paths <- function(g, target_nodes, all_target_nodes, min_flow,
                                    max_paths, min_length, log_file) {
  
  flow_vec <- E(g)$flow
  if (is.null(flow_vec)) {
    log_msg("  ERROR: graph has no $flow edge attribute", log_file)
    return(list())
  }
  
  el    <- as_edgelist(g, names = FALSE)
  paths <- list()
  iter  <- 0L
  
  ##############################################################################################
  node_flow <- rep(0.0, vcount(g))
  for (ei in seq_len(nrow(el))) {
    node_flow[el[ei, 1]] <- node_flow[el[ei, 1]] + flow_vec[ei]
    node_flow[el[ei, 2]] <- node_flow[el[ei, 2]] + flow_vec[ei]
  }
  ##############################################################################################  
  
  # FIX B: Pre-build adjacency edge index to avoid O(E) linear search
  # per edge lookup inside the path bottleneck evaluation loop.
  # edge_index[[n]] gives the indices in el[] of all edges incident to node n.
  edge_index <- vector("list", vcount(g))
  for (ei in seq_len(nrow(el))) {
    edge_index[[el[ei, 1]]] <- c(edge_index[[el[ei, 1]]], ei)
    edge_index[[el[ei, 2]]] <- c(edge_index[[el[ei, 2]]], ei)
  }
  
  # FIX I4: Initialize dijkstra_weights before the while loop.
  # Updated incrementally inside the loop (only modified edges),
  # avoiding full O(E) recomputation on every iteration.
  dijkstra_weights <- ifelse(flow_vec > min_flow, 1 / flow_vec, Inf)
  
  while (iter < max_paths) {
    
    # FIX I2: Use pre-built edge_index for target_flow computation
    # instead of O(E) which() scan per target node.
    target_flow <- sapply(target_nodes, function(t) {
      sum(flow_vec[edge_index[[t]]])
    })
    
    best_target <- target_nodes[which.max(target_flow)]
    
    if (max(target_flow, na.rm = TRUE) < min_flow) break
    
    ###   node_flow <- sapply(seq_len(vcount(g)), function(n) {
    ###      adj_e <- which(el[, 1] == n | el[, 2] == n)
    ###      sum(flow_vec[adj_e])
    ###    })
    
    # FIX A: Limit candidate_sources to top-K nodes by node_flow.
    # Previously candidate_sources could include ~1800 nodes, causing
    # shortest_paths() to expand Dijkstra to nearly all graph nodes and
    # materialize ~1800 full path objects per iteration. K_SOURCES=600
    # ensures the full intracellular domain (res. 277-548) is reachable
    # as candidate seeds while keeping runtime manageable.
    top_k_idx <- order(node_flow, decreasing = TRUE)[
      seq_len(min(K_SOURCES, vcount(g)))]
    candidate_sources <- intersect(
      which(node_flow > min_flow &
              !seq_len(vcount(g)) %in% all_target_nodes),
      top_k_idx
    )
    if (length(candidate_sources) == 0) break
    
    sp_result <- shortest_paths(
      g,
      from    = V(g)[best_target],
      to      = V(g)[candidate_sources],
      weights = dijkstra_weights,
      output  = "vpath"
    )$vpath
    
    # FIX C: Trim sp_result to paths meeting min_length before bottleneck
    # evaluation loop. Avoids iterating over short/degenerate paths that
    # would be discarded anyway, reducing loop body executions.
    sp_result_trimmed <- sp_result[
      sapply(sp_result, function(sp) length(sp) >= min_length)]
    
    best_path_nodes <- NULL
    best_path_flow  <- -Inf
    
    for (sp in sp_result_trimmed) {
      path_nodes <- as.integer(sp)
      
      # FIX B (continued): Use pre-built edge_index for O(degree) edge lookup
      # instead of O(E) which() scan over all 21,042 rows of el[].
      edge_flows <- sapply(1:(length(path_nodes) - 1), function(i) {
        e <- intersect(edge_index[[path_nodes[i]]],
                       edge_index[[path_nodes[i + 1]]])
        if (length(e) == 0) return(0)
        flow_vec[e[1]]
      })
      
      bottleneck <- min(edge_flows)
      if (bottleneck > best_path_flow) {
        best_path_flow  <- bottleneck
        best_path_nodes <- path_nodes
      }
    }
    
    if (is.null(best_path_nodes) || best_path_flow <= min_flow) {
      adj_e <- which(el[, 1] == best_target | el[, 2] == best_target)
      flow_vec[adj_e] <- pmax(0, flow_vec[adj_e] - min_flow)
      # FIX I4 (continued): update dijkstra_weights for modified edges only
      dijkstra_weights[adj_e] <- ifelse(
        flow_vec[adj_e] > min_flow, 1 / flow_vec[adj_e], Inf)
      iter <- iter + 1L
      next
    }
    
    for (i in 1:(length(best_path_nodes) - 1)) {
      e <- intersect(edge_index[[best_path_nodes[i]]],
                     edge_index[[best_path_nodes[i + 1]]])
      ###      if (length(e) > 0) flow_vec[e[1]] <- max(0, flow_vec[e[1]] - best_path_flow)
      if (length(e) > 0) {
        delta <- flow_vec[e[1]] - max(0, flow_vec[e[1]] - best_path_flow)
        flow_vec[e[1]] <- max(0, flow_vec[e[1]] - best_path_flow)
        node_flow[el[e[1], 1]] <- node_flow[el[e[1], 1]] - delta
        node_flow[el[e[1], 2]] <- node_flow[el[e[1], 2]] - delta
        # FIX I4 (continued): update dijkstra_weights for this edge only
        dijkstra_weights[e[1]] <- ifelse(
          flow_vec[e[1]] > min_flow, 1 / flow_vec[e[1]], Inf)
      }
    }
    
    paths <- append(paths, list(list(
      nodes      = rev(best_path_nodes),
      flow       = best_path_flow,
      length     = length(best_path_nodes),
      target_res = best_target
    )))
    
    iter <- iter + 1L
    
    if (iter %% 50 == 0) {
      log_msg(paste("  Decomposition iter:", iter,
                    "| paths found:", length(paths)), log_file)
    }
  }
  
  log_msg(paste("  Decomposition complete:", length(paths), "paths in",
                iter, "iterations"), log_file)
  paths
}
# --- Annotate path nodes with chain+residue labels ---
annotate_path_nodes <- function(path_nodes, map_df) {
  paste0(
    map_df$chain[match(path_nodes, map_df$node_id)],
    map_df$residue[match(path_nodes, map_df$node_id)]
  )
}

# --- Relativize chain labels in a path ---
# Source (first node) → always chain A; all others rotate consistently.
relativize_chains <- function(path_chainres) {
  chain_order <- c("A", "B", "C", "D")
  chains   <- substr(path_chainres, 1, 1)
  residues <- substring(path_chainres, 2)
  offset   <- match(chains[1], chain_order) - 1
  rel      <- chain_order[((match(chains, chain_order) - 1 - offset) %% 4) + 1]
  paste0(rel, residues)
}

# --- Path similarity (Jaccard on residue numbers, kept for reference) ---
path_similarity_resonly <- function(p1_res, p2_res) {
  inter <- length(intersect(p1_res, p2_res))
  union <- length(union(p1_res, p2_res))
  if (union == 0) return(0)
  inter / union
}

# --- Build edge dataframe from a list of annotated paths ---
build_edges_from_paths <- function(path_df, type_label) {
  edge_list <- list()
  for (i in seq_len(nrow(path_df))) {
    nodes <- str_split(path_df$path_chainres[i], "-")[[1]]
    w     <- path_df$flow[i]
    if (length(nodes) < 2) next
    for (j in 1:(length(nodes) - 1)) {
      edge_list[[length(edge_list) + 1]] <- data.frame(
        Source = nodes[j], Target = nodes[j + 1],
        weight = w, type = type_label,
        stringsAsFactors = FALSE
      )
    }
  }
  if (length(edge_list) == 0) return(data.frame())
  bind_rows(edge_list) %>%
    group_by(Source, Target, type) %>%
    summarise(weight = sum(weight), .groups = "drop")
}

# --- Build symmetric flow matrix from edge dataframe ---
build_matrix_from_edges <- function(edge_df) {
  all_nodes <- sort(unique(c(edge_df$Source, edge_df$Target)))
  mat <- matrix(0, length(all_nodes), length(all_nodes),
                dimnames = list(all_nodes, all_nodes))
  for (i in seq_len(nrow(edge_df))) {
    r <- edge_df$Source[i]; c <- edge_df$Target[i]
    mat[r, c] <- mat[r, c] + edge_df$weight[i]
    mat[c, r] <- mat[c, r] + edge_df$weight[i]
  }
  mat
}

# --- Score paths for ranking ---
# Score = flow_sum * sqrt(length): rewards long paths with high flow
score_path <- function(flow, length) flow * length^2

# --- Full pipeline: decompose → relativize → iterate → NJ cluster ---
# Encapsulates all steps for a given flow vector and target set.
# mode_label: "undirected", "fwd", or "rev" — used in output filenames.
run_full_pipeline <- function(acc_vec, target_nodes, all_target_nodes, flow_label,
                              g_ref, map_df_sys, dmat,
                              sys_base, fam_label, mode_label,
                              sys_log) {
  
  if (is.null(acc_vec) || length(target_nodes) == 0) return(NULL)
  
  # --- STEP 4: Flow decomposition ---
  E(g_ref)$flow <- acc_vec
  log_msg(paste("  [", mode_label, "] Decomposing:", flow_label), sys_log)
  
  paths_raw <- decompose_flow_to_paths(
    g_ref, target_nodes, all_target_nodes,
    MIN_FLOW_THRESHOLD, MAX_PATHS, MIN_PATH_LENGTH, sys_log
  )
  
  if (length(paths_raw) == 0) {
    log_msg(paste("  No paths found for:", flow_label), sys_log)
    return(NULL)
  }
  
  # Annotate: p$nodes are already in source→target order (reversed inside
  # decompose_flow_to_paths). Truncate at the first occurrence of the specific
  # target residue number for this path (identified via p$target_res) to remove
  # spurious intermediate target crossings allowed by the undirected graph.
  # After truncation, re-apply MIN_PATH_LENGTH and discard paths below it.
  paths_annotated <- lapply(paths_raw, function(p) {
    labels <- annotate_path_nodes(p$nodes, map_df_sys)
    # Identify the residue number of the specific target node for this path
    target_resnum <- as.character(
      map_df_sys$residue[match(p$target_res, map_df_sys$node_id)])
    # Find all positions where this target residue number appears in the path
    res_nums <- substring(labels, 2)
    target_positions <- which(res_nums == target_resnum)
    # Truncate at first occurrence if it appears before the last position
    if (length(target_positions) > 0 && target_positions[1] < length(labels)) {
      labels <- labels[1:target_positions[1]]
    }
    # Re-apply MIN_PATH_LENGTH after truncation
    if (length(labels) < MIN_PATH_LENGTH) return(NULL)
    list(path_chainres = paste(labels, collapse = "-"),
         flow = p$flow, length = length(labels), target_res = p$target_res)
  })
  
  # Remove paths that became too short after truncation
  paths_annotated <- Filter(Negate(is.null), paths_annotated)
  
  if (length(paths_annotated) == 0) {
    log_msg(paste("  No paths passed MIN_PATH_LENGTH after truncation for:",
                  flow_label), sys_log)
    return(NULL)
  }
  
  path_df <- data.frame(
    path_chainres = sapply(paths_annotated, `[[`, "path_chainres"),
    flow          = sapply(paths_annotated, `[[`, "flow"),
    length        = sapply(paths_annotated, `[[`, "length"),
    target_node   = sapply(paths_annotated, `[[`, "target_res"),
    flow_type     = flow_label,
    stringsAsFactors = FALSE
  )
  
  # --- STEP 5: Save raw paths ---
  tag_raw <- paste0(sys_base, "_", fam_label, "_", mode_label, "_", flow_label)
  
  saveRDS(path_df,
          file.path(OUTPUT_DIR, "raw",
                    paste0("paths_raw_", tag_raw, ".rds")))
  write_csv(path_df,
            file.path(OUTPUT_DIR, "raw",
                      paste0("paths_raw_", tag_raw, ".csv")))
  log_msg(paste("  Raw paths saved:", nrow(path_df)), sys_log)
  
  # --- STEP 6: Relativize chain labels ---
  log_msg("  Relativizing chain labels...", sys_log)
  
  path_df$path_chainres_rel <- sapply(path_df$path_chainres, function(p) {
    paste(relativize_chains(str_split(p, "-")[[1]]), collapse = "-")
  })
  
  # Save individual paths with path_chainres_rel before summarise collapse.
  # This preserves all individual paths (one row per path, not grouped) for
  # use in the next script's consensus computation and flexible re-clustering.
  saveRDS(path_df,
          file.path(OUTPUT_DIR, "raw",
                    paste0("paths_individual_", tag_raw, ".rds")))
  write_csv(path_df,
            file.path(OUTPUT_DIR, "raw",
                      paste0("paths_individual_", tag_raw, ".csv")))
  log_msg(paste("  Individual paths saved:", nrow(path_df)), sys_log)
  
  # Group by relativized path, compute chain coverage metrics
  rel_summary <- path_df %>%
    group_by(path_chainres_rel, flow_type) %>%
    summarise(
      n_chains_contributing = n_distinct(
        substr(sapply(str_split(path_chainres, "-"), `[`, 1), 1, 1)
      ),
      flow_sum      = sum(flow),
      flow_max      = max(flow),
      length        = first(length),
      path_chainres = path_chainres[which.max(flow)],
      .groups = "drop"
    ) %>%
    mutate(
      chain_cov_all4 = n_chains_contributing == 4,
      chain_cov_min2 = n_chains_contributing >= 2
    )
  
  saveRDS(rel_summary,
          file.path(OUTPUT_DIR, "relativized",
                    paste0("paths_rel_", tag_raw, ".rds")))
  write_csv(rel_summary,
            file.path(OUTPUT_DIR, "relativized",
                      paste0("paths_rel_", tag_raw, ".csv")))
  log_msg(paste("  Relativized paths:", nrow(rel_summary)), sys_log)
  
  # --- STEP 7: Iterative consolidation loop ---
  final_paths <- rel_summary
  
  if (RUN_ITERATIVE_CONSOLIDATION) {
    
    log_msg("  Iterative consolidation loop...", sys_log)
    current_paths <- rel_summary
    
    for (cons_iter in seq_len(MAX_CONSOLIDATION_ITER)) {
      
      # FIX C2: convergence criterion — compare path sets between iterations
      # instead of checking first chain (which is always A after relativization
      # and therefore always triggers immediate break).
      if (cons_iter > 1) {
        if (nrow(path_df_iter) == nrow(current_paths) &&
            all(sort(path_df_iter$path_chainres_rel) ==
                sort(current_paths$path_chainres_rel))) {
          log_msg(paste("  Convergence at iter:", cons_iter), sys_log)
          break
        }
      }
      
      edges_iter <- build_edges_from_paths(
        current_paths %>%
          mutate(path_chainres = path_chainres_rel) %>%
          select(path_chainres, flow = flow_sum),
        paste0(fam_label, "_iter", cons_iter)
      )
      
      if (nrow(edges_iter) == 0) break
      
      all_nodes_iter <- sort(unique(c(edges_iter$Source, edges_iter$Target)))
      node_idx_iter  <- setNames(seq_along(all_nodes_iter), all_nodes_iter)
      
      g_iter <- graph_from_data_frame(
        data.frame(from = node_idx_iter[edges_iter$Source],
                   to   = node_idx_iter[edges_iter$Target],
                   flow = edges_iter$weight),
        directed = FALSE,
        vertices = data.frame(name = seq_along(all_nodes_iter))
      )
      E(g_iter)$flow <- edges_iter$weight
      
      target_nodes_iter <- which(
        as.integer(substring(all_nodes_iter, 2)) %in%
          c(TARGET_FILTER, TARGET_GATE)
      )
      
      if (length(target_nodes_iter) == 0) break
      
      paths_iter <- decompose_flow_to_paths(
        g_iter, target_nodes_iter, target_nodes_iter,
        MIN_FLOW_THRESHOLD, MAX_PATHS, MIN_PATH_LENGTH, sys_log
      )
      
      if (length(paths_iter) == 0) break
      
      path_df_iter <- data.frame(
        path_chainres_rel = sapply(paths_iter, function(p) {
          paste(all_nodes_iter[p$nodes], collapse = "-")
        }),
        flow = sapply(paths_iter, `[[`, "flow"),
        length = sapply(paths_iter, `[[`, "length"),
        flow_type = paste0(fam_label, "_iter", cons_iter),
        stringsAsFactors = FALSE
      )
      
      path_df_iter$path_chainres_rel <- sapply(
        path_df_iter$path_chainres_rel, function(p) {
          paste(relativize_chains(str_split(p, "-")[[1]]), collapse = "-")
        }
      )
      
      # Truncate at first occurrence of target residue in each iterative path,
      # for the same reason as in the initial annotation block: the undirected
      # graph allows paths to cross target nodes before reaching best_target.
      # Also re-apply MIN_PATH_LENGTH after truncation.
      target_resnums_iter <- as.character(c(TARGET_FILTER, TARGET_GATE))
      path_df_iter <- path_df_iter[sapply(seq_len(nrow(path_df_iter)), function(ii) {
        nodes    <- str_split(path_df_iter$path_chainres_rel[ii], "-")[[1]]
        res_nums <- substring(nodes, 2)
        tpos     <- which(res_nums %in% target_resnums_iter)
        if (length(tpos) > 0 && tpos[1] < length(nodes)) {
          nodes <- nodes[1:tpos[1]]
          path_df_iter$path_chainres_rel[ii] <<- paste(nodes, collapse = "-")
          path_df_iter$length[ii]            <<- length(nodes)
        }
        length(nodes) >= MIN_PATH_LENGTH
      }), ]
      
      
      
      path_df_iter <- path_df_iter %>%
        group_by(path_chainres_rel, flow_type) %>%
        summarise(
          # FIX C3: correctly count distinct chains across all nodes of each
          # relativized path, not just the first character of the full string.
          n_chains_contributing = n_distinct(
            substr(sapply(str_split(path_chainres_rel, "-"), `[`, 1), 1, 1)
          ),
          flow_sum  = sum(flow),
          flow_max  = max(flow),
          length    = first(length),
          path_chainres = path_chainres_rel[which.max(flow)],
          .groups = "drop"
        ) %>%
        mutate(chain_cov_all4 = n_chains_contributing == 4,
               chain_cov_min2 = n_chains_contributing >= 2)
      
      saveRDS(path_df_iter,
              file.path(OUTPUT_DIR, "relativized",
                        paste0("paths_rel_iter", cons_iter, "_",
                               tag_raw, ".rds")))
      write_csv(path_df_iter,
                file.path(OUTPUT_DIR, "relativized",
                          paste0("paths_rel_iter", cons_iter, "_",
                                 tag_raw, ".csv")))
      
      log_msg(paste("  Iter", cons_iter, "paths:", nrow(path_df_iter)),
              sys_log)
      current_paths <- path_df_iter
    }
    
    final_paths <- current_paths
    log_msg(paste("  Consolidation done. Final paths:", nrow(final_paths)),
            sys_log)
  }
  
  # --- STEP 8: Final edge list and matrix ---
  edges_final <- build_edges_from_paths(
    final_paths %>%
      mutate(path_chainres = path_chainres_rel, flow = flow_sum) %>%
      select(path_chainres, flow, length, flow_type),
    paste0(fam_label, "_", mode_label)
  )
  
  write_csv(edges_final,
            file.path(OUTPUT_DIR, "edges",
                      paste0("edges_", tag_raw, ".csv")))
  
  mat_final <- build_matrix_from_edges(edges_final)
  write_csv(as.data.frame(mat_final),
            file.path(OUTPUT_DIR, "matrices",
                      paste0("matrix_", tag_raw, ".csv")))
  
  log_msg(paste("  Final edges saved:", nrow(edges_final)), sys_log)
  
  # Final scored summary
  final_summary <- final_paths %>%
    mutate(
      score        = score_path(flow_sum, length),
      path_resonly = sapply(path_chainres_rel, function(p) {
        nodes <- str_split(p, "-")[[1]]
        paste(substring(nodes, 2), collapse = " ")
      }),
      source_res = sapply(path_chainres_rel, function(p) {
        nodes <- str_split(p, "-")[[1]]
        substring(nodes[1], 2)
      }),
      target_res = sapply(path_chainres_rel, function(p) {
        nodes <- str_split(p, "-")[[1]]
        substring(nodes[length(nodes)], 2)
      }),
      system = sys_base,
      family = fam_label,
      mode   = mode_label
    ) %>%
    arrange(desc(score)) %>%
    mutate(path_id = paste0("path_", row_number())) %>%
    select(path_id, system, family, mode, path_chainres, path_chainres_rel,
           path_resonly, score, flow_sum, flow_max, length,
           source_res, target_res, n_chains_contributing,
           chain_cov_all4, chain_cov_min2)
  
  write_csv(final_summary,
            file.path(OUTPUT_DIR, "final",
                      paste0("paths_final_", tag_raw, ".csv")))
  saveRDS(final_summary,
          file.path(OUTPUT_DIR, "final",
                    paste0("paths_final_", tag_raw, ".rds")))
  
  # --- STEP 9-10: NJ clustering ---
  rivers_out <- NULL
  
  if (RUN_NJ_CLUSTERING && nrow(final_summary) >= 3) {
    
    log_msg(paste("  Computing SW distance matrix for",
                  nrow(final_summary), "paths..."), sys_log)
    
    paths_for_sw <- final_summary %>%
      arrange(desc(score)) %>%
      head(MAX_PATHS_FOR_SW)
    
    ###    dist_mat <- compute_sw_distance_matrix(
    ###      final_summary$path_chainres_rel,
    sw_result <- compute_sw_distance_matrix(
      paths_for_sw$path_chainres_rel,   
      dmat,
      SW_MATCH_SCORE, SW_NEAR_SCORE, SW_FAR_SCORE, SW_MISMATCH_SCORE,
      SW_GAP_OPEN, SW_GAP_EXTEND,
      SW_DIST_NEAR, SW_DIST_FAR, SW_DIST_MAX,
      sys_log
    )
    dist_mat   <- sw_result$dist_mat
    alignments <- sw_result$alignments
    
    write_csv(as.data.frame(dist_mat),
              file.path(OUTPUT_DIR, "trees",
                        paste0("dist_matrix_", tag_raw, ".csv")))
    
    # Save pairwise SW alignments as RDS.
    # Each element of the list is a data.frame with columns:
    #   path_i, path_j — 1-based indices into paths_for_sw
    #   pos            — alignment column index
    #   residue_i      — residue label from path_i ("-" = gap)
    #   residue_j      — residue label from path_j ("-" = gap)
    # path_ids from paths_for_sw are attached as attribute for reference.
    attr(alignments, "path_ids") <- paths_for_sw$path_id
    saveRDS(alignments,
            file.path(OUTPUT_DIR, "trees",
                      paste0("sw_alignments_", tag_raw, ".rds")))
    log_msg(paste("  SW alignments saved:", length(alignments), "pairs"),
            sys_log)
    
    log_msg("  Building NJ tree...", sys_log)
    
    nj_result <- build_nj_tree_and_cut(dist_mat, N_RIVERS, sys_log)
    
    if (!is.null(nj_result$tree)) {
      
      # Save Newick tree
      tryCatch(
        ape::write.tree(nj_result$tree,
                        file = file.path(OUTPUT_DIR, "trees",
                                         paste0("nj_tree_", tag_raw, ".nwk"))),
        error = function(e) log_msg(paste("  Newick write error:", e$message),
                                    sys_log)
      )
      
      # Save dendrogram plot for manual threshold inspection
      if (!is.null(nj_result$hclust)) {
        dendro_data <- ggdendro::dendro_data(nj_result$hclust)
        p_dendro <- ggplot() +
          geom_segment(
            data = dendro_data$segments,
            aes(x = x, y = y, xend = xend, yend = yend)
          ) +
          labs(title = paste("Hierarchical clustering:", tag_raw),
               y = "Distance", x = "") +
          theme_minimal() +
          theme(axis.text.x = element_blank())
        
        ggsave(file.path(OUTPUT_DIR, "trees",
                         paste0("dendro_", tag_raw, ".png")),
               p_dendro, width = 12, height = 6, dpi = 150)
      }
      
      # Assign river group IDs and select representative path per group
      if (!is.null(nj_result$groups)) {
        
        ###        final_summary$river_id <- nj_result$groups
        paths_for_sw$river_id <- nj_result$groups
        final_summary <- final_summary %>%
          left_join(paths_for_sw %>% select(path_id, river_id),
                    by = "path_id")
        
        # FIX I3: log number of paths outside top-50 that have river_id = NA
        n_no_river <- sum(is.na(final_summary$river_id))
        if (n_no_river > 0) {
          log_msg(paste("  NOTE:", n_no_river,
                        "paths outside top-", MAX_PATHS_FOR_SW,
                        "SW set have river_id = NA"), sys_log)
        }
        
        # Representative: path with highest score in each river group
        representatives <- final_summary %>%
          group_by(river_id) %>%
          slice_max(score, n = 1, with_ties = FALSE) %>%
          ungroup() %>%
          mutate(river_label = paste0("river_", river_id))
        
        rivers_out <- representatives
        
        write_csv(final_summary,
                  file.path(OUTPUT_DIR, "rivers",
                            paste0("all_paths_with_rivers_", tag_raw, ".csv")))
        
        write_csv(representatives,
                  file.path(OUTPUT_DIR, "rivers",
                            paste0("rivers_", tag_raw, ".csv")))
        
        saveRDS(representatives,
                file.path(OUTPUT_DIR, "rivers",
                          paste0("rivers_", tag_raw, ".rds")))
        
        log_msg(paste("  Rivers saved:", nrow(representatives),
                      "representative paths"), sys_log)
      }
    }
  }
  
  list(final_summary = final_summary, rivers = rivers_out)
}

# ==========================================================
# SECTION 7: FIND AND GROUP INPUT RDS FILES
# ==========================================================

rds_files <- list.files(INPUT_RDS_DIR, pattern = "\\.rds$",
                        full.names = TRUE)

if (length(rds_files) == 0) stop(paste("No RDS files found in:", INPUT_RDS_DIR))

log_msg(paste("RDS files found:", length(rds_files)), master_log)

parse_rds_meta <- function(filename) {
  nm <- basename(filename) %>% str_remove("\\.rds$")
  family <- case_when(
    str_detect(nm, "_intracellular$") ~ "intracellular",
    str_detect(nm, "_LBD$")          ~ "LBD",
    str_detect(nm, "_combined$")      ~ "combined",
    TRUE                              ~ "unknown"
  )
  system_name <- nm %>% str_remove("_intracellular$|_LBD$|_combined$")
  list(file = filename, system = system_name, family = family)
}

rds_meta <- lapply(rds_files, parse_rds_meta)
rds_df   <- bind_rows(lapply(rds_meta, as.data.frame,
                             stringsAsFactors = FALSE))

if (!IS_NMA) {
  rds_df <- rds_df %>%
    mutate(replica     = str_extract(system, "_rep[0-9]+"),
           system_base = str_remove(system, "_rep[0-9]+"))
} else {
  rds_df <- rds_df %>%
    mutate(replica = NA_character_, system_base = system)
}

systems_list <- rds_df %>%
  distinct(system_base, family) %>%
  split(seq(nrow(.)))

log_msg(paste("System-family combinations:", length(systems_list)), master_log)

# ==========================================================
# SECTION 8: PARALLEL PROCESSING
# ==========================================================

log_msg(paste("Starting cluster with", N_CORES, "cores (PSOCK)..."),
        master_log)

cl <- makeCluster(N_CORES)

clusterEvalQ(cl, {
  suppressPackageStartupMessages({
    library(dplyr)
    library(readr)
    library(stringr)
    library(igraph)
    library(bio3d)
    library(ape)
    library(ggplot2)
    library(ggdendro)
  })
})

clusterExport(cl, c(
  "rds_df", "mapping", "OUTPUT_DIR", "IS_NMA",
  "PDB_DIR",
  "TARGET_FILTER", "TARGET_GATE",
  "MIN_FLOW_THRESHOLD", "MIN_PATH_LENGTH", "MAX_PATHS",
  "SIMILARITY_THRESHOLD",
  "N_RIVERS", "MAX_PATHS_FOR_SW", "K_SOURCES",
  "RUN_TOTAL", "RUN_FILTER_SEP", "RUN_GATE_SEP",
  "RUN_UNDIRECTED", "RUN_DIRECTED",
  "RUN_ITERATIVE_CONSOLIDATION", "MAX_CONSOLIDATION_ITER",
  "RUN_NJ_CLUSTERING",
  "SW_MATCH_SCORE", "SW_NEAR_SCORE", "SW_FAR_SCORE", "SW_MISMATCH_SCORE",
  "SW_GAP_OPEN", "SW_GAP_EXTEND",
  "SW_DIST_NEAR", "SW_DIST_FAR", "SW_DIST_MAX",
  "decompose_flow_to_paths", "annotate_path_nodes",
  "get_structure_from_system",
  "relativize_chains", "path_similarity_resonly",
  "build_edges_from_paths", "build_matrix_from_edges",
  "build_ca_distance_matrix",
  "sw_align_score", "sw_align_traceback", "compute_sw_distance_matrix",
  "build_nj_tree_and_cut",
  "score_path", "run_full_pipeline",
  "log_msg"
))

results <- parLapply(cl, systems_list, function(sys_row) {
  
  sys_base  <- sys_row$system_base
  fam_label <- sys_row$family
  
  sys_log <- file.path(OUTPUT_DIR, "logs",
                       paste0("08b_", sys_base, "_", fam_label, "_",
                              format(Sys.time(), "%Y%m%d_%H%M%S"), ".txt"))
  
  log_msg(paste("=== Processing:", sys_base, "| family:", fam_label, "==="),
          sys_log)
  t0 <- Sys.time()
  
  # --------------------------------------------------------
  # STEP 1: LOAD AND CONSOLIDATE RDS FILES
  # --------------------------------------------------------
  
  sys_files <- rds_df %>%
    filter(system_base == sys_base, family == fam_label) %>%
    pull(file)
  
  log_msg(paste("  RDS files:", length(sys_files)), sys_log)
  
  # Initialize accumulators
  acc_total_consolidated  <- NULL
  acc_filter_consolidated <- NULL
  acc_gate_consolidated   <- NULL
  # Directed accumulators (from 08a v1.3 RDS fields)
  acc_fwd_filter_consolidated <- NULL
  acc_rev_filter_consolidated <- NULL
  acc_fwd_gate_consolidated   <- NULL
  acc_rev_gate_consolidated   <- NULL
  g_ref      <- NULL
  map_df_sys <- NULL
  
  # Also load intracellular and LBD separately to reconstruct combined fwd/rev
  # if the combined RDS doesn't have fwd/rev (08a < v1.3)
  acc_intra_fwd_filter <- NULL; acc_intra_rev_filter <- NULL
  acc_intra_fwd_gate   <- NULL; acc_intra_rev_gate   <- NULL
  acc_lbd_fwd_filter   <- NULL; acc_lbd_rev_filter   <- NULL
  acc_lbd_fwd_gate     <- NULL; acc_lbd_rev_gate     <- NULL
  
  for (f in sys_files) {
    
    obj <- tryCatch(readRDS(f), error = function(e) {
      log_msg(paste("  ERROR loading:", f, "-", e$message), sys_log)
      NULL
    })
    if (is.null(obj)) next
    
    if (is.null(g_ref)) {
      g_ref      <- obj$graph
      map_df_sys <- obj$map_df
    }
    
    # Skip this file if its flow vectors don't match g_ref edge count
    if (!is.null(obj$acc_total) && length(obj$acc_total) != ecount(g_ref)) {
      log_msg(paste("  WARNING: skipping", basename(f),
                    "— flow vector length", length(obj$acc_total),
                    "!= g_ref edges", ecount(g_ref)), sys_log)
      next
    }
    
    # Detect family from filename
    f_fam <- case_when(
      str_detect(basename(f), "_intracellular") ~ "intracellular",
      str_detect(basename(f), "_LBD")           ~ "LBD",
      str_detect(basename(f), "_combined")       ~ "combined",
      TRUE                                       ~ "unknown"
    )
    
    # Accumulate abs flows
    if (!is.null(obj$acc_total)) {
      acc_total_consolidated <- if (is.null(acc_total_consolidated))
        obj$acc_total else acc_total_consolidated + obj$acc_total
    }
    if (!is.null(obj$acc_filter)) {
      acc_filter_consolidated <- if (is.null(acc_filter_consolidated))
        obj$acc_filter else acc_filter_consolidated + obj$acc_filter
    }
    if (!is.null(obj$acc_gate)) {
      acc_gate_consolidated <- if (is.null(acc_gate_consolidated))
        obj$acc_gate else acc_gate_consolidated + obj$acc_gate
    }
    
    # Accumulate directed flows if available (08a v1.3)
    if (!is.null(obj$acc_fwd_filter)) {
      if (f_fam == "intracellular") {
        acc_intra_fwd_filter <- if (is.null(acc_intra_fwd_filter))
          obj$acc_fwd_filter else acc_intra_fwd_filter + obj$acc_fwd_filter
        acc_intra_rev_filter <- if (is.null(acc_intra_rev_filter))
          obj$acc_rev_filter else acc_intra_rev_filter + obj$acc_rev_filter
        acc_intra_fwd_gate   <- if (is.null(acc_intra_fwd_gate))
          obj$acc_fwd_gate else acc_intra_fwd_gate + obj$acc_fwd_gate
        acc_intra_rev_gate   <- if (is.null(acc_intra_rev_gate))
          obj$acc_rev_gate else acc_intra_rev_gate + obj$acc_rev_gate
      } else if (f_fam == "LBD") {
        acc_lbd_fwd_filter <- if (is.null(acc_lbd_fwd_filter))
          obj$acc_fwd_filter else acc_lbd_fwd_filter + obj$acc_fwd_filter
        acc_lbd_rev_filter <- if (is.null(acc_lbd_rev_filter))
          obj$acc_rev_filter else acc_lbd_rev_filter + obj$acc_rev_filter
        acc_lbd_fwd_gate   <- if (is.null(acc_lbd_fwd_gate))
          obj$acc_fwd_gate else acc_lbd_fwd_gate + obj$acc_fwd_gate
        acc_lbd_rev_gate   <- if (is.null(acc_lbd_rev_gate))
          obj$acc_rev_gate else acc_lbd_rev_gate + obj$acc_rev_gate
      }
    }
  }
  
  # Reconstruct combined fwd/rev from intracellular + LBD if available
  if (!is.null(acc_intra_fwd_filter) && !is.null(acc_lbd_fwd_filter)) {
    acc_fwd_filter_consolidated <- acc_intra_fwd_filter + acc_lbd_fwd_filter
    acc_rev_filter_consolidated <- acc_intra_rev_filter + acc_lbd_rev_filter
    acc_fwd_gate_consolidated   <- acc_intra_fwd_gate   + acc_lbd_fwd_gate
    acc_rev_gate_consolidated   <- acc_intra_rev_gate   + acc_lbd_rev_gate
    log_msg("  Directed combined fwd/rev reconstructed from intra + LBD",
            sys_log)
  }
  
  if (is.null(g_ref)) {
    log_msg("  ERROR: no valid RDS loaded. Skipping.", sys_log)
    return(NULL)
  }
  
  log_msg(paste("  Graph:", vcount(g_ref), "nodes,", ecount(g_ref), "edges"),
          sys_log)
  
  # --------------------------------------------------------
  # STEP 2: IDENTIFY TARGET NODES
  # --------------------------------------------------------
  
  target_filter_nodes <- map_df_sys %>%
    filter(residue == TARGET_FILTER) %>% pull(node_id) %>% as.integer()
  target_gate_nodes <- map_df_sys %>%
    filter(residue == TARGET_GATE) %>% pull(node_id) %>% as.integer()
  
  target_filter_nodes <- target_filter_nodes[
    target_filter_nodes >= 1 & target_filter_nodes <= vcount(g_ref)]
  target_gate_nodes <- target_gate_nodes[
    target_gate_nodes >= 1 & target_gate_nodes <= vcount(g_ref)]
  
  all_target_nodes <- unique(c(target_filter_nodes, target_gate_nodes))
  
  log_msg(paste("  Target FILTER nodes:", length(target_filter_nodes),
                "| GATE nodes:", length(target_gate_nodes)), sys_log)
  
  # --------------------------------------------------------
  # STEP 3: BUILD Cα DISTANCE MATRIX FOR SW ALIGNMENT
  # --------------------------------------------------------
  
  structure_name <- get_structure_from_system(sys_base)
  dmat           <- build_ca_distance_matrix(structure_name, PDB_DIR,
                                             map_df_sys, sys_log)
  
  # --------------------------------------------------------
  # STEPS 4-10: RUN PIPELINE PER FLOW TYPE AND MODE
  # --------------------------------------------------------
  
  all_results <- list()
  
  # Helper to run all flow types for a given mode
  run_mode <- function(mode_label,
                       total_vec, filter_vec, gate_vec) {
    
    if (RUN_TOTAL && !is.null(total_vec)) {
      r <- run_full_pipeline(
        total_vec, all_target_nodes, all_target_nodes, "total",
        g_ref, map_df_sys, dmat,
        sys_base, fam_label, mode_label, sys_log
      )
      if (!is.null(r)) all_results[[paste0(mode_label, "_total")]] <<- r
    }
    
    if (RUN_FILTER_SEP && !is.null(filter_vec)) {
      r <- run_full_pipeline(
        filter_vec, target_filter_nodes, all_target_nodes, "filter",
        g_ref, map_df_sys, dmat,
        sys_base, fam_label, mode_label, sys_log
      )
      if (!is.null(r)) all_results[[paste0(mode_label, "_filter")]] <<- r
    }
    
    if (RUN_GATE_SEP && !is.null(gate_vec)) {
      r <- run_full_pipeline(
        gate_vec, target_gate_nodes, all_target_nodes, "gate",
        g_ref, map_df_sys, dmat,
        sys_base, fam_label, mode_label, sys_log
      )
      if (!is.null(r)) all_results[[paste0(mode_label, "_gate")]] <<- r
    }
  }
  
  # Undirected mode
  if (RUN_UNDIRECTED) {
    log_msg("  --- UNDIRECTED MODE ---", sys_log)
    run_mode("undirected",
             acc_total_consolidated,
             acc_filter_consolidated,
             acc_gate_consolidated)
  }
  
  # Directed mode (forward: seed → target direction)
  if (RUN_DIRECTED) {
    if (!is.null(acc_fwd_filter_consolidated)) {
      log_msg("  --- DIRECTED FWD MODE ---", sys_log)
      # For directed total: sum fwd_filter + fwd_gate
      acc_fwd_total <- acc_fwd_filter_consolidated + acc_fwd_gate_consolidated
      run_mode("directed_fwd",
               acc_fwd_total,
               acc_fwd_filter_consolidated,
               acc_fwd_gate_consolidated)
      
      log_msg("  --- DIRECTED REV MODE ---", sys_log)
      acc_rev_total <- acc_rev_filter_consolidated + acc_rev_gate_consolidated
      run_mode("directed_rev",
               acc_rev_total,
               acc_rev_filter_consolidated,
               acc_rev_gate_consolidated)
    } else {
      log_msg("  WARNING: directed mode requested but fwd/rev vectors not found.",
              sys_log)
      log_msg("  Re-run 08a v1.3 to generate directed flow vectors.", sys_log)
    }
  }
  
  elapsed <- round(as.numeric(difftime(Sys.time(), t0, units = "secs")), 1)
  log_msg(paste("=== Finished:", sys_base, fam_label,
                "| elapsed:", elapsed, "s ==="), sys_log)
  
  list(system = sys_base, family = fam_label,
       n_results = length(all_results))
})

stopCluster(cl)

# ==========================================================
# SECTION 9: GLOBAL SUMMARY
# ==========================================================

completed <- Filter(Negate(is.null), results)
log_msg(paste("Completed:", length(completed), "system-family combinations"),
        master_log)

for (r in completed) {
  log_msg(paste(" ", r$system, r$family, "->", r$n_results, "analyses"),
          master_log)
}

log_msg("=== 08b_v3.1 FINISHED ===", master_log)

# ==========================================================
# SECTION 10: PARAMETER REGISTRY UPDATE
# ==========================================================

tryCatch({
  
  library(dplyr)
  
  config_file <- "config.txt"
  read_config <- function(f) {
    lines <- readLines(f)
    lines <- lines[!grepl("^#", lines) & nchar(lines) > 0]
    cfg <- list()
    for (l in lines) {
      parts <- strsplit(l, "=")[[1]]
      cfg[[trimws(parts[1])]] <- trimws(parts[2])
    }
    cfg
  }
  
  config <- tryCatch(read_config(config_file), error = function(e) NULL)
  
  if (!is.null(config)) {
    
    registry_dir  <- file.path(config$pipeline_output_dir, "registry")
    final_path    <- file.path(registry_dir, "parameters_final.csv")
    history_path  <- file.path(registry_dir, "parameters_history.csv")
    ts <- format(Sys.time(), "%Y-%m-%d %H:%M:%S")
    
    new_entries <- data.frame(
      script_number  = "08b",
      script_version = "3.1",
      script_name    = "08b_v3.1_flow_decomposition_paths",
      variable = c("INPUT_RDS_DIR", "MAP_FILE", "PDB_DIR", "OUTPUT_DIR",
                   "IS_NMA", "TARGET_FILTER", "TARGET_GATE",
                   "MIN_PATH_LENGTH", "MAX_PATHS", "MIN_FLOW_THRESHOLD",
                   "N_RIVERS", "RUN_UNDIRECTED", "RUN_DIRECTED",
                   "RUN_NJ_CLUSTERING",
                   "SW_MATCH_SCORE", "SW_NEAR_SCORE", "SW_FAR_SCORE",
                   "SW_MISMATCH_SCORE", "SW_GAP_OPEN", "SW_GAP_EXTEND",
                   "SW_DIST_NEAR", "SW_DIST_FAR", "SW_DIST_MAX",
                   "RUN_TOTAL", "RUN_FILTER_SEP", "RUN_GATE_SEP",
                   "RUN_ITERATIVE_CONSOLIDATION", "MAX_CONSOLIDATION_ITER",
                   "N_CORES", "K_SOURCES", "MAX_PATHS_FOR_SW",
                   "systems_completed"),
      value = c(INPUT_RDS_DIR, MAP_FILE, PDB_DIR, OUTPUT_DIR,
                as.character(IS_NMA),
                as.character(TARGET_FILTER), as.character(TARGET_GATE),
                as.character(MIN_PATH_LENGTH), as.character(MAX_PATHS),
                as.character(MIN_FLOW_THRESHOLD),
                as.character(N_RIVERS),
                as.character(RUN_UNDIRECTED), as.character(RUN_DIRECTED),
                as.character(RUN_NJ_CLUSTERING),
                as.character(SW_MATCH_SCORE), as.character(SW_NEAR_SCORE),
                as.character(SW_FAR_SCORE), as.character(SW_MISMATCH_SCORE),
                as.character(SW_GAP_OPEN), as.character(SW_GAP_EXTEND),
                as.character(SW_DIST_NEAR), as.character(SW_DIST_FAR),
                as.character(SW_DIST_MAX),
                as.character(RUN_TOTAL), as.character(RUN_FILTER_SEP),
                as.character(RUN_GATE_SEP),
                as.character(RUN_ITERATIVE_CONSOLIDATION),
                as.character(MAX_CONSOLIDATION_ITER),
                as.character(N_CORES), as.character(K_SOURCES),
                as.character(MAX_PATHS_FOR_SW),
                as.character(length(completed))),
      timestamp = ts,
      stringsAsFactors = FALSE
    )
    
    if (file.exists(final_path)) {
      existing <- read.csv(final_path, stringsAsFactors = FALSE,
                           colClasses = "character")
      existing <- existing[existing$script_number != "08b", ]
      updated  <- dplyr::bind_rows(existing, new_entries)
    } else {
      updated <- new_entries
    }
    write.csv(updated, final_path, row.names = FALSE)
    
    if (file.exists(history_path)) {
      history <- read.csv(history_path, stringsAsFactors = FALSE,
                          colClasses = "character")
      history <- dplyr::bind_rows(history, new_entries)
    } else {
      history <- new_entries
    }
    write.csv(history, history_path, row.names = FALSE)
    
    log_msg("Parameter registry updated.", master_log)
  }
  
}, error = function(e) {
  log_msg(paste("WARNING: registry update failed:", e$message), master_log)
})