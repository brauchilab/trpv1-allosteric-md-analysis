# ==========================================================
# 08g_v2.1_markov_transition_fundamental.R
# Author: DenyCB
#
# PURPOSE:
#   Build directed Markov transition matrices (T_fwd, T_rev) from
#   max-flow accumulated edge flows (08a output), compute the
#   absorbing Markov chain fundamental matrix N = (I - Q)^{-1}
#   for both forward and reverse directions, and produce:
#
#   1. Per-replica T matrices (fwd and rev) saved as RDS,
#      computed separately for filter and gate targets
#   2. Per-replica N matrices (fwd and rev) saved as RDS,
#      computed separately for filter and gate targets
#   3. Consolidated N matrices: mean and standard error across
#      replicas (MD only), saved as RDS and CSV,
#      separately for filter and gate targets
#   4. Edge importance matrices derived from N_mean and T_mean,
#      formatted for direct input to 08d visualization script
#      (CSV with row.names, two symmetric matrices per system
#       per target: one using fwd values, one using rev values
#       — NOT averaged together, saved separately as
#       edge_importance_fwd_filter, edge_importance_fwd_gate,
#       edge_importance_rev_filter, edge_importance_rev_gate)
#
#   The edge importance for direction d and target t is:
#     node_importance_d_t[i] = sum_j N_d_t_mean[j, i]
#       (expected visits to node i from all transient starting points)
#     edge_imp_d_t[i,j] = T_d_t_mean[i,j] x node_importance_d_t[i]
#   This is an intrinsic property of the network for each
#   direction and target, independent of user-specified seeds.
#   Seeds are used in downstream scripts (path reconstruction).
#
#   Each direction/target produces a symmetric matrix by setting:
#     E_sym[i,j] = E_sym[j,i] = (E[i,j] + E[j,i]) / 2
#   fwd and rev are NEVER averaged together.
#   filter and gate are NEVER averaged together.
#
#   Regularization (I - Q + eI)^{-1} is optionally applied to
#   avoid numerical instability when Q has eigenvalues near 1.
#
# INPUTS:
#   - RDS files from 08a (same directory as 08b input)
#
# OUTPUTS (per system, per direction fwd/rev, per target filter/gate):
#   - T matrices per replica: T_fwd_filter_<system>_rep<N>.rds
#   - N matrices per replica: N_fwd_filter_<system>_rep<N>.rds
#   - N mean matrix:  N_fwd_filter_<system>_mean.rds / .csv
#   - N stderr matrix: N_fwd_filter_<system>_stderr.rds / .csv
#   - T mean matrix:  T_fwd_filter_<system>_mean.rds
#   - Edge importance matrix (08d-compatible CSV):
#       edge_importance_fwd_filter_<system>.csv  (symmetric, fwd values)
#       edge_importance_rev_filter_<system>.csv  (symmetric, rev values)
#       edge_importance_fwd_gate_<system>.csv
#       edge_importance_rev_gate_<system>.csv
#   - logs/08g_<system>_<timestamp>.txt
#
# NOTES:
#   - Compatible with Windows 10 via PSOCK parallel cluster
#   - For NMA (IS_NMA=TRUE): no replicas; single T and N per system
#   - For MD (IS_NMA=FALSE): T and N computed per replica;
#     mean and stderr computed across replicas; T_mean is used
#     (not a single replica T) for the edge importance computation
#   - TARGET_FILTER=643 and TARGET_GATE=679 are set as separate
#     absorbing states in independent N computations per target
#   - Edge importance CSV uses write.csv with row.names=TRUE,
#     matching the format expected by 08d_v2.1
#   - fwd matrix: both [i,j] and [j,i] filled with fwd-derived value
#   - rev matrix: both [i,j] and [j,i] filled with rev-derived value
#   - fwd and rev are NEVER averaged together
#   - filter and gate are NEVER averaged together
#   - Nodes absent in some replicas are filled with NA (not zero)
#     before averaging; mean and stderr computed with na.rm=TRUE
#     over the effective n of non-NA replicas
#
# ==========================================================

suppressPackageStartupMessages({
  library(dplyr)
  library(readr)
  library(stringr)
  library(parallel)
  library(igraph)
  library(Matrix)   # for sparse matrix operations if needed
})

# ==========================================================
# SECTION 1: ANALYSIS SWITCHES
# ==========================================================

IS_NMA       <- TRUE    # TRUE = NMA (no replicas), FALSE = MD

# Which directions to process
RUN_FWD      <- TRUE    # forward: seed -> target
RUN_REV      <- TRUE    # reverse: target -> seed (hysteresis analysis)

# Regularization for matrix inversion: adds eI to (I - Q)
# to avoid numerical instability when Q has eigenvalues near 1.
# Set USE_REGULARIZATION <- FALSE to attempt plain inversion first.
USE_REGULARIZATION <- FALSE
REGULARIZATION_EPS <- 1e-6   # e value if USE_REGULARIZATION = TRUE

# ==========================================================
# SECTION 2: INPUT/OUTPUT CONFIGURATION
# ==========================================================

# --- NMA (uncomment to use) ---
#INPUT_RDS_DIR <- "C:/DinamicasMoleculares/TRPV1_pipeline/output/08_Max_Flow/08a_v1.3_allosteric_edges_NMA/rds"
#MAP_FILE      <- "C:/DinamicasMoleculares/TRPV1_pipeline/config_files/trimmed_residue_map.csv"
#OUTPUT_DIR    <- "C:/DinamicasMoleculares/TRPV1_pipeline/output/08_Max_Flow/08g_v2.1_markov_NMA"

# --- MD (comment out NMA block above and uncomment below) ---
 INPUT_RDS_DIR <- "C:/DinamicasMoleculares/TRPV1_pipeline/output/08_Max_Flow/08a_v1.3_allosteric_edges_MD/rds"
 MAP_FILE      <- "C:/DinamicasMoleculares/TRPV1_pipeline/output/residue_map_MD.csv"
 OUTPUT_DIR    <- "C:/DinamicasMoleculares/TRPV1_pipeline/output/08_Max_Flow/08g_v2.1_markov_MD"

dir.create(OUTPUT_DIR,                           recursive = TRUE, showWarnings = FALSE)
dir.create(file.path(OUTPUT_DIR, "per_rep"),     showWarnings = FALSE)
dir.create(file.path(OUTPUT_DIR, "consolidated"),showWarnings = FALSE)
dir.create(file.path(OUTPUT_DIR, "edge_importance"), showWarnings = FALSE)
dir.create(file.path(OUTPUT_DIR, "logs"),        showWarnings = FALSE)

# ==========================================================
# SECTION 3: ANALYSIS PARAMETERS
# ==========================================================

# Target residue numbers (must match 08a settings)
TARGET_FILTER <- 643
TARGET_GATE   <- 679

# Parallel cores (PSOCK cluster, Windows 10 compatible)
N_CORES <- 6

# ==========================================================
# SECTION 4: LOGGING SETUP
# ==========================================================

master_log <- file.path(
  OUTPUT_DIR, "logs",
  paste0("08g_v2.1_master_log_",
         format(Sys.time(), "%Y%m%d_%H%M%S"), ".txt")
)

log_msg <- function(msg, log_file = NULL) {
  txt <- sprintf("[%s] 08g_v2.1 | %s",
                 format(Sys.time(), "%Y-%m-%d %H:%M:%OS3"), msg)
  cat(txt, "\n")
  if (!is.null(log_file)) cat(txt, "\n", file = log_file, append = TRUE)
}

log_msg("=== STARTING 08g_v2.1 MARKOV TRANSITION + FUNDAMENTAL MATRIX ===",
        master_log)
log_msg(paste("INPUT_RDS_DIR:", INPUT_RDS_DIR), master_log)
log_msg(paste("IS_NMA:", IS_NMA,
              "| RUN_FWD:", RUN_FWD, "| RUN_REV:", RUN_REV), master_log)
log_msg(paste("USE_REGULARIZATION:", USE_REGULARIZATION,
              "| EPS:", REGULARIZATION_EPS), master_log)

# ==========================================================
# SECTION 5: LOAD RESIDUE MAPPING
# ==========================================================

log_msg("Loading residue mapping...", master_log)
mapping <- read_csv(MAP_FILE, show_col_types = FALSE)
log_msg(paste("Mapping loaded:", nrow(mapping), "rows"), master_log)

# ==========================================================
# SECTION 6: HELPER FUNCTIONS
# ==========================================================

# --- Parse RDS filename -> system_base, replica, family ---
parse_rds_filename <- function(filename) {
  nm     <- basename(filename) %>% str_remove("\\.rds$")
  family <- case_when(
    str_detect(nm, "_intracellular$") ~ "intracellular",
    str_detect(nm, "_LBD$")          ~ "LBD",
    str_detect(nm, "_combined$")      ~ "combined",
    TRUE                              ~ "unknown"
  )
  sys_part <- nm %>% str_remove("_intracellular$|_LBD$|_combined$")
  if (!IS_NMA) {
    replica     <- str_extract(sys_part, "_rep[0-9]+$")
    system_base <- str_remove(sys_part, "_rep[0-9]+$")
  } else {
    replica     <- NA_character_
    system_base <- sys_part
  }
  list(system_base = system_base, replica = replica,
       family = family, file = filename)
}

get_structure_from_system <- function(system_name) {
  str_split(system_name, "_")[[1]][1]
}

# --- Build transition matrix T from acc_fwd or acc_rev ---
# T_ij = acc[i->j] / sum_j(acc[i->j])
# Rows that sum to zero get uniform distribution (isolated nodes).
# Returns: dense numeric matrix nxn, rownames/colnames = chainresid labels
build_transition_matrix <- function(g, acc_vec, map_df, log_file) {
  
  n  <- vcount(g)
  el <- as_edgelist(g, names = FALSE)
  
  # Build label map: node_id -> chainresid
  label_map <- setNames(
    paste0(map_df$chain, map_df$residue),
    as.character(map_df$node_id)
  )
  
  labels <- label_map[as.character(seq_len(n))]
  
  # Accumulate outgoing flow per node
  T_mat <- matrix(0.0, n, n, dimnames = list(labels, labels))
  
  for (ei in seq_len(nrow(el))) {
    i <- el[ei, 1]
    j <- el[ei, 2]
    T_mat[i, j] <- T_mat[i, j] + acc_vec[ei]
  }
  
  # Normalize rows to get transition probabilities
  row_sums <- rowSums(T_mat)
  for (i in seq_len(n)) {
    if (row_sums[i] > 0) {
      T_mat[i, ] <- T_mat[i, ] / row_sums[i]
    }
    # Rows summing to zero remain as zeros (absorbing or isolated)
  }
  
  T_mat
}

# --- Compute fundamental matrix N = (I - Q)^{-1} ---
# absorbing_indices: integer vector of row/col indices for absorbing states
# (these are removed from Q before inversion)
# Returns: N matrix with rownames/colnames = transient state labels
compute_fundamental_matrix <- function(T_mat, absorbing_indices,
                                       use_reg, eps, log_file) {
  
  n_total     <- nrow(T_mat)
  transient   <- setdiff(seq_len(n_total), absorbing_indices)
  n_transient <- length(transient)
  
  log_msg(paste("  Transient states:", n_transient,
                "| Absorbing states:", length(absorbing_indices)), log_file)
  
  if (n_transient == 0) {
    log_msg("  ERROR: no transient states.", log_file)
    return(NULL)
  }
  
  # Q: submatrix of T restricted to transient states
  Q <- T_mat[transient, transient]
  
  # (I - Q), optionally regularized
  IminusQ <- diag(n_transient) - Q
  if (use_reg) {
    IminusQ <- IminusQ + diag(eps, n_transient)
    log_msg(paste("  Regularization applied: eps =", eps), log_file)
  }
  
  # Inversion
  N_mat <- tryCatch({
    solve(IminusQ)
  }, error = function(e) {
    log_msg(paste("  WARNING: plain inversion failed:", e$message,
                  "-- retrying with regularization"), log_file)
    tryCatch(
      solve(IminusQ + diag(1e-6, n_transient)),
      error = function(e2) {
        log_msg(paste("  ERROR: inversion failed even with regularization:",
                      e2$message), log_file)
        NULL
      }
    )
  })
  
  if (is.null(N_mat)) return(NULL)
  
  rownames(N_mat) <- rownames(T_mat)[transient]
  colnames(N_mat) <- colnames(T_mat)[transient]
  
  log_msg(paste("  N matrix computed:", nrow(N_mat), "x", ncol(N_mat)),
          log_file)
  N_mat
}

# --- Average a list of matrices to a common node space ---
# Expands each matrix to the union of all node labels.
# Nodes absent in a given replica are filled with NA (not zero),
# so mean and stderr are computed only over replicas where the
# node exists (na.rm = TRUE, effective n per cell).
# Returns: list(mean = matrix, stderr = matrix)
average_matrices <- function(mat_list, log_file) {
  
  if (length(mat_list) == 0) return(NULL)
  if (length(mat_list) == 1) {
    m <- mat_list[[1]]
    return(list(mean   = m,
                stderr = matrix(NA_real_, nrow(m), ncol(m),
                                dimnames = dimnames(m))))
  }
  
  all_nodes <- Reduce(union, lapply(mat_list, rownames))
  
  expand_mat <- function(m, nodes) {
    m_exp <- matrix(NA_real_, length(nodes), length(nodes),
                    dimnames = list(nodes, nodes))
    cn <- intersect(nodes, rownames(m))
    m_exp[cn, cn] <- m[cn, cn]
    m_exp
  }
  
  arr <- array(
    unlist(lapply(mat_list, expand_mat, nodes = all_nodes)),
    dim = c(length(all_nodes), length(all_nodes), length(mat_list))
  )
  
  m_mean   <- apply(arr, c(1,2), mean, na.rm = TRUE)
  m_stderr <- apply(arr, c(1,2), function(x) {
    x_valid <- x[!is.na(x)]
    if (length(x_valid) < 2) return(NA_real_)
    sd(x_valid) / sqrt(length(x_valid))
  })
  dimnames(m_mean)   <- list(all_nodes, all_nodes)
  dimnames(m_stderr) <- list(all_nodes, all_nodes)
  
  list(mean = m_mean, stderr = m_stderr)
}

# --- Compute edge importance matrix from N_mean and T_mean ---
# node_importance[i] = sum_j N_mean[j, i]
#   (expected total visits to node i from all transient starting points)
# edge_imp[i,j] = T_mean[i,j] x node_importance[i]
# The result is made symmetric within direction d:
#   E_sym[i,j] = E_sym[j,i] = (E[i,j] + E[j,i]) / 2
# fwd and rev are computed and saved separately -- never averaged together.
# filter and gate are computed and saved separately -- never averaged together.
compute_edge_importance <- function(T_mean, N_mean, log_file) {
  
  # All nodes in T_mean (including absorbing)
  all_labels <- rownames(T_mean)
  n_all      <- length(all_labels)
  
  # node_importance[i] = sum_j N_mean[j, i]
  # Sum over all transient starting points j how many times i is visited.
  # Nodes not present in N_mean (absorbing states) get importance 0.
  node_importance <- rep(0.0, n_all)
  names(node_importance) <- all_labels
  
  n_labels <- rownames(N_mean)
  idx      <- match(n_labels, all_labels)
  valid    <- !is.na(idx)
  
  col_sums <- colSums(N_mean, na.rm = TRUE)
  node_importance[idx[valid]] <- col_sums[valid]
  
  log_msg(paste("  node_importance computed for", sum(node_importance > 0),
                "nodes"), log_file)
  
  # edge_imp[i,j] = T_mean[i,j] x node_importance[i]
  # sweep multiplies row i of T_mean by node_importance[i]
  E_mat <- sweep(T_mean, 1, node_importance, "*")
  
  # Make symmetric within this direction:
  # E_sym[i,j] = E_sym[j,i] = (E[i,j] + E[j,i]) / 2
  E_sym <- (E_mat + t(E_mat)) / 2
  
  log_msg(paste("  Edge importance matrix computed:", nrow(E_sym), "x",
                ncol(E_sym)), log_file)
  E_sym
}

# ==========================================================
# SECTION 7: DISCOVER AND PARSE INPUT FILES
# ==========================================================

rds_files <- list.files(INPUT_RDS_DIR, pattern = "\\.rds$",
                        full.names = TRUE)

if (length(rds_files) == 0) stop(paste("No RDS files found in:", INPUT_RDS_DIR))

log_msg(paste("RDS files found:", length(rds_files)), master_log)

file_meta <- lapply(rds_files, parse_rds_filename)
meta_df   <- bind_rows(lapply(file_meta, as.data.frame,
                              stringsAsFactors = FALSE))

# Work only with intracellular and LBD families
meta_df <- meta_df %>% filter(family %in% c("intracellular", "LBD"))

all_systems <- unique(meta_df$system_base)
log_msg(paste("Unique system bases:", length(all_systems)), master_log)

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
    library(Matrix)
  })
})

clusterExport(cl, c(
  "meta_df", "mapping", "OUTPUT_DIR", "IS_NMA",
  "TARGET_FILTER", "TARGET_GATE",
  "RUN_FWD", "RUN_REV",
  "USE_REGULARIZATION", "REGULARIZATION_EPS",
  "get_structure_from_system",
  "build_transition_matrix",
  "compute_fundamental_matrix",
  "average_matrices",
  "compute_edge_importance",
  "log_msg"
))

results <- parLapply(cl, all_systems, function(sys_base) {
  
  sys_log <- file.path(OUTPUT_DIR, "logs",
                       paste0("08g_", sys_base, "_",
                              format(Sys.time(), "%Y%m%d_%H%M%S"), ".txt"))
  
  log_msg(paste("=== Processing:", sys_base, "==="), sys_log)
  t0 <- Sys.time()
  
  sys_files <- meta_df %>% filter(system_base == sys_base)
  replicas  <- unique(sys_files$replica)
  
  structure_name <- get_structure_from_system(sys_base)
  map_df_sys <- mapping %>%
    filter(system == structure_name) %>%
    arrange(node_id)
  
  # Identify target node labels separately for filter and gate.
  # Each target is used as the sole absorbing state in its own N computation.
  target_labels_filter <- paste0(c("A","B","C","D"), TARGET_FILTER)
  target_labels_gate   <- paste0(c("A","B","C","D"), TARGET_GATE)
  
  # ---- Storage for per-replica T and N matrices ----
  # Separate lists for filter and gate targets
  T_fwd_filter_list <- list()
  T_fwd_gate_list   <- list()
  T_rev_filter_list <- list()
  T_rev_gate_list   <- list()
  N_fwd_filter_list <- list()
  N_fwd_gate_list   <- list()
  N_rev_filter_list <- list()
  N_rev_gate_list   <- list()
  
  for (rep_tag in replicas) {
    
    rep_log_tag <- if (is.na(rep_tag)) "NMA" else rep_tag
    log_msg(paste("  Replica:", rep_log_tag), sys_log)
    
    rep_files <- sys_files %>%
      filter(is.na(replica) & is.na(rep_tag) |
               !is.na(replica) & !is.na(rep_tag) & replica == rep_tag)
    
    # Load and accumulate acc_fwd and acc_rev separately per target family
    g_ref              <- NULL
    acc_fwd_filter     <- NULL
    acc_fwd_gate       <- NULL
    acc_rev_filter     <- NULL
    acc_rev_gate       <- NULL
    
    for (ii in seq_len(nrow(rep_files))) {
      obj <- tryCatch(readRDS(rep_files$file[ii]), error = function(e) NULL)
      if (is.null(obj)) next
      if (is.null(g_ref)) g_ref <- obj$graph
      if (!is.null(obj$acc_fwd_filter)) {
        acc_fwd_filter <- if (is.null(acc_fwd_filter))
          obj$acc_fwd_filter else acc_fwd_filter + obj$acc_fwd_filter
        acc_fwd_gate <- if (is.null(acc_fwd_gate))
          obj$acc_fwd_gate else acc_fwd_gate + obj$acc_fwd_gate
        acc_rev_filter <- if (is.null(acc_rev_filter))
          obj$acc_rev_filter else acc_rev_filter + obj$acc_rev_filter
        acc_rev_gate <- if (is.null(acc_rev_gate))
          obj$acc_rev_gate else acc_rev_gate + obj$acc_rev_gate
      }
    }
    
    if (is.null(g_ref) || is.null(acc_fwd_filter)) {
      log_msg(paste("  WARNING: incomplete data for replica:", rep_log_tag,
                    "-- skipping"), sys_log)
      next
    }
    
    log_msg(paste("  Graph:", vcount(g_ref), "nodes,",
                  ecount(g_ref), "edges"), sys_log)
    
    rep_suffix <- if (is.na(rep_tag)) "" else paste0("_", rep_tag)
    
    # ---- FWD direction ----
    if (RUN_FWD) {
      
      # -- filter target --
      log_msg("  Building T_fwd_filter...", sys_log)
      T_fwd_filter <- build_transition_matrix(g_ref, acc_fwd_filter,
                                              map_df_sys, sys_log)
      
      saveRDS(T_fwd_filter,
              file.path(OUTPUT_DIR, "per_rep",
                        paste0("T_fwd_filter_", sys_base, rep_suffix, ".rds")))
      
      T_fwd_filter_list[[rep_log_tag]] <- T_fwd_filter
      
      absorbing_fwd_filter <- which(rownames(T_fwd_filter) %in% target_labels_filter)
      
      log_msg("  Computing N_fwd_filter...", sys_log)
      N_fwd_filter <- compute_fundamental_matrix(
        T_fwd_filter, absorbing_fwd_filter,
        USE_REGULARIZATION, REGULARIZATION_EPS, sys_log
      )
      
      if (!is.null(N_fwd_filter)) {
        saveRDS(N_fwd_filter,
                file.path(OUTPUT_DIR, "per_rep",
                          paste0("N_fwd_filter_", sys_base, rep_suffix, ".rds")))
        N_fwd_filter_list[[rep_log_tag]] <- N_fwd_filter
      }
      
      # -- gate target --
      log_msg("  Building T_fwd_gate...", sys_log)
      T_fwd_gate <- build_transition_matrix(g_ref, acc_fwd_gate,
                                            map_df_sys, sys_log)
      
      saveRDS(T_fwd_gate,
              file.path(OUTPUT_DIR, "per_rep",
                        paste0("T_fwd_gate_", sys_base, rep_suffix, ".rds")))
      
      T_fwd_gate_list[[rep_log_tag]] <- T_fwd_gate
      
      absorbing_fwd_gate <- which(rownames(T_fwd_gate) %in% target_labels_gate)
      
      log_msg("  Computing N_fwd_gate...", sys_log)
      N_fwd_gate <- compute_fundamental_matrix(
        T_fwd_gate, absorbing_fwd_gate,
        USE_REGULARIZATION, REGULARIZATION_EPS, sys_log
      )
      
      if (!is.null(N_fwd_gate)) {
        saveRDS(N_fwd_gate,
                file.path(OUTPUT_DIR, "per_rep",
                          paste0("N_fwd_gate_", sys_base, rep_suffix, ".rds")))
        N_fwd_gate_list[[rep_log_tag]] <- N_fwd_gate
      }
    }
    
    # ---- REV direction ----
    if (RUN_REV) {
      
      # -- filter target --
      log_msg("  Building T_rev_filter...", sys_log)
      T_rev_filter <- build_transition_matrix(g_ref, acc_rev_filter,
                                              map_df_sys, sys_log)
      
      saveRDS(T_rev_filter,
              file.path(OUTPUT_DIR, "per_rep",
                        paste0("T_rev_filter_", sys_base, rep_suffix, ".rds")))
      
      T_rev_filter_list[[rep_log_tag]] <- T_rev_filter
      
      absorbing_rev_filter <- which(rownames(T_rev_filter) %in% target_labels_filter)
      
      log_msg("  Computing N_rev_filter...", sys_log)
      N_rev_filter <- compute_fundamental_matrix(
        T_rev_filter, absorbing_rev_filter,
        USE_REGULARIZATION, REGULARIZATION_EPS, sys_log
      )
      
      if (!is.null(N_rev_filter)) {
        saveRDS(N_rev_filter,
                file.path(OUTPUT_DIR, "per_rep",
                          paste0("N_rev_filter_", sys_base, rep_suffix, ".rds")))
        N_rev_filter_list[[rep_log_tag]] <- N_rev_filter
      }
      
      # -- gate target --
      log_msg("  Building T_rev_gate...", sys_log)
      T_rev_gate <- build_transition_matrix(g_ref, acc_rev_gate,
                                            map_df_sys, sys_log)
      
      saveRDS(T_rev_gate,
              file.path(OUTPUT_DIR, "per_rep",
                        paste0("T_rev_gate_", sys_base, rep_suffix, ".rds")))
      
      T_rev_gate_list[[rep_log_tag]] <- T_rev_gate
      
      absorbing_rev_gate <- which(rownames(T_rev_gate) %in% target_labels_gate)
      
      log_msg("  Computing N_rev_gate...", sys_log)
      N_rev_gate <- compute_fundamental_matrix(
        T_rev_gate, absorbing_rev_gate,
        USE_REGULARIZATION, REGULARIZATION_EPS, sys_log
      )
      
      if (!is.null(N_rev_gate)) {
        saveRDS(N_rev_gate,
                file.path(OUTPUT_DIR, "per_rep",
                          paste0("N_rev_gate_", sys_base, rep_suffix, ".rds")))
        N_rev_gate_list[[rep_log_tag]] <- N_rev_gate
      }
    }
  }
  
  # ---- Consolidate N and T across replicas (mean + stderr) ----
  # Both T and N are averaged across replicas so that edge importance
  # uses T_mean (not a single replica T), consistent with N_mean.
  
  consolidate_and_save <- function(mat_list, direction, target, mat_type) {
    
    if (length(mat_list) == 0) return(NULL)
    
    consolidated <- average_matrices(mat_list, sys_log)
    if (is.null(consolidated)) return(NULL)
    
    m_mean   <- consolidated$mean
    m_stderr <- consolidated$stderr
    
    saveRDS(m_mean,
            file.path(OUTPUT_DIR, "consolidated",
                      paste0(mat_type, "_", direction, "_", target, "_",
                             sys_base, "_mean.rds")))
    write.csv(m_mean,
              file.path(OUTPUT_DIR, "consolidated",
                        paste0(mat_type, "_", direction, "_", target, "_",
                               sys_base, "_mean.csv")),
              row.names = TRUE)
    saveRDS(m_stderr,
            file.path(OUTPUT_DIR, "consolidated",
                      paste0(mat_type, "_", direction, "_", target, "_",
                             sys_base, "_stderr.rds")))
    write.csv(m_stderr,
              file.path(OUTPUT_DIR, "consolidated",
                        paste0(mat_type, "_", direction, "_", target, "_",
                               sys_base, "_stderr.csv")),
              row.names = TRUE)
    
    log_msg(paste(" ", mat_type, direction, target, "consolidated saved:",
                  nrow(m_mean), "x", ncol(m_mean)), sys_log)
    
    consolidated
  }
  
  T_fwd_filter_cons <- if (RUN_FWD) consolidate_and_save(T_fwd_filter_list, "fwd", "filter", "T") else NULL
  T_fwd_gate_cons   <- if (RUN_FWD) consolidate_and_save(T_fwd_gate_list,   "fwd", "gate",   "T") else NULL
  T_rev_filter_cons <- if (RUN_REV) consolidate_and_save(T_rev_filter_list, "rev", "filter", "T") else NULL
  T_rev_gate_cons   <- if (RUN_REV) consolidate_and_save(T_rev_gate_list,   "rev", "gate",   "T") else NULL
  N_fwd_filter_cons <- if (RUN_FWD) consolidate_and_save(N_fwd_filter_list, "fwd", "filter", "N") else NULL
  N_fwd_gate_cons   <- if (RUN_FWD) consolidate_and_save(N_fwd_gate_list,   "fwd", "gate",   "N") else NULL
  N_rev_filter_cons <- if (RUN_REV) consolidate_and_save(N_rev_filter_list, "rev", "filter", "N") else NULL
  N_rev_gate_cons   <- if (RUN_REV) consolidate_and_save(N_rev_gate_list,   "rev", "gate",   "N") else NULL
  
  # ---- Edge importance matrices for 08d ----
  # Uses T_mean and N_mean for each direction and target independently.
  # fwd_filter: symmetric within fwd direction, filter target only
  # fwd_gate:   symmetric within fwd direction, gate target only
  # rev_filter: symmetric within rev direction, filter target only
  # rev_gate:   symmetric within rev direction, gate target only
  # Directions never averaged together; targets never averaged together.
  
  compute_and_save_edge_imp <- function(T_consolidated, N_consolidated,
                                        direction, target) {
    
    if (is.null(T_consolidated) || is.null(N_consolidated)) return(NULL)
    
    T_mean <- T_consolidated$mean
    N_mean <- N_consolidated$mean
    
    log_msg(paste("  Computing edge importance (", direction, target, ")..."),
            sys_log)
    E_mat <- compute_edge_importance(T_mean, N_mean, sys_log)
    
    if (is.null(E_mat)) return(NULL)
    
    # Save as CSV compatible with 08d (row.names = TRUE, symmetric matrix)
    out_csv <- file.path(OUTPUT_DIR, "edge_importance",
                         paste0("matrix_", sys_base, "_", direction,
                                "_", target, ".csv"))
    write.csv(E_mat, out_csv, row.names = TRUE)
    log_msg(paste("  Edge importance saved:", basename(out_csv)), sys_log)
    
    E_mat
  }
  
  if (RUN_FWD) {
    compute_and_save_edge_imp(T_fwd_filter_cons, N_fwd_filter_cons, "fwd", "filter")
    compute_and_save_edge_imp(T_fwd_gate_cons,   N_fwd_gate_cons,   "fwd", "gate")
  }
  if (RUN_REV) {
    compute_and_save_edge_imp(T_rev_filter_cons, N_rev_filter_cons, "rev", "filter")
    compute_and_save_edge_imp(T_rev_gate_cons,   N_rev_gate_cons,   "rev", "gate")
  }
  
  elapsed <- round(as.numeric(difftime(Sys.time(), t0, units = "secs")), 1)
  log_msg(paste("=== Finished:", sys_base, "| elapsed:", elapsed, "s ==="),
          sys_log)
  
  list(system        = sys_base,
       n_reps_fwd_filter = length(N_fwd_filter_list),
       n_reps_fwd_gate   = length(N_fwd_gate_list),
       n_reps_rev_filter = length(N_rev_filter_list),
       n_reps_rev_gate   = length(N_rev_gate_list))
})

stopCluster(cl)

# ==========================================================
# SECTION 9: GLOBAL SUMMARY
# ==========================================================

completed <- Filter(Negate(is.null), results)
log_msg(paste("Completed:", length(completed), "systems"), master_log)
for (r in completed) {
  log_msg(paste(" ", r$system,
                "| fwd_filter replicas:", r$n_reps_fwd_filter,
                "| fwd_gate replicas:",   r$n_reps_fwd_gate,
                "| rev_filter replicas:", r$n_reps_rev_filter,
                "| rev_gate replicas:",   r$n_reps_rev_gate), master_log)
}
log_msg("=== 08g_v2.1 FINISHED ===", master_log)

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
    run_type <- if (IS_NMA) "NMA" else "MD"
    
    new_entries <- data.frame(
      script_number  = "08g",
      script_version = "2.1",
      script_name    = "08g_v2.1_markov_transition_fundamental",
      run_type       = run_type,
      variable = c("INPUT_RDS_DIR", "MAP_FILE", "OUTPUT_DIR",
                   "IS_NMA", "run_type",
                   "TARGET_FILTER", "TARGET_GATE",
                   "RUN_FWD", "RUN_REV",
                   "USE_REGULARIZATION", "REGULARIZATION_EPS",
                   "N_CORES",
                   "systems_completed"),
      value = c(INPUT_RDS_DIR, MAP_FILE, OUTPUT_DIR,
                as.character(IS_NMA), run_type,
                as.character(TARGET_FILTER), as.character(TARGET_GATE),
                as.character(RUN_FWD), as.character(RUN_REV),
                as.character(USE_REGULARIZATION),
                as.character(REGULARIZATION_EPS),
                as.character(N_CORES),
                as.character(length(completed))),
      timestamp = ts,
      stringsAsFactors = FALSE
    )
    
    if (file.exists(final_path)) {
      existing <- read.csv(final_path, stringsAsFactors = FALSE,
                           colClasses = "character")
      existing <- existing[!(existing$script_number == "08g" &
                               existing$run_type == run_type), ]
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
