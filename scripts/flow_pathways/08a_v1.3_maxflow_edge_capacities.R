# ==========================================================
# 08a_v1.3_maxflow_allosteric_edges.R
# Author: DenyCB
# PURPOSE:
#   Compute allosteric edge weights using network max-flow
#   between source seed residues and target residues (selectivity
#   filter and gate). Two independent seed families are used:
#
#   FAMILY 1 - Intracellular seeds:
#     Deep intracellular residues (N-terminal ankyrin region
#     and C-terminal TRP box / post-TRP region). These are
#     position-agnostic seeds designed to capture long-range
#     allosteric paths from the cytoplasmic domain.
#
#   FAMILY 2 - LBD seeds:
#     Capsaicin binding site residues and extracellular residues.
#     These are functionally motivated seeds targeting the
#     vanilloid binding pocket and upper channel region.
#
#   For each seed family, max-flow is computed for every
#   (source_node, target_node) pair and the flows are
#   accumulated across all pairs. The accumulated flow on
#   each edge represents how much allosteric capacity that
#   edge carries across all source-target combinations.
#
#   Outputs per system:
#     - Edge flow CSV (Gephi-compatible) for each family
#     - Edge flow CSV for combined (sum) of both families
#     - Flow matrix CSV (residue x residue) for each family
#     - RDS objects for downstream use in 08b (path decomposition)
#     - Log file with timing and statistics
#
# INPUTS:
#   - Network RDS files (igraph objects from 04a)
#   - Residue mapping CSV (node_id -> chain, residue, system)
#   - PDB files (for coordinate-based filters if needed)
#
# NOTES:
#   - max_flow() in igraph requires a single (source, target) pair.
#     To handle multiple chains, all nodes matching a target residue
#     number across all chains are used as individual targets,
#     and flows are accumulated.
#   - Seeds outside the resolved range of a given structure are
#     silently skipped (no error, warning logged).
#   - Compatible with Windows 10 via PSOCK parallel cluster.
#   - NMA networks (consolidated, no replicas) and MD networks
#     (per-replica or repMean) are handled by the same code;
#     select the appropriate NETWORK_DIR.
#   - accumulate_flow() now returns three vectors per family/target:
#     acc_fwd (flow in from->to edgelist direction),
#     acc_rev (flow in to->from direction),
#     acc_abs (total absolute flow, direction-agnostic).
#     All three are stored in RDS for downstream use.
#
# ==========================================================

library(igraph)
library(dplyr)
library(readr)
library(stringr)
library(parallel)

# ==========================================================
# SECTION 1: ANALYSIS SWITCHES
# ==========================================================

# Set to TRUE to run each seed family
RUN_INTRACELLULAR <- TRUE
RUN_LBD           <- TRUE

# Set to TRUE to also save the combined (summed) flow
RUN_COMBINED      <- TRUE

# ==========================================================
# SECTION 2: INPUT/OUTPUT CONFIGURATION
# ==========================================================

# --- NMA networks (uncomment to use) ---
#NETWORK_DIR <- "C:/DinamicasMoleculares/TRPV1_pipeline/output/04_Struct_Network/A_network_objects_NMA/Th0_370"
#MAP_FILE    <- "C:/DinamicasMoleculares/TRPV1_pipeline/config_files/trimmed_residue_map.csv"
#OUTPUT_DIR  <- "C:/DinamicasMoleculares/TRPV1_pipeline/output/08_Max_Flow/08a_v1.3_allosteric_edges_NMA"
#IS_NMA      <- TRUE   # TRUE = NMA (no replicas), FALSE = MD

# --- MD networks (comment out NMA block above and uncomment below) ---
 NETWORK_DIR <- "C:/DinamicasMoleculares/TRPV1_pipeline/output/04_Struct_Network/A_network_objects_md/Th0.49"
 MAP_FILE    <- "C:/DinamicasMoleculares/TRPV1_pipeline/output/residue_map_MD.csv"
 OUTPUT_DIR  <- "C:/DinamicasMoleculares/TRPV1_pipeline/output/08_Max_Flow/08a_v1.3_allosteric_edges_MD"
 IS_NMA      <- FALSE

# PDB directory (used for structure name matching)
PDB_DIR <- "C:/DinamicasMoleculares/TRPV1_pipeline/input/pdb/"

dir.create(OUTPUT_DIR, showWarnings = FALSE, recursive = TRUE)
dir.create(file.path(OUTPUT_DIR, "rds"),      showWarnings = FALSE)
dir.create(file.path(OUTPUT_DIR, "logs"),     showWarnings = FALSE)
dir.create(file.path(OUTPUT_DIR, "directed"), showWarnings = FALSE)

# ==========================================================
# SECTION 3: SEED DEFINITIONS
# ==========================================================

# --- FAMILY 1: Intracellular deep seeds ---
SEEDS_INTRACELLULAR_RANGES <- list(
  c(227, 383),   # N-terminal ankyrin repeat domain
  c(720, 766)    # C-terminal post-TRP intracellular region
)
SEEDS_INTRACELLULAR_SPECIFIC <- c()  # additional specific residues if needed

# --- FAMILY 2: LBD and extracellular seeds ---
SEEDS_LBD_RANGES <- list(
  c(533, 536),   # S4-S5 linker (cap binding)
  c(457, 467)    # S2-S3 loop (cap binding)
)
SEEDS_LBD_SPECIFIC <- c(
  510, 511, 512,  # S3 cap binding residues
  550, 557, 570   # S4 and linker cap binding residues
)

# ==========================================================
# SECTION 4: TARGET RESIDUES
# ==========================================================

TARGET_FILTER <- 643
TARGET_GATE   <- 679

# ==========================================================
# SECTION 5: ANALYSIS PARAMETERS
# ==========================================================

MIN_FLOW_THRESHOLD <- 1e-6
N_CORES <- 6

# ==========================================================
# SECTION 6: LOGGING SETUP
# ==========================================================

master_log <- file.path(
  OUTPUT_DIR, "logs",
  paste0("08a_master_log_", format(Sys.time(), "%Y%m%d_%H%M%S"), ".txt")
)

log_msg <- function(msg, log_file = NULL) {
  txt <- sprintf("[%s] %s", format(Sys.time(), "%Y-%m-%d %H:%M:%OS3"), msg)
  cat(txt, "\n")
  if (!is.null(log_file)) cat(txt, "\n", file = log_file, append = TRUE)
}

log_msg("=== STARTING 08a_v1.3 MAX-FLOW ALLOSTERIC EDGE ANALYSIS ===", master_log)
log_msg(paste("NETWORK_DIR:", NETWORK_DIR), master_log)
log_msg(paste("MAP_FILE:", MAP_FILE), master_log)
log_msg(paste("OUTPUT_DIR:", OUTPUT_DIR), master_log)
log_msg(paste("IS_NMA:", IS_NMA), master_log)
log_msg(paste("RUN_INTRACELLULAR:", RUN_INTRACELLULAR), master_log)
log_msg(paste("RUN_LBD:", RUN_LBD), master_log)
log_msg(paste("TARGET_FILTER:", TARGET_FILTER, "| TARGET_GATE:", TARGET_GATE), master_log)

# ==========================================================
# SECTION 7: LOAD RESIDUE MAPPING
# ==========================================================

log_msg("Loading residue mapping...", master_log)
mapping <- read_csv(MAP_FILE, show_col_types = FALSE)
log_msg(paste("Mapping loaded:", nrow(mapping), "rows"), master_log)

# ==========================================================
# SECTION 8: HELPER FUNCTIONS
# ==========================================================

get_full_name <- function(file) {
  basename(file) %>%
    str_remove("network_NMA_|network_MD_") %>%
    str_remove("\\.rds$")
}

get_structure_name <- function(full_name) {
  str_split(full_name, "_")[[1]][1]
}

expand_seeds <- function(ranges, specific) {
  from_ranges <- unlist(lapply(ranges, function(r) r[1]:r[2]))
  sort(unique(c(from_ranges, specific)))
}

get_nodes_for_residue <- function(map_df, residue_number) {
  map_df %>%
    filter(residue == residue_number) %>%
    pull(node_id) %>%
    as.integer()
}

extract_flow_edges <- function(g, flow_vec, min_flow) {
  el <- as_edgelist(g, names = FALSE)
  df <- data.frame(
    from_node = el[, 1],
    to_node   = el[, 2],
    flow      = as.numeric(flow_vec),
    stringsAsFactors = FALSE
  )
  df <- df[df$flow > min_flow, ]
  df
}

# --- Accumulate flow across multiple (source, target) pairs ---
# For each pair, runs max_flow() and accumulates into THREE vectors:
#   acc_fwd: sum of positive flow values (flow in from->to edgelist direction)
#   acc_rev: sum of abs(negative flow values) (flow in to->from direction)
#   acc_abs: sum of abs(all flow values) (total activity, direction-agnostic)
# The distinction matters when different source-target pairs send flow
# in opposite directions through the same edge. acc_fwd and acc_rev
# preserve this directionality; acc_abs discards it.
# Returns a named list with three numeric vectors of length ecount(g).
accumulate_flow <- function(g, source_nodes, target_nodes,
                             min_flow, log_file, system_id) {

  acc_fwd <- rep(0.0, ecount(g))
  acc_rev <- rep(0.0, ecount(g))
  acc_abs <- rep(0.0, ecount(g))
  n_pairs   <- 0L
  n_skipped <- 0L

  for (s in source_nodes) {

    if (s > vcount(g) || s < 1) {
      n_skipped <- n_skipped + 1L
      next
    }

    for (t in target_nodes) {

      if (t > vcount(g) || t < 1) {
        n_skipped <- n_skipped + 1L
        next
      }

      if (s == t) next

      result <- tryCatch({
        max_flow(g,
                 source = V(g)[s],
                 target = V(g)[t],
                 capacity = E(g)$weight)
      }, error = function(e) {
        log_msg(paste("  max_flow error s=", s, "t=", t, ":", e$message),
                log_file)
        NULL
      })

      if (is.null(result)) next

      flow_signed <- as.numeric(result$flow)

      # Positive values: flow goes in from->to direction of edgelist
      acc_fwd <- acc_fwd + pmax(flow_signed, 0)

      # Negative values: flow goes in to->from direction; store as positive
      acc_rev <- acc_rev + pmax(-flow_signed, 0)

      # Absolute value: total activity regardless of direction
      acc_abs <- acc_abs + abs(flow_signed)

      n_pairs <- n_pairs + 1L
    }
  }

  log_msg(paste("  Pairs processed:", n_pairs,
                "| Skipped (out of range):", n_skipped), log_file)

  list(acc_fwd = acc_fwd, acc_rev = acc_rev, acc_abs = acc_abs)
}

# --- Convert accumulated flow vector to edge CSV ---
# Uses acc_abs (direction-agnostic) for the non-directed output.
# Source and Target columns use chain+residue labels (e.g. A300, C426).
# Output format compatible with 05d matrix visualization script:
#   Source, Target, type, weight (normalized 0-1), weight_abs (raw flow)
flow_to_gephi <- function(g, acc_flow, map_df, min_flow, type_label) {

  el <- as_edgelist(g, names = FALSE)

  label_map <- setNames(
    paste0(map_df$chain, map_df$residue),
    map_df$node_id
  )

  df <- data.frame(
    Source     = label_map[as.character(el[, 1])],
    Target     = label_map[as.character(el[, 2])],
    type       = type_label,
    weight_abs = acc_flow,
    stringsAsFactors = FALSE
  )

  df <- df[df$weight_abs > min_flow & !is.na(df$Source) & !is.na(df$Target), ]

  max_flow_val <- max(df$weight_abs, na.rm = TRUE)
  df$weight <- if (max_flow_val > 0) df$weight_abs / max_flow_val else df$weight_abs

  df <- df[, c("Source", "Target", "type", "weight", "weight_abs")]

  df
}

# --- Convert directed flow vectors to directed edge CSV ---
# Uses real acc_fwd and acc_rev from accumulate_flow() — not abs().
# acc_fwd_filter: flow in from->to direction for filter target
# acc_rev_filter: flow in to->from direction for filter target
# acc_fwd_gate:   flow in from->to direction for gate target
# acc_rev_gate:   flow in to->from direction for gate target
# Each target produces two rows per edge (fwd and rev), saved in directed/.
flow_to_gephi_directed <- function(g, acc_fwd_filter, acc_rev_filter,
                                    acc_fwd_gate, acc_rev_gate,
                                    map_df, min_flow, type_label) {

  el <- as_edgelist(g, names = FALSE)

  label_map <- setNames(
    paste0(map_df$chain, map_df$residue),
    map_df$node_id
  )

  make_directed <- function(acc_fwd, acc_rev, target_label) {

    df_fwd <- data.frame(
      Source     = label_map[as.character(el[, 1])],
      Target     = label_map[as.character(el[, 2])],
      type       = paste0(type_label, "_to_", target_label, "_fwd"),
      weight_abs = acc_fwd,
      stringsAsFactors = FALSE
    )

    df_rev <- data.frame(
      Source     = label_map[as.character(el[, 2])],
      Target     = label_map[as.character(el[, 1])],
      type       = paste0(type_label, "_to_", target_label, "_rev"),
      weight_abs = acc_rev,
      stringsAsFactors = FALSE
    )

    df <- rbind(df_fwd, df_rev)
    df <- df[df$weight_abs > min_flow & !is.na(df$Source) & !is.na(df$Target), ]

    max_val   <- max(df$weight_abs, na.rm = TRUE)
    df$weight <- if (max_val > 0) df$weight_abs / max_val else df$weight_abs

    df <- df[, c("Source", "Target", "type", "weight", "weight_abs")]
    df
  }

  df_filter <- make_directed(acc_fwd_filter, acc_rev_filter, "filter")
  df_gate   <- make_directed(acc_fwd_gate,   acc_rev_gate,   "gate")

  rbind(df_filter, df_gate)
}

# --- Build symmetric flow matrix (chain+residue x chain+residue) ---
flow_to_matrix <- function(gephi_df, map_df) {

  df <- gephi_df %>%
    filter(!is.na(Source), !is.na(Target)) %>%
    group_by(Source, Target) %>%
    summarise(flow = sum(weight_abs), .groups = "drop")

  all_labels <- sort(unique(c(df$Source, df$Target)))
  mat        <- matrix(0, length(all_labels), length(all_labels),
                       dimnames = list(all_labels, all_labels))

  for (i in seq_len(nrow(df))) {
    r <- df$Source[i]
    c <- df$Target[i]
    mat[r, c] <- mat[r, c] + df$flow[i]
    mat[c, r] <- mat[c, r] + df$flow[i]
  }

  mat
}

# ==========================================================
# SECTION 9: FIND NETWORK FILES
# ==========================================================

network_files <- list.files(NETWORK_DIR,
                             pattern = "\\.rds$",
                             full.names = TRUE)

if (length(network_files) == 0) {
  stop(paste("No network RDS files found in:", NETWORK_DIR))
}

log_msg(paste("Network files found:", length(network_files)), master_log)

# ==========================================================
# SECTION 10: EXPAND SEEDS
# ==========================================================

seeds_intra <- expand_seeds(SEEDS_INTRACELLULAR_RANGES,
                             SEEDS_INTRACELLULAR_SPECIFIC)
seeds_lbd   <- expand_seeds(SEEDS_LBD_RANGES,
                             SEEDS_LBD_SPECIFIC)

log_msg(paste("Intracellular seeds (residue numbers):", length(seeds_intra)), master_log)
log_msg(paste("LBD seeds (residue numbers):", length(seeds_lbd)), master_log)

# ==========================================================
# SECTION 11: PARALLEL PROCESSING
# ==========================================================

log_msg(paste("Starting parallel cluster with", N_CORES, "cores (PSOCK)..."), master_log)

cl <- makeCluster(N_CORES)

clusterEvalQ(cl, {
  library(igraph)
  library(dplyr)
  library(readr)
  library(stringr)
})

clusterExport(cl, c(
  "mapping", "OUTPUT_DIR", "IS_NMA",
  "seeds_intra", "seeds_lbd",
  "TARGET_FILTER", "TARGET_GATE",
  "MIN_FLOW_THRESHOLD",
  "RUN_INTRACELLULAR", "RUN_LBD", "RUN_COMBINED",
  "get_full_name", "get_structure_name", "get_nodes_for_residue",
  "accumulate_flow", "extract_flow_edges",
  "flow_to_gephi", "flow_to_gephi_directed", "flow_to_matrix",
  "log_msg"
))

results <- parLapply(cl, network_files, function(net_file) {

  full_name  <- get_full_name(net_file)
  sys_log    <- file.path(OUTPUT_DIR, "logs",
                           paste0("08a_", full_name, "_",
                                  format(Sys.time(), "%Y%m%d_%H%M%S"), ".txt"))

  log_msg(paste("=== Processing:", full_name, "==="), sys_log)
  t0 <- Sys.time()

  obj <- tryCatch(readRDS(net_file), error = function(e) {
    log_msg(paste("ERROR loading:", net_file, "-", e$message), sys_log)
    return(NULL)
  })
  if (is.null(obj)) return(NULL)

  g <- obj$graph

  if (is.null(E(g)$weight) || length(E(g)$weight) == 0) {
    log_msg("ERROR: graph has no edge weights. Skipping.", sys_log)
    return(NULL)
  }

  log_msg(paste("Graph loaded | nodes:", vcount(g),
                "| edges:", ecount(g)), sys_log)

  structure_name <- get_structure_name(full_name)

  if (IS_NMA) {
    map_df <- mapping %>%
      filter(system == structure_name) %>%
      arrange(node_id)
  } else {
    map_df <- mapping %>%
      filter(system == structure_name) %>%
      arrange(node_id)
  }

  if (nrow(map_df) == 0) {
    log_msg(paste("WARNING: no mapping found for structure", structure_name,
                  "- skipping."), sys_log)
    return(NULL)
  }

  log_msg(paste("Mapping rows for", structure_name, ":", nrow(map_df)), sys_log)

  resolve_seeds <- function(seed_residues) {
    node_ids <- unlist(lapply(seed_residues, function(r) {
      get_nodes_for_residue(map_df, r)
    }))
    node_ids <- node_ids[node_ids >= 1 & node_ids <= vcount(g)]
    sort(unique(as.integer(node_ids)))
  }

  target_filter_nodes <- resolve_seeds(TARGET_FILTER)
  target_gate_nodes   <- resolve_seeds(TARGET_GATE)

  log_msg(paste("Target FILTER nodes (", TARGET_FILTER, "):",
                length(target_filter_nodes)), sys_log)
  log_msg(paste("Target GATE nodes (", TARGET_GATE, "):",
                length(target_gate_nodes)), sys_log)

  if (length(target_filter_nodes) == 0 && length(target_gate_nodes) == 0) {
    log_msg("WARNING: no target nodes found. Skipping.", sys_log)
    return(NULL)
  }

  family_results <- list()

  # ==========================================================
  # BLOCK A: INTRACELLULAR SEEDS
  # ==========================================================

  if (RUN_INTRACELLULAR) {

    log_msg("--- FAMILY 1: Intracellular seeds ---", sys_log)

    src_intra <- resolve_seeds(seeds_intra)
    log_msg(paste("Intracellular source nodes resolved:", length(src_intra)),
            sys_log)

    if (length(src_intra) == 0) {
      log_msg("WARNING: no intracellular source nodes in this structure.",
              sys_log)
    } else {

      # accumulate_flow returns list(acc_fwd, acc_rev, acc_abs)
      res_intra_filter <- list(acc_fwd = rep(0.0, ecount(g)),
                               acc_rev = rep(0.0, ecount(g)),
                               acc_abs = rep(0.0, ecount(g)))
      if (length(target_filter_nodes) > 0) {
        log_msg(paste("  Computing flow: intra ->", TARGET_FILTER), sys_log)
        res_intra_filter <- accumulate_flow(
          g, src_intra, target_filter_nodes,
          MIN_FLOW_THRESHOLD, sys_log, full_name
        )
      }

      res_intra_gate <- list(acc_fwd = rep(0.0, ecount(g)),
                             acc_rev = rep(0.0, ecount(g)),
                             acc_abs = rep(0.0, ecount(g)))
      if (length(target_gate_nodes) > 0) {
        log_msg(paste("  Computing flow: intra ->", TARGET_GATE), sys_log)
        res_intra_gate <- accumulate_flow(
          g, src_intra, target_gate_nodes,
          MIN_FLOW_THRESHOLD, sys_log, full_name
        )
      }

      # Use acc_abs for non-directed outputs
      acc_intra_filter <- res_intra_filter$acc_abs
      acc_intra_gate   <- res_intra_gate$acc_abs
      acc_intra_total  <- acc_intra_filter + acc_intra_gate

      gephi_intra <- flow_to_gephi(g, acc_intra_total, map_df,
                                    MIN_FLOW_THRESHOLD, "intracellular")
      mat_intra   <- flow_to_matrix(gephi_intra, map_df)

      tag <- paste0(full_name, "_intracellular")

      write_csv(gephi_intra,
                file.path(OUTPUT_DIR, paste0(tag, "_edges_gephi.csv")))
      write_csv(as.data.frame(mat_intra),
                file.path(OUTPUT_DIR, paste0(tag, "_flow_matrix.csv")))
      saveRDS(list(acc_filter     = acc_intra_filter,
                   acc_gate       = acc_intra_gate,
                   acc_total      = acc_intra_total,
                   acc_fwd_filter = res_intra_filter$acc_fwd,
                   acc_rev_filter = res_intra_filter$acc_rev,
                   acc_fwd_gate   = res_intra_gate$acc_fwd,
                   acc_rev_gate   = res_intra_gate$acc_rev,
                   gephi          = gephi_intra,
                   matrix         = mat_intra,
                   graph          = g,
                   map_df         = map_df),
              file.path(OUTPUT_DIR, "rds", paste0(tag, ".rds")))

      log_msg(paste("  Intracellular edges saved:",
                    nrow(gephi_intra), "above threshold"), sys_log)

      # Directed version uses real acc_fwd and acc_rev — not acc_abs
      gephi_intra_dir <- flow_to_gephi_directed(
        g,
        res_intra_filter$acc_fwd, res_intra_filter$acc_rev,
        res_intra_gate$acc_fwd,   res_intra_gate$acc_rev,
        map_df, MIN_FLOW_THRESHOLD, "intracellular"
      )
      write_csv(gephi_intra_dir,
                file.path(OUTPUT_DIR, "directed",
                          paste0(full_name, "_intracellular_edges_directed.csv")))
      log_msg(paste("  Intracellular directed edges saved:",
                    nrow(gephi_intra_dir), "rows"), sys_log)

      family_results$intracellular        <- acc_intra_total
      family_results$intracellular_filter <- acc_intra_filter
      family_results$intracellular_gate   <- acc_intra_gate
      family_results$intra_fwd_filter     <- res_intra_filter$acc_fwd
      family_results$intra_rev_filter     <- res_intra_filter$acc_rev
      family_results$intra_fwd_gate       <- res_intra_gate$acc_fwd
      family_results$intra_rev_gate       <- res_intra_gate$acc_rev
    }
  }

  # ==========================================================
  # BLOCK B: LBD SEEDS
  # ==========================================================

  if (RUN_LBD) {

    log_msg("--- FAMILY 2: LBD seeds ---", sys_log)

    src_lbd <- resolve_seeds(seeds_lbd)
    log_msg(paste("LBD source nodes resolved:", length(src_lbd)), sys_log)

    if (length(src_lbd) == 0) {
      log_msg("WARNING: no LBD source nodes in this structure.", sys_log)
    } else {

      res_lbd_filter <- list(acc_fwd = rep(0.0, ecount(g)),
                             acc_rev = rep(0.0, ecount(g)),
                             acc_abs = rep(0.0, ecount(g)))
      if (length(target_filter_nodes) > 0) {
        log_msg(paste("  Computing flow: LBD ->", TARGET_FILTER), sys_log)
        res_lbd_filter <- accumulate_flow(
          g, src_lbd, target_filter_nodes,
          MIN_FLOW_THRESHOLD, sys_log, full_name
        )
      }

      res_lbd_gate <- list(acc_fwd = rep(0.0, ecount(g)),
                           acc_rev = rep(0.0, ecount(g)),
                           acc_abs = rep(0.0, ecount(g)))
      if (length(target_gate_nodes) > 0) {
        log_msg(paste("  Computing flow: LBD ->", TARGET_GATE), sys_log)
        res_lbd_gate <- accumulate_flow(
          g, src_lbd, target_gate_nodes,
          MIN_FLOW_THRESHOLD, sys_log, full_name
        )
      }

      acc_lbd_filter <- res_lbd_filter$acc_abs
      acc_lbd_gate   <- res_lbd_gate$acc_abs
      acc_lbd_total  <- acc_lbd_filter + acc_lbd_gate

      gephi_lbd <- flow_to_gephi(g, acc_lbd_total, map_df,
                                   MIN_FLOW_THRESHOLD, "LBD")
      mat_lbd   <- flow_to_matrix(gephi_lbd, map_df)

      tag <- paste0(full_name, "_LBD")

      write_csv(gephi_lbd,
                file.path(OUTPUT_DIR, paste0(tag, "_edges_gephi.csv")))
      write_csv(as.data.frame(mat_lbd),
                file.path(OUTPUT_DIR, paste0(tag, "_flow_matrix.csv")))
      saveRDS(list(acc_filter     = acc_lbd_filter,
                   acc_gate       = acc_lbd_gate,
                   acc_total      = acc_lbd_total,
                   acc_fwd_filter = res_lbd_filter$acc_fwd,
                   acc_rev_filter = res_lbd_filter$acc_rev,
                   acc_fwd_gate   = res_lbd_gate$acc_fwd,
                   acc_rev_gate   = res_lbd_gate$acc_rev,
                   gephi          = gephi_lbd,
                   matrix         = mat_lbd,
                   graph          = g,
                   map_df         = map_df),
              file.path(OUTPUT_DIR, "rds", paste0(tag, ".rds")))

      log_msg(paste("  LBD edges saved:", nrow(gephi_lbd),
                    "above threshold"), sys_log)

      gephi_lbd_dir <- flow_to_gephi_directed(
        g,
        res_lbd_filter$acc_fwd, res_lbd_filter$acc_rev,
        res_lbd_gate$acc_fwd,   res_lbd_gate$acc_rev,
        map_df, MIN_FLOW_THRESHOLD, "LBD"
      )
      write_csv(gephi_lbd_dir,
                file.path(OUTPUT_DIR, "directed",
                          paste0(full_name, "_LBD_edges_directed.csv")))
      log_msg(paste("  LBD directed edges saved:",
                    nrow(gephi_lbd_dir), "rows"), sys_log)

      family_results$lbd            <- acc_lbd_total
      family_results$lbd_filter     <- acc_lbd_filter
      family_results$lbd_gate       <- acc_lbd_gate
      family_results$lbd_fwd_filter <- res_lbd_filter$acc_fwd
      family_results$lbd_rev_filter <- res_lbd_filter$acc_rev
      family_results$lbd_fwd_gate   <- res_lbd_gate$acc_fwd
      family_results$lbd_rev_gate   <- res_lbd_gate$acc_rev
    }
  }

  # ==========================================================
  # BLOCK C: COMBINED (SUM OF BOTH FAMILIES)
  # ==========================================================

  if (RUN_COMBINED &&
      !is.null(family_results$intracellular) &&
      !is.null(family_results$lbd)) {

    log_msg("--- COMBINED: intracellular + LBD ---", sys_log)

    acc_combined <- family_results$intracellular + family_results$lbd

    gephi_comb <- flow_to_gephi(g, acc_combined, map_df, MIN_FLOW_THRESHOLD,
                                 "combined")
    mat_comb   <- flow_to_matrix(gephi_comb, map_df)

    tag <- paste0(full_name, "_combined")

    write_csv(gephi_comb,
              file.path(OUTPUT_DIR, paste0(tag, "_edges_gephi.csv")))
    write_csv(as.data.frame(mat_comb),
              file.path(OUTPUT_DIR, paste0(tag, "_flow_matrix.csv")))
    saveRDS(list(acc_total = acc_combined,
                 gephi     = gephi_comb,
                 matrix    = mat_comb,
                 graph     = g,
                 map_df    = map_df),
            file.path(OUTPUT_DIR, "rds", paste0(tag, ".rds")))

    log_msg(paste("  Combined edges saved:", nrow(gephi_comb),
                  "above threshold"), sys_log)

    # Combined directed: sum fwd and rev vectors across both families
    acc_comb_fwd_filter <- family_results$intra_fwd_filter +
                           family_results$lbd_fwd_filter
    acc_comb_rev_filter <- family_results$intra_rev_filter +
                           family_results$lbd_rev_filter
    acc_comb_fwd_gate   <- family_results$intra_fwd_gate +
                           family_results$lbd_fwd_gate
    acc_comb_rev_gate   <- family_results$intra_rev_gate +
                           family_results$lbd_rev_gate

    gephi_comb_dir <- flow_to_gephi_directed(
      g,
      acc_comb_fwd_filter, acc_comb_rev_filter,
      acc_comb_fwd_gate,   acc_comb_rev_gate,
      map_df, MIN_FLOW_THRESHOLD, "combined"
    )
    write_csv(gephi_comb_dir,
              file.path(OUTPUT_DIR, "directed",
                        paste0(full_name, "_combined_edges_directed.csv")))
    log_msg(paste("  Combined directed edges saved:",
                  nrow(gephi_comb_dir), "rows"), sys_log)
  }

  elapsed <- round(as.numeric(difftime(Sys.time(), t0, units = "secs")), 1)
  log_msg(paste("=== Finished:", full_name, "| elapsed:", elapsed, "s ==="),
          sys_log)

  return(full_name)
})

stopCluster(cl)

# ==========================================================
# SECTION 12: SUMMARY
# ==========================================================

completed <- Filter(Negate(is.null), results)
log_msg(paste("Completed:", length(completed), "of",
              length(network_files), "systems"), master_log)

log_msg("=== 08a_v1.3 FINISHED ===", master_log)

# ==========================================================
# SECTION 13: PARAMETER REGISTRY UPDATE
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
      script_number  = "08a",
      script_version = "1.3",
      script_name    = "08a_v1.3_maxflow_allosteric_edges",
      variable = c("NETWORK_DIR", "MAP_FILE", "OUTPUT_DIR", "IS_NMA",
                   "TARGET_FILTER", "TARGET_GATE",
                   "seeds_intracellular_count", "seeds_LBD_count",
                   "MIN_FLOW_THRESHOLD", "N_CORES",
                   "systems_completed"),
      value    = c(NETWORK_DIR, MAP_FILE, OUTPUT_DIR, as.character(IS_NMA),
                   as.character(TARGET_FILTER), as.character(TARGET_GATE),
                   as.character(length(seeds_intra)),
                   as.character(length(seeds_lbd)),
                   as.character(MIN_FLOW_THRESHOLD),
                   as.character(N_CORES),
                   as.character(length(completed))),
      timestamp = ts,
      stringsAsFactors = FALSE
    )

    if (file.exists(final_path)) {
      existing <- read.csv(final_path, stringsAsFactors = FALSE,
                           colClasses = "character")
      existing <- existing[existing$script_number != "08a", ]
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
