# ==========================================================
# 08f_v1.3_seed_selection.R
# Author: DenyCB
#
# PURPOSE:
#   Select candidate allosteric source seeds from max-flow
#   accumulated edge flows (08a output). Seeds are identified
#   as nodes with high outgoing flow in the directed network,
#   computed as:
#     seed_score(i) = sum_j acc_fwd[i→j]   (total fwd outflow)
#   Transition probability T_ij = acc_fwd[i→j] / seed_score(i)
#   is also reported as a secondary score.
#
#   Two output CSVs are produced:
#     - seeds_WT_top100.csv   : top 100 seeds for WT systems
#     - seeds_MUT_top20.csv   : top 20 seeds for mutant systems
#   Both ordered by target, then by seed_score descending.
#   Each seed row includes cross-system equivalence counts
#   (how many chains across all replicas of each system have
#   a seed at the same residue number with score >= 10% of
#   the reported seed score).
#
# INPUTS:
#   - RDS files from 08a (intracellular and LBD per system/replica)
#   - Same format and directory structure as 08b input
#
# OUTPUTS:
#   - seeds_WT_top100.csv
#   - seeds_MUT_top20.csv
#   - logs/08f_master_log_<timestamp>.txt
#   - logs/08f_<system>_<timestamp>.txt (per system)
#
# NOTES:
#   - Compatible with Windows 10 via PSOCK parallel cluster
#   - NMA mode: one RDS per system (no _repN suffix)
#   - MD mode: multiple _repN RDS per system; scores averaged
#     across replicas before ranking
#   - Seeds are reported without intracellular/LBD distinction
#     (family label retained for reference only)
#   - Cross-system equivalence: a seed at residue R in system S
#     counts if any chain of any replica of S has seed_score >= 10%
#     of the score of the reported seed
#   - Target column uses three levels: "filter", "gate", "both"
#     when a seed scores above threshold for both targets
#
# ==========================================================

suppressPackageStartupMessages({
  library(dplyr)
  library(readr)
  library(stringr)
  library(parallel)
  library(igraph)
})

# ==========================================================
# SECTION 1: ANALYSIS SWITCHES
# ==========================================================

IS_NMA <- FALSE   # TRUE = NMA (single file per system), FALSE = MD (replicas)

# Direction to use for seed scoring
# Seeds are always from fwd (source→target direction)
USE_FWD <- TRUE   # always TRUE for seed selection; kept as explicit switch

# ==========================================================
# SECTION 2: INPUT/OUTPUT CONFIGURATION
# ==========================================================

# --- NMA (uncomment to use) ---
#INPUT_RDS_DIR <- "C:/DinamicasMoleculares/TRPV1_pipeline/output/08_Max_Flow/08a_v1.3_allosteric_edges_NMA/rds"
#MAP_FILE      <- "C:/DinamicasMoleculares/TRPV1_pipeline/config_files/trimmed_residue_map.csv"
#OUTPUT_DIR    <- "C:/DinamicasMoleculares/TRPV1_pipeline/output/08_Max_Flow/08f_v1.3_seed_selection_NMA_1"

# --- MD (comment out NMA block above and uncomment below) ---
 INPUT_RDS_DIR <- "C:/DinamicasMoleculares/TRPV1_pipeline/output/08_Max_Flow/08a_v1.3_allosteric_edges_MD/rds"
 MAP_FILE      <- "C:/DinamicasMoleculares/TRPV1_pipeline/output/residue_map_MD.csv"
 OUTPUT_DIR    <- "C:/DinamicasMoleculares/TRPV1_pipeline/output/08_Max_Flow/08f_v1.3_seed_selection_MD"

dir.create(OUTPUT_DIR,                    recursive = TRUE, showWarnings = FALSE)
dir.create(file.path(OUTPUT_DIR, "logs"), showWarnings = FALSE)

# ==========================================================
# SECTION 3: ANALYSIS PARAMETERS
# ==========================================================

# Target residue numbers (must match 08a settings)
TARGET_FILTER <- 643
TARGET_GATE   <- 679

# Top-K seeds to report per system type
TOP_K_WT  <- 500   # seeds for WT systems
TOP_K_MUT <-  20   # seeds for mutant systems

# Cross-system equivalence threshold:
# A seed at residue R in system S counts as equivalent if its
# seed_score >= EQUIV_THRESH_FRAC * (score of the reported seed)
EQUIV_THRESH_FRAC <- 0.10

# Parallel cores (PSOCK cluster, Windows 10 compatible)
N_CORES <- 6

# Residue range to exclude from seeds (too close to targets).
# Residues within this range are removed regardless of their score.
# Set to NULL to disable exclusion.
EXCLUDE_RESNUM_RANGE <- c(560, 710)   # S4-S5 linker to TRP helix




# ==========================================================
# SECTION 4: LOGGING SETUP
# ==========================================================

master_log <- file.path(
  OUTPUT_DIR, "logs",
  paste0("08f_v1.3_master_log_",
         format(Sys.time(), "%Y%m%d_%H%M%S"), ".txt")
)

log_msg <- function(msg, log_file = NULL) {
  txt <- sprintf("[%s] 08f_v1.3 | %s",
                 format(Sys.time(), "%Y-%m-%d %H:%M:%OS3"), msg)
  cat(txt, "\n")
  if (!is.null(log_file)) cat(txt, "\n", file = log_file, append = TRUE)
}

log_msg("=== STARTING 08f_v1.3 SEED SELECTION ===", master_log)
log_msg(paste("INPUT_RDS_DIR:", INPUT_RDS_DIR), master_log)
log_msg(paste("IS_NMA:", IS_NMA), master_log)
log_msg(paste("TARGET_FILTER:", TARGET_FILTER,
              "| TARGET_GATE:", TARGET_GATE), master_log)
log_msg(paste("TOP_K_WT:", TOP_K_WT, "| TOP_K_MUT:", TOP_K_MUT), master_log)
log_msg(paste("EQUIV_THRESH_FRAC:", EQUIV_THRESH_FRAC), master_log)

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
# --- Add cross-system equivalence counts (vectorized) ---
# For each seed row, counts how many chainresid instances exist in each
# other system at the same resnum with seed_score >= EQUIV_THRESH_FRAC
# × (this seed's score). Equivalence is by resID only (not by chain):
# each matching chainresid across all replicas counts independently.
add_cross_system_counts <- function(seed_df, all_system_names,
                                    all_scores, thresh_frac) {
  result <- seed_df
  
  # Expand "both" target rows to also appear as filter and gate
  # so that join works correctly for all target levels
  lookup <- all_scores %>%
    select(system, resnum, target, seed_score) %>%
    bind_rows(
      all_scores %>%
        filter(target == "both") %>%
        mutate(target = "filter") %>%
        select(system, resnum, target, seed_score),
      all_scores %>%
        filter(target == "both") %>%
        mutate(target = "gate") %>%
        select(system, resnum, target, seed_score)
    )
  
  for (sys_other in all_system_names) {
    
    col_name <- paste0("equiv_", sys_other)
    
    sys_lookup <- lookup %>%
      filter(system == sys_other) %>%
      select(resnum, target, score_other = seed_score)
    
    # Join seed_df to sys_lookup on resnum + target,
    # then apply threshold and count matches per original seed row
    counts <- seed_df %>%
      mutate(row_idx = row_number(),
             thresh  = seed_score * thresh_frac) %>%
      left_join(sys_lookup,
                by      = c("resnum", "target"),
                relationship = "many-to-many") %>%
      mutate(is_equiv = !is.na(score_other) & score_other >= thresh) %>%
      group_by(row_idx) %>%
      summarise(equiv_count = as.integer(sum(is_equiv, na.rm = TRUE)),
                .groups = "drop") %>%
      arrange(row_idx) %>%
      pull(equiv_count)
    
    result[[col_name]] <- counts
  }
  
  result
}
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

# --- Extract system name (first token of system_base) ---
get_structure_from_system <- function(system_name) {
  str_split(system_name, "_")[[1]][1]
}

# --- Compute seed scores for all nodes from acc_fwd ---
# seed_score(i) = sum_j acc_fwd[i→j]  (total outgoing fwd flow)
# max_T(i)      = max_j T_ij = max_j acc_fwd[i→j] / seed_score(i)
# Returns data.frame: node_id, chainresid, resnum, seed_score, max_T, fwd_sum
compute_seed_scores <- function(g, acc_fwd, map_df,
                                target_filter_nodes, target_gate_nodes,
                                log_file) {
  
  el <- as_edgelist(g, names = FALSE)
  n  <- vcount(g)
  
  target_nodes_all <- unique(c(target_filter_nodes, target_gate_nodes))
  
  seed_score_vec <- rep(0.0, n)
  max_T_vec      <- rep(0.0, n)
  
  for (ni in seq_len(n)) {
    out_edges <- which(el[, 1] == ni)
    if (length(out_edges) == 0) next
    total_out <- sum(acc_fwd[out_edges])
    seed_score_vec[ni] <- total_out
    if (total_out > 0) {
      max_T_vec[ni] <- max(acc_fwd[out_edges] / total_out)
    }
  }
  
  fwd_sum_vec <- seed_score_vec
  
  label_map <- setNames(
    paste0(map_df$chain, map_df$residue),
    as.character(map_df$node_id)
  )
  resnum_map <- setNames(
    map_df$residue,
    as.character(map_df$node_id)
  )
  
  df <- data.frame(
    node_id    = seq_len(n),
    chainresid = label_map[as.character(seq_len(n))],
    resnum     = resnum_map[as.character(seq_len(n))],
    seed_score = seed_score_vec,
    max_T      = max_T_vec,
    fwd_sum    = fwd_sum_vec,
    is_target  = seq_len(n) %in% target_nodes_all,
    stringsAsFactors = FALSE
  )
  
  df <- df[!df$is_target & !is.na(df$chainresid), ]
  df <- df[df$seed_score > 0, ]
  
  # Exclude residues too close to targets (e.g. pore-lining helices)
  if (!is.null(EXCLUDE_RESNUM_RANGE)) {
    df <- df[is.na(df$resnum) |
               df$resnum < EXCLUDE_RESNUM_RANGE[1] |
               df$resnum > EXCLUDE_RESNUM_RANGE[2], ]
  }
  
  df
}
# ==========================================================
# SECTION 7: DISCOVER AND PARSE INPUT FILES
# ==========================================================

rds_files <- list.files(INPUT_RDS_DIR, pattern = "\\.rds$",
                        full.names = TRUE)

if (length(rds_files) == 0) stop(paste("No RDS files found in:", INPUT_RDS_DIR))

log_msg(paste("RDS files found:", length(rds_files)), master_log)

# Parse metadata for all files
file_meta <- lapply(rds_files, function(f) {
  m <- parse_rds_filename(f)
  c(m, list(file = f))
})

meta_df <- bind_rows(lapply(file_meta, as.data.frame,
                            stringsAsFactors = FALSE))

# Work only with intracellular and LBD families (skip combined)
meta_df <- meta_df %>%
  filter(family %in% c("intracellular", "LBD"))

# List of unique system_base values
all_systems <- unique(meta_df$system_base)
log_msg(paste("Unique system bases:", length(all_systems)), master_log)

# ==========================================================
# SECTION 8: PARALLEL COMPUTATION OF SEED SCORES PER SYSTEM
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
  })
})

clusterExport(cl, c(
  "meta_df", "mapping", "OUTPUT_DIR", "IS_NMA",
  "TARGET_FILTER", "TARGET_GATE",
  "EXCLUDE_RESNUM_RANGE",          # ← agregar esta línea
  "compute_seed_scores", "get_structure_from_system",
  "log_msg"
))

# Process each system_base: load all replica files, compute scores,
# average across replicas (MD) or use directly (NMA)
system_scores <- parLapply(cl, all_systems, function(sys_base) {
  
  sys_log <- file.path(OUTPUT_DIR, "logs",
                       paste0("08f_", sys_base, "_",
                              format(Sys.time(), "%Y%m%d_%H%M%S"), ".txt"))
  
  log_msg(paste("=== Processing:", sys_base, "==="), sys_log)
  t0 <- Sys.time()
  
  sys_files <- meta_df %>%
    filter(system_base == sys_base)
  
  # Identify unique replicas (NA for NMA = single replica)
  replicas <- unique(sys_files$replica)
  
  # Accumulate scores per replica
  rep_scores_list <- list()
  
  for (rep_tag in replicas) {
    
    rep_files <- sys_files %>%
      filter(is.na(replica) & is.na(rep_tag) |
               !is.na(replica) & !is.na(rep_tag) & replica == rep_tag)
    
    if (nrow(rep_files) == 0) next
    
    g_ref      <- NULL
    map_df_sys <- NULL
    acc_fwd_filter_total <- NULL
    acc_fwd_gate_total   <- NULL
    
    for (ii in seq_len(nrow(rep_files))) {
      
      f   <- rep_files$file[ii]
      obj <- tryCatch(readRDS(f), error = function(e) {
        log_msg(paste("  ERROR loading:", basename(f), "-", e$message), sys_log)
        NULL
      })
      if (is.null(obj)) next
      
      if (is.null(g_ref)) {
        g_ref      <- obj$graph
        map_df_sys <- obj$map_df
      }
      
      # Accumulate fwd flow vectors across families within this replica
      if (!is.null(obj$acc_fwd_filter)) {
        acc_fwd_filter_total <- if (is.null(acc_fwd_filter_total))
          obj$acc_fwd_filter else acc_fwd_filter_total + obj$acc_fwd_filter
      }
      if (!is.null(obj$acc_fwd_gate)) {
        acc_fwd_gate_total <- if (is.null(acc_fwd_gate_total))
          obj$acc_fwd_gate else acc_fwd_gate_total + obj$acc_fwd_gate
      }
    }
    
    if (is.null(g_ref) || is.null(acc_fwd_filter_total)) {
      log_msg(paste("  WARNING: incomplete data for replica:", rep_tag), sys_log)
      next
    }
    
    structure_name <- get_structure_from_system(sys_base)
    map_df_use <- mapping %>%
      filter(system == structure_name) %>%
      arrange(node_id)
    
    target_filter_nodes <- map_df_use %>%
      filter(residue == TARGET_FILTER) %>% pull(node_id) %>% as.integer()
    target_gate_nodes <- map_df_use %>%
      filter(residue == TARGET_GATE) %>% pull(node_id) %>% as.integer()
    
    # Compute seed scores separately for filter and gate targets
    scores_filter <- compute_seed_scores(
      g_ref, acc_fwd_filter_total, map_df_use,
      target_filter_nodes, target_gate_nodes, sys_log
    ) %>%
      mutate(target = "filter",
             system = sys_base,
             replica = if (is.na(rep_tag)) "" else rep_tag)
    
    scores_gate <- compute_seed_scores(
      g_ref, acc_fwd_gate_total, map_df_use,
      target_filter_nodes, target_gate_nodes, sys_log
    ) %>%
      mutate(target = "gate",
             system = sys_base,
             replica = if (is.na(rep_tag)) "" else rep_tag)
    
    rep_scores_list[[length(rep_scores_list) + 1]] <-
      bind_rows(scores_filter, scores_gate)
  }
  
  if (length(rep_scores_list) == 0) {
    log_msg("  ERROR: no scores computed. Skipping.", sys_log)
    return(NULL)
  }
  
  all_rep_scores <- bind_rows(rep_scores_list)
  
  # Average seed_score and max_T across replicas for each chainresid × target
  avg_scores <- all_rep_scores %>%
    group_by(chainresid, resnum, target, system) %>%
    summarise(
      seed_score = mean(seed_score, na.rm = TRUE),
      max_T      = mean(max_T,      na.rm = TRUE),
      fwd_sum    = mean(fwd_sum,    na.rm = TRUE),
      n_replicas = n(),
      .groups = "drop"
    )
  
  # Identify seeds that score for both targets — assign "both" level
  filter_seeds <- avg_scores %>% filter(target == "filter") %>%
    select(chainresid, resnum, system, score_filter = seed_score)
  gate_seeds <- avg_scores %>% filter(target == "gate") %>%
    select(chainresid, resnum, system, score_gate = seed_score)
  
  both_seeds <- inner_join(filter_seeds, gate_seeds,
                           by = c("chainresid", "resnum", "system")) %>%
    filter(score_filter > 0 & score_gate > 0) %>%
    pull(chainresid)
  
  # Add "both" rows for seeds scoring for both targets
  both_rows <- avg_scores %>%
    filter(chainresid %in% both_seeds, target == "filter") %>%
    mutate(
      seed_score = (seed_score +
                      avg_scores$seed_score[
                        match(paste(chainresid, "gate"),
                              paste(avg_scores$chainresid, avg_scores$target))
                      ]) / 2,
      target = "both"
    )
  
  avg_scores_full <- bind_rows(avg_scores, both_rows) %>%
    arrange(target, desc(seed_score))
  
  elapsed <- round(as.numeric(difftime(Sys.time(), t0, units = "secs")), 1)
  log_msg(paste("=== Finished:", sys_base,
                "| elapsed:", elapsed, "s ==="), sys_log)
  
  avg_scores_full
})

stopCluster(cl)

# ==========================================================
# SECTION 9: ASSEMBLE GLOBAL SEED TABLE + CROSS-SYSTEM COUNTS
# ==========================================================

log_msg("Assembling global seed table...", master_log)

all_scores <- bind_rows(Filter(Negate(is.null), system_scores))

if (nrow(all_scores) == 0) stop("No seed scores computed. Check inputs.")

# Identify all systems present
all_system_names <- sort(unique(all_scores$system))
log_msg(paste("Systems with scores:", length(all_system_names)), master_log)

# For each seed row, compute cross-system equivalence counts:
# For each other system S, count how many unique chains in S
# have seed_score >= EQUIV_THRESH_FRAC * (this seed's score)
# at the same residue number (any chain, any replica).



log_msg("Computing cross-system equivalence counts...", master_log)

all_scores_with_equiv <- add_cross_system_counts(
  all_scores, all_system_names, all_scores, EQUIV_THRESH_FRAC
)

# ==========================================================
# SECTION 10: SPLIT INTO WT AND MUTANT CSVS
# ==========================================================

log_msg("Splitting into WT and mutant seed tables...", master_log)

# Identify WT systems (second token of system_base == "WT")
is_wt <- function(sys) str_split(sys, "_")[[1]][2] == "WT"

wt_systems  <- all_system_names[sapply(all_system_names, is_wt)]
mut_systems <- all_system_names[!sapply(all_system_names, is_wt)]

log_msg(paste("WT systems:", paste(wt_systems, collapse=", ")), master_log)
log_msg(paste("Mutant systems:", paste(mut_systems, collapse=", ")), master_log)

# Column order for output
base_cols  <- c("system", "target", "chainresid", "resnum",
                "seed_score", "max_T", "fwd_sum", "n_replicas")
equiv_cols <- paste0("equiv_", all_system_names)
out_cols   <- c(base_cols, equiv_cols)

# WT: top TOP_K_WT seeds per system, ordered by target then seed_score desc
seeds_wt <- all_scores_with_equiv %>%
  filter(system %in% wt_systems) %>%
  arrange(system, target, desc(seed_score)) %>%
  group_by(system, target) %>%
  slice_head(n = TOP_K_WT) %>%
  ungroup() %>%
  select(all_of(out_cols))

# Mutant: top TOP_K_MUT seeds per system
seeds_mut <- all_scores_with_equiv %>%
  filter(system %in% mut_systems) %>%
  arrange(system, target, desc(seed_score)) %>%
  group_by(system, target) %>%
  slice_head(n = TOP_K_MUT) %>%
  ungroup() %>%
  select(all_of(out_cols))

write_csv(seeds_wt,
          file.path(OUTPUT_DIR, "seeds_WT_top100.csv"))
write_csv(seeds_mut,
          file.path(OUTPUT_DIR, "seeds_MUT_top20.csv"))

log_msg(paste("WT seeds saved:", nrow(seeds_wt), "rows"), master_log)
log_msg(paste("Mutant seeds saved:", nrow(seeds_mut), "rows"), master_log)
log_msg("=== 08f_v1.3 FINISHED ===", master_log)

# ==========================================================
# SECTION 11: PARAMETER REGISTRY UPDATE
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
      script_number  = "08f",
      script_version = "1.0",
      script_name    = "08f_v1.3_seed_selection",
      run_type       = run_type,
      variable = c("INPUT_RDS_DIR", "MAP_FILE", "OUTPUT_DIR",
                   "IS_NMA", "run_type",
                   "TARGET_FILTER", "TARGET_GATE",
                   "TOP_K_WT", "TOP_K_MUT",
                   "EQUIV_THRESH_FRAC", "N_CORES",
                   "n_systems", "n_seeds_wt", "n_seeds_mut"),
      value = c(INPUT_RDS_DIR, MAP_FILE, OUTPUT_DIR,
                as.character(IS_NMA), run_type,
                as.character(TARGET_FILTER), as.character(TARGET_GATE),
                as.character(TOP_K_WT), as.character(TOP_K_MUT),
                as.character(EQUIV_THRESH_FRAC), as.character(N_CORES),
                as.character(length(all_system_names)),
                as.character(nrow(seeds_wt)),
                as.character(nrow(seeds_mut))),
      timestamp = ts,
      stringsAsFactors = FALSE
    )
    
    if (file.exists(final_path)) {
      existing <- read.csv(final_path, stringsAsFactors = FALSE,
                           colClasses = "character")
      # Keep other scripts and other run_types of 08f
      existing <- existing[!(existing$script_number == "08f" &
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