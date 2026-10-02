# =============================================================
# 03a_v1.2_compute_dccm_md.R
# Author: DenyCB
# =============================================================
# Compute Dynamic Cross-Correlation Matrices (DCCM)
# from MD trajectories
# Alignment 
# Parallel version compatible with Windows
# =============================================================

library(bio3d)
library(parallel)

timestamp <- function() format(Sys.time(), "%Y-%m-%d %H:%M:%S")

cat("[", timestamp(), "] Starting script 03a\n")

# ----------------------------
# Configuration
# ----------------------------

config_file <- "config.txt"

if (!file.exists(config_file)) stop("config.txt not found!")

config_lines <- readLines(config_file)

config <- list()

for (line in config_lines) {
  if (grepl("=", line) & !grepl("^#", line)) {

    parts <- strsplit(line, "=")[[1]]

    key <- trimws(parts[1])
    value <- trimws(parts[2])

    config[[key]] <- value
  }
}

pipeline_output_dir <- config$pipeline_output_dir
preprocessing_dir <- "C:/DinamicasMoleculares/analisis_bio3d_output/1_preprocessing"

# ----------------------------
# Load systems table
# ----------------------------

systems_file <- file.path("output","project_systems_table.csv")

if (!file.exists(systems_file))
  stop("Project systems table not found!")

systems <- read.csv(systems_file, stringsAsFactors=FALSE)

# ----------------------------
# Load residue mapping
# ----------------------------

map_file <- file.path("output","residue_structural_map.csv")

if (!file.exists(map_file))
  stop("Residue mapping table not found!")

res_map <- read.csv(map_file, stringsAsFactors=FALSE, check.names=TRUE)

# ----------------------------
# Output directory
# ----------------------------

dccm_output_dir <- file.path(pipeline_output_dir,"dccm_md")

if (!dir.exists(dccm_output_dir))
  dir.create(dccm_output_dir, recursive=TRUE)

log_dir <- file.path(dccm_output_dir, "logs")

if (!dir.exists(log_dir))
  dir.create(log_dir, recursive=TRUE)

log_file <- file.path(log_dir, paste0("03a_v1.2_compute_dccm_md_log_", format(Sys.time(), "%Y%m%d_%H%M%S"), ".txt"))

log_message <- function(...) {
  msg <- paste0("[", timestamp(), "] ", paste(..., collapse = " "))
  cat(msg, "\n")
  cat(msg, "\n", file = log_file, append = TRUE)
}

log_message("============================================================")
log_message("STARTING SCRIPT 03a_v1.2_compute_dccm_md.R")
log_message("pipeline_output_dir:", pipeline_output_dir)
log_message("preprocessing_dir:", preprocessing_dir)
log_message("dccm_output_dir:", dccm_output_dir)

# =============================================================
# PARALLEL SETUP (Windows compatible)
# =============================================================

n_cores <- 6

cat("[", timestamp(), "] Starting cluster with", n_cores, "cores\n")
log_message("Starting cluster with", n_cores, "cores")

cl <- makeCluster(n_cores)

# Load required library on each worker
clusterEvalQ(cl, library(bio3d))

systems_unique <- unique(systems[, c("structure","mutant")])

# Export objects needed by workers
clusterExport(cl,
              c("systems",
                "systems_unique",
                "timestamp",
                "log_message",
                "preprocessing_dir",
                "dccm_output_dir",
                "log_file"),
              envir=environment())

# =============================================================
# Parallel execution
# =============================================================



parLapply(cl, 1:nrow(systems_unique), function(i){

  sys <- systems_unique$structure[i]
  mut <- systems_unique$mutant[i]

  cond <- paste(sys, mut, sep="_")

  cond_dir <- file.path(preprocessing_dir, cond)

  aligned_ensemble_file <- file.path(cond_dir, paste0("aligned_ensemble_", cond, ".rds"))
  aligned_avg_file <- file.path(cond_dir, paste0("aligned_avg_", cond, ".rds"))

  log_message("------------------------------------------------------------")
  log_message("Processing:", cond)
  log_message("Aligned ensemble:", aligned_ensemble_file)
  log_message("Aligned average:", aligned_avg_file)

  if (!file.exists(aligned_ensemble_file))
    stop("Missing aligned ensemble:", aligned_ensemble_file)

  if (!file.exists(aligned_avg_file))
    stop("Missing aligned average:", aligned_avg_file)

  ens <- readRDS(aligned_ensemble_file)
  avg <- readRDS(aligned_avg_file)

  if (is.null(ens$trj))
    stop("aligned_ensemble has no trj list:", aligned_ensemble_file)

  if (is.null(avg$pdb))
    stop("aligned_avg has no selected pdb object:", aligned_avg_file)

  pdb_sel <- avg$pdb

  ca_inds <- atom.select(pdb_sel, elety="CA")

  if (length(ca_inds$atom) == 0)
    stop("No CA atoms detected in selected PDB for:", cond)

  nodes <- data.frame(
    index = seq_along(ca_inds$atom),
    chain = pdb_sel$atom$chain[ca_inds$atom],
    resno = pdb_sel$atom$resno[ca_inds$atom],
    node_label = paste0(pdb_sel$atom$chain[ca_inds$atom], pdb_sel$atom$resno[ca_inds$atom]),
    stringsAsFactors = FALSE
  )

  nodes_file <- file.path(dccm_output_dir,
                          paste0(cond,"_dccm_nodes.csv"))

  write.csv(nodes, nodes_file, row.names=FALSE)

  log_message("Nodes saved:", nodes_file)
  log_message("CA nodes:", nrow(nodes))
  log_message("Residue range detected:", paste(range(nodes$resno, na.rm=TRUE), collapse=" to "))
  log_message("Replicas detected:", paste(names(ens$trj), collapse=", "))

  dccm_list <- list()

  for (rep_name in names(ens$trj)) {

    rep <- gsub(".*_rep", "", rep_name)

    trj_fit <- ens$trj[[rep_name]]

    log_message("Processing:", sys, mut, "replica", rep)
    log_message("Frames:", nrow(trj_fit), "Columns:", ncol(trj_fit))

    xyz <- trj_fit[, ca_inds$xyz, drop=FALSE]

    if (ncol(xyz) %% 3 != 0)
      stop("XYZ coordinate vector not divisible by 3")

    log_message("CA xyz columns:", ncol(xyz))

    dccm_matrix <- dccm(xyz)

    base_name <- paste(sys, mut, paste0("rep",rep), sep="_")

    rds_file <- file.path(dccm_output_dir,
                          paste0(base_name,"_dccm.rds"))

    csv_file <- file.path(dccm_output_dir,
                          paste0(base_name,"_dccm.csv"))

    rep_nodes_file <- file.path(dccm_output_dir,
                                paste0(base_name,"_dccm_nodes.csv"))

    saveRDS(dccm_matrix, rds_file)

    write.csv(dccm_matrix, csv_file, row.names=FALSE)

    write.csv(nodes, rep_nodes_file, row.names=FALSE)

    dccm_list[[base_name]] <- dccm_matrix

    log_message("Saved replica DCCM:", rds_file)
  }

  if (length(dccm_list) == 0)
    stop("No replica DCCM matrices generated for:", cond)

  if (length(dccm_list) != 3)
    log_message("WARNING:", cond, "has", length(dccm_list), "replicas available; consensus will use available replicas.")

  dccm_mean <- Reduce("+", dccm_list) / length(dccm_list)

  base_name_mean <- paste(sys, mut, "repMean", sep="_")

  rds_file_mean <- file.path(dccm_output_dir,
                             paste0(base_name_mean,"_dccm.rds"))

  csv_file_mean <- file.path(dccm_output_dir,
                             paste0(base_name_mean,"_dccm.csv"))

  saveRDS(dccm_mean, rds_file_mean)

  write.csv(dccm_mean, csv_file_mean, row.names=FALSE)

  log_message("Saved mean DCCM:", rds_file_mean)

  if (length(dccm_list) > 1) {

    dccm_array <- simplify2array(dccm_list)

    dccm_variance <- apply(dccm_array, c(1,2), var, na.rm=TRUE)

    base_name_var <- paste(sys, mut, "repVariance", sep="_")

    rds_file_var <- file.path(dccm_output_dir,
                              paste0(base_name_var,"_dccm_variance.rds"))

    csv_file_var <- file.path(dccm_output_dir,
                              paste0(base_name_var,"_dccm_variance.csv"))

    saveRDS(dccm_variance, rds_file_var)

    write.csv(dccm_variance, csv_file_var, row.names=FALSE)

    log_message("Saved variance DCCM:", rds_file_var)
  }

  return(TRUE)

})

# ----------------------------
# Close cluster
# ----------------------------

stopCluster(cl)

log_message("============================================================")
log_message("Script 03a_v1.2 finished")
log_message("End time:", timestamp())
log_message("============================================================")

cat("[", timestamp(), "] Script 03a finished\n")


# ----------------------------------------------------------
# Parameter registry update
# ----------------------------------------------------------

tryCatch({
  
  registry_dir          <- file.path(config$pipeline_output_dir, "registry")
  registry_final_path   <- file.path(registry_dir, "parameters_final.csv")
  registry_history_path <- file.path(registry_dir, "parameters_history.csv")
  
  script_number  <- "03b"
  script_version <- "2.0"
  script_name    <- "03b_validate_nma_dccm"
  ts             <- format(Sys.time(), "%Y-%m-%d %H:%M:%S")
  
  # Variables to register
  new_entries <- data.frame(
    script_number  = script_number,
    script_version = script_version,
    script_name    = script_name,
    variable       = c("system_name",
                       "mutant_name",
                       "nma_dccm_dimensions",
                       "md_dccm_dimensions",
                       "common_residues_validation",
                       "validation_output_dir"),
    value          = c(system_name,
                       mutant_name,
                       paste(dim(dccm_nma_trim), collapse = "x"),
                       paste(dim(dccm_md_trim),  collapse = "x"),
                       as.character(length(common_labels)),
                       out_dir),
    timestamp      = ts,
    stringsAsFactors = FALSE
  )
  
  # --- Update parameters_final.csv ---
  # Read existing, remove rows from this script, append new entries
  if (file.exists(registry_final_path)) {
    existing <- read.csv(registry_final_path, stringsAsFactors = FALSE)
    # Remove any existing entries from this script
    existing <- existing[existing$script_number != script_number, ]
    updated  <- rbind(existing, new_entries)
  } else {
    updated <- new_entries
  }
  write.csv(updated, registry_final_path, row.names = FALSE)
  
  # --- Append to parameters_history.csv ---
  if (file.exists(registry_history_path)) {
    history <- read.csv(registry_history_path, stringsAsFactors = FALSE)
    history <- rbind(history, new_entries)
  } else {
    history <- new_entries
  }
  write.csv(history, registry_history_path, row.names = FALSE)
  
  cat("[", timestamp(), "] Parameter registry updated\n")
  
}, error = function(e) {
  cat("[", timestamp(), "] WARNING: registry update failed:", e$message, "\n")
})