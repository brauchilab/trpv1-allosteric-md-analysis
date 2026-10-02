# =========================================================
# Script 01 v2.1: Preprocess and Read Trajectories (Selected-PDB + Validations)
# Author: DenyCB
# Date: 2025-10-09
# =========================================================
# - Produce: aligned_ensemble_<cond>.rds
#            aligned_avg_<cond>.rds  (contiene pdb recortado, pdb_full, trj, sel)
#            average_<cond>.pdb     (selected atoms only) <- USAR por scripts 02/03
#            average_full_<cond>.pdb (backup)
#            _manifest_<cond>.csv
#            validation_summary.csv
#            preprocess_log.txt
# =========================================================

suppressPackageStartupMessages({
  library(bio3d)
  library(parallel)
  library(foreach)
  library(doParallel)
  library(tools)
})

# ---- Parameters ----
input_dir  <- "C:/DinamicasMoleculares/input/PDBs/DM"
output_dir <- "C:/DinamicasMoleculares/analisis_bio3d_output/1_preprocessing"
subsample_N <- 1   # Cambia a >1 para submuestrear (comentario en el código)
burnin_frames <- 101  # frames a descartar al inicio de cada réplica (equilibración)
cores <- 6

if(!dir.exists(output_dir)) dir.create(output_dir, recursive = TRUE)
log_file <- file.path(output_dir, "preprocess_log.txt")

log_message <- function(msg){
  ts <- format(Sys.time(), "%Y-%m-%d %H:%M:%OS6")
  cat(paste0(ts, " - ", msg, "\n"))
  cat(paste0(ts, " - ", msg, "\n"), file = log_file, append = TRUE)
}

log_message("=== Starting Preprocessing (v2.1) ===")
log_message(paste("Parameters: input_dir =", input_dir, "; output_dir =", output_dir,
                  "; subsample_N =", subsample_N, "; cores =", cores))

# Réplicas a excluir (dejar vacío para usar todas)
# Formato: vector de nombres base sin extensión
exclude_replicas <- c("7LPC_W426A_rep3")


# ---- Parallel ----
cl <- makeCluster(cores)
doParallel::registerDoParallel(cl)

# ---- discover conditions ----
pdb_files <- list.files(input_dir, pattern="\\.pdb$", full.names = TRUE)
base_names <- tools::file_path_sans_ext(basename(pdb_files))
conditions <- unique(sub("(_rep[0-9]+)$", "", base_names))
log_message(paste("Found", length(conditions), "conditions:", paste(conditions, collapse=", ")))

# ---- main loop ----
results <- foreach(cond = conditions, .packages = "bio3d") %dopar% {
  out_dir <- file.path(output_dir, cond)
  if(!dir.exists(out_dir)) dir.create(out_dir, recursive = TRUE)
  log_message(paste("Processing condition:", cond))
  
  # find replicas (pdb + dcd)
  rep_pdbs <- list.files(input_dir, pattern = paste0("^", cond, "_rep[0-9]+\\.pdb$"), full.names = TRUE)
  if(length(rep_pdbs) == 0){
    log_message(paste("No replicas found for", cond)); return(FALSE)
  }
  # Excluir réplicas marcadas como problemáticas en el parámetro exclude_replicas
  if (length(exclude_replicas) > 0) {
    rep_pdbs <- rep_pdbs[!tools::file_path_sans_ext(basename(rep_pdbs)) %in% exclude_replicas]
    if (length(rep_pdbs) == 0) {
      log_message(paste("All replicas excluded for", cond, "- skipping.")); return(FALSE)
    }
    log_message(paste("  Excluded replicas:", paste(exclude_replicas, collapse=", ")))
  }
  
  trj_list <- list()
  manifest <- data.frame()
  sel_global <- NULL
  
  # sequential per-replica processing (safer for logs)
  for(pdb_file in rep_pdbs){
    base <- tools::file_path_sans_ext(basename(pdb_file))
    dcd_file <- file.path(input_dir, paste0(base, ".dcd"))
    log_message(paste("  Replica:", base))
    tryCatch({
      pdb <- read.pdb(pdb_file)
      # selection (heavy atoms in residues 358-725)
      sel <- atom.select(pdb, elety = c("C","N","O","S","P","CA"))
      log_message(paste("    Selected atoms:", length(sel$atom), "-> (coords length)", length(sel$xyz)))
      # read trajectory and apply selection
      trj_raw <- read.dcd(dcd_file)
      log_message(paste("    Raw frames:", nrow(trj_raw)))
      trj_raw  <- trj_raw[-(1:burnin_frames), ]
      trj_sel  <- trj_raw[, sel$xyz]
      log_message(paste("    Frames after selection:", nrow(trj_sel), "Cols:", ncol(trj_sel)))
      # subsample (if >1)
      if(subsample_N > 1){
        trj_sel <- trj_sel[seq(1, nrow(trj_sel), by=subsample_N), ]
        log_message(paste("    Subsampled to", nrow(trj_sel), "frames."))
      }
      # fit trajectory to reference (selected atoms)
      trj_fit <- fit.xyz(fixed = pdb$xyz[sel$xyz], mobile = trj_sel)
      log_message("    RMS-fit completed.")
      trj_list[[base]] <- trj_fit
      manifest <- rbind(manifest, data.frame(pdb=pdb_file, dcd=dcd_file, out_dir=out_dir, stringsAsFactors = FALSE))
      sel_global <- sel  # store last sel (all reps should produce same sel)
    }, error = function(e){
      log_message(paste("    ERROR processing", base, ":", e$message))
    })
  } # end per-replica loop
  
  # Save ensemble RDS
  aligned_ensemble_file <- file.path(out_dir, paste0("aligned_ensemble_", cond, ".rds"))
  # Save ensemble as a list containing trajectories and the selection (pdb_sel) so downstream scripts find inputs$pdb_sel
  saveRDS(list(trj = trj_list, pdb_sel = sel_global), aligned_ensemble_file)
  log_message(paste("  Saved ensemble RDS:", aligned_ensemble_file))
  
  # Average across replicas (element-wise mean)
  tryCatch({
    avg_trj <- Reduce("+", trj_list) / length(trj_list)
    # full pdb template
    pdb_full <- read.pdb(rep_pdbs[1])  # backup full system for reference
    # assign averaged coords into the full pdb at selected indices
    pdb_full$xyz[sel_global$xyz] <- colMeans(avg_trj)
    # create selected-only pdb (trim)
    pdb_sel <- trim.pdb(pdb_full, sel_global)
    
    # Save selected pdb (this is the one scripts 02/03 should use)
    avg_pdb_file <- file.path(out_dir, paste0("average_", cond, ".pdb"))
    write.pdb(pdb_sel, file = avg_pdb_file)
    log_message(paste("  Saved average (selected) PDB:", avg_pdb_file))
    
    # Save full pdb as backup
    avg_full_file <- file.path(out_dir, paste0("average_full_", cond, ".pdb"))
    write.pdb(pdb_full, file = avg_full_file)
    log_message(paste("  Saved average full PDB (backup):", avg_full_file))
    
    # Save avg RDS with metadata (pdb selected, pdb_full backup, trj and sel)
    aligned_avg_file <- file.path(out_dir, paste0("aligned_avg_", cond, ".rds"))
    saveRDS(list(pdb = pdb_sel, pdb_full = pdb_full, trj = avg_trj, sel = sel_global), aligned_avg_file)
    log_message(paste("  Saved average RDS:", aligned_avg_file))
    
  }, error = function(e){
    log_message(paste("  ERROR generating average for", cond, ":", e$message))
  })
  
  # save manifest
  manifest_file <- file.path(out_dir, paste0("_manifest_", cond, ".csv"))
  write.csv(manifest, manifest_file, row.names = FALSE)
  log_message(paste("  Saved manifest:", manifest_file))
  
  return(TRUE)
} # end foreach

stopCluster(cl)
log_message("=== Finished Preprocessing; now running validations ===")

# -----------------------------
# Validation (improved)
# -----------------------------
check_rds_consistency_v2 <- function(output_dir){
  dirs <- list.dirs(output_dir, recursive = FALSE)
  summary <- data.frame()
  
  for(d in dirs){
    ens_file <- list.files(d, pattern="^aligned_ensemble_.*\\.rds$", full.names = TRUE)
    avg_file <- list.files(d, pattern="^aligned_avg_.*\\.rds$", full.names = TRUE)
    if(length(ens_file)==0 | length(avg_file)==0) next
    
    ens <- readRDS(ens_file)
    avg <- readRDS(avg_file)
    
    n_reps <- length(ens$trj)
    n_atoms <- sapply(ens$trj, ncol)
    n_frames <- sapply(ens$trj, nrow)
    atoms_match <- length(unique(n_atoms)) == 1
    frames_match <- length(unique(n_frames)) == 1
    # now compare against pdb (selected)
    pdb_atoms <- nrow(avg$pdb$atom)
    trj_atoms <- ncol(avg$trj) / 3
    pdb_consistent <- (pdb_atoms == trj_atoms)
    
    # RMSD between avg pdb coords and first frame of each replica (selected): detect if pdb equals initial frame
    vec_to_mat <- function(x) matrix(x, ncol=3, byrow=TRUE)
    pdb_mat <- vec_to_mat(avg$pdb$xyz)
    rmsd_firsts <- sapply(ens$trj, function(rep) {
      frame1_mat <- vec_to_mat(rep[1,])
      sqrt(mean(rowSums((pdb_mat - frame1_mat)^2)))
    })
    min_rmsd_first <- min(rmsd_firsts, na.rm=TRUE)
    matches_initial <- min_rmsd_first < 1e-3  # exact match threshold
    
    summary <- rbind(summary, data.frame(
      Condition = basename(d),
      Reps = n_reps,
      Atoms_Match = atoms_match,
      Frames_Match = frames_match,
      PDB_Consistent = pdb_consistent,
      PDB_matches_initial = matches_initial,
      min_rmsd_first = signif(min_rmsd_first,6),
      n_atoms = unique(n_atoms),
      n_frames = unique(n_frames)
    ))
  }
  # write summary
  write.csv(summary, file = file.path(output_dir, "validation_summary.csv"), row.names = FALSE)
  return(summary)
}

validation <- check_rds_consistency_v2(output_dir)
log_message("Validation Summary:")
capture.output(print(validation), file = log_file, append = TRUE)
log_message(paste("Validation passed for", sum(validation$Atoms_Match & validation$Frames_Match & validation$PDB_Consistent),
                  "of", nrow(validation), "conditions."))
log_message("=== Preprocessing Completed (v2.1) ===")
# End of script