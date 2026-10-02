# ================================================================
# Script 02A v5.0: NMA Generation (per-replica + consolidated)
# Author: DenyCB
# Date: 2025-04-29
# ----------------------------------------------------------------
# Inputs (from Script 01):
#   - aligned_ensemble_<cond>.rds  -> $trj (list of matrices, one per replica)
#                                  -> $pdb_sel (pdb object, heavy atoms)
#   - aligned_avg_<cond>.rds       -> $pdb (trimmed pdb, heavy atoms)
#                                     $pdb_full, $trj, $sel  [fallback]
#
# Auxiliary inputs:
#   - config.txt                   -> pipeline parameters
#   - trimmed_residue_map.csv      -> node_id / chain / residue / system
#                                     (used to validate/remap resno if needed)
#
# Outputs (under output_root/nma_data/):
#   Consolidated (one per system):
#     <system>_nma_obj.rds / .RData        <- NMA on global mean coordinates
#     <system>_enma_obj.rds                <- enma object from nma.pdbs() [convergence check]
#     <system>_rmsip_replicas.csv          <- RMSIP matrix between replicas [convergence check]
#     <system>_f_100.rds                   <- fluctuations over n_modes_full modes
#     <system>_mode1_fluct.rds / .csv
#     <system>_eigenvalues_nontrivial.csv
#     <system>_cumulative_variance.csv
#     <system>_nma_correlation_matrix.rds / .RData / .csv (if <= 2000 res)
#     <system>_NMA_trajectory_mode1.pdb
#     <system>_pdb_ca_trimmed.rds          <- Cα-only trimmed pdb used for NMA
#
#   Per-replica (under nma_data/per_rep/):
#     <system>_<rep>_nma_obj.rds / .RData
#     <system>_<rep>_f_100.rds
#     <system>_<rep>_mode1_fluct.rds / .csv
#     <system>_<rep>_eigenvalues_nontrivial.csv
#     <system>_<rep>_cumulative_variance.csv
#     <system>_<rep>_nma_correlation_matrix.rds
#     <system>_<rep>_NMA_trajectory_mode1.pdb
#     <system>_<rep>_pdb_ca_trimmed.rds
#
# ================================================================

suppressPackageStartupMessages({
  library(bio3d)
  library(parallel)
  library(foreach)
  library(doParallel)
})

# ================================================================
# SECTION 1: PARAMETERS
# ================================================================

config_file           <- "C:/DinamicasMoleculares/Pipeline_Complete/config.txt"
trimmed_map_file      <- "C:/DinamicasMoleculares/TRPV1_pipeline/config_files/trimmed_residue_map.csv"
input_preproc_dir     <- "C:/DinamicasMoleculares/analisis_bio3d_output/1_preprocessing"
output_root           <- "C:/DinamicasMoleculares/analisis_bio3d_output/2_nma_correlation"
cores                 <- 6

# Residue ranges to keep for NMA (Cα only)
# Range covers all residues resolved in common across all systems
trim_ranges           <- list(c(277, 602), c(625, 752))

# NMA fluctuation modes
# f_100 is computed over n_modes_full non-trivial modes.
# n_modes_short removed: downstream scripts select which modes to display.
n_modes_full          <- 200   # for f_100

# Mode trajectory output parameters
# cutoff_mode_trj: amplitude in sd units for mktrj.nma() visualization PDB.
# Value of 4 produces realistic-amplitude visualization (typical range 3-5).
cutoff_mode_trj       <- 4
n_trj_frames          <- 40

# ================================================================
# SECTION 2: OUTPUT DIRECTORIES
# ================================================================

dir.create(file.path(output_root, "nma_data"),          recursive = TRUE, showWarnings = FALSE)
dir.create(file.path(output_root, "nma_data", "per_rep"), recursive = TRUE, showWarnings = FALSE)
dir.create(file.path(output_root, "logs"),               recursive = TRUE, showWarnings = FALSE)

# ================================================================
# SECTION 3: LOGGING
# ================================================================

log_file_master <- file.path(
  output_root, "logs",
  paste0("02A_log_", format(Sys.time(), "%Y%m%d_%H%M%S"), ".txt")
)
log_buffer <- character()

log_message <- function(...) {
  msg <- paste0(format(Sys.time(), "%Y-%m-%d %H:%M:%OS3"), " - ", paste(..., collapse = " "))
  cat(msg, "\n")
  log_buffer <<- c(log_buffer, msg)
}

flush_log <- function(append = TRUE) {
  if (length(log_buffer) > 0) {
    cat(paste0(log_buffer, collapse = "\n"), "\n", file = log_file_master, append = append)
    log_buffer <<- character()
  }
}

log_message("=== Starting 02A v5.0: NMA Generation (per-replica + consolidated) ===")

# ================================================================
# SECTION 4: LOAD CONFIG AND TRIMMED MAP
# ================================================================

# --- Read config.txt (key = value format) ---
read_config <- function(cfg_file) {
  lines <- readLines(cfg_file)
  lines <- lines[!grepl("^\\s*#", lines) & nchar(trimws(lines)) > 0]
  cfg   <- list()
  for (ln in lines) {
    parts <- strsplit(ln, "=", fixed = TRUE)[[1]]
    if (length(parts) >= 2) {
      key       <- trimws(parts[1])
      val       <- trimws(paste(parts[-1], collapse = "="))
      cfg[[key]] <- val
    }
  }
  return(cfg)
}

cfg <- tryCatch(read_config(config_file), error = function(e) {
  log_message("WARNING: Could not read config file:", e$message, "- using script defaults.")
  list()
})
log_message("Config loaded. Keys found:", paste(names(cfg), collapse = ", "))

# --- Read trimmed residue map ---
# Columns: node_id, chain, residue, system
# Used to validate that the Cα trimming produces the expected node set,
# and to provide correct resno mapping for downstream scripts (e.g. DCCM labeling).
# If not found, validation is skipped with a WARNING (script continues normally).
trimmed_map <- tryCatch({
  df <- read.csv(trimmed_map_file, stringsAsFactors = FALSE)
  log_message("Loaded trimmed_residue_map.csv:", nrow(df), "rows,",
              "systems:", paste(unique(df$system), collapse = ", "))
  df
}, error = function(e) {
  log_message("WARNING: Could not read trimmed_residue_map.csv:", e$message,
              "- residue mapping validation will be skipped.")
  NULL
})

# ================================================================
# SECTION 5: HELPER FUNCTIONS
# ================================================================

# -- Discover condition directories from preprocessing output --
find_conditions <- function(preproc_dir) {
  dirs <- list.dirs(preproc_dir, recursive = FALSE, full.names = FALSE)
  dirs[nchar(dirs) > 0]
}

# -- Load aligned_ensemble RDS produced by Script 01 --
# Returns list with $trj (named list of matrices) and $pdb_sel (pdb object, heavy atoms)
load_ensemble <- function(cond_dir, cond_name) {
  ens_file <- file.path(cond_dir, paste0("aligned_ensemble_", cond_name, ".rds"))
  if (!file.exists(ens_file)) stop(paste("aligned_ensemble RDS not found:", ens_file))
  ens <- readRDS(ens_file)

  # Ensure $pdb_sel is a proper pdb object; if it is a 'select' index, recover from aligned_avg
  if (!inherits(ens$pdb_sel, "pdb")) {
    log_message("  pdb_sel is not a pdb object - attempting recovery from aligned_avg...")
    avg_file <- file.path(cond_dir, paste0("aligned_avg_", cond_name, ".rds"))
    if (!file.exists(avg_file)) stop(paste("aligned_avg RDS not found and pdb_sel is not a pdb:", avg_file))
    avg <- readRDS(avg_file)
    pdb_template <- if (!is.null(avg$pdb_full)) avg$pdb_full else avg$pdb
    if (is.null(pdb_template)) stop("No pdb_full/pdb found in aligned_avg for recovery.")
    sel_obj <- ens$pdb_sel
    ens$pdb_sel <- trim.pdb(pdb_template, sel = sel_obj)
    log_message("  pdb_sel recovered via trim.pdb from aligned_avg.")
  }

  return(ens)
}

# -- Extract Cα-only pdb trimmed to the specified residue ranges --
# pdb_in  : a bio3d pdb object (heavy atoms)
# ranges  : list of c(start, end) pairs, e.g. list(c(277,602), c(625,752))
# system  : system name (used to filter trimmed_map if provided)
# map_df  : optional trimmed_residue_map data.frame for validation
extract_ca_trimmed <- function(pdb_in, ranges, system = NULL, map_df = NULL) {
  # Build combined residue number vector
  resno_keep <- unlist(lapply(ranges, function(r) r[1]:r[2]))

  sel <- atom.select(pdb_in, elety = "CA", resno = resno_keep)
  if (length(sel$atom) == 0) stop("No Cα atoms found in specified residue ranges.")

  pdb_ca <- trim.pdb(pdb_in, sel = sel)

  # Log first and last 5 resno per chain for verification of atom ordering
  atom_df <- pdb_ca$atom
  for (ch in unique(atom_df$chain)) {
    resno_ch <- atom_df$resno[atom_df$chain == ch]
    log_message("  Chain", ch, "- first 5 resno:", paste(head(resno_ch, 5), collapse = ","),
                "| last 5 resno:", paste(tail(resno_ch, 5), collapse = ","),
                "| total Cα:", length(resno_ch))
  }

  # Validate against trimmed_map if available
  if (!is.null(map_df) && !is.null(system)) {
    map_sys <- map_df[map_df$system == system, ]
    if (nrow(map_sys) > 0) {
      expected_n   <- nrow(map_sys)
      obtained_n   <- nrow(pdb_ca$atom)
      expected_res <- sort(map_sys$residue)
      obtained_res <- sort(pdb_ca$atom$resno)
      if (obtained_n != expected_n) {
        log_message("  WARNING: trimmed Cα count mismatch for", system,
                    "- expected", expected_n, "got", obtained_n)
      } else if (!identical(expected_res, obtained_res)) {
        log_message("  WARNING: trimmed resno mismatch for", system,
                    "- check ranges vs trimmed_residue_map.csv")
      } else {
        log_message("  Trimmed Cα validated against trimmed_residue_map.csv:",
                    obtained_n, "atoms OK.")
      }
    }
  }

  return(pdb_ca)
}

# -- Build a mean-coordinate pdb from a trajectory matrix (nframes x 3*natoms) --
# Returns a copy of pdb_template with xyz set to the column means of trj_mat
build_mean_pdb <- function(pdb_template, trj_mat) {
  pdb_out      <- pdb_template
  pdb_out$xyz  <- as.numeric(colMeans(trj_mat, na.rm = TRUE))
  return(pdb_out)
}

# -- Build global mean-coordinate pdb averaged across ALL replicas --
# trj_list : named list of trajectory matrices (one per replica)
# pdb_template : heavy-atom pdb object to use as coordinate template
# Returns a pdb object with xyz = mean of all replica means
build_global_mean_pdb <- function(pdb_template, trj_list) {
  # Compute per-replica column means, then average across replicas
  rep_means  <- lapply(trj_list, function(trj) colMeans(trj, na.rm = TRUE))
  global_mean <- Reduce("+", rep_means) / length(rep_means)
  pdb_out      <- pdb_template
  pdb_out$xyz  <- as.numeric(global_mean)
  return(pdb_out)
}

# -- Construct a bio3d pdbs object from a list of Cα-only pdb objects --
# All pdbs must have the same number of atoms (identical trimming).
# Returns a pdbs object suitable for nma.pdbs().
build_pdbs_object <- function(pdb_list, names_vec) {
  n_structs <- length(pdb_list)
  n_atoms   <- nrow(pdb_list[[1]]$atom)

  # xyz matrix: rows = structures, cols = 3*n_atoms (x1,y1,z1,x2,...)
  xyz_mat <- do.call(rbind, lapply(pdb_list, function(p) as.numeric(p$xyz)))
  rownames(xyz_mat) <- names_vec

  # resno, resid, chain from first structure (all identical after trimming)
  ref        <- pdb_list[[1]]
  resno_vec  <- ref$atom$resno
  resid_vec  <- ref$atom$resid
  chain_vec  <- ref$atom$chain
  b_vec      <- ref$atom$b

  pdbs_obj <- list(
    xyz   = xyz_mat,
    resno = matrix(rep(resno_vec, n_structs), nrow = n_structs, byrow = TRUE),
    resid = matrix(rep(resid_vec, n_structs), nrow = n_structs, byrow = TRUE),
    chain = matrix(rep(chain_vec, n_structs), nrow = n_structs, byrow = TRUE),
    b     = matrix(rep(b_vec,     n_structs), nrow = n_structs, byrow = TRUE),
    id    = names_vec,
    call  = match.call()
  )
  class(pdbs_obj) <- "pdbs"
  return(pdbs_obj)
}

# -- Run NMA on a single Cα pdb and compute fluctuations --
run_single_nma <- function(pdb_ca, label, n_full, cutoff_trj, n_frames, outdir) {
  n_ca      <- nrow(pdb_ca$atom)
  max_modes <- (3 * n_ca) - 6
  if (max_modes < 1) stop(paste("Insufficient Cα atoms for NMA:", n_ca))

  nmodes_full  <- min(n_full,  max_modes)

  log_message("  Running nma() for", label, "| Cα atoms:", n_ca,
              "| modes (full):", nmodes_full)

  nma_obj <- tryCatch(
    nma(pdb_ca, ff = "calpha"),
    error = function(e) {
      log_message("  ERROR nma() failed for", label, ":", e$message)
      NULL
    }
  )
  if (is.null(nma_obj)) return(NULL)

  # -- Available non-trivial modes --
  av_modes <- length(nma_obj$L) - 6
  nmodes_full  <- min(nmodes_full,  av_modes)

  f_100       <- tryCatch(fluct.nma(nma_obj, mode.inds = 7:(6 + nmodes_full)),  error = function(e) NULL)
  mode1_fluct <- tryCatch(fluct.nma(nma_obj, mode.inds = 7),                    error = function(e) NULL)

  # -- Save outputs --
  save_nma_outputs(label, nma_obj, f_100, mode1_fluct,
                   pdb_ca, cutoff_trj, n_frames, outdir)

  return(list(nma_obj = nma_obj, f_100 = f_100, mode1_fluct = mode1_fluct))
}

# -- Save all NMA outputs for a given label to outdir --
save_nma_outputs <- function(label, nma_obj, f_100, mode1_fluct,
                             pdb_ca, cutoff_trj, n_frames, outdir) {
  log_message("  Saving NMA outputs for:", label)

  # Core objects
  tryCatch({
    saveRDS(nma_obj, file.path(outdir, paste0(label, "_nma_obj.rds")))
    save(nma_obj,   file = file.path(outdir, paste0(label, "_nma_obj.RData")))
    if (!is.null(f_100))       saveRDS(f_100,       file.path(outdir, paste0(label, "_f_100.rds")))
    if (!is.null(mode1_fluct)) saveRDS(mode1_fluct, file.path(outdir, paste0(label, "_mode1_fluct.rds")))
  }, error = function(e) log_message("  ERROR saving core RDS for", label, ":", e$message))

  # Eigenvalues and cumulative variance
  tryCatch({
    ev <- as.numeric(nma_obj$L)
    ev_nt  <- ev[-(1:6)]
    cumvar <- cumsum(ev_nt) / sum(ev_nt)
    write.csv(
      data.frame(mode = seq_along(ev_nt), eigenvalue = ev_nt),
      file.path(outdir, paste0(label, "_eigenvalues_nontrivial.csv")),
      row.names = FALSE
    )
    write.csv(
      data.frame(mode = seq_along(cumvar), cumulative_variance = cumvar),
      file.path(outdir, paste0(label, "_cumulative_variance.csv")),
      row.names = FALSE
    )
    # Informational log: how many modes needed to reach 80% cumulative variance
    modes_80 <- which(cumvar >= 0.80)[1]
    if (!is.na(modes_80)) {
      log_message("  INFO cumvar: modes needed for >=80% variance:", modes_80,
                  "| cumvar at mode", nmodes_full, ":", round(cumvar[min(nmodes_full, length(cumvar))], 3))
    } else {
      log_message("  INFO cumvar: 80% variance not reached within available modes.")
    }
  }, error = function(e) log_message("  ERROR saving eigenvalue CSVs for", label, ":", e$message))

  # DCCM (per-replica or single-structure; consolidated uses save_consolidated_dccm())
  tryCatch({
    cm <- dccm(nma_obj)
    saveRDS(cm,  file.path(outdir, paste0(label, "_nma_correlation_matrix.rds")))
    save(cm, file = file.path(outdir, paste0(label, "_nma_correlation_matrix.RData")))
    if (nrow(cm) <= 2000) {
      write.csv(cm, file.path(outdir, paste0(label, "_nma_correlation_matrix.csv")), row.names = FALSE)
    } else {
      log_message("  DCCM > 2000 residues - saved as RDS only.")
    }
  }, error = function(e) log_message("  ERROR computing/saving DCCM for", label, ":", e$message))

  # Mode 1 trajectory PDB
  tryCatch({
    mktrj.nma(nma_obj, mode = 7,
              file = file.path(outdir, paste0(label, "_NMA_trajectory_mode1.pdb")),
              n = n_frames, sd = cutoff_trj)
  }, error = function(e) log_message("  WARNING: mktrj.nma failed for", label, ":", e$message))

  # Mode 1 fluctuation CSV
  tryCatch({
    if (!is.null(mode1_fluct)) {
      write.csv(
        as.numeric(mode1_fluct),
        file.path(outdir, paste0(label, "_mode1_fluct.csv")),
        row.names = FALSE
      )
    }
  }, error = function(e) log_message("  ERROR saving mode1_fluct CSV for", label, ":", e$message))

  log_message("  Done saving outputs for:", label)
}

# -- Compute and save averaged DCCM across all replicas for a condition --
# enma_obj must have been computed with full=TRUE so full.nma is available.
# The averaged DCCM is the element-wise mean of the individual DCCM matrices.
save_consolidated_dccm <- function(cond, enma_obj, outdir) {
  log_message("  Computing averaged DCCM across replicas for:", cond)

  nma_list <- tryCatch(enma_obj$full.nma, error = function(e) NULL)
  if (is.null(nma_list) || length(nma_list) == 0) {
    log_message("  ERROR: enma_obj$full.nma is NULL or empty - cannot compute averaged DCCM.")
    return(invisible(NULL))
  }

  # Compute DCCM for each replica's nma object
  cm_list <- lapply(seq_along(nma_list), function(i) {
    tryCatch(
      dccm(nma_list[[i]]),
      error = function(e) {
        log_message("  WARNING: dccm() failed for replica", i, ":", e$message)
        NULL
      }
    )
  })
  cm_list <- Filter(Negate(is.null), cm_list)

  if (length(cm_list) == 0) {
    log_message("  ERROR: all per-replica DCCM computations failed for", cond)
    return(invisible(NULL))
  }
  if (length(cm_list) < length(nma_list)) {
    log_message("  WARNING: only", length(cm_list), "of", length(nma_list),
                "replica DCCMs computed - averaging over available ones.")
  }

  # Element-wise mean across replica DCCMs
  cm_avg <- Reduce("+", cm_list) / length(cm_list)

  # Save averaged DCCM
  tryCatch({
    saveRDS(cm_avg, file.path(outdir, paste0(cond, "_nma_correlation_matrix.rds")))
    save(cm_avg, file = file.path(outdir, paste0(cond, "_nma_correlation_matrix.RData")))
    if (nrow(cm_avg) <= 2000) {
      write.csv(cm_avg,
                file.path(outdir, paste0(cond, "_nma_correlation_matrix.csv")),
                row.names = FALSE)
    } else {
      log_message("  Averaged DCCM > 2000 residues - saved as RDS/RData only.")
    }
    log_message("  Saved averaged DCCM (", length(cm_list), "replicas ) for:", cond)
  }, error = function(e) {
    log_message("  ERROR saving averaged DCCM for", cond, ":", e$message)
  })

  return(invisible(cm_avg))
}

# ================================================================
# SECTION 6: MAIN PIPELINE
# ================================================================

conds <- find_conditions(input_preproc_dir)
if (length(conds) == 0) stop("No conditions found in preprocessing directory.")
log_message("Conditions found:", paste(conds, collapse = ", "))

cl <- makeCluster(cores)
doParallel::registerDoParallel(cl)

results <- foreach(
  cond = conds,
  .packages = "bio3d",
  .export = c(
    "load_ensemble", "extract_ca_trimmed", "build_mean_pdb",
    "build_global_mean_pdb", "build_pdbs_object",
    "run_single_nma", "save_nma_outputs",
    "save_consolidated_dccm",
    "log_message", "trim_ranges", "n_modes_full",
    "cutoff_mode_trj", "n_trj_frames", "output_root",
    "input_preproc_dir", "trimmed_map"
  ),
  .errorhandling = "pass"
) %dopar% {

  log_message("============================================================")
  log_message("Processing condition:", cond)

  cond_dir <- file.path(input_preproc_dir, cond)
  nma_data_dir   <- file.path(output_root, "nma_data")
  per_rep_dir    <- file.path(output_root, "nma_data", "per_rep")

  # --- 6.1 Load ensemble ---
  ens <- tryCatch(
    load_ensemble(cond_dir, cond),
    error = function(e) {
      log_message("ERROR loading ensemble for", cond, ":", e$message)
      NULL
    }
  )
  if (is.null(ens)) return(list(cond = cond, status = "error: load_ensemble failed"))

  rep_ids <- names(ens$trj)
  if (is.null(rep_ids) || length(rep_ids) == 0) {
    return(list(cond = cond, status = "error: no replica names found in ensemble$trj"))
  }
  log_message("  Replicas:", paste(rep_ids, collapse = ", "))

  # --- 6.2 Per-replica NMA ---
  rep_pdb_list   <- list()   # will accumulate Cα trimmed pdbs for convergence check
  rep_name_list  <- character()

  for (rep_id in rep_ids) {
    # Label: if rep_id already contains cond, use as-is; otherwise compose
    label <- if (grepl(cond, rep_id, fixed = TRUE)) rep_id else paste0(cond, "_", rep_id)
    log_message("  --- Replica:", label)

    trj_rep <- ens$trj[[rep_id]]  # matrix: nframes x 3*n_heavy_atoms

    # Build mean-coordinate pdb for this replica (from heavy-atom pdb template)
    pdb_mean_rep <- tryCatch(
      build_mean_pdb(ens$pdb_sel, trj_rep),
      error = function(e) {
        log_message("  ERROR build_mean_pdb for", label, ":", e$message)
        NULL
      }
    )
    if (is.null(pdb_mean_rep)) next

    # Extract Cα trimmed to the analysis ranges
    pdb_ca_rep <- tryCatch(
      extract_ca_trimmed(pdb_mean_rep, trim_ranges, system = cond, map_df = trimmed_map),
      error = function(e) {
        log_message("  ERROR extract_ca_trimmed for", label, ":", e$message)
        NULL
      }
    )
    if (is.null(pdb_ca_rep)) next

    # Save trimmed Cα pdb for reference (used by 02B for plotting context)
    saveRDS(pdb_ca_rep, file.path(per_rep_dir, paste0(label, "_pdb_ca_trimmed.rds")))
    log_message("  Saved trimmed Cα pdb:", label)

    # Run single-structure NMA
    nma_result <- run_single_nma(
      pdb_ca   = pdb_ca_rep,
      label    = label,
      n_full   = n_modes_full,
      cutoff_trj = cutoff_mode_trj,
      n_frames   = n_trj_frames,
      outdir     = per_rep_dir
    )
    if (!is.null(nma_result)) {
      rep_pdb_list[[label]]  <- pdb_ca_rep
      rep_name_list          <- c(rep_name_list, label)
      log_message("  Per-replica NMA completed:", label)
    } else {
      log_message("  WARNING: NMA failed for replica", label, "- excluded from consolidated.")
    }
  } # end per-replica loop

  # --- 6.3 Consolidated NMA on global mean coordinates ---
  # The consolidated NMA is computed on the mean coordinates averaged across
  # ALL replicas, which is the methodologically correct approach for summarizing
  # the equilibrium structure of a system sampled by multiple MD replicas.
  # This produces a single representative NMA per system for downstream analysis.

  log_message("  Building global mean pdb (averaged across all replicas) for:", cond)

  pdb_global_mean <- tryCatch(
    build_global_mean_pdb(ens$pdb_sel, ens$trj),
    error = function(e) {
      log_message("  ERROR build_global_mean_pdb for", cond, ":", e$message)
      NULL
    }
  )

  if (is.null(pdb_global_mean)) {
    return(list(cond = cond, status = "error: build_global_mean_pdb failed"))
  }

  # Extract Cα trimmed from global mean pdb
  pdb_ca_global <- tryCatch(
    extract_ca_trimmed(pdb_global_mean, trim_ranges, system = cond, map_df = trimmed_map),
    error = function(e) {
      log_message("  ERROR extract_ca_trimmed (global mean) for", cond, ":", e$message)
      NULL
    }
  )

  if (is.null(pdb_ca_global)) {
    return(list(cond = cond, status = "error: extract_ca_trimmed (global mean) failed"))
  }

  # Save trimmed Cα pdb of global mean (used by 02B for plotting context)
  saveRDS(pdb_ca_global, file.path(nma_data_dir, paste0(cond, "_pdb_ca_trimmed.rds")))
  log_message("  Saved global mean trimmed Cα pdb for:", cond)

  # Run consolidated NMA on global mean structure
  log_message("  Running consolidated NMA on global mean coordinates for:", cond)
  nma_consolidated <- run_single_nma(
    pdb_ca   = pdb_ca_global,
    label    = cond,
    n_full   = n_modes_full,
    cutoff_trj = cutoff_mode_trj,
    n_frames   = n_trj_frames,
    outdir     = nma_data_dir
  )

  if (is.null(nma_consolidated)) {
    return(list(cond = cond, status = "error: consolidated NMA failed"))
  }

  log_message("  Consolidated NMA completed for:", cond)

  # --- 6.4 Convergence check: nma.pdbs() on per-replica Cα pdbs ---
  # nma.pdbs() is used here solely to compute RMSIP between replicas,
  # which quantifies whether all replicas explored the same conformational subspace.
  # RMSIP > 0.7 indicates good convergence.
  # The enma_obj is saved for reference but is NOT the primary NMA output.

  if (length(rep_pdb_list) >= 2) {

    log_message("  Running nma.pdbs() convergence check (", length(rep_pdb_list), "replicas ) for:", cond)

    # Verify all replicas have the same atom count (required by nma.pdbs)
    atom_counts <- sapply(rep_pdb_list, function(p) nrow(p$atom))
    if (length(unique(atom_counts)) > 1) {
      log_message("  WARNING: replicas have different atom counts:", paste(atom_counts, collapse = ", "),
                  "- skipping convergence check.")
    } else {

      pdbs_obj <- tryCatch(
        build_pdbs_object(rep_pdb_list, rep_name_list),
        error = function(e) {
          log_message("  ERROR build_pdbs_object for convergence check:", cond, ":", e$message)
          NULL
        }
      )

      if (!is.null(pdbs_obj)) {
        enma_obj <- tryCatch(
          nma(pdbs_obj, ff = "calpha", fit = TRUE, full = TRUE, rm.gaps = FALSE, ncore = 1),
          error = function(e) {
            log_message("  ERROR nma.pdbs() convergence check for", cond, ":", e$message)
            NULL
          }
        )

        if (!is.null(enma_obj)) {
          # Save enma object for reference
          saveRDS(enma_obj, file.path(nma_data_dir, paste0(cond, "_enma_obj.rds")))
          log_message("  Saved enma_obj (convergence check) for:", cond)

          # Extract and save RMSIP matrix between replicas
          tryCatch({
            rmsip_mat <- enma_obj$rmsip
            if (!is.null(rmsip_mat)) {
              rmsip_df <- as.data.frame(rmsip_mat)
              colnames(rmsip_df) <- rep_name_list
              rownames(rmsip_df) <- rep_name_list
              write.csv(rmsip_df,
                        file.path(nma_data_dir, paste0(cond, "_rmsip_replicas.csv")),
                        row.names = TRUE)
              log_message("  RMSIP between replicas for", cond, ":")
              for (i in seq_len(nrow(rmsip_mat))) {
                for (j in seq_len(ncol(rmsip_mat))) {
                  if (j > i) {
                    log_message("    ", rep_name_list[i], "vs", rep_name_list[j],
                                "RMSIP =", round(rmsip_mat[i, j], 3))
                  }
                }
              }
            } else {
              log_message("  WARNING: RMSIP not found in enma_obj for", cond)
            }
          }, error = function(e) {
            log_message("  ERROR extracting RMSIP for", cond, ":", e$message)
          })

          # Save averaged DCCM from per-replica NMAs (informational complement)
          save_consolidated_dccm(cond, enma_obj, nma_data_dir)

          # Save consensus fluctuation profile (row-averaged across replicas)
          tryCatch({
            f_consolidated <- colMeans(enma_obj$fluctuations, na.rm = TRUE)
            saveRDS(f_consolidated, file.path(nma_data_dir, paste0(cond, "_f_consolidated.rds")))
            write.csv(
              data.frame(ca_index = seq_along(f_consolidated), fluct = f_consolidated),
              file.path(nma_data_dir, paste0(cond, "_f_consolidated.csv")),
              row.names = FALSE
            )
            log_message("  Saved consensus fluctuation profile (averaged across replicas):", cond)
          }, error = function(e) {
            log_message("  ERROR saving f_consolidated for", cond, ":", e$message)
          })

        } # end if enma_obj not null
      } # end if pdbs_obj not null
    } # end atom count check

  } else {
    log_message("  WARNING: fewer than 2 replicas completed NMA for", cond,
                "- skipping convergence check (need >= 2 for nma.pdbs).")
  }

  log_message("  Completed condition:", cond)
  flush_log_buffer <- function(append = TRUE) {} # no-op inside dopar; flushed in main loop
  return(list(cond = cond, status = "ok"))

} # end foreach

stopCluster(cl)

# -- Report results --
log_message("============================================================")
log_message("Pipeline finished. Results:")
for (res in results) {
  if (inherits(res, "error")) {
    log_message("  ERROR:", res$message)
  } else {
    log_message(" ", res$cond, "->", res$status)
  }
}

log_message("=== 02A v5.0 complete ===")
log_message("Data outputs saved under:", file.path(output_root, "nma_data"))
flush_log(append = FALSE)
