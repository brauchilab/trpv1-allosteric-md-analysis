# ================================================================
# Script 03A_Deformation_v2.0
# Author: DenyCB
# Version: v2.0
# ================================================================

suppressPackageStartupMessages({
  library(bio3d)
  library(parallel)
  library(ggplot2)
})

# ---- PARAMETERS (UNCHANGED) ----
input_root   <- "C:/DinamicasMoleculares/analisis_bio3d_output/2_nma_correlation/nma_data"
output_root  <- "C:/DinamicasMoleculares/analisis_bio3d_output/3_domain_analysis_modes_2_31_37,38"
#output_root  <- "C:/DinamicasMoleculares/Pipeline_Complete/output/analisis_bio3d_output/3_domain_analysis_modes_41,43"
pdb_input_folder <- "C:/DinamicasMoleculares/input/PDBs/DM"

# ============================================================
# DEFINIR RANGO DE RESIDUOS A UTILIZAR EN TODO EL PIPELINE
# ============================================================
res_range <- c(277:602, 625:752)

global_ymax <- NULL
deformation_mode_inds <- c(2:31, 37,38)
#deformation_mode_inds <- c(41,43)

# parallel
ncores <- 6 #parallel::detectCores(logical = FALSE)
if(is.na(ncores) || ncores < 1) ncores <- 1

# ensure same directories as v1.3
folders <- c("Deformation_plots","Deformation_pdb","rds","csv","logs","debug")
for (f in folders) dir.create(file.path(output_root, f), recursive = TRUE, showWarnings = FALSE)

# -------------------------
# Directorio adicional para PDB full-atom deformation
# -------------------------
out_def_dir <- file.path(output_root, "Deformation_pdb_full")
if(!dir.exists(out_def_dir)) dir.create(out_def_dir, recursive = TRUE, showWarnings = FALSE)

# log (per-run file)
log_file <- file.path(output_root, "logs", paste0("03A_deformation_log_v2.0_", format(Sys.time(), "%Y%m%d_%H%M%S"), ".txt"))
log_message <- function(...) {
  msg <- paste0(format(Sys.time(), "%Y-%m-%d %H:%M:%OS3"), " - ", paste(..., collapse = " "))
  cat(msg, "\n")
  cat(msg, "\n", file = log_file, append = TRUE)
}

log_message("=== Starting Script 03A v2.0: Deformation analysis ===")
log_message("Input root:", input_root)
log_message("Output root:", output_root)
log_message("PDB input folder:", pdb_input_folder)
log_message("Using deformation_mode_inds:", if(is.null(deformation_mode_inds)) "NULL (default)" else paste(deformation_mode_inds, collapse=","))
log_message("NOTE: deformation_mode_inds are interpreted as non-trivial mode indices matching 02C overlap plots; +6 is added internally.")

run_and_capture <- function(expr, capture_file = NULL, envir = parent.frame()) {
  # Evaluates 'expr' in 'envir', captures stdout+messages into capture_file (if given),
  # captures warnings and error, and returns list(value, warnings, error).
  zz <- NULL
  wlist <- character()
  err_msg <- NULL
  val <- NULL
  
  if (!is.null(capture_file)) {
    dir.create(dirname(capture_file), recursive = TRUE, showWarnings = FALSE)
    zz <- file(capture_file, open = "wt")
  }
  
  # Ensure sinks are always closed even on error
  if (!is.null(zz)) {
    sink(zz)
    sink(zz, type = "message")
    on.exit({
      try(sink(type = "message"), silent = TRUE)
      try(sink(), silent = TRUE)
      try(close(zz), silent = TRUE)
    }, add = TRUE)
  }
  
  # Evaluate with warning handler
  val <- tryCatch(
    withCallingHandlers(
      eval(expr, envir = envir),
      warning = function(w) {
        wlist <<- c(wlist, conditionMessage(w))
        invokeRestart("muffleWarning")
      }
    ),
    error = function(e) {
      err_msg <<- conditionMessage(e)
      NULL
    }
  )
  
  return(list(value = val, warnings = wlist, error = err_msg))
}



# map residue->atom b (unchanged)
map_residue_to_atom_b <- function(pdb, resid_vals, resid_numbers = NULL) {
  if(is.null(pdb$atom)) stop("pdb object lacks atom table")
  n_atoms <- nrow(pdb$atom)
  b_out <- rep(NA_real_, n_atoms)
  if(!is.null(resid_numbers) && length(resid_numbers) == length(resid_vals)) {
    for(i in seq_along(resid_vals)) {
      rn <- resid_numbers[i]; idxs <- which(pdb$atom$resno == rn)
      if(length(idxs)>0) b_out[idxs] <- resid_vals[i]
    }
    return(b_out)
  }
  uniq_res <- sort(unique(pdb$atom$resno))
  if(length(resid_vals) == length(uniq_res)) {
    for(i in seq_along(uniq_res)) {
      rn <- uniq_res[i]; b_out[pdb$atom$resno == rn] <- resid_vals[i]
    }
    return(b_out)
  }
  ca_sel <- atom.select(pdb, "calpha")
  if(length(resid_vals) == length(ca_sel$atom)) {
    for(i in seq_along(ca_sel$atom)) {
      ai <- ca_sel$atom[i]; rn <- pdb$atom$resno[ai]; b_out[pdb$atom$resno == rn] <- resid_vals[i]
    }
    return(b_out)
  }
  warning("map_residue_to_atom_b: No mapping coincide en longitud; devolviendo NA vector")
  return(b_out)
}

# minimal mapping heuristics using nma/def objects (unchanged pattern)
find_mapping_indices_minimal <- function(nma_obj = NULL, def_obj = NULL, expected_len = NULL, pdb_obj = NULL) {
  log_message("  [mapping] try minimal heuristics expected_len =", expected_len)
  if(!is.null(nma_obj$inds)) {
    inds <- nma_obj$inds
    if(is.numeric(inds) && !is.null(expected_len) && length(inds) == expected_len) {
      log_message("  [mapping] matched nma_obj$inds (numeric)")
      return(as.integer(inds))
    }
    if(is.list(inds)) {
      flat <- unlist(inds)
      if(!is.null(expected_len) && length(flat) == expected_len) {
        log_message("  [mapping] matched flattened nma_obj$inds list")
        return(as.integer(flat))
      }
    }
  }
  if(!is.null(def_obj)) {
    if(!is.null(def_obj$inds) && is.numeric(def_obj$inds) && !is.null(expected_len) && length(def_obj$inds) == expected_len) {
      log_message("  [mapping] used def_obj$inds")
      return(as.integer(def_obj$inds))
    }
    if(!is.null(def_obj$atom) && is.numeric(def_obj$atom) && !is.null(expected_len) && length(def_obj$atom) == expected_len) {
      log_message("  [mapping] used def_obj$atom")
      return(as.integer(def_obj$atom))
    }
  }
  # try parse rownames of def_obj$ei
  if(!is.null(def_obj) && !is.null(def_obj$ei)) {
    rn <- rownames(def_obj$ei)
    if(!is.null(rn) && length(rn) == expected_len) {
      parsed <- suppressWarnings(as.integer(gsub("[^0-9]", "", rn)))
      if(all(!is.na(parsed))) {
        if(!is.null(pdb_obj)) {
          map_inds <- integer(length(parsed))
          for(i in seq_along(parsed)) {
            cand <- which(pdb_obj$atom$resno == parsed[i])
            map_inds[i] <- if(length(cand)>0) cand[1] else NA_integer_
          }
          if(sum(is.na(map_inds)) == 0) { log_message("  [mapping] used parsed rownames(def_obj$ei)"); return(map_inds) }
        }
      }
    }
  }
  log_message("  [mapping] minimal heuristics failed")
  return(NULL)
}

# Read conditions (same pattern as before)
conds_files <- list.files(input_root, pattern = "_nma_obj\\.rds$", full.names = FALSE)
conds <- gsub("_nma_obj\\.rds$", "", conds_files)
log_message("Found", length(conds), "conditions:", paste(conds, collapse = ", "))

# processing function (almost identical to v1.9, with targeted fixes)
process_one <- function(cond) {
  tryCatch({
    log_message("------------------------------------------------------------")
    log_message("START condition:", cond)
    t0 <- Sys.time()
    nma_path <- file.path(input_root, paste0(cond, "_nma_obj.rds"))
    if(!file.exists(nma_path)) { log_message("nma_obj missing:", nma_path); return(NULL) }
    nma_obj <- readRDS(nma_path)
    log_message("Loaded nma_obj class:", paste(class(nma_obj), collapse = ","))
    if(!is.null(nma_obj$U)) log_message("nma_obj$U dims:", paste(dim(nma_obj$U), collapse = "x"))
    
    # load pdb: prefer trimmed RDS from nma_data (input_root); fallback to raw .pdb in pdb_input_folder
    pdb_obj <- NULL
    trimmed_candidate <- file.path(input_root, paste0(cond, "_pdb_ca_trimmed.rds"))
    if (file.exists(trimmed_candidate)) {
      pdb_obj <- tryCatch(readRDS(trimmed_candidate), error = function(e) { log_message("ERROR reading trimmed RDS:", e$message); NULL })
      if(!is.null(pdb_obj)) log_message("Loaded trimmed PDB RDS from", trimmed_candidate, "atoms:", if(!is.null(pdb_obj$atom)) nrow(pdb_obj$atom) else "NA")
    } else {
      pdb_candidates <- list.files(pdb_input_folder, pattern = paste0("^", cond, ".*\\.pdb$"), full.names = TRUE, ignore.case = TRUE)
      if(length(pdb_candidates) >= 1) {
        pdb_obj <- tryCatch(read.pdb(pdb_candidates[1]), error = function(e) { log_message("ERROR reading PDB:", e$message); NULL })
        if(!is.null(pdb_obj)) log_message("Loaded PDB from", pdb_candidates[1], "atoms:", if(!is.null(pdb_obj$atom)) nrow(pdb_obj$atom) else "NA")
      } else {
        log_message("No pdb candidate found for", cond)
      }
    }
    
    # call deformation.nma
    log_message("[DEFORMATION] start for", cond)
    capture_def <- file.path(output_root, "debug", paste0(cond, "_deformation_call_output.txt"))
    
    # --- FIX CRÍTICO ---
    # La llamada correcta es deformation.nma(x = nma_obj, mode.inds = ...)
    expr_def <- if (is.null(deformation_mode_inds)) {
      substitute(deformation.nma(NMA), list(NMA = as.name("nma_obj")))
    } else {
      deformation_mode_inds_real <- deformation_mode_inds + 6
      substitute(deformation.nma(NMA, mode.inds = MO),
                 list(NMA = as.name("nma_obj"), MO = deformation_mode_inds_real))
    }
    
    # --- FIN FIX ---
    
    def_run <- run_and_capture(expr_def,
                               capture_file = capture_def,
                               envir = list2env(list(nma_obj = nma_obj)))
    
    
    if(!is.null(def_run$error)) {
      log_message("[DEFORMATION] ERROR calling deformation.nma:", def_run$error)
      saveRDS(
        list(error = def_run$error,
             out = tryCatch(readLines(capture_def, warn = FALSE), error = function(e) NULL)),
        file.path(output_root, "debug", paste0(cond, "_deformation_error.rds"))
      )
    } else {
      def_obj <- def_run$value
      saveRDS(def_obj, file.path(output_root, "rds", paste0(cond, "_deformation.rds")))
      log_message("[DEFORMATION] saved rds for", cond)
      
      if(is.null(def_obj$ei)) {
        log_message("[DEFORMATION] def_obj$ei is NULL for", cond, "- saving debug and skipping")
        saveRDS(def_obj, file.path(output_root, "debug", paste0(cond, "_def_obj_no_ei.rds")))
      } else {
        ei_mat <- def_obj$ei
        energy_tot <- rowSums(ei_mat, na.rm = TRUE)
        n_energy <- length(energy_tot)
        if(!is.null(pdb_obj)) {
          n_atoms <- nrow(pdb_obj$atom)
          uniq_resnos <- sort(unique(pdb_obj$atom$resno))
          n_residues_unique <- length(uniq_resnos)
          ca_sel <- atom.select(pdb_obj, "calpha")
          n_ca <- length(ca_sel$atom)
        } else {
          n_atoms <- NA; n_residues_unique <- NA; n_ca <- NA
        }
        log_message("[DEFORMATION] lengths -> energy_tot:", n_energy,
                    "n_atoms:", n_atoms, "n_residues_unique:", n_residues_unique,
                    "n_ca:", n_ca)
        
        df_res <- NULL
        used_method <- NA_character_
        
        # Method 1: exact match to unique residues (original v1.3 path)
        if(!is.na(n_residues_unique) && n_energy == n_residues_unique) {
          log_message("[DEFORMATION] mapping method: exact unique-residue (lengths equal)")
          df_res <- data.frame(resno = uniq_resnos, deformation_energy = as.numeric(energy_tot))
          used_method <- "unique_res_exact"
        }
        
        # Method 2: energy per CA
        if(is.null(df_res) && !is.na(n_ca) && n_energy == n_ca && !is.null(pdb_obj)) {
          log_message("[DEFORMATION] mapping method: direct CA (energy vector matches CA count)")
          ca_inds <- ca_sel$atom
          ca_resnos <- pdb_obj$atom$resno[ca_inds]
          df_res <- data.frame(resno = ca_resnos, deformation_energy = as.numeric(energy_tot))
          # collapse duplicates (if multiple CA per resno) by mean
          df_res <- aggregate(deformation_energy ~ resno, data = df_res, FUN = mean)
          used_method <- "direct_CA"
        }
        
        # Method 3: use nma_obj$inds or def_obj$inds if present and matching
        if(is.null(df_res)) {
          mapping_inds <- find_mapping_indices_minimal(
            nma_obj = nma_obj,
            def_obj = def_obj,
            expected_len = n_energy,
            pdb_obj = pdb_obj
          )
          if(!is.null(mapping_inds) && !is.null(pdb_obj) && length(mapping_inds) == n_energy) {
            log_message("[DEFORMATION] mapping method: nma/def indices mapping found")
            df_atom <- data.frame(resno = pdb_obj$atom$resno[mapping_inds],
                                  energy = energy_tot)
            df_res <- aggregate(energy ~ resno, data = df_atom, FUN = mean, na.rm = TRUE)
            colnames(df_res) <- c("resno","deformation_energy")
            used_method <- "nma_def_inds"
          }
        }
        
        # Method 4: group-resample fallback REMOVED.
        # If no method succeeded, stop this condition with an explicit error.
        # This prevents silent production of incorrect output.
        # The outer tryCatch will catch this stop() and log it as a critical error.
        if(is.null(df_res)) {
          debug_info <- list(
            energy_tot_len    = n_energy,
            n_atoms           = n_atoms,
            n_residues_unique = n_residues_unique,
            n_ca              = n_ca,
            sample_energy     = head(energy_tot, 80)
          )
          saveRDS(debug_info,
                  file.path(output_root, "debug",
                            paste0(cond, "_deformation_length_mismatch_debug_v2.0.rds")))
          stop(paste0(
            "CRITICAL: no mapping method succeeded for ", cond,
            ". energy_tot length=", n_energy,
            ", n_residues_unique=", n_residues_unique,
            ", n_ca=", n_ca,
            ". Debug info saved. Stopping condition to prevent incorrect output."
          ))
        }

        # final sanity checks: if many NAs or all zeros -> mark WARN
        na_frac <- sum(is.na(df_res$deformation_energy))/nrow(df_res)
        zero_frac <- sum(df_res$deformation_energy == 0, na.rm = TRUE)/nrow(df_res)
        warn_flag <- (na_frac > 0.5) || (zero_frac > 0.9)
        out_prefix <- if(warn_flag) "WARN_" else ""
        
        csv_out <- file.path(output_root, "csv", paste0(out_prefix, cond, "_deformation_residuewise.csv"))
        png_out <- file.path(output_root, "Deformation_plots", paste0(out_prefix, cond, "_deformation_residuewise.png"))
        pdb_out <- file.path(output_root, "Deformation_pdb", paste0(out_prefix, cond, "_deformation.pdb"))
        
        write.csv(df_res, csv_out, row.names = FALSE)
        log_message("[DEFORMATION] wrote CSV:", csv_out, "mapping_used:", used_method, "na_frac:", round(na_frac,3), "zero_frac:", round(zero_frac,3))
        
        # plot (resno on x)
        png(png_out, width = 1400, height = 700, res = 150)
        par(mar = c(6,6,4,2))
        plot(df_res$resno, df_res$deformation_energy, type = "h",
             main = paste0(cond, " – Deformation per residue [", used_method, "]"),
             xlab = "Residue (resno)", ylab = "Deformation energy (a.u.)")
        dev.off()
        log_message("[DEFORMATION] saved PNG:", png_out)
        
        # map to pdb atoms (b-factors) and write pdb
        if(!is.null(pdb_obj)) {
          bvals <- map_residue_to_atom_b(pdb_obj, df_res$deformation_energy, df_res$resno)
          if(!all(is.na(bvals))) {
            tryCatch({
              write.pdb(pdb_obj, file = pdb_out, b = bvals)
              log_message("[DEFORMATION] wrote deformation-colored PDB:", pdb_out)
            }, error = function(e) {
              log_message("[DEFORMATION] ERROR writing PDB:", e$message)
              saveRDS(list(pdb_head = head(pdb_obj$atom,50), head_b = head(bvals,100)), file.path(output_root, "debug", paste0(cond, "_deformation_pdb_write_error_v2.0.rds")))
            })
          } else {
            log_message("[DEFORMATION] mapping produced all NA b-values -> saved CSV/PNG but skipping PDB write")
          }

          
          # ============================================================
          # FULL-ATOM deformation PDB (todos los átomos)
          # Respeta res_range definido globalmente
          # ============================================================
          
          # --- identificar PDB original según prefijo ---
          pdb_prefix <- substr(cond, 1, 4)
          pdb_name_map <- c("7LP9"="7LP9.pdb", "7LPB"="7LPB.pdb", "7LPC"="7LPC.pdb")
          original_pdb_path <- if (pdb_prefix %in% names(pdb_name_map))
            file.path("D:/EstructurasTODO/TRPV/TRPV1/2021", pdb_name_map[pdb_prefix]) else NA
          
          if (!is.na(original_pdb_path) && file.exists(original_pdb_path)) {
            pdb_full <- tryCatch(read.pdb(original_pdb_path), error=function(e){
              log_message("[DEFORMATION] ERROR reading full PDB:", e$message); NULL
            })
            if (!is.null(pdb_full))
              log_message("[DEFORMATION] Loaded original full-atom PDB:", original_pdb_path)
          } else {
            pdb_full <- NULL
            log_message("[DEFORMATION] ERROR: cannot find full-atom PDB for prefix:", pdb_prefix)
          }
          
          if (!is.null(df_res) && !is.null(pdb_full)) {
            
            # --------------------------
            # 1) Filtrar df_res al rango
            # --------------------------
            df_res_ranged <- df_res[df_res$resno %in% res_range, ]

            if (nrow(df_res_ranged) > 0) {

              # ----------------------------------------------
              # 2) Filtrar full PDB a res_range (sub-PDB real)
              # ----------------------------------------------
              sel_full <- atom.select(pdb_full, resno = res_range)
              pdb_full_sub <- trim.pdb(pdb_full, sel_full)
              
              log_message("[DEFORMATION] Full-atom trimmed to res_range. Atoms:",
                          nrow(pdb_full_sub$atom))
              
              # -------------------------------------------------------
              # 3) Intento directo: asignar b por resno coincidente
              # -------------------------------------------------------
              res_in_sub <- as.character(pdb_full_sub$atom$resno)
              common <- intersect(as.character(df_res_ranged$resno), res_in_sub)
              
              if (length(common) > 0) {
                log_message("[DEFORMATION] Direct resno mapping available. Matches:", length(common))
                
                def_map <- setNames(df_res_ranged$deformation_energy,
                                    as.character(df_res_ranged$resno))
                bfac <- def_map[as.character(pdb_full_sub$atom$resno)]
                bfac[is.na(bfac)] <- 0
                pdb_full_sub$atom$b <- as.numeric(bfac)
                
              } else {
                # -----------------------------------------------------
                # 4) Fallback: mapear por orden de CA (muy robusto)
                # -----------------------------------------------------
                log_message("[DEFORMATION] No direct matches — using ordered CA fallback.")
                
                # CA en trimmed original (pdb_obj) pero filtrado por rango
                ca_trim <- which(pdb_obj$atom$elety == "CA" & pdb_obj$atom$resno %in% res_range)
                ca_trim_resno <- pdb_obj$atom$resno[ca_trim]
                
                # CA en full-atom sub-PDB
                ca_full <- which(pdb_full_sub$atom$elety == "CA")
                ca_full_resno <- pdb_full_sub$atom$resno[ca_full]
                
                if (length(ca_trim_resno) == 0 || length(ca_full_resno) < length(ca_trim_resno)) {
                  log_message("[DEFORMATION] FATAL: CA fallback impossible. Assigning b=0.")
                  pdb_full_sub$atom$b <- 0
                } else {
                  map_len <- length(ca_trim_resno)
                  mapped_full_res <- ca_full_resno[seq_len(map_len)]
                  
                  # mapping trimmed_resno -> full_resno
                  map_tf <- setNames(as.character(mapped_full_res),
                                     as.character(ca_trim_resno))
                  
                  # construir nuevo df_res usando resno_full
                  df_res_ranged$resno_full <- map_tf[as.character(df_res_ranged$resno)]
                  
                  def_map_full <- setNames(df_res_ranged$deformation_energy,
                                           df_res_ranged$resno_full)
                  
                  bfac <- def_map_full[as.character(pdb_full_sub$atom$resno)]
                  bfac[is.na(bfac)] <- 0
                  pdb_full_sub$atom$b <- as.numeric(bfac)
                }
              }
              
              # -----------------------------------------
              # 5) Guardar PDB full-atom FINAL
              # -----------------------------------------
              out_def_pdb_full <- file.path(out_def_dir,
                                            paste0(cond, "_deformation_fullatom_ranged.pdb"))
              
              tryCatch({
                write.pdb(pdb_full_sub, file = out_def_pdb_full)
                log_message("[DEFORMATION] Saved RANGED FULL-ATOM PDB:", out_def_pdb_full)
              }, error=function(e){
                log_message("[DEFORMATION] ERROR writing ranged full-atom PDB:", e$message)
              })

            } else {
              log_message("[DEFORMATION] WARNING: df_res has no residues in res_range -> skipping full-atom PDB")
            }

          } else {
            log_message("[DEFORMATION] Skipped full-atom block: df_res or pdb_full missing.")
          }
          
          
        } else {
          log_message("[DEFORMATION] pdb_obj missing -> CSV/PNG written but cannot write CA-PDB or full-atom PDB")
        }
      } # end def_obj$ei
    } # end def_run success
    log_message("FINISHED condition:", cond, "elapsed(s):",
                round(as.numeric(difftime(Sys.time(), t0, units = "secs")), 2))
    return(TRUE)
  }, error = function(e) {
    log_message("ERROR processing", cond, ":", e$message)
    saveRDS(list(error = e$message),
            file.path(output_root, "debug",
                      paste0(cond, "_deformation_uncaught_error_v2.0.rds")))
    return(FALSE)
  })
}

# run parallel or sequential as before
if(ncores > 1 && length(conds) > 1) {
  ncores_use <- min(length(conds), ncores)
  log_message("Running in parallel with PSOCK cluster. Workers:", ncores_use)
  cl <- makeCluster(ncores_use)
  clusterEvalQ(cl, {
    library(bio3d)
    library(parallel)
    library(ggplot2)
  })
  clusterExport(cl, ls(), envir = environment())
  parLapply(cl, conds, process_one)
  stopCluster(cl)
} else {
  log_message("Running sequentially.")
  lapply(conds, process_one)
}

# ============================================================
# Replot deformation PNGs using common global Y scale
# ============================================================
replot_deformation_common_y <- function() {
  csv_files <- list.files(file.path(output_root, "csv"),
                          pattern = "_deformation_residuewise\\.csv$",
                          full.names = TRUE)
  
  if(length(csv_files) == 0) {
    log_message("[REPLOT] No deformation CSVs found. Skipping common Y-scale replot.")
    return(invisible(FALSE))
  }
  
  global_ymax <- 0
  for(cf in csv_files) {
    df_tmp <- tryCatch(read.csv(cf), error = function(e) NULL)
    if(!is.null(df_tmp) && "deformation_energy" %in% names(df_tmp)) {
      vals <- df_tmp$deformation_energy
      vals <- vals[is.finite(vals)]
      if(length(vals) > 0) global_ymax <- max(global_ymax, max(vals, na.rm = TRUE), na.rm = TRUE)
    }
  }
  
  if(!is.finite(global_ymax) || global_ymax <= 0) {
    log_message("[REPLOT] Invalid global_ymax. Skipping common Y-scale replot.")
    return(invisible(FALSE))
  }
  
  log_message("[REPLOT] Global deformation ymax from CSVs:", global_ymax)
  
  for(cf in csv_files) {
    df_res <- tryCatch(read.csv(cf), error = function(e) NULL)
    if(is.null(df_res) || !("resno" %in% names(df_res)) || !("deformation_energy" %in% names(df_res))) next
    
    cond_name <- basename(cf)
    cond_name <- gsub("_deformation_residuewise\\.csv$", "", cond_name)
    
    png_out <- file.path(output_root, "Deformation_plots",
                         paste0(cond_name, "_deformation_residuewise.png"))
    
    png(png_out, width = 1400, height = 700, res = 150)
    par(mar = c(6,6,4,2))
    plot(df_res$resno, df_res$deformation_energy, type = "h",
         ylim = c(0, global_ymax),
         main = paste0(cond_name, " – Deformation per residue"),
         xlab = "Residue (resno)", ylab = "Deformation energy (a.u.)")
    dev.off()
    
    log_message("[REPLOT] saved common-scale PNG:", png_out)
  }
  
  return(invisible(TRUE))
}

replot_deformation_common_y()

log_message("=== Script 03A v2.0 complete ===")
