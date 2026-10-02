###############################################################################
# 07a_v2.1_build_distance_arrays.R
# Author: DenyCB
# PURPOSE:
# Build residue-residue distance arrays from MD trajectories of TRPV1 systems.
# Uses Cbeta atoms (except Glycine -> Calpha).
#
# INPUT FOLDER:
# C:/DinamicasMoleculares/TRPV1_pipeline/input/md
#
# OUTPUT FOLDER:
# C:/DinamicasMoleculares/TRPV1_pipeline/output/MD_distances/07a_build_arrays_v1.1
#
# OUTPUTS PER REPLICA:
#   *_mean_dist.rds          [res x res] mean distances in Angstrom
#   *_var_dist.rds           [res x res] variance of distances
#   *_contact4_freq.rds      [res x res] contact frequency <4A
#   *_contact8_freq.rds      [res x res] contact frequency <8A
#   *_contact12_freq.rds     [res x res] contact frequency <12A
#   *_meta.rds               residue labels / metadata
#
# FEATURES:
# - Automatic pdb+dcd pairing
# - Logs with timestamps
# - Optional frame trimming
# - Optional fitting/alignment
# - Parallel-ready (serial default for robustness)
# - Reduced RAM mode
#
###############################################################################

suppressPackageStartupMessages({
  library(bio3d)
  library(parallel)
})

###############################################################################
# USER SETTINGS
###############################################################################

INPUT_DIR  <- "C:/DinamicasMoleculares/TRPV1_pipeline/input/md"
OUTPUT_DIR <- "C:/DinamicasMoleculares/TRPV1_pipeline/output/MD_distances/07a_build_arrays_v2.2"

# Number of cores to use (recommended max 6 by user)
N_CORES <- 1

# Frames to use (NULL = all)
FRAME_RANGE <- 101:501

# Use every nth frame (1 = all, 2 = every second frame, etc.)
FRAME_STRIDE <- 1

# Residue range to keep
RES_RANGE <- c(300:602, 625:752)

# Optional structural fitting before extracting coords
DO_FIT <- FALSE

# Save contact arrays
SAVE_CONTACT_4  <- FALSE
SAVE_CONTACT_8  <- FALSE
SAVE_CONTACT_12 <- FALSE

# Overwrite existing outputs
OVERWRITE <- FALSE

# Parallel execution
USE_PARALLEL <- FALSE

###############################################################################
# CREATE OUTPUT FOLDER
###############################################################################

if (!dir.exists(OUTPUT_DIR)) dir.create(OUTPUT_DIR, recursive = TRUE)

###############################################################################
# LOGGING
###############################################################################

LOG_FILE <- file.path(
  OUTPUT_DIR,
  paste0("07a_v2.1_log_", format(Sys.time(), "%Y%m%d_%H%M%S"), ".txt")
)

log_msg <- function(...) {
  txt <- paste0("[", format(Sys.time(), "%Y-%m-%d %H:%M:%S"), "] ",
                paste(..., collapse=""))
  cat(txt, "\n")
  cat(txt, "\n", file = LOG_FILE, append = TRUE)
}

###############################################################################
# START
###############################################################################

log_msg("============================================================")
log_msg("STARTING SCRIPT 07a_v2.1_build_distance_arrays.R")
log_msg("Input folder: ", INPUT_DIR)
log_msg("Output folder: ", OUTPUT_DIR)
log_msg("Cores: ", N_CORES)
log_msg("DO_FIT: ", DO_FIT)
log_msg("FRAME_STRIDE: ", FRAME_STRIDE)

###############################################################################
# HELPERS
###############################################################################

select_rep_atoms <- function(pdb) {
  
  atm <- pdb$atom
  prot <- which(atm$type == "ATOM")
  atm <- atm[prot, ]
  
  key <- paste(atm$chain, atm$resno, sep="_")
  selected <- integer(0)
  ukeys <- unique(key)
  
  for (k in ukeys) {
    
    idx <- which(key == k)
    sub <- atm[idx, ]
    
    cb <- idx[sub$elety == "CB"]
    ca <- idx[sub$elety == "CA"]
    
    if (length(cb) > 0) {
      selected <- c(selected, cb[1])
    } else if (length(ca) > 0) {
      selected <- c(selected, ca[1])
    }
  }
  
  selected
}

make_labels <- function(atom_df) {
  paste0(atom_df$chain, atom_df$resno)
}

xyz_to_distmat <- function(v) {
  xyz <- matrix(v, ncol = 3, byrow = TRUE)
  as.matrix(dist(xyz))
}

###############################################################################
# FIND FILES
###############################################################################

pdb_files <- list.files(INPUT_DIR, pattern="\\.pdb$", full.names=TRUE)
dcd_files <- list.files(INPUT_DIR, pattern="\\.dcd$", full.names=TRUE)

if (length(pdb_files) == 0) stop("No PDB files found.")
if (length(dcd_files) == 0) stop("No DCD files found.")

base_pdb <- tools::file_path_sans_ext(basename(pdb_files))
base_dcd <- tools::file_path_sans_ext(basename(dcd_files))

common <- intersect(base_pdb, base_dcd)

log_msg("Detected systems with PDB+DCD pairs: ", length(common))

###############################################################################
# PROCESS ONE SYSTEM
###############################################################################

process_system <- function(sysname) {
  
  log_msg("------------------------------------------------------------")
  log_msg("Processing: ", sysname)
  
  out_file <- file.path(OUTPUT_DIR, paste0(sysname, "_mean_dist.rds"))
  
  if (file.exists(out_file) && !OVERWRITE) {
    log_msg("Skipping existing file.")
    return(NULL)
  }
  
  pdb_path <- file.path(INPUT_DIR, paste0(sysname, ".pdb"))
  dcd_path <- file.path(INPUT_DIR, paste0(sysname, ".dcd"))
  
  pdb <- read.pdb(pdb_path)
  trj <- read.dcd(dcd_path)
  
  nframes_total <- nrow(trj)
  
  log_msg("Frames in trajectory: ", nframes_total)
  
  frames <- if (is.null(FRAME_RANGE)) seq_len(nframes_total) else FRAME_RANGE
  frames <- frames[seq(1, length(frames), by = FRAME_STRIDE)]
  
  trj <- trj[frames, , drop = FALSE]
  
  log_msg("Frames used: ", length(frames))
  
  sel_atoms <- select_rep_atoms(pdb)
  
  atom_meta <- pdb$atom[sel_atoms, ]
  
  keep <- atom_meta$resno %in% RES_RANGE
  
  sel_atoms <- sel_atoms[keep]
  atom_meta <- atom_meta[keep, ]
  
  inds <- atom2xyz(sel_atoms)
  labels <- make_labels(atom_meta)
  
  nres <- length(sel_atoms)
  nfr  <- nrow(trj)
  
  log_msg("Representative residues selected: ", nres)
  
  if (DO_FIT) {
    log_msg("Applying fitting/alignment...")
    trj <- fit.xyz(
      fixed = pdb$xyz,
      mobile = trj,
      fixed.inds = inds,
      mobile.inds = inds
    )
  }
  
  mean_mat <- matrix(0, nres, nres)
  m2_mat   <- matrix(0, nres, nres)
  
  if (SAVE_CONTACT_4)  c4  <- matrix(0, nres, nres)
  if (SAVE_CONTACT_8)  c8  <- matrix(0, nres, nres)
  if (SAVE_CONTACT_12) c12 <- matrix(0, nres, nres)
  
  for (f in seq_len(nfr)) {
    
    if (f %% 25 == 0 || f == 1 || f == nfr)
      log_msg("Frame ", f, "/", nfr)
    
    xyz <- trj[f, inds]
    dm <- xyz_to_distmat(xyz)
    
    delta <- dm - mean_mat
    mean_mat <- mean_mat + delta / f
    delta2 <- dm - mean_mat
    m2_mat <- m2_mat + delta * delta2
    
    if (SAVE_CONTACT_4)  c4  <- c4  + (dm < 4)
    if (SAVE_CONTACT_8)  c8  <- c8  + (dm < 8)
    if (SAVE_CONTACT_12) c12 <- c12 + (dm < 12)
  }
  
  var_mat <- m2_mat / max(1, nfr - 1)
  
  mean_mat[lower.tri(mean_mat)] <- NA
  var_mat[lower.tri(var_mat)]   <- NA
  
  if (SAVE_CONTACT_4)  c4[lower.tri(c4)]   <- NA
  if (SAVE_CONTACT_8)  c8[lower.tri(c8)]   <- NA
  if (SAVE_CONTACT_12) c12[lower.tri(c12)] <- NA
  
  if (SAVE_CONTACT_4)  c4  <- c4 / nfr
  if (SAVE_CONTACT_8)  c8  <- c8 / nfr
  if (SAVE_CONTACT_12) c12 <- c12 / nfr
  
  dimnames(mean_mat) <- list(labels, labels)
  dimnames(var_mat)  <- list(labels, labels)
  
  saveRDS(mean_mat, file.path(OUTPUT_DIR, paste0(sysname, "_mean_dist.rds")))
  saveRDS(var_mat,  file.path(OUTPUT_DIR, paste0(sysname, "_var_dist.rds")))
  
  if (SAVE_CONTACT_4)
    saveRDS(c4, file.path(OUTPUT_DIR, paste0(sysname, "_contact4_freq.rds")))
  
  if (SAVE_CONTACT_8)
    saveRDS(c8, file.path(OUTPUT_DIR, paste0(sysname, "_contact8_freq.rds")))
  
  if (SAVE_CONTACT_12)
    saveRDS(c12, file.path(OUTPUT_DIR, paste0(sysname, "_contact12_freq.rds")))
  
  meta <- list(
    system = sysname,
    labels = labels,
    atom_meta = atom_meta,
    frames_used = frames,
    nres = nres,
    nframes = nfr,
    atom_definition = "CB except residues lacking CB -> CA"
  )
  
  saveRDS(meta, file.path(OUTPUT_DIR, paste0(sysname, "_meta.rds")))
  
  log_msg("Saved outputs for ", sysname)
}

###############################################################################
# RUN ALL
###############################################################################

if (!USE_PARALLEL) {
  
  for (nm in common) {
    tryCatch(
      process_system(nm),
      error = function(e) {
        log_msg("ERROR in ", nm, ": ", e$message)
      }
    )
  }
  
} else {
  
  cl <- makeCluster(N_CORES)
  
  clusterEvalQ(cl, library(bio3d))
  
  clusterExport(
    cl,
    varlist = c(
      "common","process_system","INPUT_DIR","OUTPUT_DIR","FRAME_RANGE",
      "FRAME_STRIDE","RES_RANGE","DO_FIT","SAVE_CONTACT_4",
      "SAVE_CONTACT_8","SAVE_CONTACT_12","OVERWRITE",
      "select_rep_atoms","make_labels","xyz_to_distmat","log_msg"
    ),
    envir = environment()
  )
  
  parLapply(cl, common, process_system)
  
  stopCluster(cl)
}

###############################################################################
# END
###############################################################################

log_msg("============================================================")
log_msg("SCRIPT FINISHED SUCCESSFULLY")
log_msg("End time: ", format(Sys.time(), "%Y-%m-%d %H:%M:%S"))
log_msg("============================================================")