###############################################################################
# 07c_v2.1_make_distograms_from_saved_matrices.R
# Author: DenyCB
# PURPOSE:
# Combine the 3 replicas of each MD system from saved matrices and generate
# distograms for:
#   1) Mean distance
#   2) Normalized variance to distance: var / mean^2
#
# Additionally generates:
#   3) Delta distograms: (Cap - Apo) and (Heat - Apo) per mutant,
#      for both mean distance and normalized variance matrices
#   4) Delta-delta distograms: [mutant delta] - [WT delta],
#      for Cap-Apo and Heat-Apo comparisons, for both matrix types
#
# INPUT FOLDER:
# C:/DinamicasMoleculares/TRPV1_pipeline/output/MD_distances/07a_build_arrays_v2.2
#
# OUTPUT FOLDER:
# C:/DinamicasMoleculares/TRPV1_pipeline/output/MD_distances/07c_distograms_v2.1
#
# REQUIRED INPUTS PER REPLICA:
#   *_mean_dist.rds
#   *_var_dist.rds
#   *_meta.rds
#
# FEATURES:
# - Automatic grouping of replicas by system
# - Pooled combination of replica means and variances
# - Delta matrices: Cap-Apo and Heat-Apo per mutant (mean and normvar)
# - Delta-delta matrices: mutant delta vs WT delta (mean and normvar)
# - Loop over multiple zoom ranges (ZOOM_RANGES), one output folder per range
# - Parallelization over systems within each range (parallel, Win10 compatible)
# - Logs with timestamps
# - PNG outputs
###############################################################################

suppressPackageStartupMessages({
  library(viridisLite)
  library(fields)
  library(parallel)
})

###############################################################################
# USER SETTINGS
###############################################################################

INPUT_DIR  <- "C:/DinamicasMoleculares_Prev/TRPV1_pipeline/output/MD_distances/07a_build_arrays_v2.2"
OUTPUT_DIR <- "C:/DinamicasMoleculares/TRPV1_pipeline/output/07_MD_distances/07c_distograms_v2.1_2"

# Multiple zoom ranges to process in a single run.
# Each entry generates its own subfolder inside OUTPUT_DIR.
# Format: list(X = c("start_label", "end_label"), Y = c("start_label", "end_label"))
ZOOM_RANGES <- list(
  list(X = c("A300","D752"), Y = c("A300","D752")),
  list(X = c("A300","A752"), Y = c("D300","D752")),
  list(X = c("B300","C752"), Y = c("B300","C752")),
  list(X = c("A300","A752"), Y = c("A300","A752")),
  list(X = c("A300","A752"), Y = c("B300","C752")),
  list(X = c("B300","C752"), Y = c("D300","D752")),
  list(X = c("D300","D752"), Y = c("D300","D752")),
  list(X = c("A400","A590"), Y = c("A400","A590")),
  list(X = c("A400","A560"), Y = c("D560","D690"))
)

# PNG settings
PNG_W   <- 3600
PNG_H   <- 3200
PNG_DPI <- 600

RUN_PER_SYSTEM  <- TRUE
RUN_DELTA       <- TRUE
RUN_DELTADELTA  <- TRUE


# Number of colors in heatmaps
N_COLS <- 256

# Small epsilon to avoid division by zero in normalized variance
EPSILON <- 1e-12

# Path to residue map CSV (used to enforce canonical chain order ABCD)
RESIDUE_MAP <- "C:/DinamicasMoleculares/TRPV1_pipeline/output/residue_map_MD.csv"

# Overwrite existing outputs
OVERWRITE <- TRUE

# Number of cores for parallelization (per-system loop within each range)
N_CORES <- 6

# ===============================================================
# residuos especiales (idénticos al original)
# ===============================================================

special_res1 <- c("A426","A697","A441","B426","B697","B441","C426","C697","C441","D426","D697","D441","A563","B563","C563","D563")

special_res2 <- c("A679","A643","B679","B643","C679","C643","D679","D643")

special_res3 <- c("D593","A543","B543","A593")#"A434", "A449" )#"A505","B505","C505","D505","A409","B409","C409","D409","A690","B690","C690","D690")

special_res4 <- c("B300","C300","D300")#,"A593")

###############################################################################
# CREATE OUTPUT FOLDER
###############################################################################

if (!dir.exists(OUTPUT_DIR)) dir.create(OUTPUT_DIR, recursive = TRUE)

###############################################################################
# LOGGING
###############################################################################

LOG_FILE <- file.path(
  OUTPUT_DIR,
  paste0("07c_v2.1_log_", format(Sys.time(), "%Y%m%d_%H%M%S"), ".txt")
)

log_msg <- function(...) {
  txt <- paste0("[", format(Sys.time(), "%Y-%m-%d %H:%M:%S"), "] ", paste(..., collapse=""))
  cat(txt, "\n")
  cat(txt, "\n", file = LOG_FILE, append = TRUE)
}

###############################################################################
# START
###############################################################################

log_msg("============================================================")
log_msg("STARTING SCRIPT 07c_v2.1_make_distograms_from_saved_matrices.R")
log_msg("Input folder: ", INPUT_DIR)
log_msg("Output folder: ", OUTPUT_DIR)

###############################################################################
# HELPERS
###############################################################################

# Keep matrix labels between start and end labels (inclusive)
subset_nodes <- function(node_names, range_vec) {
  
  if (is.null(range_vec) || length(range_vec) != 2) {
    return(seq_along(node_names))
  }
  
  i1 <- match(range_vec[1], node_names)
  i2 <- match(range_vec[2], node_names)
  
  if (is.na(i1) || is.na(i2)) {
    stop(paste0("Range labels not found: ", paste(range_vec, collapse = " to ")))
  }
  
  seq(min(i1, i2), max(i1, i2))
}

# Rebuild full symmetric matrix from upper-triangle+diagonal matrix
upper_to_full <- function(mat) {
  full <- mat
  full[lower.tri(full)] <- t(full)[lower.tri(full)]
  full
}

# Pooled combination of means and variances across replicas
combine_replica_stats <- function(mean_list, var_list, n_list) {
  
  if (length(mean_list) == 0) stop("Empty mean_list.")
  if (length(mean_list) != length(var_list) || length(mean_list) != length(n_list)) {
    stop("mean_list, var_list and n_list must have same length.")
  }
  
  total_n <- sum(n_list)
  
  # weighted mean across replicas
  mean_num <- mean_list[[1]] * n_list[1]
  if (length(mean_list) > 1) {
    for (i in 2:length(mean_list)) {
      mean_num <- mean_num + mean_list[[i]] * n_list[i]
    }
  }
  mean_pool <- mean_num / total_n
  
  # pooled variance across replicas
  ss_total <- matrix(0, nrow(mean_pool), ncol(mean_pool))
  
  for (i in seq_along(mean_list)) {
    mu_i  <- mean_list[[i]]
    var_i <- var_list[[i]]
    n_i   <- n_list[i]
    
    ss_within  <- (n_i - 1) * var_i
    ss_between <- n_i * (mu_i - mean_pool)^2
    
    ss_total <- ss_total + ss_within + ss_between
  }
  
  var_pool <- ss_total / max(1, total_n - 1)
  
  list(mean = mean_pool, var = var_pool, n = total_n)
}

# Simple classifier for pretty labels
classify_system <- function(sysname) {
  
  parts <- strsplit(sysname, "_")[[1]]
  
  pdb <- parts[1]
  mut <- parts[2]
  
  cond <- ifelse(pdb == "7LP9", "Apo",
                 ifelse(pdb == "7LPB", "Cap", "Heat"))
  
  list(pdb = pdb, mut = mut, cond = cond)
}

# Plot one heatmap
plot_heatmap <- function(mat, labels_x, labels_y, outfile, main_txt, legend_txt, invert_col = FALSE, zlim_fixed = NULL) {
  
  png(outfile, width = PNG_W, height = PNG_H, res = PNG_DPI)
  
  old_par <- par(no.readonly = TRUE)
  on.exit({
    par(old_par)
    dev.off()
  }, add = TRUE)
  
  if (!is.null(zlim_fixed)) {
    zlim <- zlim_fixed
  } else {
    zlim <- range(mat, na.rm = TRUE)
    if (!all(is.finite(zlim))) zlim <- c(0, 1)
    if (diff(zlim) == 0) zlim <- zlim + c(-0.5, 0.5)
  }
  
  nx <- length(labels_x)
  ny <- length(labels_y)
  
  STEP_TICKS <- 10   # <-- lo puedes mover a USER SETTINGS si quieres
  
  tick_x <- seq(1, nx, by = STEP_TICKS)
  tick_y <- seq(1, ny, by = STEP_TICKS)
  
  par(mar = c(4, 4, 4, 4))
  
  fields::image.plot(
    x = seq_len(ncol(mat)),
    y = seq_len(nrow(mat)),
    z = t(mat[nrow(mat):1, , drop = FALSE]),
    col = if (invert_col) rev(viridis(N_COLS)) else colorRampPalette(c("white", "blue", "red"))(N_COLS),
    axes = FALSE,
    xlab = "",
    ylab = "",
    main = main_txt,
    zlim = zlim,
    legend.lab = legend_txt,
    legend.line = 2.5
  )
  # ===============================================================
  # GRID (líneas cada STEP_TICKS)
  # ===============================================================
  
  STEP_TICKS <- 20
  
  nx <- ncol(mat)
  ny <- nrow(mat)
  
  #  abline(v = seq(1, nx, by = STEP_TICKS), col = "grey85", lwd = 1)
  # abline(h = ny - seq(1, ny, by = STEP_TICKS) + 1, col = "grey85", lwd = 1)
  
  
  # ===============================================================
  # RESIDUOS ESPECIALES
  # ===============================================================
  
  for(v in list(
    list(special_res1,"orange"),
    list(special_res2,"magenta"),
    list(special_res3,"darkgray"),
    list(special_res4,"black")
  )){
    
    idx_x <- match(v[[1]], colnames(mat))
    idx_x <- idx_x[!is.na(idx_x)]
    
    idx_y <- match(v[[1]], rownames(mat))
    idx_y <- idx_y[!is.na(idx_y)]
    
    if(length(idx_x) > 0){
      abline(
        v   = idx_x,
        col = v[[2]],
        lwd = 0.5,
        lty = 2
      )
    }
    
    if(length(idx_y) > 0){
      abline(
        h   = ny - idx_y + 1,
        col = v[[2]],
        lwd = 0.5,
        lty = 2
      )
    }
  }
  
  # abline(v = seq(1, nx, by = STEP_TICKS), col = "grey85", lwd = 1)
  # abline(h = ny - seq(1, ny, by = STEP_TICKS) + 1, col = "grey85", lwd = 1)
  
  axis(1, at = tick_x, labels = labels_x[tick_x], las = 2, cex.axis = 0.4)
  axis(2, at = tick_y, labels = rev(labels_y)[tick_y], las = 2, cex.axis = 0.4)
  
  box()
}

# Plot one delta/delta-delta heatmap (red = negative, gray = zero, blue = positive)
plot_heatmap_delta <- function(mat, labels_x, labels_y, outfile, main_txt, legend_txt) {
  
  png(outfile, width = PNG_W, height = PNG_H, res = PNG_DPI)
  
  old_par <- par(no.readonly = TRUE)
  on.exit({
    par(old_par)
    dev.off()
  }, add = TRUE)
  
  max_abs <- max(abs(mat), na.rm = TRUE)
  if (!is.finite(max_abs) || max_abs == 0) max_abs <- 1
  zlim <- c(-max_abs, max_abs)
  
  cols_delta <- colorRampPalette(c("red", "white", "blue"))(N_COLS)
  
  nx <- length(labels_x)
  ny <- length(labels_y)
  
  STEP_TICKS <- 10
  
  tick_x <- seq(1, nx, by = STEP_TICKS)
  tick_y <- seq(1, ny, by = STEP_TICKS)
  
  par(mar = c(4, 4, 4, 4))
  
  fields::image.plot(
    x = seq_len(ncol(mat)),
    y = seq_len(nrow(mat)),
    z = t(mat[nrow(mat):1, , drop = FALSE]),
    col = cols_delta,
    axes = FALSE,
    xlab = "",
    ylab = "",
    main = main_txt,
    zlim = zlim,
    legend.lab = legend_txt,
    legend.line = 2.5
  )
  
  # ===============================================================
  # GRID (líneas cada STEP_TICKS)
  # ===============================================================
  
  STEP_TICKS <- 20
  
  nx <- ncol(mat)
  ny <- nrow(mat)
  
  #  abline(v = seq(1, nx, by = STEP_TICKS), col = "grey85", lwd = 1)
  # abline(h = ny - seq(1, ny, by = STEP_TICKS) + 1, col = "grey85", lwd = 1)
  
  
  # ===============================================================
  # RESIDUOS ESPECIALES
  # ===============================================================
  
  for(v in list(
    list(special_res1,"orange"),
    list(special_res2,"magenta"),
    list(special_res3,"darkgray"),
    list(special_res4,"black")
  )){
    
    idx_x <- match(v[[1]], colnames(mat))
    idx_x <- idx_x[!is.na(idx_x)]
    
    idx_y <- match(v[[1]], rownames(mat))
    idx_y <- idx_y[!is.na(idx_y)]
    
    if(length(idx_x) > 0){
      abline(
        v   = idx_x,
        col = v[[2]],
        lwd = 0.5,
        lty = 2
      )
    }
    
    if(length(idx_y) > 0){
      abline(
        h   = ny - idx_y + 1,
        col = v[[2]],
        lwd = 0.5,
        lty = 2
      )
    }
  }
  
  # abline(v = seq(1, nx, by = STEP_TICKS), col = "grey85", lwd = 1)
  # abline(h = ny - seq(1, ny, by = STEP_TICKS) + 1, col = "grey85", lwd = 1)
  
  axis(1, at = tick_x, labels = labels_x[tick_x], las = 2, cex.axis = 0.4)
  axis(2, at = tick_y, labels = rev(labels_y)[tick_y], las = 2, cex.axis = 0.4)
  
  box()
}

###############################################################################
# FIND FILES
###############################################################################

mean_files <- list.files(INPUT_DIR, pattern = "_mean_dist\\.rds$", full.names = TRUE)
var_files  <- list.files(INPUT_DIR, pattern = "_var_dist\\.rds$",  full.names = TRUE)
meta_files <- list.files(INPUT_DIR, pattern = "_meta\\.rds$",      full.names = TRUE)

if (length(mean_files) == 0) stop("No mean_dist files found.")
if (length(var_files)  == 0) stop("No var_dist files found.")
if (length(meta_files) == 0) stop("No meta files found.")

mean_names <- basename(mean_files)
rep_ids <- sub("_mean_dist\\.rds$", "", mean_names)

# Base system name without replica suffix
system_ids <- sub("_rep[0-9]+$", "", rep_ids)
systems <- unique(system_ids)

log_msg("Detected replica mean matrices: ", length(rep_ids))
log_msg("Detected unique systems: ", length(systems))

###############################################################################
# COMBINE REPLICAS PER SYSTEM (done once, outside zoom loop)
###############################################################################

log_msg("------------------------------------------------------------")
log_msg("Combining replicas for all systems...")

combined_systems <- list()

for (sys in systems) {
  
  log_msg("Processing system: ", sys)
  
  rep_mean_files <- list.files(
    INPUT_DIR,
    pattern = paste0("^", sys, "_rep[0-9]+_mean_dist\\.rds$"),
    full.names = TRUE
  )
  
  rep_var_files <- list.files(
    INPUT_DIR,
    pattern = paste0("^", sys, "_rep[0-9]+_var_dist\\.rds$"),
    full.names = TRUE
  )
  
  rep_meta_files <- list.files(
    INPUT_DIR,
    pattern = paste0("^", sys, "_rep[0-9]+_meta\\.rds$"),
    full.names = TRUE
  )
  
  if (length(rep_mean_files) == 0 || length(rep_var_files) == 0 || length(rep_meta_files) == 0) {
    log_msg("Missing replica files. Skipping system.")
    next
  }
  
  rep_mean_files <- sort(rep_mean_files)
  rep_var_files  <- sort(rep_var_files)
  rep_meta_files <- sort(rep_meta_files)
  
  log_msg("Replicas detected: ", length(rep_mean_files))
  
  mean_list <- list()
  var_list  <- list()
  n_list    <- numeric(0)
  labels_ref <- NULL
  atom_meta_ref <- NULL
  
  for (i in seq_along(rep_mean_files)) {
    
    mean_i <- readRDS(rep_mean_files[i])
    var_i  <- readRDS(rep_var_files[i])
    meta_i <- readRDS(rep_meta_files[i])
    
    if (is.null(labels_ref)) {
      labels_ref <- meta_i$labels
      atom_meta_ref <- meta_i$atom_meta
    } else {
      if (!identical(labels_ref, meta_i$labels)) {
        stop(paste0("Label mismatch among replicas for system: ", sys))
      }
    }
    
    mean_list[[i]] <- mean_i
    var_list[[i]]  <- var_i
    n_list[i]      <- meta_i$nframes
  }
  
  comb <- combine_replica_stats(mean_list, var_list, n_list)
  
  mean_full <- upper_to_full(comb$mean)
  var_full  <- upper_to_full(comb$var)
  
  normvar_full <- var_full / (mean_full^2 + EPSILON)
  
  rownames(mean_full)    <- labels_ref
  colnames(mean_full)    <- labels_ref
  rownames(var_full)     <- labels_ref
  colnames(var_full)     <- labels_ref
  rownames(normvar_full) <- labels_ref
  colnames(normvar_full) <- labels_ref
  
  # Enforce canonical chain order (ABCD) using residue_map_MD.csv
  rmap      <- read.csv(RESIDUE_MAP, stringsAsFactors = FALSE)
  pdb_id    <- info$pdb
  rmap_sys  <- rmap[rmap$system == pdb_id, ]
  rmap_sys  <- rmap_sys[order(rmap_sys$chain, rmap_sys$residue), ]
  canon_labels <- paste0(rmap_sys$chain, rmap_sys$residue)
  ord       <- canon_labels[canon_labels %in% labels_ref]
  
  if (length(ord) == length(labels_ref)) {
    mean_full    <- mean_full[ord, ord, drop = FALSE]
    var_full     <- var_full[ord, ord, drop = FALSE]
    normvar_full <- normvar_full[ord, ord, drop = FALSE]
    labels_ref   <- ord
  } else {
    log_msg("WARNING: residue_map order incomplete for ", sys, " — keeping original order")
  }
  
  combined_systems[[sys]] <- list(
    mean    = mean_full,
    normvar = normvar_full,
    var     = var_full,
    labels  = labels_ref,
    atom_meta = atom_meta_ref,
    n_list  = n_list,
    info    = classify_system(sys)
  )
  
  log_msg("Replicas combined for: ", sys)
}

log_msg("All systems combined. Total: ", length(combined_systems))

###############################################################################
# COMPUTE DELTA MATRICES (done once, outside zoom loop)
# Delta = condition - Apo, for each mutant and each matrix type (mean, normvar)
# Comparisons: Cap-Apo and Heat-Apo only
###############################################################################

log_msg("------------------------------------------------------------")
log_msg("Computing delta matrices...")

# Index systems by mutant and condition
sys_index <- list()
for (sys in names(combined_systems)) {
  info <- combined_systems[[sys]]$info
  mut  <- info$mut
  cond <- info$cond
  if (is.null(sys_index[[mut]])) sys_index[[mut]] <- list()
  sys_index[[mut]][[cond]] <- sys
}

delta_mats <- list()   # named: e.g. "7LPB_WT_delta_Cap_mean"

for (mut in names(sys_index)) {
  
  apo_sys <- sys_index[[mut]][["Apo"]]
  
  if (is.null(apo_sys)) {
    log_msg("No Apo system found for mutant: ", mut, " — skipping deltas")
    next
  }
  
  for (cond in c("Cap", "Heat")) {
    
    cond_sys <- sys_index[[mut]][[cond]]
    
    if (is.null(cond_sys)) {
      log_msg("No ", cond, " system found for mutant: ", mut, " — skipping delta")
      next
    }
    
    for (mat_type in c("mean", "normvar")) {
      
      mat_apo  <- combined_systems[[apo_sys]][[mat_type]]
      mat_cond <- combined_systems[[cond_sys]][[mat_type]]
      
      delta <- mat_cond - mat_apo
      
      delta_name <- paste0(mut, "_delta_", cond, "_", mat_type)
      delta_mats[[delta_name]] <- delta
      
      log_msg("Delta computed: ", delta_name)
    }
  }
}

log_msg("Delta matrices computed: ", length(delta_mats))

###############################################################################
# COMPUTE DELTA-DELTA MATRICES (done once, outside zoom loop)
# Delta-delta = mutant delta - WT delta, for Cap-Apo and Heat-Apo
###############################################################################

log_msg("------------------------------------------------------------")
log_msg("Computing delta-delta matrices...")

dd_mats <- list()   # named: e.g. "W426A_dd_Cap_mean"

mutants_non_wt <- setdiff(names(sys_index), "WT")

for (mut in mutants_non_wt) {
  
  for (cond in c("Cap", "Heat")) {
    
    for (mat_type in c("mean", "normvar")) {
      
      wt_delta_name  <- paste0("WT_delta_",  cond, "_", mat_type)
      mut_delta_name <- paste0(mut, "_delta_", cond, "_", mat_type)
      
      wt_delta  <- delta_mats[[wt_delta_name]]
      mut_delta <- delta_mats[[mut_delta_name]]
      
      if (is.null(wt_delta) || is.null(mut_delta)) {
        log_msg("Missing delta for dd: ", mut, " ", cond, " ", mat_type, " — skipping")
        next
      }
      
      dd <- mut_delta - wt_delta
      
      dd_name <- paste0(mut, "_dd_", cond, "_", mat_type)
      dd_mats[[dd_name]] <- dd
      
      log_msg("Delta-delta computed: ", dd_name)
    }
  }
}

log_msg("Delta-delta matrices computed: ", length(dd_mats))

###############################################################################
# ZOOM LOOP
###############################################################################

log_msg("============================================================")
log_msg("Starting zoom range loop. Total ranges: ", length(ZOOM_RANGES))

for (zoom_idx in seq_along(ZOOM_RANGES)) {
  
  X_RANGE <- ZOOM_RANGES[[zoom_idx]]$X
  Y_RANGE <- ZOOM_RANGES[[zoom_idx]]$Y
  
  zoom_tag <- paste0(
    X_RANGE[1], "_", X_RANGE[2],
    "_vs_",
    Y_RANGE[1], "_", Y_RANGE[2]
  )
  
  log_msg("------------------------------------------------------------")
  log_msg("ZOOM ", zoom_idx, " of ", length(ZOOM_RANGES), " | tag: ", zoom_tag)
  
  zoom_out <- file.path(OUTPUT_DIR, zoom_tag)
  
  dir_distograms  <- file.path(zoom_out, "distograms")
  dir_delta       <- file.path(zoom_out, "delta")
  dir_deltadelta  <- file.path(zoom_out, "delta_delta")
  
  for (d in c(dir_distograms, dir_delta, dir_deltadelta)) {
    if (!dir.exists(d)) dir.create(d, recursive = TRUE)
  }
  
  # ============================================================
  # PARALLELIZED PER-SYSTEM DISTOGRAMS + RDS SAVES
  # ============================================================
  
  if (RUN_PER_SYSTEM) {
  
  log_msg("Generating per-system distograms (parallel, ", N_CORES, " cores)...")
  
    ZLIM_MEAN <- c(
      0,
      quantile(unlist(lapply(combined_systems, function(s) s$mean)), 0.95, na.rm = TRUE)
    )
    ZLIM_NORMVAR <- c(
      0,
      quantile(unlist(lapply(combined_systems, function(s) s$normvar)), 0.99997, na.rm = TRUE)
    )
  
  cl <- makeCluster(N_CORES)
  
  clusterExport(cl, varlist = c(
    "combined_systems","info",
    "X_RANGE", "Y_RANGE",
    "zoom_out", "dir_distograms",
    "PNG_W", "PNG_H", "PNG_DPI",
    "N_COLS", "EPSILON", "OVERWRITE",
    "special_res1", "special_res2", "special_res3", "special_res4",
    "subset_nodes", "plot_heatmap", "classify_system",
    "log_msg", "LOG_FILE", "ZLIM_MEAN", "ZLIM_NORMVAR"
  ), envir = environment())
  
  clusterEvalQ(cl, {
    suppressPackageStartupMessages({
      library(viridisLite)
      library(fields)
    })
  })
  
  parLapply(cl, names(combined_systems), function(sys) {
    
    info <- combined_systems[[sys]]$info
    
    out_mean_png <- file.path(dir_distograms, paste0(sys, "_mean_distogram.png"))
    out_norm_png <- file.path(dir_distograms, paste0(sys, "_normvar_distogram.png"))
    out_mean_rds <- file.path(dir_distograms, paste0(sys, "_combined_mean_dist.rds"))
    out_var_rds  <- file.path(dir_distograms, paste0(sys, "_combined_var_dist.rds"))
    out_nv_rds   <- file.path(dir_distograms, paste0(sys, "_combined_normvar_dist.rds"))
    out_meta_rds <- file.path(dir_distograms, paste0(sys, "_combined_meta.rds"))
    
    if (file.exists(out_mean_png) && file.exists(out_norm_png) && !OVERWRITE) {
      return(invisible(NULL))
    }
    
    mean_full    <- combined_systems[[sys]]$mean
    normvar_full <- combined_systems[[sys]]$normvar
    var_full     <- combined_systems[[sys]]$var
    labels_ref   <- combined_systems[[sys]]$labels
    atom_meta_ref <- combined_systems[[sys]]$atom_meta
    n_list       <- combined_systems[[sys]]$n_list
    
    # Zoom selection
    nodes_x <- subset_nodes(colnames(mean_full), X_RANGE)
    nodes_y <- subset_nodes(rownames(mean_full), Y_RANGE)
    
    mean_z    <- mean_full[nodes_y, nodes_x, drop = FALSE]
    normvar_z <- normvar_full[nodes_y, nodes_x, drop = FALSE]
    
    labels_x <- colnames(mean_full)[nodes_x]
    labels_y <- rownames(mean_full)[nodes_y]
    
    # Save combined matrices
    saveRDS(mean_full,    out_mean_rds)
    saveRDS(var_full,     out_var_rds)
    saveRDS(normvar_full, out_nv_rds)
    
    saveRDS(
      list(
        system = sys,
        pdb = info$pdb,
        mutation = info$mut,
        condition = info$cond,
        labels = labels_ref,
        atom_meta = atom_meta_ref,
        replicas_combined = length(n_list),
        total_frames = sum(n_list),
        x_range = X_RANGE,
        y_range = Y_RANGE,
        normalized_variance_definition = "var / mean^2"
      ),
      out_meta_rds
    )
    
    plot_heatmap(
      mat = mean_z,
      labels_x = labels_x,
      labels_y = labels_y,
      outfile = out_mean_png,
      main_txt = paste0(sys, " | Mean distance"),
      legend_txt = "Distance (Å)",
      invert_col = TRUE,
      zlim_fixed = ZLIM_MEAN
    )
    
    plot_heatmap(
      mat = normvar_z,
      labels_x = labels_x,
      labels_y = labels_y,
      outfile = out_norm_png,
      main_txt = paste0(sys, " | Normalized variance"),
      legend_txt = "Var / Mean^2",
      invert_col = FALSE,
      zlim_fixed = ZLIM_NORMVAR
    )
    
    return(invisible(NULL))
  })
  
  stopCluster(cl)
  
  log_msg("Per-system distograms done for zoom: ", zoom_tag)
  }
  
  # ============================================================
  # PARALLELIZED DELTA DISTOGRAMS
  # ============================================================
  
  if (RUN_DELTA) {
  
  log_msg("Generating delta distograms (parallel, ", N_CORES, " cores)...")
  
  cl <- makeCluster(N_CORES)
  
  clusterExport(cl, varlist = c(
    "delta_mats",
    "X_RANGE", "Y_RANGE",
    "dir_delta",
    "PNG_W", "PNG_H", "PNG_DPI",
    "N_COLS", "OVERWRITE",
    "special_res1", "special_res2", "special_res3", "special_res4",
    "subset_nodes", "plot_heatmap_delta",
    "log_msg", "LOG_FILE"
  ), envir = environment())
  
  clusterEvalQ(cl, {
    suppressPackageStartupMessages({
      library(viridisLite)
      library(fields)
    })
  })
  
  parLapply(cl, names(delta_mats), function(delta_name) {
    
    out_png <- file.path(dir_delta, paste0(delta_name, "_distogram.png"))
    out_rds <- file.path(dir_delta, paste0(delta_name, ".rds"))
    
    if (file.exists(out_png) && !OVERWRITE) return(invisible(NULL))
    
    mat_full <- delta_mats[[delta_name]]
    
    nodes_x <- subset_nodes(colnames(mat_full), X_RANGE)
    nodes_y <- subset_nodes(rownames(mat_full), Y_RANGE)
    
    mat_z    <- mat_full[nodes_y, nodes_x, drop = FALSE]
    labels_x <- colnames(mat_full)[nodes_x]
    labels_y <- rownames(mat_full)[nodes_y]
    
    saveRDS(mat_full, out_rds)
    
    # legend label: detect matrix type from name
    leg_txt <- if (grepl("_mean$", delta_name)) "Delta Distance (Å)" else "Delta Var / Mean^2"
    
    plot_heatmap_delta(
      mat = mat_z,
      labels_x = labels_x,
      labels_y = labels_y,
      outfile = out_png,
      main_txt = paste0(delta_name, " | Delta distogram"),
      legend_txt = leg_txt
    )
    
    return(invisible(NULL))
  })
  
  stopCluster(cl)
  
  log_msg("Delta distograms done for zoom: ", zoom_tag)
  }
  
  # ============================================================
  # PARALLELIZED DELTA-DELTA DISTOGRAMS
  # ============================================================
  
  if (RUN_DELTADELTA) {
  
  log_msg("Generating delta-delta distograms (parallel, ", N_CORES, " cores)...")
  
  cl <- makeCluster(N_CORES)
  
  clusterExport(cl, varlist = c(
    "dd_mats",
    "X_RANGE", "Y_RANGE",
    "dir_deltadelta",
    "PNG_W", "PNG_H", "PNG_DPI",
    "N_COLS", "OVERWRITE",
    "special_res1", "special_res2", "special_res3", "special_res4",
    "subset_nodes", "plot_heatmap_delta",
    "log_msg", "LOG_FILE"
  ), envir = environment())
  
  clusterEvalQ(cl, {
    suppressPackageStartupMessages({
      library(viridisLite)
      library(fields)
    })
  })
  
  parLapply(cl, names(dd_mats), function(dd_name) {
    
    out_png <- file.path(dir_deltadelta, paste0(dd_name, "_distogram.png"))
    out_rds <- file.path(dir_deltadelta, paste0(dd_name, ".rds"))
    
    if (file.exists(out_png) && !OVERWRITE) return(invisible(NULL))
    
    mat_full <- dd_mats[[dd_name]]
    
    nodes_x <- subset_nodes(colnames(mat_full), X_RANGE)
    nodes_y <- subset_nodes(rownames(mat_full), Y_RANGE)
    
    mat_z    <- mat_full[nodes_y, nodes_x, drop = FALSE]
    labels_x <- colnames(mat_full)[nodes_x]
    labels_y <- rownames(mat_full)[nodes_y]
    
    saveRDS(mat_full, out_rds)
    
    # legend label: detect matrix type from name
    leg_txt <- if (grepl("_mean$", dd_name)) "Delta-Delta Distance (Å)" else "Delta-Delta Var / Mean^2"
    
    plot_heatmap_delta(
      mat = mat_z,
      labels_x = labels_x,
      labels_y = labels_y,
      outfile = out_png,
      main_txt = paste0(dd_name, " | Delta-delta distogram"),
      legend_txt = leg_txt
    )
    
    return(invisible(NULL))
  })
  
  stopCluster(cl)
  
  log_msg("Delta-delta distograms done for zoom: ", zoom_tag)
  }
  log_msg("ZOOM ", zoom_idx, " COMPLETE: ", zoom_tag)
}

###############################################################################
# END
###############################################################################

log_msg("============================================================")
log_msg("SCRIPT FINISHED SUCCESSFULLY")
log_msg("End time: ", format(Sys.time(), "%Y-%m-%d %H:%M:%S"))
log_msg("============================================================")
