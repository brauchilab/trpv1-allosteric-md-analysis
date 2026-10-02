# ================================================================
# Script 02C_DIAG v3.0: Diagnostic NMA Overlap Analysis
# Author: DenyCB
# Date: 2025-05-05
# ----------------------------------------------------------------
# Based directly on 02C v3.0. Changes vs 02C v3.0:
#   1. output_root -> overlap_diag/
#   2. Additional parameters: restricted_ranges, n_ca_per_chain,
#      chain_labels, chain_colors
#   3. Additional subdirectories for diagnostic outputs
#   4. New function: make_perchain_panel()
#   5. New function: compute_overlap_sub()
#   6. compute_dv(): corrected atom-level intersect (not xyz-level)
#   7. foreach: additionally computes ov_restricted and ov_chains
#   8. Assembly Section 6: only desde_apo; adds restricted and
#      per-chain figure sets; no inter-mutant figures
#   9. compute_common_ymax(): guards for empty list and non-df items
#
# PURPOSE:
#   Diagnostic companion to 02C v3.0. Runs two additional overlap
#   analyses alongside the standard analysis:
#
#   ANALYSIS 1 - Standard (desde_apo only, full range, all chains):
#     Reproduced for direct comparison with diagnostic analyses.
#
#   ANALYSIS 2 - Restricted residue range (400-602 + 625-715):
#     Tests whether distributed overlap is driven by flexible
#     distal regions outside the TMD+TRP helix core.
#
#   ANALYSIS 3 - Per-chain (full range, position-based):
#     Overlap calculated independently for each chain
#     (atoms 1-454=chain1, 455-908=chain2, etc.).
#     NOTE: arbitrary assignment - diagnostic only.
#
# INPUTS/OUTPUTS: see overlap_diag/ subfolders
# ================================================================

suppressPackageStartupMessages({
  library(bio3d)
  library(ggplot2)
  library(gridExtra)
  library(parallel)
  library(foreach)
  library(doParallel)
})

# ================================================================
# SECTION 1: PARAMETERS — edit here to configure the analysis
# ================================================================

# --- Input: NMA data directory (consolidated outputs from 02A v4.0) ---
nma_data_dir <- "C:/DinamicasMoleculares/analisis_bio3d_output/2_nma_correlation/nma_data"

# --- Input: residue mapping (for robust coord matching by resid) ---
trimmed_map_file <- "C:/DinamicasMoleculares/TRPV1_pipeline/config_files/trimmed_residue_map.csv"

# --- Output root (diagnostic, separate from main 02C output) ---
output_root <- "C:/DinamicasMoleculares/analisis_bio3d_output/2_nma_correlation/overlap_diag"

# --- Parallel cores (PSOCK cluster, Windows 10 compatible) ---
cores <- 6

# --- Number of modes to include in overlap calculation ---
# bio3d::overlap() uses the first nmodes non-trivial modes (modes 7 to 7+nmodes-1)
nmodes_overlap <- 200

# --- X-axis range for acotado plots (low-frequency modes detail) ---
# Mode index 1 = first non-trivial mode (mode 7 in NMA numbering)
plot_mode_range <- c(1, 45)

# --- X-axis range for extendido plots (full range to detect high-frequency contributions) ---
plot_mode_range_extended <- c(1, 500)

# --- System identifiers (must match directory/file naming from 01/02A) ---
# PDB base codes: 7LP9 = apo, 7LPB = cap (capsaicin), 7LPC = heat
pdb_apo  <- "7LP9"
pdb_cap  <- "7LPB"
pdb_heat <- "7LPC"

# --- Mutants to analyze (must match naming convention in file names) ---
mutants <- c("WT", "W426A", "W697A", "Y441A")

# --- Figure dimensions ---
fig_width_px  <- 1600   # pixels wide per figure
fig_height_px <- 400    # pixels tall per panel (total height = n_panels * fig_height_px)
fig_res       <- 150    # resolution (dpi)

# --- ANALYSIS 2: Restricted residue range ---
# Cα in these resno ranges (all 4 chains) used for dv and nma_obj$U subsetting.
restricted_ranges <- list(c(400, 602), c(625, 715))

# --- ANALYSIS 3: Per-chain parameters ---
# Position-based chain assignment: atoms 1:n_ca_per_chain = chain1, etc.
# Verify n_ca_per_chain = nrow(pdb_ca$atom) / 4 before running.
n_ca_per_chain <- 454
chain_labels   <- c("Chain1", "Chain2", "Chain3", "Chain4")
chain_colors   <- c("Chain1" = "#2166AC",
                    "Chain2" = "#D6604D",
                    "Chain3" = "#4DAC26",
                    "Chain4" = "#8E44AD")

# ================================================================
# SECTION 2: OUTPUT DIRECTORIES AND LOGGING
# ================================================================

subdirs <- c(
  "figures_standard/acotado",
  "figures_standard/extendido",
  "figures_restricted/acotado",
  "figures_restricted/extendido",
  "figures_perchain/acotado",
  "figures_perchain/extendido",
  "logs"
)
for (sd in subdirs) {
  dir.create(file.path(output_root, sd), recursive = TRUE, showWarnings = FALSE)
}

log_file <- file.path(
  output_root, "logs",
  paste0("02C_DIAG_log_", format(Sys.time(), "%Y%m%d_%H%M%S"), ".txt")
)

# Buffered log (minimize disk I/O, same pattern as 02A v4.0)
log_buffer <- character()

log_message <- function(...) {
  msg <- paste0(format(Sys.time(), "%Y-%m-%d %H:%M:%OS3"), " - ",
                paste(..., collapse = " "))
  cat(msg, "\n")
  log_buffer <<- c(log_buffer, msg)
}

flush_log <- function(append = TRUE) {
  if (length(log_buffer) > 0) {
    cat(paste0(log_buffer, collapse = "\n"), "\n",
        file = log_file, append = append)
    log_buffer <<- character()
  }
}

log_message("=== Starting NMA Overlap Diagnostic Analysis (02C_DIAG v3.0) ===")
log_message("Parameters:")
log_message("  nma_data_dir         =", nma_data_dir)
log_message("  output_root          =", output_root)
log_message("  nmodes_overlap       =", nmodes_overlap)
log_message("  plot_mode_range      =", paste(plot_mode_range, collapse = "-"))
log_message("  plot_mode_range_ext  =", paste(plot_mode_range_extended, collapse = "-"))
log_message("  pdb_apo/cap/heat     =", pdb_apo, "/", pdb_cap, "/", pdb_heat)
log_message("  mutants              =", paste(mutants, collapse = ", "))
log_message("  restricted_ranges    =",
            paste(sapply(restricted_ranges, function(r) paste(r, collapse="-")),
                  collapse=" + "))
log_message("  n_ca_per_chain       =", n_ca_per_chain)

# ================================================================
# SECTION 3: HELPER FUNCTIONS
# ================================================================

# --- Build system label from pdb code and mutant name ---
# Example: make_label("7LP9", "W426A") -> "7LP9_W426A"
make_label <- function(pdb_code, mutant) {
  paste0(pdb_code, "_", mutant)
}

# --- Load NMA object and trimmed Cα pdb for a given system label ---
# Returns list($nma_obj, $pdb_ca) or NULL if files missing.
load_nma_system <- function(label, nma_dir) {
  nma_file <- file.path(nma_dir, paste0(label, "_nma_obj.rds"))
  pdb_file <- file.path(nma_dir, paste0(label, "_pdb_ca_trimmed.rds"))

  if (!file.exists(nma_file)) {
    log_message("  MISSING nma_obj:", nma_file)
    return(NULL)
  }
  if (!file.exists(pdb_file)) {
    log_message("  MISSING pdb_ca_trimmed:", pdb_file)
    return(NULL)
  }

  list(
    nma_obj = readRDS(nma_file),
    pdb_ca  = readRDS(pdb_file)
  )
}

# --- Extract Cα coordinate vector (length 3N) from a pdb_ca object ---
# Returns named numeric vector with names like "A_300_x", "A_300_y", etc.
# for robust matching by chain+resid rather than position.
get_ca_coords <- function(pdb_ca) {
  atom  <- pdb_ca$atom
  xyz   <- as.numeric(pdb_ca$xyz)
  n_ca  <- nrow(atom)

  # Build atom labels: chain_resno for each Cα
  labels_x <- paste0(atom$chain, "_", atom$resno, "_x")
  labels_y <- paste0(atom$chain, "_", atom$resno, "_y")
  labels_z <- paste0(atom$chain, "_", atom$resno, "_z")

  # xyz vector is stored as x1,y1,z1, x2,y2,z2, ...
  all_labels <- as.vector(rbind(labels_x, labels_y, labels_z))
  names(xyz) <- all_labels
  return(xyz)
}

# --- Compute displacement vector dv between two systems ---
# Matches atoms by chain+resno label (atom-level, not xyz-level) to handle
# potential ordering differences. Aligns B onto A via fit.xyz().
# Returns numeric vector of length 3*n_common_atoms.
# dv = aligned_coords_B - coords_A  (direction: A -> B)
compute_dv <- function(pdb_ca_A, pdb_ca_B, label_A, label_B) {
  coords_A <- get_ca_coords(pdb_ca_A)
  coords_B <- get_ca_coords(pdb_ca_B)

  # Build atom-level labels (strip _x/_y/_z suffix) for matching
  atom_labels_A <- unique(sub("_(x|y|z)$", "", names(coords_A)))
  atom_labels_B <- unique(sub("_(x|y|z)$", "", names(coords_B)))
  common_atoms  <- intersect(atom_labels_A, atom_labels_B)

  if (length(common_atoms) == 0) {
    stop(paste("No common Cα atoms found between", label_A, "and", label_B))
  }
  if (length(common_atoms) < length(atom_labels_A)) {
    log_message("  WARNING: only", length(common_atoms), "of", length(atom_labels_A),
                "atoms matched between", label_A, "and", label_B,
                "- using common subset.")
  }

  # Expand atom labels to xyz coordinate labels (x,y,z interleaved)
  common_xyz <- as.vector(rbind(
    paste0(common_atoms, "_x"),
    paste0(common_atoms, "_y"),
    paste0(common_atoms, "_z")
  ))

  coords_A_vec <- as.numeric(coords_A[common_xyz])
  coords_B_vec <- as.numeric(coords_B[common_xyz])

  # Alineamiento de B sobre A usando Cα del rango trimmed como referencia.
  # fit.xyz() minimiza RMSD por rotación/traslación (superposición de Procrustes).
  # fixed  = coordenadas de A (referencia)
  # mobile = coordenadas de B (a rotar/trasladar)
  # Ambos vectores deben estar en orden x1,y1,z1, x2,y2,z2,...
  coords_B_aligned <- tryCatch(
    as.numeric(fit.xyz(fixed = coords_A_vec, mobile = coords_B_vec)),
    error = function(e) {
      log_message("  WARNING: fit.xyz() failed for", label_A, "->", label_B,
                  ":", e$message, "- using unaligned coords.")
      coords_B_vec
    }
  )

  dv <- coords_B_aligned - coords_A_vec
  return(as.numeric(dv))  # length = 3 * length(common_atoms)
}

# --- Compute overlap scores using bio3d::overlap() ---
# nma_obj : NMA object (class "nma") used as the conformational reference
# dv      : displacement vector of length 3N
# nmodes  : number of non-trivial modes to include
# label   : string label for logging
# Returns data.frame with columns: mode_index, overlap, overlap_cum
compute_overlap_sub <- function(nma_obj, U_sub, dv_sub, nmodes, label) {
  if (length(dv_sub) != nrow(U_sub)) {
    log_message("  ERROR dimension mismatch for", label,
                ": dv_sub =", length(dv_sub), ", U_sub rows =", nrow(U_sub))
    return(NULL)
  }
  nma_sub      <- nma_obj
  nma_sub$U    <- U_sub
  nma_sub$xyz  <- dv_sub  # bio3d::overlap() uses xyz for dimension check
  nmodes_use   <- min(nmodes, ncol(U_sub) - 6)
  log_message("  Computing overlap (sub) for:", label,
              "| dv_sub length:", length(dv_sub),
              "| modes used:", nmodes_use)
  ov <- tryCatch(
    bio3d::overlap(nma_sub, dv = dv_sub, nmodes = nmodes_use),
    error = function(e) {
      log_message("  ERROR bio3d::overlap() for", label, ":", e$message)
      NULL
    }
  )
  if (is.null(ov)) return(NULL)
  data.frame(
    mode_index   = seq_len(nmodes_use),
    overlap      = as.numeric(ov$overlap),
    overlap_cum  = as.numeric(ov$overlap.cum),
    stringsAsFactors = FALSE
  )
}

# --- Compute overlap with subsetted U matrix and dv (for restricted/per-chain) ---
# nma_obj : full NMA object (metadata; U is replaced internally)
# U_sub   : subsetted nma_obj$U (rows = 3*n_subset_atoms, cols = all modes)
# dv_sub  : displacement vector, length = nrow(U_sub)
# nmodes  : number of non-trivial modes to include
# label   : string label for logging
# Returns data.frame with columns: mode_index, overlap, overlap_cum
compute_overlap_sub <- function(nma_obj, U_sub, dv_sub, nmodes, label) {
  if (length(dv_sub) != nrow(U_sub)) {
    log_message("  ERROR dimension mismatch for", label,
                ": dv_sub =", length(dv_sub), ", U_sub rows =", nrow(U_sub))
    return(NULL)
  }
  nma_sub   <- nma_obj
  nma_sub$U <- U_sub
  nmodes_use <- min(nmodes, ncol(U_sub) - 6)
  log_message("  Computing overlap (sub) for:", label,
              "| dv_sub length:", length(dv_sub),
              "| modes used:", nmodes_use)
  ov <- tryCatch(
    bio3d::overlap(nma_sub, dv = dv_sub, nmodes = nmodes_use),
    error = function(e) {
      log_message("  ERROR bio3d::overlap() for", label, ":", e$message)
      NULL
    }
  )
  if (is.null(ov)) return(NULL)
  data.frame(
    mode_index   = seq_len(nmodes_use),
    overlap      = as.numeric(ov$overlap),
    overlap_cum  = as.numeric(ov$overlap.cum),
    stringsAsFactors = FALSE
  )
}

# --- Build a per-chain ggplot2 panel (Analysis 3) ---
# df_chains  : list of 4 data.frames from compute_overlap_sub() (NULLs allowed)
# title      : panel title string
# mode_range : c(min, max) for x-axis
# y_max_bars : common upper limit for bars Y axis (left), shared across figure
# Returns a ggplot with 4 colored bar sets (side-by-side) and 4 cumulative lines.
make_perchain_panel <- function(df_chains, chain_labels, chain_colors,
                                title, mode_range = c(1, 45),
                                y_max_bars = NULL) {
  combined <- do.call(rbind, lapply(seq_along(df_chains), function(i) {
    df <- df_chains[[i]]
    if (is.null(df)) return(NULL)
    df$chain <- chain_labels[i]
    df
  }))
  if (is.null(combined) || nrow(combined) == 0) {
    return(ggplot() + labs(title = paste(title, "— NO DATA")) +
             theme_void() + theme(plot.title = element_text(size = 11, color = "red")))
  }
  combined <- combined[combined$mode_index >= mode_range[1] &
                       combined$mode_index <= mode_range[2], ]
  combined$chain <- factor(combined$chain, levels = chain_labels)
  if (is.null(y_max_bars)) y_max_bars <- max(combined$overlap, na.rm = TRUE) * 1.05
  if (!is.finite(y_max_bars) || y_max_bars <= 0) y_max_bars <- 0.01
  scale_factor <- y_max_bars
  ggplot(combined, aes(x = mode_index)) +
    geom_col(aes(y = overlap, fill = chain, color = chain),
             position = position_dodge(width = 0.9), alpha = 0.75, width = 0.85) +
    geom_line(aes(y = overlap_cum * scale_factor, color = chain, group = chain),
              linewidth = 0.7) +
    geom_point(aes(y = overlap_cum * scale_factor, color = chain, group = chain),
               size = 0.9) +
    scale_fill_manual(values  = chain_colors, name = "Chain") +
    scale_color_manual(values = chain_colors, name = "Chain") +
    scale_x_continuous(breaks = seq(mode_range[1], mode_range[2], by = 5),
                       limits = mode_range) +
    scale_y_continuous(
      name     = "Squared overlap per chain (bars)",
      limits   = c(0, y_max_bars),
      sec.axis = sec_axis(transform = ~ . / scale_factor,
                          name = "Cumulative overlap per chain (line)")
    ) +
    labs(title = title, x = "Mode index (non-trivial; 1 = lowest frequency)") +
    theme_minimal(base_size = 12) +
    theme(plot.title = element_text(size = 13, face = "bold"),
          axis.title = element_text(size = 11),
          axis.text  = element_text(size = 9),
          axis.title.y.right = element_text(color = "grey30", size = 11),
          axis.text.y.right  = element_text(color = "grey30", size = 9),
          legend.position    = "right",
          panel.grid.minor   = element_blank())
}

# --- Build a single ggplot2 panel for one overlap result ---
# df_overlap  : data.frame from compute_overlap()
# title       : panel title string
# mode_range  : c(min, max) for x-axis
# y_max_bars  : common upper limit for left (bars) Y axis across all panels in figure.
#               Set to max(overlap) * 1.05 computed over all panels of the figure.
# Returns a ggplot object with dual Y axes:
#   Left  axis: squared overlap per mode (bars), scale 0 to y_max_bars
#   Right axis: cumulative overlap (line), scale 0 to 1
make_overlap_panel <- function(df_overlap, title, mode_range = c(1, 45),
                               y_max_bars = NULL) {
  df_plot <- df_overlap[
    df_overlap$mode_index >= mode_range[1] &
    df_overlap$mode_index <= mode_range[2], ]

  # If y_max_bars not provided, use panel-local max * 1.05
  if (is.null(y_max_bars)) {
    y_max_bars <- max(df_plot$overlap, na.rm = TRUE) * 1.05
  }
  # Avoid y_max_bars = 0 (edge case with all-zero overlaps)
  if (!is.finite(y_max_bars) || y_max_bars <= 0) y_max_bars <- 0.01

  # Scale factor to map cumulative overlap (0-1) onto bars axis (0-y_max_bars)
  # Right axis shows 0-1; bars axis shows 0-y_max_bars.
  # ggplot sec_axis requires a transformation: right = left / scale_factor
  scale_factor <- y_max_bars  # right axis value = left axis value / y_max_bars

  ggplot(df_plot, aes(x = mode_index)) +
    # Bars: individual squared overlap per mode (left axis)
    geom_col(aes(y = overlap), fill = "steelblue", alpha = 0.75, width = 0.8) +
    # Line + points: cumulative overlap, scaled to left axis units
    geom_line(aes(y = overlap_cum * scale_factor), color = "firebrick", linewidth = 0.8) +
    geom_point(aes(y = overlap_cum * scale_factor), color = "firebrick", size = 1.2) +
    scale_x_continuous(
      breaks = seq(mode_range[1], mode_range[2], by = 5),
      limits = mode_range
    ) +
    scale_y_continuous(
      name   = "Squared overlap (bars)",
      limits = c(0, y_max_bars),
      sec.axis = sec_axis(
        transform = ~ . / scale_factor,
        name      = "Cumulative overlap (line)"
      )
    ) +
    labs(
      title = title,
      x     = "Mode index (non-trivial; 1 = lowest frequency)"
    ) +
    theme_minimal(base_size = 12) +
    theme(
      plot.title        = element_text(size = 13, face = "bold"),
      axis.title        = element_text(size = 11),
      axis.text         = element_text(size = 9),
      axis.title.y.right = element_text(color = "firebrick", size = 11),
      axis.text.y.right  = element_text(color = "firebrick", size = 9),
      panel.grid.minor  = element_blank()
    )
}

# --- Compute common y_max_bars across a list of overlap data.frames ---
# panels_data : list of data.frames (from compute_overlap), may contain NULLs
# mode_range  : c(min, max) to restrict to the plotted range
# Returns scalar: max(overlap) * 1.05 across all non-NULL panels
compute_common_ymax <- function(panels_data, mode_range) {
  if (length(panels_data) == 0) return(0.01)
  max_vals <- sapply(panels_data, function(df) {
    if (is.null(df) || !is.data.frame(df)) return(NA_real_)
    df_sub <- df[df$mode_index >= mode_range[1] & df$mode_index <= mode_range[2], ]
    if (nrow(df_sub) == 0) return(NA_real_)
    max(df_sub$overlap, na.rm = TRUE)
  })
  global_max <- max(as.numeric(max_vals), na.rm = TRUE)
  if (!is.finite(global_max) || global_max <= 0) return(0.01)
  return(global_max * 1.05)
}

# --- Save a list of ggplot panels as a stacked multi-panel PNG ---
# panels     : named list of ggplot objects (top to bottom)
# filepath   : full output path including filename
# panel_h_px : height in pixels per panel
# width_px   : total width in pixels
save_multipanel <- function(panels, filepath, panel_h_px, width_px, res_dpi) {
  n_panels    <- length(panels)
  total_h_px  <- n_panels * panel_h_px

  tryCatch({
    png(filepath, width = width_px, height = total_h_px, res = res_dpi)
    gridExtra::grid.arrange(grobs = panels, ncol = 1, nrow = n_panels)
    dev.off()
    log_message("  Saved figure:", filepath)
  }, error = function(e) {
    log_message("  ERROR saving figure", filepath, ":", e$message)
    try(dev.off(), silent = TRUE)
  })
}

# ================================================================
# SECTION 4: DEFINE ALL COMPARISONS
# ================================================================
# Each comparison is defined by:
#   label_A : system label for state A (reference / "from" state)
#   label_B : system label for state B (target / "to" state)
#   dv_description : human-readable description of the transition
#   comparison_type: "apo_vs_stim" or "inter_mutant"
#   mutant  : which mutant this comparison belongs to
#   stim    : "cap" or "heat" (for apo_vs_stim comparisons)

build_comparisons <- function(pdb_apo, pdb_cap, pdb_heat, mutants) {
  comps <- list()

  for (mut in mutants) {
    # apo -> cap
    comps[[length(comps) + 1]] <- list(
      label_A          = make_label(pdb_apo, mut),
      label_B          = make_label(pdb_cap, mut),
      dv_description   = paste0(mut, " apo->cap (", pdb_apo, " -> ", pdb_cap, ")"),
      comparison_type  = "apo_vs_stim",
      mutant           = mut,
      stim             = "cap"
    )
    # apo -> heat
    comps[[length(comps) + 1]] <- list(
      label_A          = make_label(pdb_apo, mut),
      label_B          = make_label(pdb_heat, mut),
      dv_description   = paste0(mut, " apo->heat (", pdb_apo, " -> ", pdb_heat, ")"),
      comparison_type  = "apo_vs_stim",
      mutant           = mut,
      stim             = "heat"
    )
    # inter-mutant: WT vs this mutant, within each condition
    if (mut != "WT") {
      for (pdb_code in c(pdb_apo, pdb_cap, pdb_heat)) {
        stim_label <- switch(pdb_code,
          "7LP9" = "apo", "7LPB" = "cap", "7LPC" = "heat", pdb_code)
        comps[[length(comps) + 1]] <- list(
          label_A         = make_label(pdb_code, "WT"),
          label_B         = make_label(pdb_code, mut),
          dv_description  = paste0(stim_label, " WT->", mut,
                                   " (", pdb_code, ")"),
          comparison_type = "inter_mutant",
          mutant          = mut,
          stim            = stim_label
        )
      }
    }
  }
  return(comps)
}

all_comparisons <- build_comparisons(pdb_apo, pdb_cap, pdb_heat, mutants)
log_message("Total comparisons defined:", length(all_comparisons))

# ================================================================
# SECTION 5: PARALLEL OVERLAP COMPUTATION
# ================================================================
# Each worker computes overlaps for one comparison (both perspectives)
# and saves per-comparison CSVs. Returns result list for figure assembly.

log_message("Setting up parallel cluster (", cores, "cores, PSOCK)...")
cl <- makeCluster(cores)
doParallel::registerDoParallel(cl)

overlap_results <- foreach(
  comp       = all_comparisons,
  .packages  = c("bio3d"),
  .export    = c(
    "nma_data_dir", "output_root", "nmodes_overlap", "plot_mode_range",
    "log_message", "load_nma_system", "get_ca_coords",
    "compute_dv", "compute_overlap", "compute_overlap_sub", "make_overlap_panel",
    "restricted_ranges", "n_ca_per_chain", "chain_labels"
  ),
  .errorhandling = "pass"
) %dopar% {

  label_A <- comp$label_A
  label_B <- comp$label_B
  desc    <- comp$dv_description

  log_message("------------------------------------------------------------")
  log_message("Processing comparison:", desc)

  # --- Load both systems ---
  sys_A <- load_nma_system(label_A, nma_data_dir)
  sys_B <- load_nma_system(label_B, nma_data_dir)

  if (is.null(sys_A) || is.null(sys_B)) {
    log_message("  SKIPPING: one or both systems could not be loaded.")
    return(list(comp = comp, status = "skipped: missing files",
                ov_from_A = NULL, ov_from_B = NULL,
                ov_restricted = NULL, ov_chains = NULL))
  }

  # --- Compute displacement vector dv = coords_B - coords_A ---
  dv <- tryCatch(
    compute_dv(sys_A$pdb_ca, sys_B$pdb_ca, label_A, label_B),
    error = function(e) {
      log_message("  ERROR computing dv:", e$message)
      NULL
    }
  )
  if (is.null(dv)) {
    return(list(comp = comp, status = "error: dv computation failed",
                ov_from_A = NULL, ov_from_B = NULL,
                ov_restricted = NULL, ov_chains = NULL))
  }

  log_message("  dv computed: length =", length(dv),
              "| norm =", round(sqrt(sum(dv^2)), 3), "Angstrom")

  # --- ANALYSIS 1: Standard overlap (full range, both perspectives) ---
  ov_from_A <- compute_overlap(
    nma_obj = sys_A$nma_obj,
    dv      = dv,
    nmodes  = nmodes_overlap,
    label   = paste0(desc, " [desde ", label_A, "]")
  )

  ov_from_B <- compute_overlap(
    nma_obj = sys_B$nma_obj,
    dv      = dv,
    nmodes  = nmodes_overlap,
    label   = paste0(desc, " [desde ", label_B, "]")
  )

  # --- ANALYSIS 2: Restricted residue range ---
  # Subset dv and nma_obj$U to atom positions matching restricted_ranges resno.
  # Atom index i in pdb_ca$atom corresponds to rows 3i-2, 3i-1, 3i in U and dv.
  ov_restricted <- tryCatch({
    atom_df    <- sys_A$pdb_ca$atom
    resno_keep <- unlist(lapply(restricted_ranges, function(r) r[1]:r[2]))
    keep_atoms <- which(atom_df$resno %in% resno_keep)
    if (length(keep_atoms) == 0) stop("No atoms in restricted range.")
    log_message("  Analysis 2: restricted Cα =", length(keep_atoms))
    xyz_idx <- as.vector(rbind(keep_atoms*3-2, keep_atoms*3-1, keep_atoms*3))
    dv_sub  <- dv[xyz_idx]
    U_sub   <- sys_A$nma_obj$U[xyz_idx, , drop = FALSE]
    compute_overlap_sub(sys_A$nma_obj, U_sub, dv_sub, nmodes_overlap,
                        paste0(desc, " [restricted]"))
  }, error = function(e) {
    log_message("  ERROR Analysis 2:", e$message); NULL
  })

  # --- ANALYSIS 3: Per-chain (position-based, full range) ---
  # Atoms 1:n_ca_per_chain = chain1, etc. Returns list of 4 data.frames.
  n_ca_total <- nrow(sys_A$pdb_ca$atom)
  ov_chains  <- lapply(seq_len(4), function(ch_idx) {
    tryCatch({
      atom_start <- (ch_idx - 1) * n_ca_per_chain + 1
      atom_end   <- min(ch_idx * n_ca_per_chain, n_ca_total)
      ch_atoms   <- atom_start:atom_end
      xyz_idx    <- as.vector(rbind(ch_atoms*3-2, ch_atoms*3-1, ch_atoms*3))
      valid      <- xyz_idx[xyz_idx <= length(dv)]
      if (length(valid) == 0) return(NULL)
      dv_ch <- dv[valid]
      U_ch  <- sys_A$nma_obj$U[valid, , drop = FALSE]
      log_message("  Analysis 3 chain", ch_idx, ": atoms", atom_start, "-", atom_end)
      compute_overlap_sub(sys_A$nma_obj, U_ch, dv_ch, nmodes_overlap,
                          paste0(desc, " [chain", ch_idx, "]"))
    }, error = function(e) {
      log_message("  ERROR Analysis 3 chain", ch_idx, ":", e$message); NULL
    })
  })

  log_message("  Completed:", desc)
  return(list(
    comp          = comp,
    status        = "ok",
    ov_from_A     = ov_from_A,
    ov_from_B     = ov_from_B,
    ov_restricted = ov_restricted,
    ov_chains     = ov_chains,
    label_A       = label_A,
    label_B       = label_B
  ))
}

stopCluster(cl)
log_message("Parallel computation complete.")

# Filter successful results
ok_results <- Filter(function(x) {
  !inherits(x, "error") && !is.null(x) && x$status == "ok"
}, overlap_results)

log_message("Successful comparisons:", length(ok_results), "of",
            length(all_comparisons))
flush_log()

# ================================================================
# SECTION 6: FIGURE ASSEMBLY
# ================================================================
# Only desde_apo figures. Three figure types per mutant:
#   standard (Analysis 1), restricted (Analysis 2), per-chain (Analysis 3).
# Each produced in acotado and extendido versions.
# Pattern identical to 02C v3.0 Section 6A — panels_ov_* are plain
# lists of data.frames (or NULLs), built with c() to safely include NULLs.

log_message("============================================================")
log_message("Assembling diagnostic figures...")

# Helper: find result for a specific comparison by label_A and label_B
find_result <- function(results, lA, lB) {
  for (r in results) {
    if (!is.null(r$label_A) && r$label_A == lA &&
        !is.null(r$label_B) && r$label_B == lB) return(r)
  }
  return(NULL)
}

# --- 6A: Per-mutant figures (4 panels each): apo->stim, desde_apo only ---
# Structure (4 panels, top to bottom):
#   WT   apo->cap
#   WT   apo->heat
#   MUT  apo->cap
#   MUT  apo->heat

for (mut in mutants[mutants != "WT"]) {

  # Initialise collection lists for this mutant
  panels_ov_std  <- list()   # Analysis 1: ov_from_A (data.frame or NULL)
  panels_ov_rest <- list()   # Analysis 2: ov_restricted (data.frame or NULL)
  panels_ov_ch   <- list()   # Analysis 3: ov_chains (list of 4 df/NULL, or NULL)
  panels_meta    <- list()

  for (wt_or_mut in c("WT", mut)) {
    for (stim in c("cap", "heat")) {

      pdb_stim <- if (stim == "cap") pdb_cap else pdb_heat
      lA <- make_label(pdb_apo,  wt_or_mut)
      lB <- make_label(pdb_stim, wt_or_mut)

      res <- find_result(ok_results, lA, lB)

      panel_title <- paste0(
        wt_or_mut, "  apo -> ", stim,
        "  [NMA desde ", lA, "]"
      )

      # Use c(list, list(x)) to safely append even when x is NULL
      panels_ov_std  <- c(panels_ov_std,  list(if (!is.null(res)) res$ov_from_A     else NULL))
      panels_ov_rest <- c(panels_ov_rest, list(if (!is.null(res)) res$ov_restricted  else NULL))
      panels_ov_ch   <- c(panels_ov_ch,   list(if (!is.null(res)) res$ov_chains      else NULL))
      panels_meta    <- c(panels_meta,    list(list(title = panel_title, res = res)))
    }
  }

  for (version in c("acotado", "extendido")) {
    mode_range_use <- if (version == "acotado") plot_mode_range else plot_mode_range_extended

    # --- Analysis 1: Standard (desde_apo, full range) ---
    y_max_std <- compute_common_ymax(panels_ov_std, mode_range_use)
    panels_std <- lapply(seq_along(panels_meta), function(idx) {
      meta  <- panels_meta[[idx]]
      ov_df <- panels_ov_std[[idx]]
      if (is.null(meta$res) || is.null(ov_df)) {
        ggplot() + labs(title = paste(meta$title, "— NO DATA")) +
          theme_void() + theme(plot.title = element_text(size = 11, color = "red"))
      } else {
        make_overlap_panel(ov_df, meta$title, mode_range_use, y_max_std)
      }
    })
    save_multipanel(panels_std,
                    file.path(output_root, "figures_standard", version,
                              paste0("WT_vs_", mut, "_standard_", version, ".png")),
                    fig_height_px, fig_width_px, fig_res)

    # --- Analysis 2: Restricted range ---
    y_max_rest <- compute_common_ymax(panels_ov_rest, mode_range_use)
    panels_rest <- lapply(seq_along(panels_meta), function(idx) {
      meta  <- panels_meta[[idx]]
      ov_df <- panels_ov_rest[[idx]]
      title_r <- paste0(meta$title, " [restricted 400-715]")
      if (is.null(meta$res) || is.null(ov_df)) {
        ggplot() + labs(title = paste(title_r, "— NO DATA")) +
          theme_void() + theme(plot.title = element_text(size = 11, color = "red"))
      } else {
        make_overlap_panel(ov_df, title_r, mode_range_use, y_max_rest)
      }
    })
    save_multipanel(panels_rest,
                    file.path(output_root, "figures_restricted", version,
                              paste0("WT_vs_", mut, "_restricted_", version, ".png")),
                    fig_height_px, fig_width_px, fig_res)

    # --- Analysis 3: Per-chain ---
    # panels_ov_ch is a list of 4 elements, each is either a list of 4 df/NULLs or NULL.
    # Flatten one level to get a plain list of df/NULLs for compute_common_ymax.
    all_chain_dfs <- do.call(c, lapply(panels_ov_ch, function(ch_list) {
      if (is.null(ch_list)) return(list(NULL, NULL, NULL, NULL))
      ch_list  # already list of 4 data.frames or NULLs
    }))
    y_max_ch <- compute_common_ymax(all_chain_dfs, mode_range_use)

    panels_ch <- lapply(seq_along(panels_meta), function(idx) {
      meta       <- panels_meta[[idx]]
      ov_ch_list <- panels_ov_ch[[idx]]
      title_c    <- paste0(meta$title, " [per-chain diagnostic]")
      if (is.null(meta$res) || is.null(ov_ch_list)) {
        ggplot() + labs(title = paste(title_c, "— NO DATA")) +
          theme_void() + theme(plot.title = element_text(size = 11, color = "red"))
      } else {
        make_perchain_panel(ov_ch_list, chain_labels, chain_colors,
                            title_c, mode_range_use, y_max_ch)
      }
    })
    save_multipanel(panels_ch,
                    file.path(output_root, "figures_perchain", version,
                              paste0("WT_vs_", mut, "_perchain_", version, ".png")),
                    fig_height_px, fig_width_px, fig_res)
  }
}

# ================================================================
# SECTION 7: COMPLETION
# ================================================================

log_message("============================================================")
log_message("=== 02C_DIAG v3.0 complete ===")
log_message("Outputs:")
log_message("  Standard figures    ->", file.path(output_root, "figures_standard"))
log_message("  Restricted figures  ->", file.path(output_root, "figures_restricted"))
log_message("  Per-chain figures   ->", file.path(output_root, "figures_perchain"))
log_message("NOTE: Per-chain assignment is position-based (arbitrary).")
log_message("NOTE: Restricted range uses residues 400-602 and 625-715.")
log_message("KNOWN LIMITATION: chain symmetry not corrected (see header).")
flush_log(append = FALSE)

# End of script 02C_DIAG v3.0
