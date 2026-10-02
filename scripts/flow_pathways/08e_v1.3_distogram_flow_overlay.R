# ===============================================================
# 08e_v1.3_overlay_distogram_flow.R
# SUPERPOSICIÓN DE DISTOGRAMAS (07c) Y MATRICES DE FLUJO (08d)
# Author: DenyCB
# ---------------------------------------------------------------
# OBJETIVO:
# Combina en un solo gráfico:
#   - Fondo: distograma de distancias medias (07c, rev(viridis))
#     o delta/delta-delta de distancias (rojo-gris-azul)
#   - Superposición: flujo de información (08d), representado
#     como cruces (pch=3) coloreadas por valor de flujo (per system)
#     o por categoría aparece/desaparece/mantiene (delta/dd)
#
# ANÁLISIS GENERADOS:
#   1) Per system: distancia media + flujo absoluto (cruces blanco→rojo→azul)
#   2) Delta (Cap-Apo, Heat-Apo): delta distancia + flujo aparece/desaparece
#   3) Delta-delta: dd distancia + flujo aparece/desaparece
#   4) Cap-Heat versión 1: dist(7LPB) - dist(7LPC), directo
#   5) Cap-Heat versión 2: (Cap-Apo) - (Heat-Apo), desde deltas
#      En Cap-Heat: triángulo inferior = flujos Cap-Apo,
#                   triángulo superior = flujos Heat-Apo
#
# INPUTS:
#   - Distancias: *_combined_mean_dist.rds y *_combined_normvar_dist.rds
#     desde la carpeta de rango completo del 07c (BASE_DIST_DIR)
#   - Flujos: matrices CSV desde DIR_MATRICES (carpeta matrices/ del 08d)
#
# OUTPUTS:
#   - PNGs en subcarpetas por zoom dentro de OUTPUT_DIR
#
# FILOSOFÍA:
#   - Lee inputs preprocesados, no recalcula nada
#   - Loop sobre ZOOM_RANGES idéntico al 07c y 08d
#   - Paralelización por sistema dentro de cada zoom (parallel, Win10)


#OJO, buscar dist_z = dist_z, flow_z = NULL, reemplazarlo por dist_z = dist_z, flow_z = flow_z,
# ===============================================================

suppressPackageStartupMessages({
  library(viridisLite)
  library(fields)
  library(parallel)
  library(data.table)
})

# ===============================================================
# 🔧 SETTINGS
# ===============================================================

# --- Carpeta base del 07c con rango completo (fuente de .rds de distancias) ---
BASE_DIST_DIR <- "C:/DinamicasMoleculares/TRPV1_pipeline/output/07_MD_distances/07c_distograms_v2.0_1/A300_D752_vs_A300_D752/distograms"

# --- Carpeta de matrices de flujo del 08d ---
DIR_MATRICES <- "C:/DinamicasMoleculares/TRPV1_pipeline/output/08_Max_Flow/08b_v3.1_flow_paths_MD_min7/matrices"

# --- Output ---
OUTPUT_DIR <- "C:/DinamicasMoleculares/TRPV1_pipeline/output/08_Max_Flow/08e_v1.3_overlay_MD_soloDistogram"

# --- Suffix para nombres de archivos ---
SUFFIX_SYSTEM <- "_MD"


RUN_PER_SYSTEM  <- TRUE
RUN_DELTA       <- FALSE
RUN_DELTADELTA  <- FALSE
RUN_CAP_HEAT    <- FALSE



# --- Full canvas expansion ---
# Debe coincidir con lo usado en el 08d al generar las matrices de flujo
EXPAND_TO_FULL_CANVAS <- TRUE
CANVAS_MAP_FILE_MD  <- "C:/DinamicasMoleculares/TRPV1_pipeline/output/residue_map_MD.csv"
CANVAS_MAP_FILE_NMA <- "C:/DinamicasMoleculares/TRPV1_pipeline/config_files/trimmed_residue_map.csv"
CANVAS_GAP_RESNUM   <- 603:624

# --- Zoom ranges (idénticos al 07c y 08d) ---
ZOOM_RANGES <- list(
  list(X = c("A430","A560"), Y = c("A430","A560")),
  list(X = c("A300","D752"), Y = c("A300","D752")),
  list(X = c("A300","A752"), Y = c("D300","D752")),
  #  list(X = c("B300","C752"), Y = c("B300","C752")),
  list(X = c("A300","A752"), Y = c("A300","A752")),
  #  list(X = c("A300","A752"), Y = c("B300","C752")),
  #  list(X = c("B300","C752"), Y = c("D300","D752")),
  #  list(X = c("D300","D752"), Y = c("D300","D752")),
  list(X = c("A400","A590"), Y = c("D560","D690")),
  list(X = c("A400","A560"), Y = c("A400","A560"))
)

# --- PNG settings ---
PNG_W   <- 3600
PNG_H   <- 3200
PNG_DPI <- 600

# --- Colores ---
N_COLS <- 256

# --- Tamaño y tipo de marca para flujos superpuestos ---
# pch=3 es cruz (+). Para cambiar a pixel sólido: pch=15
# Para ajustar grosor de la cruz: lwd dentro de points() — ver función draw_overlay
FLOW_PCH <- 3
FLOW_CEX <- 0.2   # tamaño de la cruz; aumentar para hacerla más visible

# --- Colores de flujo para delta/dd ---
# Cap-Apo: aparece = verde, desaparece = cyan
COL_CAP_APARECE    <- "green"
COL_CAP_DESAPARECE <- "cyan"
# Heat-Apo: aparece = darkgreen, desaparece = "mediumorchid" (lila)
COL_HEAT_APARECE    <- "darkgreen"
COL_HEAT_DESAPARECE <- "mediumorchid"
# Mantiene (presente en ambos): gris
COL_MANTIENE <- "gray70"

# --- Umbral de flujo para considerar un edge como presente ---
# (mismo criterio que make_delta del 08d: 0.05% del máximo)
FLOW_THRESH_FRAC <- 0.0005

# --- Epsilon distancias ---
EPSILON <- 1e-12

# --- Overwrite ---
OVERWRITE <- TRUE

# --- Cores ---
N_CORES <- 6

# ===============================================================
# residuos especiales
# ===============================================================

special_res1 <- c("A426","A697","A441","B426","B697","B441","C426","C697","C441","D426","D697","D441","A563","B563","C563","D563")
special_res2 <- c("A679","A643","B679","B643","C679","C643","D679","D643")
special_res3 <- c("B300","C300","D300")
special_res4 <- c("D593","A543","B543","A593")

# ===============================================================
# CREATE OUTPUT FOLDER + LOGGING
# ===============================================================

if (!dir.exists(OUTPUT_DIR)) dir.create(OUTPUT_DIR, recursive = TRUE)

LOG_FILE <- file.path(
  OUTPUT_DIR,
  paste0("08e_v1.3_log_", format(Sys.time(), "%Y%m%d_%H%M%S"), ".txt")
)

log_msg <- function(...) {
  txt <- paste0("[", format(Sys.time(), "%Y-%m-%d %H:%M:%S"), "] ", paste(..., collapse = ""))
  cat(txt, "\n")
  cat(txt, "\n", file = LOG_FILE, append = TRUE)
}

# ===============================================================
# HELPERS — RANGO DE NODOS
# ===============================================================

subset_nodes <- function(node_names, range_vec) {
  if (is.null(range_vec) || length(range_vec) != 2) return(seq_along(node_names))
  i1 <- match(range_vec[1], node_names)
  i2 <- match(range_vec[2], node_names)
  if (is.na(i1) || is.na(i2)) stop(paste0("Range labels not found: ", paste(range_vec, collapse = " to ")))
  seq(min(i1, i2), max(i1, i2))
}

# ===============================================================
# HELPERS — LECTURA DE MATRICES DE FLUJO (idéntico al 08d)
# ===============================================================

read_matrix_csv <- function(f) {
  x <- read.csv(f, row.names = NULL, check.names = FALSE)
  if (colnames(x)[1] == "") {
    m <- as.matrix(x[, -1])
    rownames(m) <- x[, 1]
  } else {
    m <- as.matrix(x)
    rownames(m) <- colnames(m)
  }
  storage.mode(m) <- "numeric"
  m
}

get_system_from_matrix_filename <- function(nm) {
  nm2 <- sub("^matrix_", "", nm)
  parts <- strsplit(nm2, "_")[[1]]
  paste(parts[1], parts[2], sep = "_")
}

# ===============================================================
# HELPERS — CLASIFICACIÓN DE SISTEMA (idéntico al 07c)
# ===============================================================

classify_system <- function(sysname) {
  parts <- strsplit(sysname, "_")[[1]]
  pdb <- parts[1]
  mut <- parts[2]
  cond <- ifelse(pdb == "7LP9", "Apo",
                 ifelse(pdb == "7LPB", "Cap", "Heat"))
  list(pdb = pdb, mut = mut, cond = cond)
}

# ===============================================================
# HELPERS — RESIDUOS ESPECIALES (líneas sobre el gráfico)
# ===============================================================

draw_special_res <- function(mat) {
  ny <- nrow(mat)
  for (v in list(
    list(special_res1, "orange"),
    list(special_res2, "magenta"),
    list(special_res3, "darkgray"),
    list(special_res4, "black")
  )) {
    idx_x <- match(v[[1]], colnames(mat))
    idx_x <- idx_x[!is.na(idx_x)]
    idx_y <- match(v[[1]], rownames(mat))
    idx_y <- idx_y[!is.na(idx_y)]
    if (length(idx_x) > 0) abline(v = idx_x, col = v[[2]], lwd = 0.5, lty = 2)
    if (length(idx_y) > 0) abline(h = ny - idx_y + 1, col = v[[2]], lwd = 0.5, lty = 2)
  }
}

# ===============================================================
# HELPERS — SUPERPOSICIÓN DE FLUJO PER SYSTEM
# Cruces coloreadas blanco→rojo→azul según valor de flujo.
# Para cambiar símbolo: modificar FLOW_PCH en settings.
# Para cambiar grosor de la cruz: modificar lwd en points() abajo.
# ===============================================================

draw_flow_persys <- function(flow_z, max_flow) {
  ny <- nrow(flow_z)
  cols_flow <- colorRampPalette(c("white", "red", "blue"))(100)
  for (col_i in seq_len(ncol(flow_z))) {
    for (row_i in seq_len(nrow(flow_z))) {
      val <- flow_z[row_i, col_i]
      if (is.na(val) || val <= 0) next
      col_idx <- max(1, min(100, round(val / max_flow * 99) + 1))
      points(
        x   = col_i,
        y   = ny - row_i + 1,
        pch = FLOW_PCH,
        cex = FLOW_CEX,
        col = cols_flow[col_idx],
        lwd = 0.5   # grosor de la cruz; aumentar para hacerla más gruesa
      )
    }
  }
}

# ===============================================================
# HELPERS — SUPERPOSICIÓN DE FLUJO DELTA (aparece/desaparece/mantiene)
# col_aparece / col_desaparece: colores para esta comparación
# triangle: "lower" dibuja solo triángulo inferior (col <= row en imagen)
#           "upper" dibuja solo triángulo superior
#           NULL dibuja toda la matriz
# ===============================================================

draw_flow_delta <- function(flow_z, col_aparece, col_desaparece, triangle = NULL) {
  ny  <- nrow(flow_z)
  ncx <- ncol(flow_z)
  for (col_i in seq_len(ncx)) {
    for (row_i in seq_len(ny)) {
      val <- flow_z[row_i, col_i]
      if (is.na(val)) next
      
      # posición en imagen (y invertido)
      img_y <- ny - row_i + 1
      
      # filtro de triángulo
      if (!is.null(triangle)) {
        if (triangle == "lower" && col_i > (ny + 1 - img_y)) next
        if (triangle == "upper" && col_i <= (ny + 1 - img_y)) next
      }
      
      if (val == 999) {
        col_use <- col_aparece
      } else if (val == -999) {
        col_use <- col_desaparece
      } else {
        col_use <- COL_MANTIENE
      }
      
      points(
        x   = col_i,
        y   = img_y,
        pch = FLOW_PCH,
        cex = FLOW_CEX,
        col = col_use,
        lwd = 0.5   # grosor de la cruz; aumentar para hacerla más gruesa
      )
    }
  }
}

# ===============================================================
# HELPERS — FUNCIÓN GENÉRICA DE FONDO (distancia o delta)
# ===============================================================

open_bg_distogram <- function(mat, outfile, main_txt, legend_txt,
                              zlim = NULL, invert_col = FALSE,
                              is_delta = FALSE) {
  png(outfile, width = PNG_W, height = PNG_H, res = PNG_DPI)
  
  old_par <- par(no.readonly = TRUE)
  on.exit({ par(old_par); dev.off() }, add = TRUE)
  
  if (is.null(zlim)) {
    zlim <- range(mat, na.rm = TRUE)
    if (!all(is.finite(zlim))) zlim <- c(0, 1)
    if (diff(zlim) == 0) zlim <- zlim + c(-0.5, 0.5)
    if (is_delta) {
      mx <- max(abs(zlim))
      zlim <- c(-mx, mx)
    }
  }
  
  nx <- ncol(mat)
  ny <- nrow(mat)
  
  STEP_TICKS <- 10
  tick_x <- seq(1, nx, by = STEP_TICKS)
  tick_y <- seq(1, ny, by = STEP_TICKS)
  
  labels_x <- colnames(mat)
  labels_y <- rownames(mat)
  
  par(mar = c(4, 4, 4, 4))
  
  if (is_delta) {
    cols_bg <- colorRampPalette(c("red", "grey85", "blue"))(N_COLS)
  } else {
    cols_bg <- if (invert_col) rev(viridis(N_COLS)) else viridis(N_COLS)
  }
  
  fields::image.plot(
    x = seq_len(nx),
    y = seq_len(ny),
    z = t(mat[ny:1, , drop = FALSE]),
    col = cols_bg,
    axes = FALSE,
    xlab = "",
    ylab = "",
    main = main_txt,
    zlim = zlim,
    legend.lab = legend_txt,
    legend.line = 2.5
  )
  
  # STEP_TICKS redeclarado igual que en 07c original
  STEP_TICKS <- 20
  nx <- ncol(mat)
  ny <- nrow(mat)
  
  draw_special_res(mat)
  
  axis(1, at = tick_x, labels = labels_x[tick_x], las = 2, cex.axis = 0.4)
  axis(2, at = tick_y, labels = rev(labels_y)[tick_y], las = 2, cex.axis = 0.4)
  
  box()
  # NOTE: dev.off() es llamado por on.exit — el llamador dibuja encima ANTES de salir
  # Para poder dibujar encima del fondo, esta función NO llama dev.off() explícitamente.
  # on.exit() lo hará al retornar al llamador.
  # El llamador debe dibujar sus puntos DENTRO de esta función o usar environment tricks.
  # Por eso usamos el patrón: abrir png, dibujar fondo, dibujar overlay, cerrar.
  # Esta función retorna invisible el old_par; el cierre lo maneja on.exit.
  invisible(old_par)
}

# ===============================================================
# FUNCIÓN COMPLETA DE OVERLAY: fondo + puntos + cierre
# ===============================================================

plot_overlay <- function(dist_z, flow_z, outfile, main_txt, legend_txt,
                         mode = "persys",
                         max_flow = NULL,
                         flow_z2 = NULL,
                         col_aparece = COL_CAP_APARECE,
                         col_desaparece = COL_CAP_DESAPARECE,
                         col_aparece2 = COL_HEAT_APARECE,
                         col_desaparece2 = COL_HEAT_DESAPARECE,
                         is_delta = FALSE,
                         invert_col = FALSE) {
  
  png(outfile, width = PNG_W, height = PNG_H, res = PNG_DPI)
  
  old_par <- par(no.readonly = TRUE)
  on.exit({ par(old_par); dev.off() }, add = TRUE)
  
  # --- fondo ---
  zlim <- range(dist_z, na.rm = TRUE)
  if (!all(is.finite(zlim))) zlim <- c(0, 1)
  if (diff(zlim) == 0) zlim <- zlim + c(-0.5, 0.5)
  if (is_delta) {
    mx <- max(abs(zlim))
    zlim <- c(-mx, mx)
  }
  
  nx <- ncol(dist_z)
  ny <- nrow(dist_z)
  
  STEP_TICKS <- 10
  tick_x <- seq(1, nx, by = STEP_TICKS)
  tick_y <- seq(1, ny, by = STEP_TICKS)
  
  labels_x <- colnames(dist_z)
  labels_y <- rownames(dist_z)
  
  par(mar = c(4, 4, 4, 4))
  
  if (is_delta) {
    cols_bg <- colorRampPalette(c("red", "grey85", "blue"))(N_COLS)
  } else {
    cols_bg <- if (invert_col) rev(viridis(N_COLS)) else viridis(N_COLS)
  }
  
  fields::image.plot(
    x = seq_len(nx),
    y = seq_len(ny),
    z = t(dist_z[ny:1, , drop = FALSE]),
    col = cols_bg,
    axes = FALSE,
    xlab = "",
    ylab = "",
    main = main_txt,
    zlim = zlim,
    legend.lab = legend_txt,
    legend.line = 2.5
  )
  
  STEP_TICKS <- 20
  nx <- ncol(dist_z)
  ny <- nrow(dist_z)
  
  draw_special_res(dist_z)
  
  # --- overlay de flujo ---
  if (!is.null(flow_z)) {
    if (mode == "persys") {
      if (is.null(max_flow) || max_flow == 0) max_flow <- max(flow_z, na.rm = TRUE)
      if (is.na(max_flow) || max_flow == 0) max_flow <- 1
      draw_flow_persys(flow_z, max_flow)
    } else if (mode == "delta") {
      # un solo set de flujos, toda la matriz
      draw_flow_delta(flow_z, col_aparece, col_desaparece, triangle = NULL)
    } else if (mode == "cap_heat") {
      # dos sets: flow_z = Cap-Apo (triángulo inferior imagen: col >= img_y)
      #           flow_z2 = Heat-Apo (triángulo superior imagen: col <= img_y)
      draw_flow_delta(flow_z,  col_aparece,  col_desaparece,  triangle = "lower")
      if (!is.null(flow_z2)) {
        draw_flow_delta(flow_z2, col_aparece2, col_desaparece2, triangle = "upper")
      }
    }
  }
  
  # Ticks y tamaño de letra proporcionales al tamaño del rango graficado.
  # Para rangos grandes (~450 nodos): cada 25 residuos, cex.axis reducido.
  # Para rangos pequeños (~160 nodos): cada 10 residuos, cex.axis normal.
  # Función continua interpolada entre los dos extremos de referencia.
  n_range      <- max(nx, ny)
  step_ticks   <- max(10, round(10 + (n_range - 160) / (450 - 160) * (25 - 10)))
  cex_ax       <- max(0.28, 0.4 - (n_range - 160) / (450 - 160) * (0.4 - 0.28))
  tick_x_final <- seq(1, nx, by = step_ticks)
  tick_y_final <- seq(1, ny, by = step_ticks)
  
  # ticks menores: un tick entre cada par de labels mayores (step_ticks / 2)
  step_minor   <- max(1, round(step_ticks / 2))
  minor_x      <- seq(1, nx, by = step_minor)
  minor_y      <- seq(1, ny, by = step_minor)
  # quitar de los menores los que ya son mayores (para no solapar)
  minor_x      <- setdiff(minor_x, tick_x_final)
  minor_y      <- setdiff(minor_y, tick_y_final)
  
  axis(1, at = tick_x_final, labels = labels_x[tick_x_final], las = 2, cex.axis = cex_ax)
  axis(2, at = tick_y_final, labels = rev(labels_y)[tick_y_final], las = 2, cex.axis = cex_ax)
  axis(1, at = minor_x, labels = FALSE, tcl = -0.25)
  axis(2, at = minor_y, labels = FALSE, tcl = -0.25)
  
  box()
}

# ===============================================================
# HELPERS — MAKE_DELTA DE FLUJO (idéntico al 08d v2.0)
# ===============================================================

make_delta_flow <- function(a, b) {
  all_nodes <- sort(unique(c(rownames(a), rownames(b))))
  expand_mat <- function(m, nodes) {
    m_exp <- matrix(0, length(nodes), length(nodes), dimnames = list(nodes, nodes))
    common <- intersect(nodes, rownames(m))
    m_exp[common, common] <- m[common, common]
    m_exp
  }
  a_exp <- expand_mat(a, all_nodes)
  b_exp <- expand_mat(b, all_nodes)
  a_exp[a_exp == 0] <- NA
  b_exp[b_exp == 0] <- NA
  max_combined <- max(c(a_exp, b_exp), na.rm = TRUE)
  thresh <- FLOW_THRESH_FRAC * max_combined
  na_both <- (is.na(a_exp) | a_exp < thresh) & (is.na(b_exp) | b_exp < thresh)
  na_ref  <- (is.na(a_exp) | a_exp < thresh) & !is.na(b_exp) & b_exp >= thresh
  na_cmp  <- !is.na(a_exp) & a_exp >= thresh & (is.na(b_exp) | b_exp < thresh)
  a_exp[is.na(a_exp) | a_exp < thresh] <- 0
  b_exp[is.na(b_exp) | b_exp < thresh] <- 0
  out <- b_exp - a_exp
  out[na_ref]  <-  999
  out[na_cmp]  <- -999
  out[na_both] <- NA
  diag(out) <- NA
  out
}

# ===============================================================
# START
# ===============================================================

log_msg("============================================================")
log_msg("STARTING SCRIPT 08e_v1.3_overlay_distogram_flow.R")
log_msg("BASE_DIST_DIR: ", BASE_DIST_DIR)
log_msg("DIR_MATRICES: ", DIR_MATRICES)
log_msg("OUTPUT_DIR: ", OUTPUT_DIR)

# ===============================================================
# LEER MATRICES DE DISTANCIA (desde BASE_DIST_DIR)
# ===============================================================

log_msg("------------------------------------------------------------")
log_msg("Reading distance matrices from BASE_DIST_DIR...")

dist_mean_files <- list.files(BASE_DIST_DIR, pattern = "_combined_mean_dist\\.rds$", full.names = TRUE)

if (length(dist_mean_files) == 0) stop("No _combined_mean_dist.rds files found in BASE_DIST_DIR.")

dist_systems <- list()

for (f in dist_mean_files) {
  bn  <- basename(f)
  sys <- sub("_combined_mean_dist\\.rds$", "", bn)
  mean_full    <- readRDS(f)
  normvar_file <- file.path(BASE_DIST_DIR, paste0(sys, "_combined_normvar_dist.rds"))
  normvar_full <- if (file.exists(normvar_file)) readRDS(normvar_file) else NULL
  dist_systems[[sys]] <- list(mean = mean_full, normvar = normvar_full,
                              info = classify_system(sys))
  log_msg("  Loaded distance matrix: ", sys,
          " | nodes: ", nrow(mean_full))
}

log_msg("Distance systems loaded: ", length(dist_systems))

# ===============================================================
# LEER Y SUMAR MATRICES DE FLUJO (idéntico al 08d)
# ===============================================================

log_msg("------------------------------------------------------------")
log_msg("Reading and summing flow matrices...")

ff_all <- list.files(DIR_MATRICES, pattern = "\\.csv$", full.names = TRUE)
if (length(ff_all) == 0) stop("No CSV files found in DIR_MATRICES.")

system_names_all <- sapply(
  tools::file_path_sans_ext(basename(ff_all)),
  get_system_from_matrix_filename
)

unique_systems_flow <- unique(system_names_all)
log_msg("Unique flow systems: ", length(unique_systems_flow))

flow_mats <- list()

for (sys in unique_systems_flow) {
  idx <- which(system_names_all == sys)
  mat_sum <- NULL
  for (i in idx) {
    m <- tryCatch(read_matrix_csv(ff_all[i]), error = function(e) NULL)
    if (is.null(m)) next
    if (is.null(mat_sum)) {
      mat_sum <- m
    } else {
      all_nodes <- sort(unique(c(rownames(mat_sum), rownames(m))))
      expand_to <- function(x, nodes) {
        x_exp <- matrix(0, length(nodes), length(nodes), dimnames = list(nodes, nodes))
        cn <- intersect(nodes, rownames(x))
        x_exp[cn, cn] <- x[cn, cn]
        x_exp
      }
      mat_sum <- expand_to(mat_sum, all_nodes) + expand_to(m, all_nodes)
    }
  }
  if (!is.null(mat_sum)) {
    flow_mats[[sys]] <- mat_sum
    log_msg("  Flow matrix summed: ", sys, " | nodes: ", nrow(mat_sum))
  }
}

# ===============================================================
# FULL CANVAS EXPANSION DE FLUJOS (idéntico al 08d)
# ===============================================================

global_nodes_flow <- sort(unique(unlist(lapply(flow_mats, rownames))))

expand_to_canvas <- function(m, nodes) {
  m_exp <- matrix(0, length(nodes), length(nodes), dimnames = list(nodes, nodes))
  cn <- intersect(nodes, rownames(m))
  m_exp[cn, cn] <- m[cn, cn]
  m_exp
}

flow_mats_canvas <- lapply(flow_mats, expand_to_canvas, nodes = global_nodes_flow)

if (EXPAND_TO_FULL_CANVAS) {
  canvas_map_file <- if (SUFFIX_SYSTEM == "_NMA") CANVAS_MAP_FILE_NMA else CANVAS_MAP_FILE_MD
  log_msg("EXPAND_TO_FULL_CANVAS = TRUE | map: ", canvas_map_file)
  rmap_canvas   <- read.csv(canvas_map_file, stringsAsFactors = FALSE)
  canvas_labels <- paste0(rmap_canvas$chain, rmap_canvas$residue)
  gap_labels_canvas <- as.vector(outer(c("A","B","C","D"), CANVAS_GAP_RESNUM, paste0))
  canvas_labels <- canvas_labels[!(canvas_labels %in% gap_labels_canvas)]
  canvas_chain  <- substr(canvas_labels, 1, 1)
  canvas_resnum <- as.numeric(substring(canvas_labels, 2))
  canvas_labels <- canvas_labels[order(canvas_chain, canvas_resnum)]
  log_msg("  Full canvas nodes (after gap removal): ", length(canvas_labels))
  expand_to_full <- function(m, nodes) {
    m_exp <- matrix(0, length(nodes), length(nodes), dimnames = list(nodes, nodes))
    cn <- intersect(nodes, rownames(m))
    m_exp[cn, cn] <- m[cn, cn]
    m_exp
  }
  flow_mats_canvas <- lapply(flow_mats_canvas, expand_to_full, nodes = canvas_labels)
  log_msg("  flow_mats_canvas expanded to full canvas.")
}

# ===============================================================
# ÍNDICE DE SISTEMAS POR MUTANTE Y CONDICIÓN
# ===============================================================

sys_index_dist <- list()
for (sys in names(dist_systems)) {
  info <- dist_systems[[sys]]$info
  mut  <- info$mut
  cond <- info$cond
  if (is.null(sys_index_dist[[mut]])) sys_index_dist[[mut]] <- list()
  sys_index_dist[[mut]][[cond]] <- sys
}

sys_index_flow <- list()
for (sys in names(flow_mats_canvas)) {
  info <- classify_system(sys)
  mut  <- info$mut
  cond <- info$cond
  if (is.null(sys_index_flow[[mut]])) sys_index_flow[[mut]] <- list()
  sys_index_flow[[mut]][[cond]] <- sys
}

# ===============================================================
# CALCULAR DELTAS DE DISTANCIA (Cap-Apo, Heat-Apo, Cap-Heat v1, v2)
# ===============================================================

log_msg("------------------------------------------------------------")
log_msg("Computing distance deltas...")

delta_dist_mats <- list()  # nombre: e.g. "WT_delta_Cap_mean"
capHeat_v1_mats <- list()  # nombre: e.g. "WT_capHeat_v1_mean"
capHeat_v2_mats <- list()  # nombre: e.g. "WT_capHeat_v2_mean"

for (mut in names(sys_index_dist)) {
  apo_sys  <- sys_index_dist[[mut]][["Apo"]]
  cap_sys  <- sys_index_dist[[mut]][["Cap"]]
  heat_sys <- sys_index_dist[[mut]][["Heat"]]
  
  for (mat_type in c("mean", "normvar")) {
    mat_apo  <- if (!is.null(apo_sys))  dist_systems[[apo_sys]][[mat_type]]  else NULL
    mat_cap  <- if (!is.null(cap_sys))  dist_systems[[cap_sys]][[mat_type]]  else NULL
    mat_heat <- if (!is.null(heat_sys)) dist_systems[[heat_sys]][[mat_type]] else NULL
    
    # Cap-Apo
    if (!is.null(mat_apo) && !is.null(mat_cap)) {
      delta_dist_mats[[paste0(mut, "_delta_Cap_", mat_type)]] <- mat_cap - mat_apo
      log_msg("  Delta Cap-Apo computed: ", mut, " ", mat_type)
    }
    
    # Heat-Apo
    if (!is.null(mat_apo) && !is.null(mat_heat)) {
      delta_dist_mats[[paste0(mut, "_delta_Heat_", mat_type)]] <- mat_heat - mat_apo
      log_msg("  Delta Heat-Apo computed: ", mut, " ", mat_type)
    }
    
    # Cap-Heat v1: directo
    if (!is.null(mat_cap) && !is.null(mat_heat)) {
      capHeat_v1_mats[[paste0(mut, "_capHeat_v1_", mat_type)]] <- mat_cap - mat_heat
      log_msg("  Cap-Heat v1 computed: ", mut, " ", mat_type)
    }
    
    # Cap-Heat v2: (Cap-Apo) - (Heat-Apo)
    cap_apo_name  <- paste0(mut, "_delta_Cap_",  mat_type)
    heat_apo_name <- paste0(mut, "_delta_Heat_", mat_type)
    if (!is.null(delta_dist_mats[[cap_apo_name]]) && !is.null(delta_dist_mats[[heat_apo_name]])) {
      capHeat_v2_mats[[paste0(mut, "_capHeat_v2_", mat_type)]] <-
        delta_dist_mats[[cap_apo_name]] - delta_dist_mats[[heat_apo_name]]
      log_msg("  Cap-Heat v2 computed: ", mut, " ", mat_type)
    }
  }
}

# ===============================================================
# CALCULAR DELTAS DE FLUJO (Cap-Apo, Heat-Apo por mutante)
# ===============================================================

log_msg("------------------------------------------------------------")
log_msg("Computing flow deltas...")

delta_flow_mats <- list()  # nombre: e.g. "WT_flow_delta_Cap"

for (mut in names(sys_index_flow)) {
  apo_sys  <- sys_index_flow[[mut]][["Apo"]]
  cap_sys  <- sys_index_flow[[mut]][["Cap"]]
  heat_sys <- sys_index_flow[[mut]][["Heat"]]
  
  if (!is.null(apo_sys) && !is.null(cap_sys)) {
    delta_flow_mats[[paste0(mut, "_flow_delta_Cap")]] <-
      make_delta_flow(flow_mats_canvas[[apo_sys]], flow_mats_canvas[[cap_sys]])
    log_msg("  Flow delta Cap-Apo: ", mut)
  }
  
  if (!is.null(apo_sys) && !is.null(heat_sys)) {
    delta_flow_mats[[paste0(mut, "_flow_delta_Heat")]] <-
      make_delta_flow(flow_mats_canvas[[apo_sys]], flow_mats_canvas[[heat_sys]])
    log_msg("  Flow delta Heat-Apo: ", mut)
  }
}

# ===============================================================
# CALCULAR DELTA-DELTA DE DISTANCIA Y FLUJO
# ===============================================================

log_msg("------------------------------------------------------------")
log_msg("Computing delta-delta matrices...")

dd_dist_mats <- list()   # nombre: e.g. "W426A_dd_Cap_mean"
dd_flow_mats <- list()   # nombre: e.g. "W426A_dd_flow_Cap"

mutants_non_wt <- setdiff(names(sys_index_dist), "WT")

for (mut in mutants_non_wt) {
  for (cond in c("Cap", "Heat")) {
    for (mat_type in c("mean", "normvar")) {
      wt_name  <- paste0("WT_delta_",  cond, "_", mat_type)
      mut_name <- paste0(mut, "_delta_", cond, "_", mat_type)
      wt_d  <- delta_dist_mats[[wt_name]]
      mut_d <- delta_dist_mats[[mut_name]]
      if (!is.null(wt_d) && !is.null(mut_d)) {
        dd_dist_mats[[paste0(mut, "_dd_", cond, "_", mat_type)]] <- mut_d - wt_d
        log_msg("  DD dist computed: ", mut, " ", cond, " ", mat_type)
      }
    }
    wt_flow_name  <- paste0("WT_flow_delta_",  cond)
    mut_flow_name <- paste0(mut, "_flow_delta_", cond)
    wt_f  <- delta_flow_mats[[wt_flow_name]]
    mut_f <- delta_flow_mats[[mut_flow_name]]
    if (!is.null(wt_f) && !is.null(mut_f)) {
      # dd de flujo: mismo make_delta_flow aplicado sobre los deltas
      # (conserva lógica de aparece/desaparece/mantiene)
      dd_flow_mats[[paste0(mut, "_dd_flow_", cond)]] <-
        make_delta_flow(wt_f, mut_f)
      log_msg("  DD flow computed: ", mut, " ", cond)
    }
  }
}

# ===============================================================
# FUNCIÓN DE ALINEACIÓN DE MATRICES DIST Y FLOW AL MISMO ESPACIO
# ===============================================================

align_to_dist <- function(flow_mat, dist_mat) {
  # Expande flow_mat para que tenga exactamente los mismos nodos que dist_mat
  nodes <- rownames(dist_mat)
  m_exp <- matrix(0, length(nodes), length(nodes), dimnames = list(nodes, nodes))
  cn <- intersect(nodes, rownames(flow_mat))
  m_exp[cn, cn] <- flow_mat[cn, cn]
  m_exp
}

align_delta_to_dist <- function(flow_mat, dist_mat) {
  # Igual que align_to_dist pero preserva NA (para deltas de flujo)
  nodes <- rownames(dist_mat)
  m_exp <- matrix(NA, length(nodes), length(nodes), dimnames = list(nodes, nodes))
  cn <- intersect(nodes, rownames(flow_mat))
  m_exp[cn, cn] <- flow_mat[cn, cn]
  m_exp
}

# ===============================================================
# ZOOM LOOP
# ===============================================================

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
  
  zoom_out       <- file.path(OUTPUT_DIR, zoom_tag)
  dir_persys     <- file.path(zoom_out, "per_system")
  dir_delta      <- file.path(zoom_out, "delta")
  dir_dd         <- file.path(zoom_out, "delta_delta")
  dir_capHeat    <- file.path(zoom_out, "cap_heat")
  
  for (d in c(dir_persys, dir_delta, dir_dd, dir_capHeat)) {
    if (!dir.exists(d)) dir.create(d, recursive = TRUE)
  }
  
  # ==========================================================
  # PER SYSTEM (paralelo)
  # ==========================================================
  
  if (RUN_PER_SYSTEM) {
  
  log_msg("Generating per-system overlays (parallel, ", N_CORES, " cores)...")
  
  GLOBAL_MAX_FLOW_PER <- max(unlist(lapply(flow_mats_canvas, function(m) m[m > 0])),
                             na.rm = TRUE)
  
  cl <- makeCluster(N_CORES)
  clusterExport(cl, varlist = c(
    "dist_systems", "flow_mats_canvas", "sys_index_flow",
    "X_RANGE", "Y_RANGE", "dir_persys",
    "PNG_W", "PNG_H", "PNG_DPI", "N_COLS", "OVERWRITE",
    "FLOW_PCH", "FLOW_CEX", "GLOBAL_MAX_FLOW_PER",
    "special_res1", "special_res2", "special_res3", "special_res4",
    "subset_nodes", "plot_overlay", "draw_special_res",
    "draw_flow_persys", "draw_flow_delta",
    "COL_CAP_APARECE", "COL_CAP_DESAPARECE",
    "COL_HEAT_APARECE", "COL_HEAT_DESAPARECE", "COL_MANTIENE",
    "align_to_dist", "align_delta_to_dist",
    "classify_system", "log_msg", "LOG_FILE"
  ), envir = environment())
  clusterEvalQ(cl, {
    suppressPackageStartupMessages({ library(viridisLite); library(fields) })
  })
  
  parLapply(cl, names(dist_systems), function(sys) {
    info     <- dist_systems[[sys]]$info
    mut      <- info$mut
    cond     <- info$cond
    flow_sys <- sys_index_flow[[mut]][[cond]]
    
    dist_mat <- dist_systems[[sys]]$mean
    nodes_x  <- subset_nodes(colnames(dist_mat), X_RANGE)
    nodes_y  <- subset_nodes(rownames(dist_mat), Y_RANGE)
    dist_z   <- dist_mat[nodes_y, nodes_x, drop = FALSE]
    
    flow_z <- NULL
    if (!is.null(flow_sys) && !is.null(flow_mats_canvas[[flow_sys]])) {
      flow_full <- align_to_dist(flow_mats_canvas[[flow_sys]], dist_mat)
      flow_z    <- flow_full[nodes_y, nodes_x, drop = FALSE]
    }
    
    outfile <- file.path(dir_persys, paste0(sys, "_overlay_persys.png"))
    if (!file.exists(outfile) || OVERWRITE) {
      plot_overlay(
        dist_z = dist_z, flow_z = flow_z,
        outfile = outfile,
        main_txt = paste0(sys, " | Distance + Flow overlay"),
        legend_txt = "Distance (Å)",
        mode = "persys",
        max_flow = GLOBAL_MAX_FLOW_PER,
        invert_col = TRUE,
        is_delta = FALSE
      )
    }
    
    normvar_full <- dist_systems[[sys]]$normvar
    if (!is.null(normvar_full)) {
      normvar_z    <- normvar_full[nodes_y, nodes_x, drop = FALSE]
      outfile_nv   <- file.path(dir_persys, paste0(sys, "_overlay_persys_normvar.png"))
      if (!file.exists(outfile_nv) || OVERWRITE) {
        plot_overlay(
          dist_z = normvar_z, flow_z = flow_z,
          outfile = outfile_nv,
          main_txt = paste0(sys, " | Normvar + Flow overlay"),
          legend_txt = "Var / Mean^2",
          mode = "persys",
          max_flow = GLOBAL_MAX_FLOW_PER,
          invert_col = FALSE,
          is_delta = FALSE
        )
      }
    }
    return(invisible(NULL))
  })
  
  stopCluster(cl)
  log_msg("Per-system overlays done for zoom: ", zoom_tag)
  }
  
  # ==========================================================
  # DELTA OVERLAYS (Cap-Apo, Heat-Apo) — paralelo
  # ==========================================================
  
  if (RUN_DELTA) {
  
  log_msg("Generating delta overlays (parallel, ", N_CORES, " cores)...")
  
  cl <- makeCluster(N_CORES)
  clusterExport(cl, varlist = c(
    "delta_dist_mats", "delta_flow_mats", "dist_systems",
    "X_RANGE", "Y_RANGE", "dir_delta",
    "PNG_W", "PNG_H", "PNG_DPI", "N_COLS", "OVERWRITE",
    "FLOW_PCH", "FLOW_CEX",
    "special_res1", "special_res2", "special_res3", "special_res4",
    "subset_nodes", "plot_overlay", "draw_special_res",
    "draw_flow_persys", "draw_flow_delta",
    "COL_CAP_APARECE", "COL_CAP_DESAPARECE",
    "COL_HEAT_APARECE", "COL_HEAT_DESAPARECE", "COL_MANTIENE",
    "align_to_dist", "align_delta_to_dist",
    "classify_system", "log_msg", "LOG_FILE"
  ), envir = environment())
  clusterEvalQ(cl, {
    suppressPackageStartupMessages({ library(viridisLite); library(fields) })
  })
  
  delta_keys <- names(delta_dist_mats)
  
  parLapply(cl, delta_keys, function(dname) {
    # dname e.g. "WT_delta_Cap_mean"
    parts    <- strsplit(dname, "_")[[1]]
    mut      <- parts[1]
    cond     <- parts[3]   # Cap or Heat
    mat_type <- parts[4]   # mean or normvar
    
    dist_delta_full <- delta_dist_mats[[dname]]
    
    # referencia para alinear: dist_systems del sistema Apo del mutante
    # usamos rownames del delta (ya tienen el espacio correcto)
    nodes_x <- subset_nodes(colnames(dist_delta_full), X_RANGE)
    nodes_y <- subset_nodes(rownames(dist_delta_full), Y_RANGE)
    dist_z  <- dist_delta_full[nodes_y, nodes_x, drop = FALSE]
    
    flow_delta_name <- paste0(mut, "_flow_delta_", cond)
    flow_z <- NULL
    if (!is.null(delta_flow_mats[[flow_delta_name]])) {
      flow_full <- align_delta_to_dist(delta_flow_mats[[flow_delta_name]], dist_delta_full)
      flow_z    <- flow_full[nodes_y, nodes_x, drop = FALSE]
    }
    
    col_ap  <- if (cond == "Cap") COL_CAP_APARECE    else COL_HEAT_APARECE
    col_des <- if (cond == "Cap") COL_CAP_DESAPARECE else COL_HEAT_DESAPARECE
    
    leg_txt <- if (mat_type == "mean") "Delta Distance (Å)" else "Delta Var / Mean^2"
    
    outfile <- file.path(dir_delta, paste0(dname, "_overlay.png"))
    if (file.exists(outfile) && !OVERWRITE) return(invisible(NULL))
    
    plot_overlay(
      dist_z = dist_z, flow_z =NULL,# flow_z,
      outfile = outfile,
      main_txt = paste0(dname, " | Delta overlay"),
      legend_txt = leg_txt,
      mode = "delta",
      col_aparece = col_ap, col_desaparece = col_des,
      is_delta = TRUE, invert_col = FALSE
    )
    return(invisible(NULL))
  })
  
  stopCluster(cl)
  log_msg("Delta overlays done for zoom: ", zoom_tag)
  }
  
  
  # ==========================================================
  # DELTA-DELTA OVERLAYS — paralelo
  # ==========================================================
  
  if (RUN_DELTADELTA) {
  log_msg("Generating delta-delta overlays (parallel, ", N_CORES, " cores)...")
  
  cl <- makeCluster(N_CORES)
  clusterExport(cl, varlist = c(
    "dd_dist_mats", "dd_flow_mats",
    "X_RANGE", "Y_RANGE", "dir_dd",
    "PNG_W", "PNG_H", "PNG_DPI", "N_COLS", "OVERWRITE",
    "FLOW_PCH", "FLOW_CEX",
    "special_res1", "special_res2", "special_res3", "special_res4",
    "subset_nodes", "plot_overlay", "draw_special_res",
    "draw_flow_persys", "draw_flow_delta",
    "COL_CAP_APARECE", "COL_CAP_DESAPARECE",
    "COL_HEAT_APARECE", "COL_HEAT_DESAPARECE", "COL_MANTIENE",
    "align_to_dist", "align_delta_to_dist",
    "classify_system", "log_msg", "LOG_FILE"
  ), envir = environment())
  clusterEvalQ(cl, {
    suppressPackageStartupMessages({ library(viridisLite); library(fields) })
  })
  
  dd_keys <- names(dd_dist_mats)
  
  parLapply(cl, dd_keys, function(dname) {
    # dname e.g. "W426A_dd_Cap_mean"
    parts    <- strsplit(dname, "_")[[1]]
    mut      <- parts[1]
    cond     <- parts[3]
    mat_type <- parts[4]
    
    dist_dd_full <- dd_dist_mats[[dname]]
    
    nodes_x <- subset_nodes(colnames(dist_dd_full), X_RANGE)
    nodes_y <- subset_nodes(rownames(dist_dd_full), Y_RANGE)
    dist_z  <- dist_dd_full[nodes_y, nodes_x, drop = FALSE]
    
    flow_dd_name <- paste0(mut, "_dd_flow_", cond)
    flow_z <- NULL
    if (!is.null(dd_flow_mats[[flow_dd_name]])) {
      flow_full <- align_delta_to_dist(dd_flow_mats[[flow_dd_name]], dist_dd_full)
      flow_z    <- flow_full[nodes_y, nodes_x, drop = FALSE]
    }
    
    col_ap  <- if (cond == "Cap") COL_CAP_APARECE    else COL_HEAT_APARECE
    col_des <- if (cond == "Cap") COL_CAP_DESAPARECE else COL_HEAT_DESAPARECE
    
    leg_txt <- if (mat_type == "mean") "DD Distance (Å)" else "DD Var / Mean^2"
    
    outfile <- file.path(dir_dd, paste0(dname, "_overlay.png"))
    if (file.exists(outfile) && !OVERWRITE) return(invisible(NULL))
    
    plot_overlay(
      dist_z = dist_z, flow_z =NULL,# flow_z,
      outfile = outfile,
      main_txt = paste0(dname, " | Delta-delta overlay"),
      legend_txt = leg_txt,
      mode = "delta",
      col_aparece = col_ap, col_desaparece = col_des,
      is_delta = TRUE, invert_col = FALSE
    )
    return(invisible(NULL))
  })
  
  stopCluster(cl)
  log_msg("Delta-delta overlays done for zoom: ", zoom_tag)
  }
  
  
  # ==========================================================
  # CAP-HEAT OVERLAYS (v1 y v2) con triángulos — paralelo
  # ==========================================================
  
  if (RUN_CAP_HEAT) {
  log_msg("Generating Cap-Heat overlays (parallel, ", N_CORES, " cores)...")
  
  cl <- makeCluster(N_CORES)
  clusterExport(cl, varlist = c(
    "capHeat_v1_mats", "capHeat_v2_mats",
    "delta_flow_mats",
    "X_RANGE", "Y_RANGE", "dir_capHeat",
    "PNG_W", "PNG_H", "PNG_DPI", "N_COLS", "OVERWRITE",
    "FLOW_PCH", "FLOW_CEX",
    "special_res1", "special_res2", "special_res3", "special_res4",
    "subset_nodes", "plot_overlay", "draw_special_res",
    "draw_flow_persys", "draw_flow_delta",
    "COL_CAP_APARECE", "COL_CAP_DESAPARECE",
    "COL_HEAT_APARECE", "COL_HEAT_DESAPARECE", "COL_MANTIENE",
    "align_to_dist", "align_delta_to_dist",
    "classify_system", "log_msg", "LOG_FILE"
  ), envir = environment())
  clusterEvalQ(cl, {
    suppressPackageStartupMessages({ library(viridisLite); library(fields) })
  })
  
  # Construir lista de trabajos: una entrada por mutante x mat_type x version
  ch_jobs <- list()
  for (mut in names(sys_index_dist)) {
    for (mat_type in c("mean", "normvar")) {
      for (ver in c("v1", "v2")) {
        key <- paste0(mut, "_capHeat_", ver, "_", mat_type)
        mat_list <- if (ver == "v1") capHeat_v1_mats else capHeat_v2_mats
        if (!is.null(mat_list[[key]])) {
          ch_jobs[[length(ch_jobs) + 1]] <- list(
            key = key, mut = mut, mat_type = mat_type, ver = ver
          )
        }
      }
    }
  }
  
  parLapply(cl, ch_jobs, function(job) {
    key      <- job$key
    mut      <- job$mut
    mat_type <- job$mat_type
    ver      <- job$ver
    
    mat_list    <- if (ver == "v1") capHeat_v1_mats else capHeat_v2_mats
    dist_ch_full <- mat_list[[key]]
    
    nodes_x <- subset_nodes(colnames(dist_ch_full), X_RANGE)
    nodes_y <- subset_nodes(rownames(dist_ch_full), Y_RANGE)
    dist_z  <- dist_ch_full[nodes_y, nodes_x, drop = FALSE]
    
    # flujo Cap-Apo para triángulo inferior
    cap_flow_name  <- paste0(mut, "_flow_delta_Cap")
    heat_flow_name <- paste0(mut, "_flow_delta_Heat")
    
    flow_z  <- NULL
    flow_z2 <- NULL
    
    if (!is.null(delta_flow_mats[[cap_flow_name]])) {
      flow_full <- align_delta_to_dist(delta_flow_mats[[cap_flow_name]], dist_ch_full)
      flow_z    <- flow_full[nodes_y, nodes_x, drop = FALSE]
    }
    if (!is.null(delta_flow_mats[[heat_flow_name]])) {
      flow_full2 <- align_delta_to_dist(delta_flow_mats[[heat_flow_name]], dist_ch_full)
      flow_z2    <- flow_full2[nodes_y, nodes_x, drop = FALSE]
    }
    
    leg_txt <- if (mat_type == "mean") "Cap-Heat Distance (Å)" else "Cap-Heat Var / Mean^2"
    
    outfile <- file.path(dir_capHeat, paste0(key, "_overlay.png"))
    if (file.exists(outfile) && !OVERWRITE) return(invisible(NULL))
    
    plot_overlay(
      dist_z = dist_z, flow_z =NULL,# flow_z,
      outfile = outfile,
      main_txt = paste0(key, " | Cap-Heat overlay (", ver, ")"),
      legend_txt = leg_txt,
      mode = "cap_heat",
      flow_z2 = flow_z2,
      col_aparece    = COL_CAP_APARECE,
      col_desaparece = COL_CAP_DESAPARECE,
      col_aparece2    = COL_HEAT_APARECE,
      col_desaparece2 = COL_HEAT_DESAPARECE,
      is_delta = TRUE, invert_col = FALSE
    )
    return(invisible(NULL))
  })
  
  stopCluster(cl)
  log_msg("Cap-Heat overlays done for zoom: ", zoom_tag)
  
  log_msg("ZOOM ", zoom_idx, " COMPLETE: ", zoom_tag)
}
}
# ===============================================================
# END
# ===============================================================

log_msg("============================================================")
log_msg("SCRIPT FINISHED SUCCESSFULLY")
log_msg("End time: ", format(Sys.time(), "%Y-%m-%d %H:%M:%S"))
log_msg("============================================================")
