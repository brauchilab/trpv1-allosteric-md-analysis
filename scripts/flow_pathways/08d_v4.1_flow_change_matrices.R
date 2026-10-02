# ===============================================================
# 08d_v4.1_matrix_visualization_intercambio_de_chains_7LPC.R
# VISUALIZACIÓN DE MATRICES DE FLUJO DESDE CSV DEL 08b
# Author: DenyCB
# ---------------------------------------------------------------
# OBJETIVO:
# Lee las matrices CSV generadas por 08b_v4.1 (una por sistema-
# familia-modo-flow_type), las suma por sistema para obtener una
# matriz de flujo total por edge, y genera gráficos de:
#   - Matrices por sistema (per system)
#   - Matrices delta (condición vs apo)
#   - Matrices delta-delta (mutante vs WT, en espacio de deltas)
#
# ===============================================================

suppressPackageStartupMessages({
  library(data.table)
})

# ===============================================================
# 🔧 SETTINGS
# ===============================================================

# --- NMA (uncomment to use) ---
#BASE_OUTPUT_DIR <- "C:/DinamicasMoleculares/TRPV1_pipeline/output/08_Max_Flow/08b_v3.1_flow_paths_NMA_score_not_sqrt_min9/"
SUFFIX_SYSTEM <- "_MD"

# --- MD (comment out NMA block above and uncomment below) ---
BASE_OUTPUT_DIR <- "C:/DinamicasMoleculares/TRPV1_pipeline/output/08_Max_Flow/08b_v3.1_flow_paths_MD_min7"
###BASE_OUTPUT_DIR <- "C:/DinamicasMoleculares/TRPV1_pipeline/output/08_Max_Flow/08g_v2.1_markov_NMA"
###OJO CAMBIOS EN 262 y 869 a 876
#BASE_OUTPUT_DIR <- "C:/DinamicasMoleculares/TRPV1_pipeline/output/08_Max_Flow/08h_v1.4_paths_MD"


RUN_PER_SYSTEM  <- TRUE
RUN_COMPARATIVE <- TRUE
RUN_DELTADELTA  <- TRUE
# Full canvas expansion settings
# Set TRUE to expand all matrices to the full residue map canvas
# (including nodes with zero flow), FALSE to keep only nodes with flow
EXPAND_TO_FULL_CANVAS <- FALSE

# Residue map files for full canvas expansion (select via SUFFIX_SYSTEM)
CANVAS_MAP_FILE_MD  <- "C:/DinamicasMoleculares/TRPV1_pipeline/output/residue_map_MD.csv"
CANVAS_MAP_FILE_NMA <- "C:/DinamicasMoleculares/TRPV1_pipeline/config_files/trimmed_residue_map.csv"

# Gap residues to exclude explicitly from full canvas (all chains)
CANVAS_GAP_RESNUM <- 603:624

# ---------------------------------------------------------------
# zoom ranges
# ---------------------------------------------------------------
# Lista de pares (X_RANGE, Y_RANGE) a procesar en una sola ejecución.
# Cada par genera su propia carpeta de output con su zoom_tag.
# Ejemplos de rango:
#   list(c("A277","D752"), c("A277","D752"))   # full matrix
#   list(c("A440","A555"), c("D575","D670"))   # zoom LBD vs pore
# ---------------------------------------------------------------

ZOOM_RANGES <- list(
#  list(X = c("A400","A590"), Y = c("D560","D690")),
  list(X = c("A500","A580"), Y = c("A500","A580"))#,
#  list(X = c("A300","A752"), Y = c("A300","A752")),
#  list(X = c("A300","A752"), Y = c("B300","C752")),
#  list(X = c("B300","C752"), Y = c("B300","C752")),
#  list(X = c("D300","D752"), Y = c("D300","D752")),
#  list(X = c("A300","A752"), Y = c("D300","D752")),
#  list(X = c("B300","C752"), Y = c("D300","D752")),
#  list(X = c("A300","D752"), Y = c("A300","D752"))
)

# ---------------------------------------------------------------
# tamaño imagen
# ---------------------------------------------------------------

PNG_WIDTH  <- 9000
PNG_HEIGHT <- 9000
PNG_RES    <- 300

# ---------------------------------------------------------------
# barra de color
# ---------------------------------------------------------------
# ancho = distancia horizontal
# xleft/xright mueven barra
# ---------------------------------------------------------------

LEGEND_XLEFT_OFFSET  <- 2
LEGEND_XRIGHT_OFFSET <- 3

# ===============================================================
# 🔔 LOGGING
# ===============================================================

stamp <- function(){
  format(Sys.time(), "[%Y-%m-%d %H:%M:%S]")
}

msg <- function(...){
  cat(stamp(), ..., "\n")
}

# ===============================================================
# 📁 INPUTS
# ===============================================================

# Matrices del 08b están en la subcarpeta "matrices/"
DIR_MATRICES <- file.path(BASE_OUTPUT_DIR, "matrices/")#/Method2/")#"7LP9_W426A/")#/edge_importance")#

# ===============================================================
# 📁 OUTPUT BASE (subcarpetas por zoom se crean dentro del loop)
# ===============================================================

OUT_BASE <- file.path(BASE_OUTPUT_DIR, paste0("08d_viz_v4.1_2", SUFFIX_SYSTEM,"FLIP"))#(BASE_OUTPUT_DIR, paste0("08d_viz_v4.1", SUFFIX_SYSTEM))
dir.create(OUT_BASE, recursive = TRUE, showWarnings = FALSE)

# ===============================================================
# 🔎 HELPERS
# ===============================================================

# --- Leer matriz CSV del 08b ---
# Las matrices del 08b se guardan con write_csv (sin rownames).
# Como son simétricas, los rownames se reconstruyen desde los colnames:
# el header ya contiene los labels correctos (ej. A297, A308, ...),
# y como colnames == rownames en una matriz simétrica, se asignan
# directamente sin descartar ninguna columna de datos.
read_matrix_csv <- function(f){
  x <- read.csv(f, row.names = NULL, check.names = FALSE)
  # Detectar si la primera columna es de rownames (header vacío)
  # CSVs del 08b: sin rownames explícitos, colnames == rownames
  # CSVs del 08d: primera columna sin nombre contiene los rownames
  if(colnames(x)[1] == ""){
    m <- as.matrix(x[, -1])
    rownames(m) <- x[, 1]
  } else {
    m <- as.matrix(x)
    rownames(m) <- colnames(m)
  }
  storage.mode(m) <- "numeric"
  m
}

# --- Extraer nombre de sistema desde nombre de archivo del 08b ---
# Formato: matrix_<system>_<family>_<mode>_<flowtype>.csv
# Ejemplo: matrix_7LP9_WT_combined_undirected_filter.csv
# Sistema = primeras dos partes tras "matrix_": "7LP9_WT"
get_system_from_matrix_filename <- function(nm){
  # Eliminar prefijo "matrix_"
  nm2 <- sub("^matrix_", "", nm)
  # Las primeras dos partes separadas por "_" son pdb+mutante
  parts <- strsplit(nm2, "_")[[1]]
  paste(parts[1], parts[2], sep = "_")
}

range_nodes <- function(nodes, lim){
  
  parse_node <- function(x){
    chain <- substr(x,1,1)
    num   <- as.numeric(substring(x,2))
    list(chain=chain,num=num)
  }
  
  p1 <- parse_node(lim[1])
  p2 <- parse_node(lim[2])
  
  chain_nodes <- nodes[
    substr(nodes,1,1) %in% c(p1$chain,p2$chain)
  ]
  
  nums <- as.numeric(substring(chain_nodes,2))
  
  # inicio = primer nodo >= solicitado
  i1 <- which(
    substr(chain_nodes,1,1)==p1$chain &
      nums >= p1$num
  )[1]
  
  # fin = último nodo <= solicitado
  i2 <- tail(which(
    substr(chain_nodes,1,1)==p2$chain &
      nums <= p2$num
  ),1)
  
  if(length(i1)==0 || length(i2)==0){
    stop("No se pudo resolver rango: ",
         paste(lim, collapse=" - "))
  }
  
  a <- match(chain_nodes[i1], nodes)
  b <- match(chain_nodes[i2], nodes)
  
  nodes[seq(a,b)]
}

# ===============================================================
# residuos especiales (idénticos al original)
# ===============================================================

special_res1 <- c("A426","A697","A441","B426","B697","B441","C426","C697","C441","D426","D697","D441","A563","B563","C563","D563")

special_res2 <- c("A679","A643","B679","B643","C679","C643","D679","D643")

special_res3 <- c("A574","A510")#"D591","A547","A445","B547","A591")

special_res4 <- c("D300", "B300", "C300")

# ===============================================================
# 🎨 DRAW PER SYSTEM (idéntico al 05d + zoom)
# ===============================================================

draw_per_system <- function(z, outfile, main_txt){
  
  cell_px <- 10   # tamaño de cada celda en píxeles
  
  img_w <- min(18000,max(4000, ncol(z) * cell_px +2200))
  img_h <- min(18000,max(4000, nrow(z) * cell_px +2200))
  
  png(outfile, width=img_w, height=img_h, res=PNG_RES)
  
  par(mar = c(14,14,6,6))
  
  cols_use <- colorRampPalette(
    c("white", "red", "green", "blue")
  )(100)
  
  image(
    1:ncol(z),
    1:nrow(z),
    t(z[nrow(z):1, ]),
    col = cols_use,
    ###    zlim = c(0, GLOBAL_MAX_PER),
    zlim = c(0, GLOBAL_MAX_PER),
    axes = FALSE,
    xlab = "",
    ylab = "",
    main = main_txt,
    useRaster = TRUE
  )
  
  axis(1, at=1:ncol(z), labels=colnames(z), las=2, cex.axis=0.75)
  axis(2, at=1:nrow(z), labels=rev(rownames(z)), las=2, cex.axis=0.75)
  axis(3, at=1:ncol(z), labels=colnames(z), las=2, cex.axis=0.75)
  axis(4, at=1:nrow(z), labels=rev(rownames(z)), las=2, cex.axis=.75)
  
  chain_cols <- c(
    A="red",
    B="blue",
    C="green3",
    D="orange"
  )
  
  chains_x <- substr(colnames(z),1,1)
  chains_y <- substr(rownames(z),1,1)
  
  rect(
    xleft=(1:ncol(z))-0.5,
    ybottom=nrow(z)+0.1,
    xright=(1:ncol(z))+0.5,
    ytop=nrow(z)+1.1,
    col=chain_cols[chains_x],
    border=NA,
    xpd=TRUE
  )
  
  rect(
    xleft=-0,
    ybottom=(1:nrow(z))-0.5,
    xright=-1,
    ytop=(1:nrow(z))+0.5,
    col=rev(chain_cols[chains_y]),
    border=NA,
    xpd=TRUE
  )
  
  
  for(v in list(
    list(special_res1,"orange"),
    list(special_res2,"magenta"),
    list(special_res3,"darkgray"),
    list(special_res4,"black")
  )){
    
    # posiciones reales en eje X
    idx_x <- match(v[[1]], colnames(z))
    idx_x <- idx_x[!is.na(idx_x)]
    
    # posiciones reales en eje Y
    idx_y <- match(v[[1]], rownames(z))
    idx_y <- idx_y[!is.na(idx_y)]
    
    # líneas verticales
    if(length(idx_x) > 0){
      abline(
        v   = idx_x,
        col = v[[2]],
        lwd = 3,
        lty = 2
      )
    }
    
    # líneas horizontales
    if(length(idx_y) > 0){
      abline(
        h   = nrow(z) - idx_y + 1,
        col = v[[2]],
        lwd = 3,
        lty = 2
      )
    }
  }
  
  par(xpd=TRUE)
  
  yseq <- seq(0,1,length.out=100)
  
  yy <- seq(1, nrow(z), length.out=200)
  
  for(i in 1:199){
    rect(
      xleft  = ncol(z)+2,
      ybottom= yy[i],
      xright = ncol(z)+3,
      ytop   = yy[i+1],
      col    = cols_use[i],
      border = cols_use[i]
    )
  }
  
  box()
  dev.off()
}

# ===============================================================
# 🎨 DRAW COMPARATIVE (idéntico al 05d + zoom)
# ===============================================================

draw_comp <- function(z, outfile, main_txt){
  
  cell_px <- 10   # tamaño de cada celda en píxeles
  
  img_w <- min(18000, max(4000, ncol(z) * cell_px +2200))
  img_h <- min(18000, max(4000, nrow(z) * cell_px +2200))
  
  png(outfile, width=img_w, height=img_h, res=PNG_RES)
  
  par(mar = c(14,14,6,6))
  
  max_abs <- GLOBAL_MAX_COMP
  
  cols_use <- c(
    "magenta",
    colorRampPalette(c("red","gray85"))(46),
    rep("gray85",6),
    colorRampPalette(c("gray85","blue"))(46),
    "cyan"
  )
  
  z_plot <- z
  z_plot[z == -999] <- -max_abs * 1.05
  z_plot[z ==  999] <-  max_abs * 1.05
  
  image(
    1:ncol(z),
    1:nrow(z),
    t(z_plot[nrow(z):1, ]),
    col = cols_use,
    zlim = c(-max_abs*1.05, max_abs*1.05),
    axes = FALSE,
    xlab = "",
    ylab = "",
    main = "",
    useRaster = FALSE
  )
  title(main = main_txt, line = 5)
  
  axis(1, at=1:ncol(z), labels=colnames(z), las=2, cex.axis=1)
  axis(2, at=1:nrow(z), labels=rev(rownames(z)), las=2, cex.axis=1)
  axis(3, at=1:ncol(z), labels=colnames(z), las=2, cex.axis=1)
  axis(4, at=1:nrow(z), labels=rev(rownames(z)), las=2, cex.axis=1)
  
  chain_cols <- c(
    A="red",
    B="blue",
    C="green3",
    D="orange"
  )
  
  chains_x <- substr(colnames(z),1,1)
  chains_y <- substr(rownames(z),1,1)
  
  # rect(
  #   xleft=(1:ncol(z))-0.5,
  #   ybottom=nrow(z)+0.1,
  #   xright=(1:ncol(z))+0.5,
  #   ytop=nrow(z)+1.1,
  #   col=chain_cols[chains_x],
  #   border=NA,
  #   xpd=TRUE
  # )
  # 
  # rect(
  #   xleft=-0,
  #   ybottom=(1:nrow(z))-0.5,
  #   xright=-1,
  #   ytop=(1:nrow(z))+0.5,
  #   col=rev(chain_cols[chains_y]),
  #   border=NA,
  #   xpd=TRUE
  # )
  
  
  
  for(v in list(
    list(special_res1,"orange"),
    list(special_res2,"magenta"),
    list(special_res3,"darkgray"),
    list(special_res4,"black")
  )){
    
    # posiciones reales en eje X
    idx_x <- match(v[[1]], colnames(z))
    idx_x <- idx_x[!is.na(idx_x)]
    
    # posiciones reales en eje Y
    idx_y <- match(v[[1]], rownames(z))
    idx_y <- idx_y[!is.na(idx_y)]
    
    # líneas verticales
    if(length(idx_x) > 0){
      abline(
        v   = idx_x,
        col = v[[2]],
        lwd = 3,
        lty = 2
      )
    }
    
    # líneas horizontales
    if(length(idx_y) > 0){
      abline(
        h   = nrow(z) - idx_y + 1,
        col = v[[2]],
        lwd = 3,
        lty = 2
      )
    }
  }
  
 # par(xpd=TRUE)
  
#  yseq <- seq(0,1,length.out=100)
  
#  yy <- seq(1, nrow(z), length.out=200)
  
#  for(i in 1:199){
#    rect(
#      xleft  = ncol(z)+2,
#      ybottom= yy[i],
#      xright = ncol(z)+3,
#      ytop   = yy[i+1],
#      col    = cols_use[i],
#      border = cols_use[i]
#    )
#  }
  
  box()
  dev.off()
}

# ===============================================================
# 🔧 DELTA HELPERS
# ===============================================================

# --- Calcula delta entre dos matrices sumadas por sistema ---
# a = matriz sistema de referencia (apo o WT)
# b = matriz sistema de comparación (cap, heat o mutante)
# Las matrices son de flujo absoluto (no normalizadas).
# Lógica de sentinels:
#   - edge presente solo en b (aparece):  delta = +999
#   - edge presente solo en a (desaparece): delta = -999
#   - edge presente en ambos: delta = (b - a) normalizado a [-1, 1]
#   - edge ausente en ambos: NA
# "Presente" = valor > 0 en la matriz sumada.
make_delta <- function(a, b){
  
  all_nodes <- sort(unique(c(rownames(a), rownames(b))))
  
  expand_mat <- function(m, nodes){
    m_exp <- matrix(0, length(nodes), length(nodes),
                    dimnames = list(nodes, nodes))
    common <- intersect(nodes, rownames(m))
    m_exp[common, common] <- m[common, common]
    m_exp
  }
  
  a_exp <- expand_mat(a, all_nodes)
  b_exp <- expand_mat(b, all_nodes)
  
  # Ceros exactos → NA (ausencia genuina, igual que el original)
  a_exp[a_exp == 0] <- NA
  b_exp[b_exp == 0] <- NA
  
  # Umbral: 0.05% del máximo entre ambas matrices combinadas.
  # Edges por debajo son ruido y se tratan como ausentes.
  max_combined <- max(c(a_exp, b_exp), na.rm = TRUE)
  thresh <- 0.0005 * max_combined
  
  na_both <- (is.na(a_exp) | a_exp < thresh) &
    (is.na(b_exp) | b_exp < thresh)
  na_ref  <- (is.na(a_exp) | a_exp < thresh) &
    !is.na(b_exp) & b_exp >= thresh
  na_cmp  <- !is.na(a_exp) & a_exp >= thresh &
    (is.na(b_exp) | b_exp < thresh)
  
  a_exp[is.na(a_exp) | a_exp < thresh] <- 0
  b_exp[is.na(b_exp) | b_exp < thresh] <- 0
  
  out <- b_exp - a_exp
  
  out[na_ref] <-  999
  out[na_cmp] <- -999
  out[na_both] <- NA
  diag(out) <- NA
  out
}

# --- Calcula delta-delta desde dos matrices delta ---
# Idéntico a make_dd del 05d pero operando sobre deltas normalizados.
# wt_delta  = delta del WT (referencia)
# mut_delta = delta del mutante
# ΔΔ = mut_delta - wt_delta, con lógica de sentinels.
make_dd <- function(wt_delta, mut_delta){
  
  # Alinear dimensiones
  all_nodes <- sort(unique(c(rownames(wt_delta), rownames(mut_delta))))
  
  expand_mat <- function(m, nodes){
    m_exp <- matrix(NA, length(nodes), length(nodes),
                    dimnames = list(nodes, nodes))
    common <- intersect(nodes, rownames(m))
    m_exp[common, common] <- m[common, common]
    m_exp
  }
  
  a <- expand_mat(wt_delta,  all_nodes)  # WT
  b <- expand_mat(mut_delta, all_nodes)  # mutante
  
  # Clasificar cada celda en: aparece (999), desaparece (-999), no cambia (0), ausente (NA)
  cat_a <- ifelse(is.na(a),   NA,
                  ifelse(a ==  999,  1,    # aparece en WT
                         ifelse(a == -999, -1,    # desaparece en WT
                                0)))  # no cambia en WT
  
  cat_b <- ifelse(is.na(b),   NA,
                  ifelse(b ==  999,  1,    # aparece en mut
                         ifelse(b == -999, -1,    # desaparece en mut
                                0)))  # no cambia en mut
  
  # Matriz de output: mismas dimensiones
  out <- matrix(NA, length(all_nodes), length(all_nodes),
                dimnames = list(all_nodes, all_nodes))
  
  # Aplicar tabla lógica:
  # mut\wt   Aparece(1)  NoCambia(0)  Desaparece(-1)
  # Aparece(1)    0           2            999
  # NoCambia(0)  -2           0              1
  # Desaparece(-1) -999      -1              0
  
  both_present <- !is.na(cat_a) & !is.na(cat_b)
  
  out[both_present & cat_b ==  1 & cat_a ==  1] <-    0   # Gris
  out[both_present & cat_b ==  1 & cat_a ==  0] <-    2   # Blue
  out[both_present & cat_b ==  1 & cat_a == -1] <-  999   # Green
  out[both_present & cat_b ==  0 & cat_a ==  1] <-   -2   # Red
  out[both_present & cat_b ==  0 & cat_a ==  0] <-    0   # Gris
  out[both_present & cat_b ==  0 & cat_a == -1] <-    1   # Blue
  out[both_present & cat_b == -1 & cat_a ==  1] <- -999   # Cyan
  out[both_present & cat_b == -1 & cat_a ==  0] <-   -1   # Red
  out[both_present & cat_b == -1 & cat_a == -1] <-    0   # Gris
  
  diag(out) <- NA
  out
}

# ===============================================================
# 🚀 START
# ===============================================================

msg("SCRIPT START — 08d_v4.1")
msg("BASE_OUTPUT_DIR:", BASE_OUTPUT_DIR)
msg("SUFFIX_SYSTEM:", SUFFIX_SYSTEM)
msg("DIR_MATRICES:", DIR_MATRICES)

# ===============================================================
# LEER Y SUMAR MATRICES POR SISTEMA
# ===============================================================
# Para cada sistema (ej. 7LP9_WT), se suman todas las matrices
# disponibles (todas las familias, modos y flow_types) en una
# sola matriz de flujo total. Esta matriz sumada es la que se
# usa para graficar per system y para calcular deltas.
# ===============================================================

msg("Reading and summing matrices by system...")

ff_all <- list.files(DIR_MATRICES, pattern = "\\.csv$", full.names = TRUE)

if(length(ff_all) == 0){
  stop("No CSV files found in: ", DIR_MATRICES)
}

msg(paste("  Matrix files found:", length(ff_all)))

# Agrupar por sistema
system_names_all <- sapply(
  tools::file_path_sans_ext(basename(ff_all)),
  get_system_from_matrix_filename
)

unique_systems <- unique(system_names_all)
msg(paste("  Unique systems:", length(unique_systems)))

# Sumar todas las matrices del mismo sistema
system_mats <- list()

for(sys in unique_systems){
  
  msg(paste("  Summing matrices for system:", sys))
  
  idx <- which(system_names_all == sys)
  
  mat_sum <- NULL
  
  for(i in idx){
    
    m <- tryCatch(
      read_matrix_csv(ff_all[i]),
      error = function(e){
        msg(paste("  WARNING: could not read", basename(ff_all[i]),
                  "-", e$message))
        NULL
      }
    )
    
    if(is.null(m)) next
    
    if(is.null(mat_sum)){
      mat_sum <- m
    } else {
      # Alinear y sumar (en caso de que haya diferencias de nodos entre matrices)
      all_nodes <- sort(unique(c(rownames(mat_sum), rownames(m))))
      expand_to <- function(x, nodes){
        x_exp <- matrix(0, length(nodes), length(nodes),
                        dimnames = list(nodes, nodes))
        cn <- intersect(nodes, rownames(x))
        x_exp[cn, cn] <- x[cn, cn]
        x_exp
      }
      mat_sum <- expand_to(mat_sum, all_nodes) + expand_to(m, all_nodes)
    }
  }
  
  if(!is.null(mat_sum)){
    
    # Fix 7LPC chain label swap: C↔D relabeling to match ABCD convention
    # In 7LPC the PDB order is CABD, causing C and D labels to be swapped
    # relative to the functional ABCD assignment used in all other systems.
 #   if(grepl("^7LPC_", sys)){
 #     rn <- rownames(mat_sum)
#      rn_fixed <- ifelse(substr(rn,1,1) == "C", paste0("D", substring(rn,2)),
#                         ifelse(substr(rn,1,1) == "D", paste0("C", substring(rn,2)),
#                                rn))
#      rownames(mat_sum) <- rn_fixed
#      colnames(mat_sum) <- rn_fixed
#      msg(paste("    -> 7LPC C<->D chain relabeling applied"))
#    }
    
    system_mats[[sys]] <- mat_sum
    msg(paste("    -> matrix size:", nrow(mat_sum), "x", ncol(mat_sum),
              "| max value:", round(max(mat_sum, na.rm = TRUE), 2)))
  }
}

msg(paste("Systems loaded:", length(system_mats)))

# ===============================================================
# MATRIZ CANVAS GLOBAL
# ===============================================================
# Construye el espacio de nodos global: la unión de todos los nodos
# presentes en al menos una de las matrices sumadas por sistema.
# Todas las matrices graficadas se expandirán a este espacio común
# (rellenando con 0 donde no hay dato), garantizando que todos los
# gráficos tengan exactamente el mismo tamaño y que un pixel que
# aparece en cualquier sistema aparezca en todos.
# ===============================================================

msg("Building global canvas node space...")

global_nodes <- sort(unique(unlist(lapply(system_mats, rownames))))
msg(paste("  Global canvas nodes:", length(global_nodes)))

expand_to_canvas <- function(m, nodes){
  m_exp <- matrix(0, length(nodes), length(nodes),
                  dimnames = list(nodes, nodes))
  cn <- intersect(nodes, rownames(m))
  m_exp[cn, cn] <- m[cn, cn]
  m_exp
}

# Expandir todas las matrices sumadas al espacio global
# Expandir todas las matrices sumadas al espacio global
system_mats_canvas <- lapply(system_mats, expand_to_canvas, nodes = global_nodes)

# ===============================================================
# FULL CANVAS EXPANSION (optional, controlled by EXPAND_TO_FULL_CANVAS)
# ===============================================================
# If TRUE, expands all matrices to the full residue space defined
# by the residue map (MD or NMA), excluding gap residues explicitly.
# This allows overlaying flow matrices with distogram outputs.
# ===============================================================

if (EXPAND_TO_FULL_CANVAS) {
  
  canvas_map_file <- if (SUFFIX_SYSTEM == "_NMA") CANVAS_MAP_FILE_NMA else CANVAS_MAP_FILE_MD
  
  msg("EXPAND_TO_FULL_CANVAS = TRUE")
  msg("Using canvas map: ", canvas_map_file)
  
  rmap_canvas     <- read.csv(canvas_map_file, stringsAsFactors = FALSE)
  canvas_labels   <- paste0(rmap_canvas$chain, rmap_canvas$residue)
  
  # Exclude gap residues from all chains
  gap_labels_canvas <- as.vector(outer(c("A","B","C","D"), CANVAS_GAP_RESNUM, paste0))
  canvas_labels     <- canvas_labels[!(canvas_labels %in% gap_labels_canvas)]
  
  # Sort canonically by chain then residue number
  canvas_chain  <- substr(canvas_labels, 1, 1)
  canvas_resnum <- as.numeric(substring(canvas_labels, 2))
  canvas_labels <- canvas_labels[order(canvas_chain, canvas_resnum)]
  
  msg(paste("  Full canvas nodes (after gap removal):", length(canvas_labels)))
  
  expand_to_full <- function(m, nodes) {
    m_exp <- matrix(0, length(nodes), length(nodes),
                    dimnames = list(nodes, nodes))
    cn <- intersect(nodes, rownames(m))
    m_exp[cn, cn] <- m[cn, cn]
    m_exp
  }
  
  system_mats_canvas <- lapply(system_mats_canvas, expand_to_full, nodes = canvas_labels)
  
  msg("  system_mats_canvas expanded to full canvas.")
}

# ===============================================================
# CALCULAR DELTAS EN MEMORIA (una sola vez, fuera del loop de rangos)
# ===============================================================

all_per_mats  <- system_mats
all_comp_mats <- list()

if(RUN_COMPARATIVE || RUN_DELTADELTA){
  
  msg("Computing delta matrices...")
  
  # Mapeo PDB → condición
  pdb_condition <- list(
    "7LP9" = "Apo",
    "7LPB" = "Cap",
    "7LPC" = "Heat"
  )
  
  get_pdb  <- function(sys) strsplit(sys, "_")[[1]][1]
  get_mut  <- function(sys) strsplit(sys, "_")[[1]][2]
  
  mutants    <- unique(sapply(names(system_mats), get_mut))
  conditions <- c("Cap", "Heat")
  
  for(mut in mutants){
    for(cond in conditions){
      
      apo_sys  <- names(system_mats)[
        sapply(names(system_mats), function(s)
          get_mut(s) == mut && pdb_condition[[get_pdb(s)]] == "Apo")
      ]
      cond_sys <- names(system_mats)[
        sapply(names(system_mats), function(s)
          get_mut(s) == mut && pdb_condition[[get_pdb(s)]] == cond)
      ]
      
      if(length(apo_sys) == 0 || length(cond_sys) == 0){
        msg(paste("  Missing system for:", mut, cond, "-> skipped"))
        next
      }
      
      # Use canvas-expanded matrices for delta calculation
      apo_mat  <- system_mats_canvas[[apo_sys[1]]]
      cond_mat <- system_mats_canvas[[cond_sys[1]]]
      
      delta_name <- paste0(mut, "_Apo_vs_", cond, SUFFIX_SYSTEM)
      msg(paste("  Computing delta:", delta_name))
      
      delta <- make_delta(apo_mat, cond_mat)
      all_comp_mats[[delta_name]] <- delta
    }
  }
  
  msg(paste("  Delta matrices computed:", length(all_comp_mats)))
}

# ===============================================================
# LOOP DE RANGOS DE ZOOM
# ===============================================================
# Itera sobre todos los pares definidos en ZOOM_RANGES.
# Para cada rango genera su propia carpeta de output y
# corre los bloques per system, comparative y delta-delta.
# ===============================================================

all_vals_per <- unlist(lapply(system_mats_canvas, function(m) m[m > 0]))
GLOBAL_MIN_PER <- quantile(all_vals_per, probs = 0.05, na.rm = TRUE)
GLOBAL_MAX_PER <- quantile(all_vals_per, probs = 0.95, na.rm = TRUE)

GLOBAL_MAX_COMP <- quantile(
  abs(unlist(lapply(all_comp_mats, function(m){
    m[m != 999 & m != -999 & !is.na(m)]
  }))),
  probs = 0.995,
  na.rm = TRUE,
  names = FALSE,
  type = 7
)

for(zoom_idx in seq_along(ZOOM_RANGES)){
  
  X_RANGE <- ZOOM_RANGES[[zoom_idx]]$X
  Y_RANGE <- ZOOM_RANGES[[zoom_idx]]$Y
  
  # ==========================================================
  # PRECOMPUTE ZOOM NODE SUBSETS (same for all matrices)
  # ==========================================================
  
  canvas_nodes_x <- colnames(system_mats_canvas[[1]])
  canvas_nodes_y <- rownames(system_mats_canvas[[1]])
  
  xnodes_zoom <- range_nodes(canvas_nodes_x, X_RANGE)
  ynodes_zoom <- range_nodes(canvas_nodes_y, Y_RANGE)
  
  zoom_tag <- paste0(
    X_RANGE[1], "_", X_RANGE[2],
    "vs",
    Y_RANGE[1], "_", Y_RANGE[2]
  )
  
  msg(paste("=== ZOOM", zoom_idx, "of", length(ZOOM_RANGES),
            "| range:", zoom_tag, "==="))
  
  OUT_PER <- file.path(
    OUT_BASE,
    paste0("per_system", SUFFIX_SYSTEM, "_", zoom_tag)
  )
  
  OUT_COMP <- file.path(
    OUT_BASE,
    paste0("comparative", SUFFIX_SYSTEM, "_", zoom_tag)
  )
  
  dir.create(OUT_PER,  recursive = TRUE, showWarnings = FALSE)
  dir.create(OUT_COMP, recursive = TRUE, showWarnings = FALSE)
  
  # ===============================================================
  # PER SYSTEM
  # ===============================================================
  
  if(RUN_PER_SYSTEM){
    
    msg("START PER SYSTEM")
    
    #    GLOBAL_MAX_PER <- max(
    #      unlist(lapply(system_mats_canvas, function(m) m[m > 0])),
    #      na.rm = TRUE
    #    )
    
    #    all_vals_per <- unlist(lapply(system_mats_canvas, function(m) m[m > 0]))
    #    GLOBAL_MIN_PER <- quantile(all_vals_per, probs = 0.05, na.rm = TRUE)
    #    GLOBAL_MAX_PER <- quantile(all_vals_per, probs = 0.95, na.rm = TRUE)
    
    msg(paste("  GLOBAL_MAX_PER:", round(GLOBAL_MAX_PER, 4)))
    
    for(sys in names(all_per_mats)){
      
      msg("Processing:", sys)
      
      z <- system_mats_canvas[[sys]]
      
      zz <- z[ynodes_zoom, xnodes_zoom, drop=FALSE]
      
      outpng <- file.path(
        OUT_PER,
        paste0(
          sys, SUFFIX_SYSTEM,
          "_",
          zoom_tag,
          ".png"
        )
      )
      
      draw_per_system(
        zz,
        outpng,
        paste("Flow Matrix -", sys, SUFFIX_SYSTEM)
      )
      
      # Guardar también la matriz sumada como CSV para referencia
      ###      write.csv(
      ###        zz,
      ###        file.path(
      ###          OUT_PER,
      ###          paste0("matrix_", sys, SUFFIX_SYSTEM, "_", zoom_tag, ".csv")
      ###        ),
      ###        row.names = TRUE
      ###      )
    }
    
    msg("END PER SYSTEM")
  }
  
  # ===============================================================
  # COMPARATIVE (DELTA)
  # ===============================================================
  
  if(RUN_COMPARATIVE){
    
    msg("START COMPARATIVE (DELTA)")
    
    if(length(all_comp_mats) == 0){
      msg("  WARNING: no delta matrices computed — check system names")
    } else {
      
      
      
      msg(paste("  GLOBAL_MAX_COMP:", round(GLOBAL_MAX_COMP, 4)))
      
      for(nm in names(all_comp_mats)){
        
        msg("Processing:", nm)
        
        z <- all_comp_mats[[nm]]
        
        zz <- z[ynodes_zoom, xnodes_zoom, drop=FALSE]
        
        outpng <- file.path(
          OUT_COMP,
          paste0(nm, "_", zoom_tag, ".png")
        )
        
        draw_comp(
          zz,
          outpng,
          paste("Delta Flow Matrix -", nm)
        )
        
        # Guardar CSV del delta
        write.csv(
          zz,
          file.path(
            OUT_COMP,
            paste0("matrix_", nm, "_", zoom_tag, ".csv")
          ),
          row.names = TRUE
        )
      }
    }
    
    msg("END COMPARATIVE")
  }
  
  # ===============================================================
  # DELTA-DELTA MATRICES
  # ===============================================================
  # Usa matrices delta ya calculadas en memoria (all_comp_mats)
  # y calcula:
  #
  # ΔΔ = Δ(mutante) - Δ(WT)
  #
  # Manteniendo lógica de aparición/desaparición idéntica al 05d.
  # ===============================================================
  
  if(RUN_DELTADELTA){
    
    msg("START DELTA-DELTA")
    
    OUT_DD <- file.path(
      OUT_BASE,
      paste0("delta_delta", SUFFIX_SYSTEM, "_", zoom_tag)
    )
    
    dir.create(OUT_DD, recursive=TRUE, showWarnings=FALSE)
    
    # ------------------------------------------------------------
    # comparaciones definidas
    # ΔΔ = Δ(mutante, misma condición) - Δ(WT, misma condición)
    # ------------------------------------------------------------
    
    dd_list <- list(
      
      c(
        paste0("WT_Apo_vs_Cap", SUFFIX_SYSTEM),
        paste0("W426A_Apo_vs_Cap", SUFFIX_SYSTEM),
        paste0("DD_W426A_vs_WT_ApoCap", SUFFIX_SYSTEM)
      ),
      
      c(
        paste0("WT_Apo_vs_Heat", SUFFIX_SYSTEM),
        paste0("W426A_Apo_vs_Heat", SUFFIX_SYSTEM),
        paste0("DD_W426A_vs_WT_ApoHeat", SUFFIX_SYSTEM)
      ),
      
      c(
        paste0("WT_Apo_vs_Cap", SUFFIX_SYSTEM),
        paste0("W697A_Apo_vs_Cap", SUFFIX_SYSTEM),
        paste0("DD_W697A_vs_WT_ApoCap", SUFFIX_SYSTEM)
      ),
      
      c(
        paste0("WT_Apo_vs_Heat", SUFFIX_SYSTEM),
        paste0("W697A_Apo_vs_Heat", SUFFIX_SYSTEM),
        paste0("DD_W697A_vs_WT_ApoHeat", SUFFIX_SYSTEM)
      ),
      
      c(
        paste0("WT_Apo_vs_Cap", SUFFIX_SYSTEM),
        paste0("Y441A_Apo_vs_Cap", SUFFIX_SYSTEM),
        paste0("DD_Y441A_vs_WT_ApoCap", SUFFIX_SYSTEM)
      ),
      
      c(
        paste0("WT_Apo_vs_Heat", SUFFIX_SYSTEM),
        paste0("Y441A_Apo_vs_Heat", SUFFIX_SYSTEM),
        paste0("DD_Y441A_vs_WT_ApoHeat", SUFFIX_SYSTEM)
      )
    )
    
    # ------------------------------------------------------------
    # loop
    # ------------------------------------------------------------
    all_dd_mats <- list()
    
    for(v in dd_list){
      
      wt_name  <- v[1]
      mut_name <- v[2]
      out_name <- v[3]
      
      msg("Processing:", out_name)
      
      wt  <- all_comp_mats[[wt_name]]
      mut <- all_comp_mats[[mut_name]]
      
      if(is.null(wt) || is.null(mut)){
        msg("Missing input -> skipped")
        next
      }
      
      z <- make_dd(wt, mut)
      all_dd_mats[[out_name]] <- z
      
      # zoom
      zz <- z[ynodes_zoom, xnodes_zoom, drop=FALSE]
      
      # guardar csv
      write.csv(
        zz,
        file.path(
          OUT_DD,
          paste0(
            "matrix_",
            out_name,
            "_",
            zoom_tag,
            ".csv"
          )
        ),
        row.names=TRUE
      )
      
    }
    
    if(length(all_dd_mats) > 0){
      
      GLOBAL_MAX_DD <- quantile(
        abs(unlist(lapply(all_comp_mats, function(m){
          m[m != 999 & m != -999 & !is.na(m)]
        }))),
        probs = 0.995,
        na.rm = TRUE,
        names = FALSE,
        type = 7
      )
      
      GLOBAL_MAX_COMP <- GLOBAL_MAX_DD
      
      msg(paste("  GLOBAL_MAX_DD (triangles):", round(GLOBAL_MAX_DD, 4)))
      
      for(v in dd_list){
        
        wt_name  <- v[1]
        mut_name <- v[2]
        out_name <- v[3]
        
        msg("Processing triangle plot:", out_name)
        
        wt_delta  <- all_comp_mats[[wt_name]]
        mut_delta <- all_comp_mats[[mut_name]]
        
        if(is.null(wt_delta) || is.null(mut_delta)){
          msg("Missing input -> skipped")
          next
        }
        
        # zoom
        wt_zz  <- wt_delta[ynodes_zoom,  xnodes_zoom, drop=FALSE]
        mut_zz <- mut_delta[ynodes_zoom, xnodes_zoom, drop=FALSE]
        
        # construir matriz híbrida: triángulo inferior = mutante, superior = WT
        # diagonal definida por: col > (ny + 1 - img_row) → superior (WT)
        #                         col <= (ny + 1 - img_row) → inferior (mut)
        # En términos de índices de matriz (no imagen): diagonal \
        # superior derecha (i < j) → WT; inferior izquierda (i >= j) → mut
        ny  <- nrow(wt_zz)
        ncx <- ncol(wt_zz)
        
        hybrid <- matrix(NA, ny, ncx,
                         dimnames = list(rownames(wt_zz), colnames(wt_zz)))
        
        for(i in seq_len(ny)){
          for(j in seq_len(ncx)){
            if(i < j){
              hybrid[i, j] <- wt_zz[i, j]   # triángulo superior → WT
            } else {
              hybrid[i, j] <- mut_zz[i, j]  # triángulo inferior → mut
            }
          }
        }
        
        draw_comp(
          hybrid,
          file.path(
            OUT_DD,
            paste0("matrix_", out_name, "_", zoom_tag, ".png")
          ),
          paste("WT (upper) vs", sub("_Apo_vs_.*","",mut_name),
                "(lower) —", out_name)
        )
      }
    }
    
    msg("END DELTA-DELTA")
  }
  
  msg(paste("=== ZOOM", zoom_idx, "DONE ==="))
  
} # end ZOOM_RANGES loop

#msg("SCRIPT END — 08d_v4.1")
#    msg("END DELTA-DELTA")
#  }
  
#  msg(paste("=== ZOOM", zoom_idx, "DONE ==="))
  
#} # end ZOOM_RANGES loop

msg("SCRIPT END — 08d_v4.1")
