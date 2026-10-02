# ==============================================================
# 08m_v1.0_Nfield_to_pdb_bfactor.R
# Author: DenyCB
# ==============================================================
#
# PROPOSITO
# ---------
#   Proyectar el "campo de visitacion" derivado de la matriz
#   fundamental de estados absorbentes N = (I - Q)^-1 (producida
#   por 08g) sobre la columna B-factor de los PDB originales
#   (7LP9 / 7LPB / 7LPC), de modo que CADA residuo de la estructura
#   reciba un valor escalar, INCLUIDOS los residuos que nunca
#   aparecieron en ningun path/cluster reconstruido por 08h.
#
#   Este es justamente el objetivo: a diferencia de los analisis
#   por-path (que solo contienen residuos visitados por algun path),
#   el campo de visitacion de N asigna un valor a TODO residuo de la
#   red, permitiendo colorear la estructura completa por participacion
#   en la comunicacion alosterica.
#
# DEFINICION DEL CAMPO (identica a compute_node_importance de 08h)
# ----------------------------------------------------------------
#   Para cada residuo destino j:
#
#       campo(j) = sum_seed [ (w_seed / W) * N[seed, j] ]
#
#   donde:
#     - seed recorre el conjunto SEEDS que esta presente en las FILAS
#       (rownames) de N,
#     - w_seed es el peso del seed,
#     - W = sum de los pesos de los seeds presentes (normalizacion),
#     - N[seed, j] = numero esperado de visitas al residuo j partiendo
#       del residuo seed, antes de la absorcion en el target.
#
#   IMPORTANTE: N NO es simetrica. Las FILAS son el origen (i) y las
#   COLUMNAS el destino (j). El campo vive sobre las COLUMNAS j (el
#   residuo "visitado"). En los .rds, rownames == colnames y en el
#   mismo orden (verificado: identical(rownames, colnames) == TRUE),
#   por lo que el mapeo etiqueta -> residuo es directo (p.ej. "A412"
#   = cadena A, residuo 412), sin necesidad de residue_map_MD.csv.
#
# FUENTE DE DATOS
# ---------------
#   Se leen los archivos .rds (NO los .csv). El .csv consolidado de
#   08g presenta un encabezado de columnas corrupto (etiquetas escritas
#   como "0") y filas extra ajenas a la matriz NxN; el .rds conserva
#   los dimnames de filas y columnas intactos, por lo que es la unica
#   fuente fiable para mapear columna -> residuo.
#
#   Nombre de archivo esperado:
#       N_{dir}_{target}_{sistema}_mean.rds
#   ej.:
#       N_fwd_filter_7LP9_W426A_mean.rds
#   en la carpeta:
#       {INPUT_08G_DIR}/consolidated/
#
# SALIDAS (controladas por switches en la SECCION 1)
# ---------------------------------------------------
#   OUTPUT_DIR/per_combination/{sys}/{sys}_{dir}_{target}.pdb
#       Un PDB por combinacion (fwd_filter, fwd_gate, rev_filter,
#       rev_gate). Campo crudo de cada N individual.
#
#   OUTPUT_DIR/filter_gate_combined/{sys}/{sys}_{dir}_combined.pdb
#       filter + gate SUMADOS dentro de la misma direccion. Significado:
#       participacion del residuo en llevar senal hacia el dominio de
#       poro (filtro O compuerta) bajo esa direccion de senal. filter
#       (G643) y gate (I679) son ambos elementos del poro -> su suma es
#       biologicamente interpretable. NO se mezclan fwd con rev.
#
#   OUTPUT_DIR/consolidated/{sys}/{sys}_consolidated.pdb
#       Suma de TODAS las combinaciones disponibles del sistema (replica
#       la logica de E_acc en 08h). ADVERTENCIA: mezcla fwd y rev, que
#       son configuraciones seed/target distintas (rev es la direccion
#       de construccion target->seed). Se ofrece por consistencia con
#       08h, pero interpretar con cautela.
#
#   OUTPUT_DIR/delta/...  (si GENERATE_DELTA = TRUE)
#       Mapas de diferencia mutante - WT (campo normalizado), pintados
#       sobre el PDB compartido (la mutacion es in silico => mismas
#       coordenadas que el WT). Permite ver residuos "mas rojos/azules"
#       (gano/perdio participacion) igual que en las matrices de flujo.
#       Nota: el delta cruzado entre condiciones (cap-apo, heat-apo)
#       NO se automatiza porque involucra estructuras distintas (7LPB,
#       7LPC vs 7LP9); pintar ese delta exige elegir sobre que estructura
#       hacerlo y queda como configuracion manual documentada abajo.
#
#   OUTPUT_DIR/logs/
#       Log con timestamp (nombre de script, version y variables clave).
#
# CONTRATO DE ESCRITURA DEL PDB
# -----------------------------
#   - Se copia el PDB fuente y se reescribe SOLO las columnas 61-66
#     (B-factor) de los registros ATOM de proteina en las cadenas A-D.
#   - Residuo presente en el PDB pero ausente de N -> B-factor 0.00.
#   - Registros HETATM (capsaicina, lipidos, aguas) -> B-factor 0.00.
#   - Todo lo demas (coordenadas, ocupancia, cadenas distintas de A-D,
#     headers, conectividades) se conserva intacto.
#
# DEPENDENCIAS
# ------------
#   Solo R base + paquete 'parallel' (incluido en R base). Sin CRAN.
# ==============================================================


# ==============================================================
# SECCION 1: CONFIGURACION
#   Modifica esta seccion antes de ejecutar.
# ==============================================================

VERSION <- "08m_v1.0"

# --- Directorios de entrada ---
# Carpeta que contiene la subcarpeta consolidated/ con los .rds de N.
INPUT_08G_DIR <- "C:/DinamicasMoleculares/TRPV1_pipeline/output/08_Max_Flow/08g_v2.1_markov_MD"

# Carpeta con los PDB fuente y su mapeo prefijo -> nombre de archivo.
PDB_DIR <- "C:/DinamicasMoleculares/TRPV1_pipeline/input/pdb"
PDB_FILES <- c(
  "7LP9" = "7LP9",   # apo
  "7LPB" = "7LPB",   # capsaicina
  "7LPC" = "7LPC"    # calor (48C)
)

# --- Directorio de salida ---
OUTPUT_DIR <- "C:/DinamicasMoleculares/TRPV1_pipeline/output/08_Max_Flow/08m_v1.0_Nfield_pdb"

# --- Combinaciones a procesar ---
# Direcciones y targets de N a leer. Deben coincidir con los .rds de 08g.
RUN_FWD <- TRUE     # direccion forward (seed -> target)
RUN_REV <- TRUE     # direccion reverse (target -> seed)
TARGETS <- c("filter", "gate")   # nombres tal como aparecen en los .rds

# --- Switches de salida (que PDB generar) ---
GENERATE_PER_COMBINATION    <- TRUE   # 1 PDB por (dir,target): hasta 4 por sistema
GENERATE_FILTER_GATE_COMBINED <- TRUE # filter+gate sumados por direccion: hasta 2 por sistema
GENERATE_CONSOLIDATED       <- TRUE   # suma de TODAS las combinaciones: 1 por sistema
GENERATE_DELTA              <- TRUE   # mapas mutante - WT (mismo prefijo PDB)

# --- Identificacion del WT para los deltas automaticos ---
# Sufijo que identifica al sistema WT dentro de cada prefijo PDB.
# Ej.: "7LP9_WT" es el WT de "7LP9_W426A", "7LP9_W697A", etc.
WT_TAG <- "WT"

# (OPCIONAL/MANUAL) Delta cruzado entre condiciones (estructuras distintas).
# Cada entrada: c(sistema_A, sistema_B, prefijo_PDB_sobre_el_que_pintar).
# El campo escrito sera campo(A) - campo(B). Dejar vacio para no generarlos.
# Ej.: list(c("7LPB_WT","7LP9_WT","7LPB"))  # cap - apo, pintado sobre 7LPB
DELTA_MANUAL_PAIRS <- list()

# --- Normalizacion del B-factor ---
#   "global_per_label" : por cada etiqueta de salida (p.ej. "fwd_filter",
#                        "consolidated"), se escala 0..99.99 usando el
#                        maximo global de esa etiqueta a traves de todos
#                        los sistemas. => mismos tipos comparables y
#                        RESTABLES entre sistemas (necesario para delta).
#   "raw"              : valores crudos (se truncan a 999.99 con aviso si
#                        exceden el ancho de la columna B-factor).
NORMALIZATION_MODE <- "global_per_label"

# Tope de seguridad para el ancho de la columna B-factor (cols 61-66).
# %6.2f admite hasta 999.99 (positivos) y hasta -99.99 (negativos).
BFACTOR_ABS_CAP_POS <- 99.99    # escala de campos absolutos -> [0, 99.99]
BFACTOR_ABS_CAP_DELTA <- 99.99  # escala de deltas -> [-99.99, 99.99]

# --- Tratamiento de atomos en el PDB ---
PROTEIN_CHAINS  <- c("A", "B", "C", "D")  # cadenas que reciben el campo
SET_HETATM_ZERO <- TRUE                   # poner B-factor de HETATM en 0.00

# --- Procesamiento paralelo (Windows 10 compatible, PSOCK) ---
N_CORES <- 6


# ==============================================================
# SECCION 2: SEEDS
#   Conjunto identico al usado por 08h (compute_node_importance).
#   Vector nombrado: nombre = "cadena+residuo", valor = peso.
# ==============================================================

SEEDS <- c(
  "C289"=1.0, "B289"=1.0,
  "C296"=1.0, "B296"=1.0,
  "A297"=1.0, "C297"=1.0, "B297"=1.0, "D297"=1.0,
  "C298"=1.0, "A301"=1.0,
  "B337"=1.0, "D337"=1.0, "A337"=1.0,
  "A338"=1.0, "D338"=1.0, "B338"=1.0,
  "B341"=1.0, "C341"=1.0, "D341"=1.0,
  "D344"=1.0, "B344"=1.0,
  "C345"=1.0, "B345"=1.0, "A345"=1.0,
  "D346"=1.0, "B346"=1.0, "C346"=1.0, "A346"=1.0,
  "B347"=1.0, "A347"=1.0, "C347"=1.0, "D347"=1.0,
  "C349"=1.0, "A349"=1.0, "B349"=1.0, "D349"=1.0,
  "C350"=1.0, "A350"=1.0, "D350"=1.0, "B350"=1.0,
  "A351"=1.0, "A352"=1.0, "C352"=1.0, "D354"=1.0,
  "D378"=1.0, "D382"=1.0,
  "A387"=1.0, "C387"=1.0, "A388"=1.0,
  "D395"=1.0, "C395"=1.0, "C396"=1.0,
  "D399"=1.0, "C399"=1.0, "A399"=1.0, "B399"=1.0,
  "A400"=1.0, "B400"=1.0, "D400"=1.0, "C400"=1.0,
  "A408"=1.0, "C408"=1.0, "B408"=1.0, "D408"=1.0,
  "C409"=1.0, "B409"=1.0, "D409"=1.0, "A409"=1.0,
  "B410"=1.0, "C410"=1.0, "A410"=1.0, "D410"=1.0,
  "C411"=1.0, "A411"=1.0, "D411"=1.0, "B411"=1.0,
  "D412"=1.0, "B412"=1.0, "A412"=1.0, "C412"=1.0,
  "C413"=1.0, "A413"=1.0, "B413"=1.0, "D413"=1.0,
  "A414"=1.0, "D414"=1.0, "C414"=1.0, "B414"=1.0,
  "D418"=1.0, "A418"=1.0, "B418"=1.0, "C418"=1.0,
  "D419"=1.0, "A419"=1.0, "B419"=1.0, "C419"=1.0,
  "B421"=1.0, "C426"=1.0, "B430"=1.0, "A430"=1.0,
  "D745"=1.0, "A750"=1.0, "D746"=1.0, "A749"=1.0,
  "D747"=1.0, "B746"=1.0, "A746"=1.0, "C745"=1.0,
  "B731"=1.0, "D726"=1.0, "A747"=1.0, "C726"=1.0,
  "B729"=1.0, "C746"=1.0, "B730"=1.0, "D744"=1.0,
  "B747"=1.0, "A730"=1.0, "A731"=1.0, "D728"=1.0,
  "A551"=1.0, "B551"=1.0, "C551"=1.0, "D551"=1.0,
  "A550"=1.0, "B550"=1.0, "C550"=1.0, "D550"=1.0,
  "A570"=1.0, "B570"=1.0, "C570"=1.0, "D570"=1.0,
  "A571"=1.0, "B571"=1.0, "C571"=1.0, "D571"=1.0,
  "A510"=1.0, "B510"=1.0, "C510"=1.0, "D510"=1.0,
  "A511"=1.0, "B511"=1.0, "C511"=1.0, "D511"=1.0,
  "A512"=1.0, "B512"=1.0, "C512"=1.0, "D512"=1.0,
  "A516"=1.0, "B516"=1.0, "C516"=1.0, "D516"=1.0,
  "A462"=1.0, "C461"=1.0, "B461"=1.0, "D461"=1.0,
  "A457"=1.0, "B457"=1.0, "D462"=1.0, "A455"=1.0,
  "A461"=1.0, "B455"=1.0, "B460"=1.0, "C460"=1.0,
  "D460"=1.0, "A456"=1.0, "A458"=1.0, "C457"=1.0,
  "B536"=1.0, "A536"=1.0
)


# ==============================================================
# SECCION 3: LOGGING (consola + archivo, con timestamp)
# ==============================================================

# Rutas derivadas
CONSOLIDATED_DIR <- file.path(INPUT_08G_DIR, "consolidated")
LOGS_DIR         <- file.path(OUTPUT_DIR, "logs")
dir.create(OUTPUT_DIR, recursive = TRUE, showWarnings = FALSE)
dir.create(LOGS_DIR,   recursive = TRUE, showWarnings = FALSE)

.TIMESTAMP <- format(Sys.time(), "%Y%m%d_%H%M%S")
.LOG_FILE  <- file.path(LOGS_DIR, paste0(VERSION, "_log_", .TIMESTAMP, ".txt"))

# log_msg: imprime en consola y anexa al archivo de log, con timestamp.
log_msg <- function(msg) {
  stamp <- format(Sys.time(), "%Y-%m-%d %H:%M:%S")
  line  <- sprintf("[%s] %s | %s", stamp, VERSION, msg)
  cat(line, "\n", sep = "")
  try(cat(line, "\n", file = .LOG_FILE, sep = "", append = TRUE), silent = TRUE)
}


# ==============================================================
# SECCION 4: HELPERS — DESCUBRIMIENTO Y CALCULO DEL CAMPO
# ==============================================================

# discover_tasks: lista los .rds disponibles para las direcciones y
# targets solicitados, y extrae el nombre de sistema de cada uno.
# Devuelve un data.frame con (system, direction, target, path).
discover_tasks <- function(consolidated_dir, directions, targets) {
  files <- list.files(consolidated_dir, pattern = "^N_.*_mean\\.rds$",
                       full.names = TRUE)
  rows <- list()
  for (f in files) {
    base <- basename(f)
    # Patron: N_{dir}_{target}_{sys}_mean.rds  (sys puede tener "_")
    m <- regmatches(base, regexec(
      "^N_(fwd|rev)_(filter|gate)_(.+)_mean\\.rds$", base))[[1]]
    if (length(m) != 4) next
    dir_i <- m[2]; tgt_i <- m[3]; sys_i <- m[4]
    if (!(dir_i %in% directions)) next
    if (!(tgt_i %in% targets))    next
    rows[[length(rows) + 1]] <- data.frame(
      system = sys_i, direction = dir_i, target = tgt_i,
      path = f, stringsAsFactors = FALSE
    )
  }
  if (length(rows) == 0) return(NULL)
  do.call(rbind, rows)
}

# compute_base_field: lee una matriz N (.rds) y calcula el campo de
# visitacion por residuo, replicando exactamente compute_node_importance
# de 08h. Esta funcion se ejecuta en procesos hijos (parallel), por lo
# que es autocontenida y solo usa R base.
#
#   Entrada : task = lista con $system,$direction,$target,$path ; seeds.
#   Salida  : lista con el campo (vector nombrado por residuo) y metadatos.
compute_base_field <- function(task, seeds) {
  out <- list(system = task$system, direction = task$direction,
              target = task$target, field = NULL,
              n_valid_seeds = 0L, status = "ok", msg = "")
  N <- tryCatch(readRDS(task$path),
                error = function(e) { out$status <<- "error"
                                      out$msg <<- conditionMessage(e); NULL })
  if (is.null(N)) return(out)

  rn <- rownames(N); cn <- colnames(N)
  if (is.null(rn) || is.null(cn)) {
    out$status <- "error"; out$msg <- "matriz sin dimnames"; return(out)
  }

  # Seeds presentes en las FILAS (origen i) de N
  seed_names <- names(seeds)
  valid      <- intersect(seed_names, rn)
  out$n_valid_seeds <- length(valid)
  if (length(valid) == 0) {
    out$status <- "warn"; out$msg <- "ningun seed presente en rownames"
    out$field  <- setNames(rep(0.0, length(cn)), cn)
    return(out)
  }

  # Pesos normalizados (igual que total_w en 08h)
  w <- as.numeric(seeds[valid])
  w <- w / sum(w)

  # Submatriz de filas-seed; campo(j) = sum_seed w_seed * N[seed, j]
  # (w como fila 1xS) %*% (N_sub SxC) = 1xC  -> vector sobre columnas j
  N_sub <- N[valid, , drop = FALSE]
  field <- as.numeric(w %*% N_sub)
  names(field) <- cn

  out$field <- field
  out
}


# ==============================================================
# SECCION 5: HELPERS — COMBINACION Y NORMALIZACION DE CAMPOS
# ==============================================================

# sum_fields: suma una lista de vectores nombrados (por residuo),
# alineando por nombre (union); valores ausentes cuentan como 0.
sum_fields <- function(field_list) {
  field_list <- Filter(Negate(is.null), field_list)
  if (length(field_list) == 0) return(NULL)
  all_names <- Reduce(union, lapply(field_list, names))
  acc <- setNames(rep(0.0, length(all_names)), all_names)
  for (fv in field_list) {
    acc[names(fv)] <- acc[names(fv)] + fv
  }
  acc
}

# scale_to_cap: escala un vector a [0, cap] por su maximo (o por un
# maximo externo provisto). Si max <= 0, devuelve ceros.
scale_to_cap <- function(vec, cap, ext_max = NULL) {
  mx <- if (is.null(ext_max)) max(vec, na.rm = TRUE) else ext_max
  if (!is.finite(mx) || mx <= 0) return(vec * 0.0)
  vec * (cap / mx)
}

# scale_symmetric: escala un vector con signo a [-cap, cap] por su
# maximo valor absoluto (o por un maximo externo). Preserva el 0.
scale_symmetric <- function(vec, cap, ext_absmax = NULL) {
  amx <- if (is.null(ext_absmax)) max(abs(vec), na.rm = TRUE) else ext_absmax
  if (!is.finite(amx) || amx <= 0) return(vec * 0.0)
  vec * (cap / amx)
}


# ==============================================================
# SECCION 6: HELPERS — ESCRITURA DEL PDB CON B-FACTOR
# ==============================================================

# write_pdb_bfactor: copia el PDB fuente y reescribe SOLO las columnas
# 61-66 (B-factor) de los ATOM de proteina en PROTEIN_CHAINS, segun el
# campo provisto (nombrado por "cadena+residuo", p.ej. "A412"). HETATM
# se ponen en 0.00 si SET_HETATM_ZERO. Todo lo demas se conserva.
#
#   field_named : vector nombrado (residuo -> valor B ya escalado)
#   default_val : valor para residuos ausentes del campo (0.00)
write_pdb_bfactor <- function(src_pdb_path, out_pdb_path, field_named,
                              protein_chains, set_hetatm_zero,
                              default_val = 0.0) {
  if (!file.exists(src_pdb_path)) {
    log_msg(sprintf("  ERROR: PDB fuente no encontrado: %s", src_pdb_path))
    return(FALSE)
  }
  lines <- readLines(src_pdb_path, warn = FALSE)
  out   <- character(length(lines))

  fmt_b <- function(v) sprintf("%6.2f", v)  # ancho exacto cols 61-66

  for (i in seq_along(lines)) {
    line <- lines[i]
    rec  <- substr(line, 1, 6)

    is_atom   <- startsWith(line, "ATOM")
    is_hetatm <- startsWith(line, "HETATM")

    if (!is_atom && !is_hetatm) { out[i] <- line; next }

    # Asegurar largo minimo de 66 columnas para poder reescribir B-factor
    if (nchar(line) < 66) line <- formatC(line, width = 66, flag = "-")

    bval <- NA_real_

    if (is_atom) {
      chain  <- trimws(substr(line, 22, 22))
      resnum <- suppressWarnings(as.integer(trimws(substr(line, 23, 26))))
      if (chain %in% protein_chains && !is.na(resnum)) {
        key <- paste0(chain, resnum)
        if (!is.null(field_named) && key %in% names(field_named)) {
          bval <- as.numeric(field_named[[key]])
        } else {
          bval <- default_val   # residuo proteico A-D ausente de N -> 0
        }
      }
      # ATOM en cadenas distintas de A-D: se deja intacto (bval = NA)
    } else if (is_hetatm && set_hetatm_zero) {
      bval <- 0.0
    }

    if (is.na(bval)) {
      out[i] <- line                       # no se toca
    } else {
      out[i] <- paste0(substr(line, 1, 60), fmt_b(bval),
                       substr(line, 67, nchar(line)))
    }
  }

  dir.create(dirname(out_pdb_path), recursive = TRUE, showWarnings = FALSE)
  writeLines(out, out_pdb_path)
  TRUE
}


# ==============================================================
# SECCION 7: MAIN
# ==============================================================

main <- function() {

  t0 <- Sys.time()
  log_msg("============================================================")
  log_msg(sprintf("INICIO %s", VERSION))
  log_msg("============================================================")
  log_msg(sprintf("INPUT_08G_DIR      : %s", INPUT_08G_DIR))
  log_msg(sprintf("CONSOLIDATED_DIR   : %s", CONSOLIDATED_DIR))
  log_msg(sprintf("PDB_DIR            : %s", PDB_DIR))
  log_msg(sprintf("OUTPUT_DIR         : %s", OUTPUT_DIR))
  log_msg(sprintf("RUN_FWD/RUN_REV    : %s / %s", RUN_FWD, RUN_REV))
  log_msg(sprintf("TARGETS            : %s", paste(TARGETS, collapse = ", ")))
  log_msg(sprintf("PER_COMBINATION    : %s", GENERATE_PER_COMBINATION))
  log_msg(sprintf("FILTER_GATE_COMB   : %s", GENERATE_FILTER_GATE_COMBINED))
  log_msg(sprintf("CONSOLIDATED       : %s", GENERATE_CONSOLIDATED))
  log_msg(sprintf("DELTA (mut - WT)   : %s  (WT_TAG=%s)", GENERATE_DELTA, WT_TAG))
  log_msg(sprintf("NORMALIZATION_MODE : %s", NORMALIZATION_MODE))
  log_msg(sprintf("PROTEIN_CHAINS     : %s", paste(PROTEIN_CHAINS, collapse = "")))
  log_msg(sprintf("SET_HETATM_ZERO    : %s", SET_HETATM_ZERO))
  log_msg(sprintf("N_CORES            : %d", N_CORES))
  log_msg(sprintf("N seeds definidos  : %d", length(SEEDS)))

  # ---- Direcciones efectivas ----
  directions <- c()
  if (RUN_FWD) directions <- c(directions, "fwd")
  if (RUN_REV) directions <- c(directions, "rev")
  if (length(directions) == 0) { log_msg("ERROR: sin direcciones."); return(invisible()) }

  # ---- Descubrir tareas (system x direction x target con .rds) ----
  tasks_df <- discover_tasks(CONSOLIDATED_DIR, directions, TARGETS)
  if (is.null(tasks_df)) { log_msg("ERROR: no se hallaron .rds de N."); return(invisible()) }
  systems <- sort(unique(tasks_df$system))
  log_msg(sprintf("Sistemas detectados: %d -> %s",
                  length(systems), paste(systems, collapse = ", ")))
  log_msg(sprintf("Tareas (matrices N) a procesar: %d", nrow(tasks_df)))

  # ---- Calculo paralelo de los campos base (uno por matriz N) ----
  task_list <- lapply(seq_len(nrow(tasks_df)), function(k) as.list(tasks_df[k, ]))

  use_parallel <- (N_CORES > 1 && length(task_list) > 1)
  if (use_parallel) {
    log_msg(sprintf("Calculando campos base en paralelo (%d nucleos, PSOCK)...", N_CORES))
    cl <- parallel::makeCluster(N_CORES)   # PSOCK: compatible Windows 10
    on.exit(try(parallel::stopCluster(cl), silent = TRUE), add = TRUE)
    parallel::clusterExport(cl, varlist = c("compute_base_field", "SEEDS"),
                            envir = globalenv())
    base_results <- parallel::parLapply(cl, task_list,
                                        function(t) compute_base_field(t, SEEDS))
    parallel::stopCluster(cl)
  } else {
    log_msg("Calculando campos base en modo secuencial...")
    base_results <- lapply(task_list, function(t) compute_base_field(t, SEEDS))
  }

  # ---- Indexar campos base por sistema y reportar estado ----
  # base_fields[[system]][["{dir}_{target}"]] = vector nombrado
  base_fields <- list()
  for (r in base_results) {
    label <- paste0(r$direction, "_", r$target)
    if (r$status == "error") {
      log_msg(sprintf("  ERROR %s/%s: %s", r$system, label, r$msg)); next
    }
    if (r$status == "warn") {
      log_msg(sprintf("  AVISO %s/%s: %s", r$system, label, r$msg))
    }
    if (is.null(base_fields[[r$system]])) base_fields[[r$system]] <- list()
    base_fields[[r$system]][[label]] <- r$field
    log_msg(sprintf("  Campo OK: %s/%s  (%d seeds validos, %d residuos)",
                    r$system, label, r$n_valid_seeds, length(r$field)))
  }

  # ---- Construir items de salida (sistema, etiqueta, campo crudo, prefijo) ----
  # Cada item: list(system, label, field, prefix)
  out_items <- list()
  add_item <- function(system, label, field) {
    if (is.null(field)) return(invisible())
    prefix <- strsplit(system, "_", fixed = TRUE)[[1]][1]
    out_items[[length(out_items) + 1]] <<- list(
      system = system, label = label, field = field, prefix = prefix
    )
  }

  for (sys in names(base_fields)) {
    bf <- base_fields[[sys]]

    # (a) Por combinacion
    if (GENERATE_PER_COMBINATION) {
      for (label in names(bf)) add_item(sys, label, bf[[label]])
    }

    # (b) filter + gate combinados, por direccion (significado: poro)
    if (GENERATE_FILTER_GATE_COMBINED) {
      for (dir_i in directions) {
        f_lab <- paste0(dir_i, "_filter")
        g_lab <- paste0(dir_i, "_gate")
        if (!is.null(bf[[f_lab]]) || !is.null(bf[[g_lab]])) {
          comb <- sum_fields(list(bf[[f_lab]], bf[[g_lab]]))
          add_item(sys, paste0(dir_i, "_combined"), comb)
        }
      }
    }

    # (c) Consolidado: suma de TODAS las combinaciones del sistema
    if (GENERATE_CONSOLIDATED) {
      cons <- sum_fields(bf)
      add_item(sys, "consolidated", cons)
    }
  }

  if (length(out_items) == 0) { log_msg("ERROR: sin items de salida."); return(invisible()) }

  # ---- Normalizacion ----
  # En "global_per_label": un maximo por etiqueta, comun a todos los
  # sistemas -> mismos tipos comparables y restables (necesario p/delta).
  if (NORMALIZATION_MODE == "global_per_label") {
    label_max <- list()
    for (it in out_items) {
      mx <- max(it$field, na.rm = TRUE)
      if (is.null(label_max[[it$label]]) || mx > label_max[[it$label]])
        label_max[[it$label]] <- mx
    }
    for (lab in names(label_max))
      log_msg(sprintf("  Norma global [%s]: max=%.6g -> escala a [0,%.2f]",
                      lab, label_max[[lab]], BFACTOR_ABS_CAP_POS))
  }

  # ---- Escritura de PDB para cada item (campos absolutos) ----
  log_msg("Escribiendo PDB de campos absolutos...")
  # Guardamos los campos NORMALIZADOS por (system,label) para el delta
  norm_fields <- list()

  n_written <- 0
  for (it in out_items) {
    if (NORMALIZATION_MODE == "global_per_label") {
      fld <- scale_to_cap(it$field, BFACTOR_ABS_CAP_POS,
                          ext_max = label_max[[it$label]])
    } else { # "raw"
      fld <- it$field
      if (max(fld, na.rm = TRUE) > 999.99) {
        log_msg(sprintf("  AVISO %s/%s: campo crudo excede 999.99; se trunca.",
                        it$system, it$label))
        fld[fld > 999.99] <- 999.99
      }
    }
    # Registrar para delta
    if (is.null(norm_fields[[it$system]])) norm_fields[[it$system]] <- list()
    norm_fields[[it$system]][[it$label]] <- fld

    src <- file.path(PDB_DIR, paste0(PDB_FILES[[it$prefix]], ".pdb"))
    out_sub <- if (it$label == "consolidated") "consolidated"
               else if (grepl("_combined$", it$label)) "filter_gate_combined"
               else "per_combination"
    out_pdb <- file.path(OUTPUT_DIR, out_sub, it$system,
                         paste0(it$system, "_", it$label, ".pdb"))

    ok <- write_pdb_bfactor(src, out_pdb, fld, PROTEIN_CHAINS,
                            SET_HETATM_ZERO, default_val = 0.0)
    if (ok) { n_written <- n_written + 1
              log_msg(sprintf("  PDB OK: %s", basename(out_pdb))) }
  }
  log_msg(sprintf("PDB de campos absolutos escritos: %d", n_written))

  # ---- Delta mutante - WT (automatico, mismo prefijo PDB) ----
  if (GENERATE_DELTA) {
    log_msg("Generando mapas delta (mutante - WT)...")
    n_delta <- 0
    # Agrupar sistemas por prefijo y localizar el WT de cada prefijo
    prefixes <- unique(sapply(names(norm_fields),
                              function(s) strsplit(s, "_", fixed = TRUE)[[1]][1]))
    for (pref in prefixes) {
      sys_in_pref <- names(norm_fields)[
        sapply(names(norm_fields),
               function(s) strsplit(s, "_", fixed = TRUE)[[1]][1] == pref)]
      wt_sys <- paste0(pref, "_", WT_TAG)
      if (!(wt_sys %in% sys_in_pref)) {
        log_msg(sprintf("  AVISO: sin WT (%s) para prefijo %s; delta omitido.",
                        wt_sys, pref)); next
      }
      mut_systems <- setdiff(sys_in_pref, wt_sys)
      for (mut in mut_systems) {
        common_labels <- intersect(names(norm_fields[[mut]]),
                                   names(norm_fields[[wt_sys]]))
        for (lab in common_labels) {
          d <- sum_fields(list(norm_fields[[mut]][[lab]],
                               -1.0 * norm_fields[[wt_sys]][[lab]]))
          # Escala simetrica para coloreo divergente (rojo/azul)
          d <- scale_symmetric(d, BFACTOR_ABS_CAP_DELTA)
          src <- file.path(PDB_DIR, paste0(PDB_FILES[[pref]], ".pdb"))
          out_pdb <- file.path(OUTPUT_DIR, "delta",
                               paste0(mut, "_minus_", WT_TAG),
                               paste0(mut, "_minus_", wt_sys, "_", lab, ".pdb"))
          ok <- write_pdb_bfactor(src, out_pdb, d, PROTEIN_CHAINS,
                                  SET_HETATM_ZERO, default_val = 0.0)
          if (ok) { n_delta <- n_delta + 1
                    log_msg(sprintf("  DELTA OK: %s", basename(out_pdb))) }
        }
      }
    }
    log_msg(sprintf("Mapas delta (mut - WT) escritos: %d", n_delta))

    # ---- Delta cruzado manual (estructuras distintas) ----
    if (length(DELTA_MANUAL_PAIRS) > 0) {
      log_msg("Generando deltas manuales (entre condiciones)...")
      for (pr in DELTA_MANUAL_PAIRS) {
        sysA <- pr[1]; sysB <- pr[2]; paint_pref <- pr[3]
        if (is.null(norm_fields[[sysA]]) || is.null(norm_fields[[sysB]])) {
          log_msg(sprintf("  AVISO: par manual (%s,%s) sin campos; omitido.",
                          sysA, sysB)); next
        }
        common_labels <- intersect(names(norm_fields[[sysA]]),
                                   names(norm_fields[[sysB]]))
        for (lab in common_labels) {
          d <- sum_fields(list(norm_fields[[sysA]][[lab]],
                               -1.0 * norm_fields[[sysB]][[lab]]))
          d <- scale_symmetric(d, BFACTOR_ABS_CAP_DELTA)
          src <- file.path(PDB_DIR, paste0(PDB_FILES[[paint_pref]], ".pdb"))
          out_pdb <- file.path(OUTPUT_DIR, "delta", "manual",
                               paste0(sysA, "_minus_", sysB, "_", lab,
                                      "_on_", paint_pref, ".pdb"))
          ok <- write_pdb_bfactor(src, out_pdb, d, PROTEIN_CHAINS,
                                  SET_HETATM_ZERO, default_val = 0.0)
          if (ok) log_msg(sprintf("  DELTA MANUAL OK: %s", basename(out_pdb)))
        }
      }
    }
  }

  # ---- Resumen final ----
  t1 <- Sys.time()
  elapsed <- as.numeric(difftime(t1, t0, units = "secs"))
  log_msg("============================================================")
  log_msg(sprintf("FIN %s — %.1f s", VERSION, elapsed))
  log_msg(sprintf("Sistemas: %d | Items absolutos: %d",
                  length(base_fields), length(out_items)))
  log_msg("============================================================")
}


# ==============================================================
# PUNTO DE ENTRADA
# ==============================================================

if (sys.nframe() == 0) {
  main()
}
