###############################################################
# 06a_v2.1 HELIX AXES FROM MD TRAYECTORY
# Autor: DenyCB
# Fecha: 2026-03-10
###############################################################

library(bio3d)
library(parallel)
library(future.apply)
library(ggplot2)
library(dplyr)

###############################################################
# 0. Configuración inicial
###############################################################

project_root <- "C:/DinamicasMoleculares/TRPV1_pipeline"

output_dir <- file.path(project_root,"output","angles_helics")
angles_dir <- file.path(output_dir,"angles")
if(!dir.exists(angles_dir)) dir.create(angles_dir, recursive = TRUE)

num_cores <- 6

helix_preS1_set <- TRUE
helix_s1_set    <- TRUE
helix_s2_set    <- TRUE
helix_s3_set    <- TRUE
helix_s4_set    <- TRUE
helix_linker_set<- TRUE
helix_s5_set    <- TRUE
helix_s6_set    <- TRUE
helix_trph_set  <- TRUE
compute_pore_axis_set <- TRUE   # activar/desactivar cálculo del eje del poro


helices_to_calc <- c()
if(helix_preS1_set) helices_to_calc <- c(helices_to_calc, "preS1")
if(helix_s1_set)    helices_to_calc <- c(helices_to_calc, "S1")
if(helix_s2_set)    helices_to_calc <- c(helices_to_calc, "S2")
if(helix_s3_set)    helices_to_calc <- c(helices_to_calc, "S3")
if(helix_s4_set)    helices_to_calc <- c(helices_to_calc, "S4")
if(helix_linker_set)helices_to_calc <- c(helices_to_calc, "Linker")
if(helix_s5_set)    helices_to_calc <- c(helices_to_calc, "S5")
if(helix_s6_set)    helices_to_calc <- c(helices_to_calc, "S6")
if(helix_trph_set)  helices_to_calc <- c(helices_to_calc, "TRPh")

cat("Hélices a calcular:", paste(helices_to_calc, collapse=", "), "\n")

###############################################################
# Leer config.txt
###############################################################

config_file <- file.path(project_root, "config.txt")
config_lines <- readLines(config_file)

helix_residues <- list()

for(h in helices_to_calc){
  
  pattern <- paste0("helix_", h, "=")
  line <- grep(pattern, config_lines, value = TRUE)
  
  range_text <- sub(pattern, "", line)
  bounds <- as.numeric(unlist(strsplit(range_text, "-")))
  
  helix_residues[[h]] <- seq(bounds[1], bounds[2])
}

cat("Residuos por hélice cargados correctamente\n")

###############################################################
# BLOQUE OPCIONAL: CALCULO DEL EJE DEL PORO
###############################################################



if(compute_pore_axis_set){
  
  cat("Cálculo del eje del poro ACTIVADO\n")
  
  ###############################################################
  # Leer rango del gate desde config.txt
  ###############################################################
  
  gate_line <- grep("gate_residue_range", config_lines, value = TRUE)
  
  if(length(gate_line) == 0){
    stop("No se encontró 'gate_residue_range' en config.txt")
  }
  
  # eliminar todo antes del número
  gate_range_text <- gsub(".*=", "", gate_line)
  
  # eliminar espacios
  gate_range_text <- trimws(gate_range_text)
  
  bounds <- as.numeric(unlist(strsplit(gate_range_text, "-")))
  
  if(any(is.na(bounds))){
    stop("No se pudo interpretar el rango del gate en config.txt")
  }
  
  gate_residues <- seq(bounds[1], bounds[2])
  
  # ========================
  # Definir filtro de selectividad
  # ========================
  
  selectivity_filter_residues <- 643:648
  
  # ========================
  # Residuos totales del poro
  # ========================
  
  pore_residues <- sort(unique(
    c(gate_residues, selectivity_filter_residues)
  ))
  
  cat("Residuos del poro:", paste(pore_residues, collapse=" "), "\n")
  
}

###############################################################
# Leer tablas del proyecto
###############################################################

systems <- read.csv(
  file.path(output_dir,"project_systems_table.csv"),
  stringsAsFactors = FALSE
)

###############################################################
# FUNCION PCA HELICE POR FRAME
###############################################################

compute_axis_frame <- function(xyz_frame){
  
  xyz <- matrix(xyz_frame, ncol=3, byrow=TRUE)
  
  pca <- prcomp(xyz)
  
  axis <- pca$rotation[,1]
  
  return(axis)
}

###############################################################
# FUNCION PCA EJE DEL PORO POR FRAME
###############################################################

compute_pore_axis <- function(pdb, traj, pore_residues){
  
  cat("  Calculando eje del poro...\n")
  
  sel <- atom.select(
    pdb,
    elety="CA",
    resno=pore_residues
  )
  
  if(length(sel$xyz)==0){
    
    stop("No se encontraron CA para el poro")
    
  }
  
  xyz_traj <- traj[, sel$xyz]
  
  compute_axis <- function(frame_xyz){
    
    xyz <- matrix(frame_xyz, ncol=3, byrow=TRUE)
    
    pca <- prcomp(xyz)
    
    axis <- pca$rotation[,1]
    
    return(axis)
    
  }
  
  pore_axes <- t(apply(xyz_traj,1,compute_axis))
  
  return(pore_axes)
  
}


###############################################################
# FUNCION PROCESAR HELICE POR CADENA
###############################################################

process_chain_helix <- function(chain, helix_name, pdb, traj){
  
  residues <- helix_residues[[helix_name]]
  
  sel <- atom.select(
    pdb,
    chain=chain,
    elety="CA",
    resno=residues
  )
  
  if(length(sel$xyz)==0){
    
    warning(paste("No se encontraron CA para",chain,helix_name))
    return(NULL)
    
  }
  
  xyz_traj <- traj[, sel$xyz]
  
  axes <- t(apply(xyz_traj,1,compute_axis_frame))
  
  return(axes)
}

###############################################################
# PROCESAR SISTEMA COMPLETO
###############################################################

process_system <- function(system_row){
  
  system_id <- paste(
    system_row$structure,
    system_row$mutant,
    paste0("rep", system_row$replica),
    sep="_"
  )
  
  cat("\n[",Sys.time(),"] Procesando sistema:",system_id,"\n")
  
  pdb_file <- file.path(project_root, system_row$pdb)
  dcd_file <- file.path(project_root, system_row$dcd)
  
  pdb <- read.pdb(pdb_file)
  traj <- read.dcd(dcd_file)
  
  chains <- unique(pdb$atom$chain)
  
  results <- list()
  
  for(chain in chains){
    
    for(h in helices_to_calc){
      
      key <- paste(chain,h,sep="_")
      
      cat("  Calculando",key,"\n")
      
      axes <- process_chain_helix(chain,h,pdb,traj)
      
      results[[key]] <- axes
      
    }
    
  }
  
  ###############################################################
  # Calcular eje del poro (opcional)
  ###############################################################
  
  if(compute_pore_axis_set){
    
    cat("  Calculando eje del poro\n")
    
    pore_axes <- compute_pore_axis(pdb, traj, pore_residues)
    
    results[["PORE_AXIS"]] <- pore_axes
    
  }
  
  out_file <- file.path(
    angles_dir,
    paste0("axes_",system_id,".rds")
  )
  
  saveRDS(results,out_file)
  
  cat("[",Sys.time(),"] Sistema procesado:",system_id,"\n")
  
  return(results)
}
###############################################################
# PROCESAR TODOS LOS SISTEMAS
###############################################################

plan(multisession, workers=num_cores)

all_results <- future_lapply(
  1:nrow(systems),
  function(i) process_system(systems[i,]),
  future.seed = TRUE
)

cat("\n[",Sys.time(),"] Todos los sistemas procesados\n")