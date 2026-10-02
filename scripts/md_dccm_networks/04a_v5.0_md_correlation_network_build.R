###############################################################
# TRPV1 STRUCTURAL DYNAMICS PIPELINE
# Script: 04a_v5.0_build_cna_network_MD.R
# Author: DenyCB
###############################################################

library(bio3d)
library(igraph)
library(parallel)

###############################################################
# 1 LOAD CONFIGURATION FILE
###############################################################

start_time <- Sys.time()
log_file <- paste0("04a_v5.0_log_", format(start_time,"%Y%m%d_%H%M%S"), ".txt")

config_path <- file.path("config.txt")

if(!file.exists(config_path)){
  stop("ERROR: config.txt not found")
}

read_config <- function(config_file){
  
  lines <- readLines(config_file)
  lines <- lines[!grepl("^#", lines)]
  lines <- lines[nchar(lines) > 0]
  
  config <- list()
  
  for(line in lines){
    parts <- strsplit(line,"=")[[1]]
    key <- trimws(parts[1])
    value <- trimws(parts[2])
    config[[key]] <- value
  }
  
  return(config)
}

config <- read_config(config_path)

###############################################################
# 2 DEFINE DIRECTORIES
###############################################################

project_root <- config$project_root
output_dir <- config$pipeline_output_dir

dccm_md_dir <- file.path(config$pipeline_output_dir, "dccm_md")

network_dir_base <- file.path(output_dir,"network_objects_md")

if(!dir.exists(network_dir_base)){
  dir.create(network_dir_base, recursive = TRUE)
}

###############################################################
# 🔧 NUEVO: DEFINE THRESHOLD SWEEP
###############################################################

thresholds <- seq(0.47, 0.49, by = 0.01)

###############################################################
# 3 LOAD SYSTEM TABLE
###############################################################

systems_file <- file.path(output_dir,"project_systems_table.csv")

if(!file.exists(systems_file)){
  stop("ERROR: project_systems_table.csv not found")
}

systems <- read.csv(systems_file)

###############################################################
# 4 LOAD RESIDUE MAPPING
###############################################################

mapping_file <- file.path(output_dir,"residue_structural_map.csv")

if(!file.exists(mapping_file)){
  stop("ERROR: residue_structural_map.csv not found")
}

residue_map <- read.csv(mapping_file)

###############################################################
# 5 READ NETWORK PARAMETERS
###############################################################

dccm_threshold <- as.numeric(config$dccm_threshold)
contact_cutoff <- as.numeric(config$contact_cutoff)

# Distance cutoff for Cα-Cα contacts (Angstroms).
# Edges are only retained if the two residues have Cα within this distance.
# This removes spurious long-range correlations that lack physical contact.
ca_distance_cutoff <- 12.0

###############################################################
# 6 BUILD NETWORKS
###############################################################

for(th in thresholds){
  
  cat("====================================\n")
  cat("Running threshold:", th, "\n")
  cat("====================================\n")
  
  th_label <- paste0("Th", format(th, nsmall=2))
  
  
  network_dir <- file.path(network_dir_base, th_label)
  
  if(!dir.exists(network_dir)){
    dir.create(network_dir, recursive = TRUE)
  }
  
  cl <- makeCluster(6)
  
  clusterEvalQ(cl, {
    library(bio3d)
    library(igraph)
  })
  
  clusterExport(cl, c(
    "systems",
    "output_dir",
    "project_root",
    "contact_cutoff",
    "ca_distance_cutoff",
    "network_dir",
    "dccm_md_dir",
    "th"
  ))
  
  parLapply(cl, 1:nrow(systems), function(i){
    
    system_id <- paste(
      systems$structure[i],
      systems$mutant[i],
      paste0("rep", systems$replica[i]),
      sep="_"
    )
    
    cat("Building network for:", system_id,"\n")
    
    corr_file <- file.path(
      dccm_md_dir,
      paste0(
        systems$structure[i],
        "_",
        systems$mutant[i],
        "_rep",
        systems$replica[i],
        "_dccm.rds"
      )
    )
    
    if(!file.exists(corr_file)){
      stop(paste("Missing correlation file:", corr_file))
    }
    
    corr_matrix <- readRDS(corr_file)
    
    ###############################################################
    # LOAD PDB USED FOR MD
    ###############################################################
    
    pdb_file <- file.path(
      project_root,
      "TRPV1_pipeline",
      systems$pdb[i]
    )
    
    if(!file.exists(pdb_file)){
      stop(paste("Missing PDB file:", pdb_file))
    }
    
    pdb <- read.pdb(pdb_file)
    
    ca <- which(pdb$atom$elety == "CA")
    
    ###############################################################
    # CONSISTENCY CHECK BETWEEN CORRELATION MATRIX AND STRUCTURE
    ###############################################################
    
    if(nrow(corr_matrix) != length(ca)){
      stop(paste("Mismatch between DCCM size and CA atoms for:", system_id))
    }

    ###############################################################
    # BUILD Cα DISTANCE MATRIX AND CONTACT MASK
    # Only residue pairs with Cα-Cα distance <= ca_distance_cutoff
    # are eligible for edges, regardless of their correlation value.
    # This removes long-range correlations without physical contact.
    # pdb$xyz is a flat vector (x1,y1,z1,x2,y2,z2,...) indexed by atom;
    # correct Cα xyz indices are 3*ca-2, 3*ca-1, 3*ca.
    ###############################################################

    ca_xyz_idx <- as.vector(rbind(ca*3-2, ca*3-1, ca*3))
    ca_xyz <- matrix(as.numeric(pdb$xyz[ca_xyz_idx]), ncol = 3, byrow = TRUE)

    dist_mat <- as.matrix(dist(ca_xyz))

    contact_mask <- dist_mat <= ca_distance_cutoff
    diag(contact_mask) <- FALSE

    ###############################################################
    
    adj <- abs(corr_matrix)
    
    adj[adj < th] <- 0

    # Apply contact mask: zero out edges between residues
    # whose Cα atoms are farther apart than ca_distance_cutoff
    adj[!contact_mask] <- 0

    diag(adj) <- 0
    
    g <- graph_from_adjacency_matrix(
      adj,
      mode = "undirected",
      weighted = TRUE
    )
    
    comm <- cluster_walktrap(g, steps=3)
    
    membership_vec <- membership(comm)
    
    net <- list(
      graph = g,
      communities = comm,
      membership = membership_vec
    )
    
    out_file <- file.path(
      network_dir,
      paste0("network_MD_",system_id,".rds")
    )
    
    saveRDS(net, out_file)
    
  })
  
  stopCluster(cl)
  
  cat("Threshold", th, "completed\n")
}

cat("All networks constructed successfully\n")

writeLines(c(
  "SCRIPT 04a_v5.0 NETWORK ANALYSIS MD",
  paste("Start:", format(start_time)),
  paste("End:", format(Sys.time()))
), con = file.path(network_dir_base, log_file))

###############################################################
# PARAMETER REGISTRY UPDATE
###############################################################

tryCatch({

  library(dplyr)

  registry_dir          <- file.path(config$pipeline_output_dir, "registry")
  registry_final_path   <- file.path(registry_dir, "parameters_final.csv")
  registry_history_path <- file.path(registry_dir, "parameters_history.csv")

  script_number  <- "04a_MD"
  script_version <- "5.0"
  script_name    <- "04a_v5.0_build_cna_network_MD"
  ts             <- format(Sys.time(), "%Y-%m-%d %H:%M:%S")

  new_entries <- data.frame(
    script_number  = script_number,
    script_version = script_version,
    script_name    = script_name,
    variable = c("dccm_md_dir",
                 "network_dir_base",
                 "thresholds",
                 "ca_distance_cutoff",
                 "contact_cutoff",
                 "systems_processed",
                 "walktrap_steps"),
    value    = c(dccm_md_dir,
                 network_dir_base,
                 paste(range(thresholds), collapse = " to "),
                 as.character(ca_distance_cutoff),
                 as.character(contact_cutoff),
                 as.character(nrow(systems)),
                 "3"),
    timestamp = ts,
    stringsAsFactors = FALSE
  )

  if(file.exists(registry_final_path)){
    existing <- read.csv(registry_final_path, stringsAsFactors = FALSE,
                         colClasses = "character")
    existing <- existing[existing$script_number != script_number, ]
    updated  <- dplyr::bind_rows(existing, new_entries)
  } else {
    updated <- new_entries
  }
  write.csv(updated, registry_final_path, row.names = FALSE)

  if(file.exists(registry_history_path)){
    history <- read.csv(registry_history_path, stringsAsFactors = FALSE,
                        colClasses = "character")
    history <- dplyr::bind_rows(history, new_entries)
  } else {
    history <- new_entries
  }
  write.csv(history, registry_history_path, row.names = FALSE)

  cat("Parameter registry updated\n")

}, error = function(e){
  cat("WARNING: registry update failed:", e$message, "\n")
})
