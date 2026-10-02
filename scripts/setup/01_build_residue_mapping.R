###############################################################
# TRPV1 STRUCTURAL DYNAMICS PIPELINE
#
# Script: 01_build_residue_mapping.R
#
# Purpose:
# Automatically generate a structural residue mapping table
# linking Bio3D trimmed indices with original PDB residues
# and structural annotations.
#
# Output:
# output/residue_structural_map.csv
#
###############################################################

cat("\n============================================\n")
cat("TRPV1 PIPELINE - BUILD RESIDUE MAPPING\n")
cat("============================================\n\n")

###############################################################
# 1 — LOAD LIBRARIES
###############################################################

suppressMessages(library(bio3d))
suppressMessages(library(dplyr))
suppressMessages(library(stringr))

###############################################################
# 2 — READ CONFIG FILE
###############################################################

read_config <- function(config_file) {
  
  lines <- readLines(config_file)
  
  lines <- lines[!grepl("^#", lines)]
  lines <- lines[nchar(lines) > 0]
  
  config <- list()
  
  for (line in lines) {
    
    parts <- strsplit(line, "=")[[1]]
    
    key <- trimws(parts[1])
    value <- trimws(parts[2])
    
    config[[key]] <- value
  }
  
  return(config)
}

config_path <- file.path(getwd(), "config.txt")

if (!file.exists(config_path)) {
  stop("config.txt not found.")
}

config <- read_config(config_path)

###############################################################
# 3 — LOCATE PROJECT DIRECTORIES
###############################################################

input_md_dir <- file.path(getwd(), "input", "md")
output_dir <- file.path(getwd(), "output")

if (!dir.exists(input_md_dir)) {
  stop("input/md directory not found.")
}

###############################################################
# 4 — LOAD MASTER SYSTEM TABLE
###############################################################

systems_table_path <- file.path(output_dir, "project_systems_table.csv")

if (!file.exists(systems_table_path)) {
  stop("project_systems_table.csv not found. Run script 00 first.")
}

systems <- read.csv(systems_table_path)

###############################################################
# 5 — SELECT A REFERENCE PDB
###############################################################

# All systems share the same residue numbering.
# Use the first PDB as structural reference.

# -------------------------------------------------------------
# Build correct PDB path assuming working directory is project root
# -------------------------------------------------------------

reference_pdb_path <- gsub("^\\.\\./", "", systems$pdb[1])

if (!file.exists(reference_pdb_path)) {
  
  stop(
    paste(
      "Reference PDB file not found:",
      reference_pdb_path
    )
  )
  
}

cat("Loading reference PDB:\n", reference_pdb_path, "\n")

pdb <- read.pdb(reference_pdb_path)

###############################################################
# 6 — EXTRACT C-ALPHA RESIDUES
###############################################################

ca <- atom.select(pdb, "calpha")

residue_numbers <- pdb$atom$resno[ca$atom]
chains <- pdb$atom$chain[ca$atom]

trimmed_index <- seq_along(residue_numbers)

residue_table <- data.frame(
  trimmed_index = trimmed_index,
  chain = chains,
  residue_number = residue_numbers
)

###############################################################
# 7 — PARSE HELIX DEFINITIONS FROM CONFIG
###############################################################

helix_entries <- names(config)[grepl("^helix_", names(config))]

residue_table$helix <- NA

for (helix_name in helix_entries) {
  
  helix_label <- gsub("helix_", "", helix_name)
  
  range <- config[[helix_name]]
  
  bounds <- as.numeric(strsplit(range, "-")[[1]])
  
  start <- bounds[1]
  end <- bounds[2]
  
  idx <- residue_table$residue_number >= start &
    residue_table$residue_number <= end
  
  residue_table$helix[idx] <- helix_label
}

###############################################################
# 8 — FLAG MUTATION RESIDUES
###############################################################

mutation_residues <- as.numeric(strsplit(
  config$mutation_residues, ","
)[[1]])

residue_table$is_mutation_site <-
  residue_table$residue_number %in% mutation_residues

###############################################################
# 9 — FLAG GATE RESIDUES
###############################################################

gate_range <- config$gate_residue_range

gate_bounds <- as.numeric(strsplit(gate_range, "-")[[1]])

residue_table$is_gate_residue <-
  residue_table$residue_number >= gate_bounds[1] &
  residue_table$residue_number <= gate_bounds[2]
###############################################################
# 9B — FLAG SELECTIVITY FILTER RESIDUES
###############################################################

sf_range <- config$selectivity_filter_residue_range
sf_bounds <- as.numeric(strsplit(sf_range, "-")[[1]])

residue_table$is_selectivity_filter <-
  residue_table$residue_number >= sf_bounds[1] &
  residue_table$residue_number <= sf_bounds[2]

###############################################################
# 10 — WRITE OUTPUT
###############################################################

output_file <- file.path(output_dir, "residue_structural_map.csv")

write.csv(residue_table, output_file, row.names = FALSE)

cat("\nResidue mapping table generated:\n")
cat(output_file, "\n\n")