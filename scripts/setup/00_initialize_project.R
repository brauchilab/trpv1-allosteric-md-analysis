# ==========================================================
# TRPV1 PIPELINE
# Script: 00_initialize_project.R
# Purpose:
#   Initialize project environment and detect MD systems
# ==========================================================

rm(list = ls())

cat("Initializing TRPV1 pipeline...\n")

# ----------------------------------------------------------
# 1 Load required libraries
# ----------------------------------------------------------

library(stringr)
library(dplyr)
library(readr)

# ----------------------------------------------------------
# 2 Read config file
# ----------------------------------------------------------

read_config <- function(config_file){
  
  lines <- readLines(config_file)
  
  lines <- lines[!grepl("^#", lines)]
  lines <- lines[nchar(lines) > 0]
  
  config <- list()
  
  for(line in lines){
    
    parts <- strsplit(line, "=")[[1]]
    
    key <- trimws(parts[1])
    value <- trimws(parts[2])
    
    config[[key]] <- value
  }
  
  return(config)
}

config <- read_config("config.txt")

cat("Config file loaded\n")

# ----------------------------------------------------------
# 3 Validate directories
# ----------------------------------------------------------

required_dirs <- c(
  config$pipeline_output_dir,
  config$log_dir,
  config$nma_analysis_dir
)

for(dir in required_dirs){
  
  if(!dir.exists(dir)){
    stop(paste("Directory does not exist:", dir))
  }
  
}

cat("Directory validation completed\n")

# ----------------------------------------------------------
# 4 Detect MD trajectory files
# ----------------------------------------------------------

md_dir <- "input/md"

dcd_files <- list.files(
  md_dir,
  pattern = "\\.dcd$",
  full.names = TRUE
)

pdb_files <- list.files(
  md_dir,
  pattern = "\\.pdb$",
  full.names = TRUE
)

if(length(dcd_files) == 0){
  stop("No DCD files detected")
}

cat(length(dcd_files), "trajectory files detected\n")

# ----------------------------------------------------------
# 5 Parse system names
# Expected format:
# structure_mutant_repX.dcd
# Example:
# 7LP9_W426A_rep1.dcd
# ----------------------------------------------------------

parse_md_filename <- function(file){
  
  name <- basename(file)
  name <- gsub(".dcd","",name)
  
  parts <- strsplit(name,"_")[[1]]
  
  structure <- parts[1]
  mutant <- parts[2]
  
  replica <- str_extract(parts[3],"[0-9]+")
  
  return(data.frame(
    structure = structure,
    mutant = mutant,
    replica = as.integer(replica),
    dcd = file
  ))
  
}

systems <- do.call(
  rbind,
  lapply(dcd_files, parse_md_filename)
)

# ----------------------------------------------------------
# 6 Match PDB files
# ----------------------------------------------------------

systems$pdb <- NA

for(i in 1:nrow(systems)){
  
  expected_name <- paste0(
    systems$structure[i], "_",
    systems$mutant[i], "_rep",
    systems$replica[i], ".pdb"
  )
  
  match <- pdb_files[
    grepl(expected_name, pdb_files)
  ]
  
  if(length(match) == 1){
    systems$pdb[i] <- match
  } else {
    warning(paste("PDB not found for", expected_name))
  }
  
}

# ----------------------------------------------------------
# 7 Sort table
# ----------------------------------------------------------

systems <- systems %>%
  arrange(structure, mutant, replica)

# ----------------------------------------------------------
# 7.5 Check for missing replicas
# ----------------------------------------------------------

cat("\nChecking for missing MD replicas...\n")

replica_check <- systems %>%
  group_by(structure, mutant) %>%
  summarise(
    replicas_present = paste(sort(replica), collapse = ","),
    min_rep = min(replica),
    max_rep = max(replica),
    n_rep = n(),
    .groups = "drop"
  )

for(i in 1:nrow(replica_check)){
  
  expected <- seq(
    replica_check$min_rep[i],
    replica_check$max_rep[i]
  )
  
  observed <- systems$replica[
    systems$structure == replica_check$structure[i] &
      systems$mutant == replica_check$mutant[i]
  ]
  
  missing <- setdiff(expected, observed)
  
  if(length(missing) > 0){
    
    warning(
      paste(
        "Missing replicas for",
        replica_check$structure[i],
        replica_check$mutant[i],
        ":",
        paste(missing, collapse = ",")
      )
    )
    
  }
  
}

cat("Replica check completed\n")

# ----------------------------------------------------------
# 8 Save master systems table
# ----------------------------------------------------------

output_file <- file.path(
  config$pipeline_output_dir,
  "project_systems_table.csv"
)

write_csv(systems, output_file)

cat("Systems table saved to:\n")
cat(output_file, "\n")

# ----------------------------------------------------------
# 9 Summary report
# ----------------------------------------------------------

cat("\nProject summary:\n")

cat("Structures detected:\n")
print(unique(systems$structure))

cat("\nMutants detected:\n")
print(unique(systems$mutant))

cat("\nReplicas detected:\n")
print(unique(systems$replica))

cat("\nTotal systems:", nrow(systems), "\n")

cat("\nInitialization complete\n")

# ----------------------------------------------------------
# 10 Initialize parameter registry files
# ----------------------------------------------------------

registry_dir <- file.path(config$pipeline_output_dir, "registry")
dir.create(registry_dir, recursive = TRUE, showWarnings = FALSE)

registry_final_path   <- file.path(registry_dir, "parameters_final.csv")
registry_history_path <- file.path(registry_dir, "parameters_history.csv")

registry_cols <- c("script_number", "script_version", "script_name",
                   "variable", "value", "timestamp")

# Create empty files with headers if they don't exist yet
if (!file.exists(registry_final_path)) {
  write.csv(
    data.frame(matrix(ncol = length(registry_cols), nrow = 0,
                      dimnames = list(NULL, registry_cols))),
    registry_final_path, row.names = FALSE
  )
  cat("Registry file created:", registry_final_path, "\n")
}

if (!file.exists(registry_history_path)) {
  write.csv(
    data.frame(matrix(ncol = length(registry_cols), nrow = 0,
                      dimnames = list(NULL, registry_cols))),
    registry_history_path, row.names = FALSE
  )
  cat("History file created:", registry_history_path, "\n")
}

cat("Registry initialized\n")