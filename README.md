# TRPV1 allosteric communication analysis: analysis code

Code used to analyse the molecular dynamics (MD) simulations and elastic-network models of TRPV1 reported in:



Release `v0.1-submission` contains the code as it was run for the manuscript submission. It was developed interactively; paths are hard-coded for the original Windows workstation and some scripts need minor adaptation to run elsewhere (see [Running the code](#running-the-code) and [Known limitations](#known-limitations)). A cleaned release with a DOI will accompany the accepted paper.

The repository is meant to let readers (i) audit how each figure was obtained and (ii) reproduce the analysis on their own trajectories.

---

## Contents

- [Systems and data](#systems-and-data)
- [Pipeline overview](#pipeline-overview)
- [Repository layout](#repository-layout)
- [Script reference](#script-reference)
- [Figure map](#figure-map)
- [Residue ranges](#residue-ranges)
- [Support files](#support-files)
- [Requirements](#requirements)
- [Running the code](#running-the-code)
- [Numbering of the scripts](#numbering-of-the-scripts)
- [Known limitations](#known-limitations)
- [License, citation and contact](#license-citation-and-contact)

---

## Systems and data

| Item | Description |
|---|---|
| Starting structures | Cryo-EM structures PDB 7LP9 (apo), 7LPB (capsaicin-bound) and 7LPC (elevated temperature) |
| Variants | WT, W426A, W697A, Y441A (introduced in silico) |
| Replicates | 3 per system, named `rep1`-`rep3` |
| Frames | 501 frames per replicate saved every 0.1 ns (50 ns); the first 101 frames of each replicate were discarded as equilibration |
| File naming | `<structure>_<variant>_rep<N>.pdb` / `.dcd`, e.g. `7LPB_W426A_rep2.dcd` |

The trajectories are **not** included in this repository because of their size. [tbd]. The starting structures are public (PDB IDs above); the chain order of 7LPC differs from the other two structures (see [Known limitations](#known-limitations)).

Chains are labelled A-D. In the path analyses (08b onwards) chain labels are additionally *relativised* (the chain where a path starts becomes "A") so that paths starting in different subunits can be compared.

---

## Pipeline overview

```
 input/md  (PDB + DCD per system)
     |
     +--> 00 initialize --> project_systems_table.csv
     |
     +--> 01 preprocess (align, average) --> aligned_ensemble / aligned_avg (.rds)
     |        |
     |        +--> 02A NMA ------> 02C mode overlap, 03A deformation energy
     |        |
     |        +--> 03a DCCM from MD --> 04a correlation networks
     |                                        |
     |                                        +--> 08a max-flow edge capacities
     |                                                 |
     |                                                 +--> 08b flow decomposition --> 08d change matrices
     |                                                 |                         \--> 08e distogram + flow overlay
     |                                                 +--> 08f seed selection
     |                                                 +--> 08g Markov / fundamental matrix
     |                                                          |
     |                                                          +--> 08h path reconstruction --> 08k paths to PDB
     |                                                          \--> 08m visitation field to PDB
     |
     +--> 06a helix axes --> 06b inter-helical angles --> 06d distributions
     +--> 06e helix bending (Python)
     +--> 07a distance mean/variance --> 07c distograms
     \--> 07b distance time series (protein pairs; capsaicin pairs)
```

---

## Repository layout

```
scripts/
  setup/                   00_initialize_project.R
			   01_build_residue_mapping.R
                           01_v2.1b_preprocess_aligned_ensemble.R
  nma/                     02A_v5.0_nma_compute.R
                           02C_v3.0_nma_mode_overlap_diagnostic.R
                           03A_v2.0_nma_deformation_energy.R
  md_dccm_networks/        03a_v1.2_md_dccm_per_replica.R
                           04a_v5.0_md_correlation_network_build.R
  helix_geometry/          06a_v2.1_helix_axes_from_trajectory.R
                           06b_v3_helix_interhelical_angles.R
                           06d_helix_angle_distributions.R
                           06e_v1.0_helix_bending_analysis.py
  distances/               07a_v2.1_distance_matrices_mean_variance.R
                           07b_v2.1_distance_timeseries_protein_pairs.R
                           07b_v3.1_distance_timeseries_capsaicin_pairs.R
                           07c_v2.1_distograms_from_matrices.R
  flow_pathways/           08a_v1.3_maxflow_edge_capacities.R
                           08b_v3.1_flow_decomposition_paths.R
                           08d_v4.1_flow_change_matrices.R
                           08e_v1.3_distogram_flow_overlay.R
                           08f_v1.3_seed_selection.R
                           08g_v2.1_markov_fundamental_matrix.R
                           08h_v1.4_path_reconstruction.py
                           08k_v3.3_path_clusters_to_pdb.py
                           08m_v1.0_visitation_field_to_pdb_bfactor.R
config.txt
config_files/              trimmed_residue_map.csv
                           trimmed_residue_map_300_752.csv
                           residue_map_MD.csv
README.md
LICENSE
```

The numbers in the file names follow the order in which the analyses were developed and are kept unchanged ([details](#numbering-of-the-scripts)). Folder names group the scripts by analysis block.

---

## Script reference

Output folders are those created by the scripts under their configured output directory.

### setup

| Script | What it does | Main inputs | Main outputs |
|---|---|---|---|
| `00_initialize_project.R` | Finds the `.dcd`/`.pdb` pairs in `input/md`, parses structure, variant and replicate from the file names, flags missing PDBs or replicates, and creates the parameter-registry files | `config.txt`, `input/md/*.dcd`, `input/md/*.pdb` | `project_systems_table.csv`, `registry/parameters_final.csv`, `registry/parameters_history.csv` |
| `01_v2.1b_preprocess_aligned_ensemble.R` | Per system: discards the first 101 frames, selects backbone atoms (N, CA, C, O), least-squares fits each replicate on its starting structure, and averages coordinates across replicates | PDB/DCD pairs of each system | `1_preprocessing/<system>/aligned_ensemble_<system>.rds`, `aligned_avg_<system>.rds`, average PDBs, `validation_summary.csv` |

### nma

| Script | What it does | Main inputs | Main outputs |
|---|---|---|---|
| `02A_v5.0_nma_compute.R` | Cα elastic-network NMA (Bio3D) per replicate and on the replicate-averaged structure (consolidated); fluctuation profiles (first 200 non-trivial modes), DCCM from NMA, inter-replicate convergence check | `aligned_ensemble`/`aligned_avg` (`01`), `config_files/trimmed_residue_map.csv` | `nma_data/` (NMA objects, fluctuations, DCCM, eigenvalues, trimmed Cα PDB); `nma_data/per_rep/` |
| `02C_v3.0_nma_mode_overlap_diagnostic.R` | Overlap (Bio3D `overlap`) between NMA modes of different conditions; used to choose the low-frequency modes for the deformation analysis | `nma_data/` (`02A`), `trimmed_residue_map.csv` | `overlap_diag/` figures |
| `03A_v2.0_nma_deformation_energy.R` | Per-residue local deformation energy (`deformation.nma`) over modes 2-31, 37 and 38; maps it onto Cα and all-atom PDB B-factors | `nma_data/` (`02A`), starting PDB files | `csv/`, `Deformation_plots/`, `Deformation_pdb/`, `Deformation_pdb_full/` |

### md_dccm_networks

| Script | What it does | Main inputs | Main outputs |
|---|---|---|---|
| `03a_v1.2_md_dccm_per_replica.R` | Dynamic cross-correlation matrix (DCCM) of Cα displacements per replicate from the aligned trajectories; mean and variance across replicates | `aligned_ensemble`/`aligned_avg` (`01`), `project_systems_table.csv` (`00`) | `dccm_md/<system>_rep<N>_dccm.rds` and replicate mean/variance matrices |
| `04a_v5.0_md_correlation_network_build.R` | Per replicate, builds a weighted undirected network: edge between residues i, j if \|C_ij\| >= threshold (scan 0.47-0.49; 0.49 reported) and Cα-Cα distance <= 12 Å in the starting structure; edge weight = \|C_ij\|; Walktrap communities | DCCMs (`03a`), starting PDB, `project_systems_table.csv` | `network_objects_md/Th0.47/`, `Th0.48/`, `Th0.49/` (`network_MD_<system>.rds`) |

### helix_geometry

| Script | What it does | Main inputs | Main outputs |
|---|---|---|---|
| `06a_v2.1_helix_axes_from_trajectory.R` | Axis of each helix per frame and chain (first principal component of its Cα atoms) and of the pore (selectivity filter + lower gate residues) | `config.txt` (helix and gate ranges), PDB/DCD, `project_systems_table.csv` | `angles/axes_<system>.rds` |
| `06b_v3_helix_interhelical_angles.R` | Angles between helix axes (S4-Linker, Linker-S5, Linker-S6, S6-TRPh) and S4 tilt relative to the pore axis; 5-frame running mean | axes (`06a`) | `angle_csv/`, `angle_summary/`, `graphics2/` |
| `06d_helix_angle_distributions.R` | Violin plots and histograms of the angles (data from 10 ns onwards) by variant and condition | `angle_csv/` (`06b`) | `distributions2/violins/`, `distributions2/histograms/` |
| `06e_v1.0_helix_bending_analysis.py` | Helix bending (Python/MDAnalysis): local angle per residue (5-Cα window) and global angle per helix (triangle and half-vector methods) | PDB/DCD of each replicate | `06e_helix_analysis/` (CSV, plots, summary PDF) |

### distances

| Script | What it does | Main inputs | Main outputs |
|---|---|---|---|
| `07a_v2.1_distance_matrices_mean_variance.R` | Residue-residue distances (Cβ; Cα for Gly) along frames 101-501, accumulating mean and variance with a single-pass algorithm | PDB/DCD of each replicate | `<system>_mean_dist.rds`, `<system>_var_dist.rds`, `<system>_meta.rds` |
| `07b_v2.1_distance_timeseries_protein_pairs.R` | Time series, distributions and variance of the distance between selected residue pairs (all four subunits) for each variant and condition | PDB/DCD of each replicate | `cache_xyz/`, per-pair figures, `distributions/` |
| `07b_v3.1_distance_timeseries_capsaicin_pairs.R` | Same, for atom-atom distances between protein residues and capsaicin (7LPB systems) | PDB/DCD of 7LPB replicates | per-pair figures, `distributions/` |
| `07c_v2.1_distograms_from_matrices.R` | Combines the three replicates (frame-weighted mean, pooled variance), computes normalized variance and draws distograms, differences between conditions (delta) and between variants (delta-delta) | matrices (`07a`), `residue_map_MD.csv` | `<range>/distograms`, `delta`, `delta_delta` (PNG) |

### flow_pathways

| Script | What it does | Main inputs | Main outputs |
|---|---|---|---|
| `08a_v1.3_maxflow_edge_capacities.R` | Maximum flow (igraph `max_flow`, capacity = \|C_ij\|) between source residues (intracellular and ligand-binding families) and the targets residue 643 (selectivity filter) and 679 (gate) in the four subunits; accumulates per-edge flow over all pairs | networks (`04a`), `residue_map_MD.csv` | `08a_.../rds/`, edge and flow-matrix CSV |
| `08b_v3.1_flow_decomposition_paths.R` | Decomposes the accumulated flow into paths, relativises chains, re-aggregates, aligns the top paths (Smith-Waterman, distance-aware scores) and groups them into path families | `08a` RDS, `residue_map_MD.csv`, starting PDB | `raw/`, `relativized/`, `edges/`, `matrices/`, `final/`, `rivers/`, `trees/` |
| `08d_v4.1_flow_change_matrices.R` | Flux-change matrices between conditions and between variants, with zoom ranges | `08b` matrices | `08d_viz_.../` PNG/CSV by range |
| `08e_v1.3_distogram_flow_overlay.R` | Overlays the flow on the mean-distance distogram | `07c` distograms, `08b` matrices, `residue_map_MD.csv` | PNG by zoom range |
| `08f_v1.3_seed_selection.R` | Ranks residues as candidate sources by outgoing flow towards each target | `08a` RDS, `residue_map_MD.csv` | `seeds_WT_top*.csv`, `seeds_MUT_top*.csv` |
| `08g_v2.1_markov_fundamental_matrix.R` | Row-normalised transition matrix from accumulated flow; target nodes absorbing; fundamental matrix N = (I-Q)^-1; mean and standard error across replicates; edge importance | `08a` RDS, `residue_map_MD.csv` | `per_rep/`, `consolidated/` (RDS), `edge_importance/` |
| `08h_v1.4_path_reconstruction.py` | k = 5 shortest simple paths (Yen) per source seed on the transition graph; path alignment (Needleman-Wunsch, Cβ-distance score), clustering, medoid and consensus path per family | `08g` matrices, `residue_map_MD.csv`, starting PDB | `method1_*/`, `method2_*/` (clusters, alignments, dendrograms), Gephi tables |
| `08k_v3.3_path_clusters_to_pdb.py` | Converts representative paths to poly-glycine PDB models and VMD/PyMOL scripts for rendering | `08h` cluster tables, starting PDB | `<system>/{fwd,rev}/*.pdb`, `.tcl`, `.pml` |
| `08m_v1.0_visitation_field_to_pdb_bfactor.R` | Maps the expected-visitation field of N onto PDB B-factors | `08g` consolidated RDS, starting PDB | `per_combination/`, `consolidated/`, `delta/` PDB files |

The seed residues used in `08h` and `08m` were written into those scripts after inspecting the output of `08f`; they are not read from its tables.

---

## Figure map

| Manuscript item | Scripts |
|---|---|
| Fig. 2f (inter-helical geometry along MD) | `06a` -> `06b` -> `06d` |
| Fig. 3c, Supp. Fig. 3g (normalized-variance matrices) | `07a` -> `07c` |
| Fig. 3e,f and Supp. Fig. 3c,f (NMA deformability) | `01` -> `02A` -> `03A`; mode selection from `02C`|
| Supp. Fig. 3b (capsaicin-protein distances) | `07b_v3.1_distance_timeseries_capsaicin_pairs.R` |
| Fig. 4 (flux-change matrices) | `01` -> `03a` -> `04a` -> `08a` -> `08b` -> `08d` |
| Fig. 5 (pathways) | `08a` -> `08b` -> `08g` -> `08h` -> `08k` (seeds from `08f`) -> `08m`; structural rendering in VMD |

---

## Residue ranges

The three starting structures cover different residue ranges (`config_files/residue_map_MD.csv`):

| Structure | Residues per chain | Range | Missing segment |
|---|---|---|---|
| 7LP9 | 654 | 113-766 | none |
| 7LPB | 618 | 113-752 | 603-624 |
| 7LPC | 454 | 277-752 | 603-624 |

To compare all systems over the same minimal set of residues, the analyses are restricted to a common subset:

- **NMA** (`02A`, `02C`, `03A`): C-alpha atoms of residues 277-602 and 625-752 in the four chains (454 residues per chain, 1816 nodes; `trimmed_residue_map.csv`).
- **Distance matrices** (`07a`): residues 300-602 and 625-752 (431 residues per chain, 1724 nodes), as set in the script.
- **DCCM, networks and flow** (`03a`, `04a`, `08a`-`08m`): .

`config_files/` contains the residue maps (see [Support files](#support-files)).

---

## Support files

| File | Content | Read by |
|---|---|---|
| `config.txt` (repository root) | Directories, network parameters (0.49 and 12 Å) and the helix and gate ranges used by `06a` | `00`, `03a`, `04a`, `06a`, and the parameter-registry blocks of `08a`, `08b`, `08f`, `08g` |
| `config_files/trimmed_residue_map.csv` | NMA residue map: `node_id, chain, residue, system` (`system` = structure ID), residues 277-602 and 625-752 | `02A`, `02C` |
| `config_files/trimmed_residue_map_300_752.csv` | Same columns, residues 300-602 and 625-752 |
| `config_files/residue_map_MD.csv` | Mapping of all residues present in each PDB (`node_id, chain, residue, system`; `system` = structure ID) used to translate network nodes to chain/residue and to map results onto PDB files (see [Residue ranges](#residue-ranges)) | `07c`, `08a`-`08h`, `08m` |

`config.txt` uses the format `parameter = value`; the helix lines (`name=range`, without spaces) must keep that form because `06a` reads them with a plain text search.

---

## Requirements

**R** (RStudio 2025.09.1) with: `bio3d` 2.4-5, `igraph` 2.1.4, `stringr`, `dplyr`, `readr`, `tidyr`, `purrr`, `ggplot2`, `gridExtra`, `zoo`, `mclust`, `fields`, `viridisLite`, `data.table`, `ape`, `ggdendro`, `Matrix`, `future.apply`, `foreach`, `doParallel` (and base `parallel`, `tools`).

**Python** 3.12.7 with `numpy`, `pandas`, `networkx` 3.3, `scipy` 1.13.1, `scikit-learn` 1.5.1, `matplotlib`, `MDAnalysis` 2.10.0.

**Rendering:** VMD 1.9.2 (scripts `.tcl`/`.pml` are also provided for PyMOL); figures were assembled in Inkscape 0.92.4.

Most scripts run in parallel; the number of workers is set at the top of each script (usually 6).

---

## Running the code

1. Place the trajectories and structures as `<structure>_<variant>_rep<N>.pdb/.dcd` in `input/md/` (and in the folder read by `01`, see its header).
2. Edit `config.txt` (kept at the root, where the scripts look for it) and the path variables at the top of each script (they point to `C:/DinamicasMoleculares/...`).
3. Run in numerical order: `00`, `01`, then the blocks you need (see [Pipeline overview](#pipeline-overview)). Scripts are written to be run from RStudio with the project directory as working directory.
4. Each script writes its own log; the parameter registry in `registry/` records parameters of the scripts that write to it.

The output folders created by one script are read by the next one with the names given in [Script reference](#script-reference); if you change a name, update the matching path in the downstream script.

---

## Numbering of the scripts

The numbers reflect the order in which the analyses were developed, and are kept so that the original run order can be traced. Numbers 02, 02B, 02D, 03b, 03B, 04b-04f, 05, 06c, 08c, 08i, 08j, 08l and 08n are absent: they correspond to exploratory analyses, earlier versions of included scripts, or analyses that are not part of the manuscript, and were not necessary to audit or reproduce its results. Where two scripts share a number (07b), they differ in function and are named accordingly.

**Versions in the file names.** The version in each name (e.g. `08a_v1.3_...`) is the development version of the script that produced the results; the same label appears in the script header. `00_initialize_project.R` and `06d_helix_angle_distributions.R` carry no version label. Four scripts received minimal edits for this release that do not change what they compute, described in their headers: `03a` (a definition moved before the parallel export and a registry block copied from another script removed), `07b_v3.1_distance_timeseries_capsaicin_pairs.R` (an undefined object removed from the parallel export), and `06b`, `06d` (text of log messages).

---

## Known limitations

This release is the code as run for the submission. It has the following known limitations; they will be addressed in the release that accompanies the accepted paper.

- **Paths and sessions.** Paths are absolute and Windows-specific; several scripts were developed interactively in RStudio.
- **Different residue sets.** The residue set differs between analysis blocks (NMA: 277-752; distance matrices: 300-752; see [Residue ranges](#residue-ranges)); comparisons across systems are made on a common set of residues.
- **Chain order of 7LPC.** The chain order in the 7LPC structure differs from the others. The relabelling block in `08d_v4.1_flow_change_matrices.R` is commented out. Path analyses assume the tetramer is approximately C4-symmetric when relativised chain labels are compared against distances computed on absolute chains.
- **Direction labels in `08a`.** igraph stores the flow of undirected edges with a sign that follows vertex order (lower to higher vertex ID). The `fwd`/`rev` outputs of `08a`, and the quantities derived from them in `08f`, `08g`, `08h` and `08k`, follow that convention and are not a physical source-to-target direction. As stated in the manuscript, the maximum flow between a source and a target is direction-independent.
- **Path length used in `08b`.** The released script sets a minimum path length of 8, while `08d`, `08e` and the downstream steps read an output folder named `min7`.
- **Parameters read but not applied.** `04a` reads `dccm_threshold` and `contact_cutoff` from `config.txt` but applies a threshold scan (0.47-0.49) and a Cα-Cα cut-off of 12 Å defined inside the script; `config.txt` records the values reported in the manuscript.
- **`08h` path method.** The script implements two path-reconstruction methods; the manuscript describes the transition-matrix method (method 1). 

---

## License, citation and contact

- License: Apache License 2.0 (see `LICENSE`)
- Citation: [tbd]
- Contact corresponding author:sbrauchi@uach.cl
