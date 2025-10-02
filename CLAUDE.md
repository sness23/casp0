# CLAUDE.md

This file provides guidance to Claude Code (claude.ai/code) when working with code in this repository.

## Overview

This repository contains CASP (Critical Assessment of protein Structure Prediction) data for protein-ligand binding predictions. The data is organized by CASP edition (casp15, casp16) and includes experimental structures, predictions, results, and pharmaceutical ligand datasets.

## Repository Structure

### CASP15
- `casp15/ligand.csv` - Results data with metrics including target_id, submission_file, lddt_pli, rmsd, lddt_bs, and chain mappings
- `casp15/predictions/ligands/` - Compressed prediction files (.tar.gz format) organized by target ID (e.g., H1171, T1170, T1127)

### CASP16
- `casp16/predictions/ligands/` - Compressed prediction files (.tar.gz, .tgz)
  - Regular targets: R1262, R1264, etc.
  - Pharmaceutical targets: `pharma/` subdirectory with L1000, L2000, L3000, L4000 series
- `casp16/results/ligands/` - Evaluation results
- `casp16/pharma_ligands/pharma_ligands/` - Pharmaceutical datasets:
  - **L1000 (Chymase)**: 17 protein-ligand complexes, both Pose (P) and Affinity (A) prediction
  - **L2000 (Cathepsin G)**: 2 protein-ligand complexes, only Pose (P) prediction
  - **L3000 (Autotaxin)**: 219 protein-ligand complexes (93 P+A, 96 P only, 30 A only)
  - **L4000 (Mpro)**: 25 protein-ligand complexes, only Pose (P) prediction (homodimer, some covalently bound)

## Data Formats

### Pharmaceutical Ligand Datasets

Each pharmaceutical dataset (L1000-L4000) contains:
- `*_README.txt` - Dataset description, tasks, and special considerations
- `*_exper_struct.tar.gz` - Experimental structures (per-ligand .tgz archives)
- `*_exper_affinity.csv` - Affinity data with columns: Target ID, STRUCID, ligand_smiles, Compound ID, IC50 values, Task (PA/P/A), binding_affinity, rank
- `*.SMILES.tar.gz` - SMILES strings for ligands (directory with per-ligand .tsv files)

### Prediction Files
- Archived as .tar.gz or .tgz files
- Named by target ID (e.g., H1171.tar.gz, L1000_L.tgz, L3000.LG0_T.tgz)
- Affinity predictions: .AFFNTY files (e.g., L1000.AFFNTY, L3000.AFFNTY)

### Results Data (casp15/ligand.csv)
Key columns:
- `target_id`, `submission_file`, `model_num`, `pose_num`
- `ref_lig`, `ref_lig_compound`, `ref_lig_formula`
- `mdl_lig`, `mdl_lig_name`, `mdl_lig_formula`
- Metrics: `lddt_pli`, `rmsd`, `lddt_pli_rmsd`, `lddt_pli_symmetry`, `lddt_pli_n_contacts`
- `chain_mapping` (JSON format)
- Binding site metrics: `lddt_bs`, `bs_bb_rmsd`, `bs_num_res`, `bs_num_overlap_res`

## Important Notes

### Pharmaceutical Datasets
- **Stereochemistry**: Known and provided for all ligands
- **Protonation**: May not be optimal and needs determination
- **L3000 (Autotaxin)**:
  - N-glycosylation on Asn497
  - Two zinc ions in active site
  - 22 complexes have multiple ligands (< 4.5 Å from main ligand)
- **L4000 (Mpro)**:
  - Homodimer structure
  - 20 complexes have multiple ligands
  - Some ligands covalently bound to cysteine (L4003, L4013, L4019, L4023)
  - L4007 has S-Hydroperoxycysteine modification on Cys145

## Working with Archives

Extract prediction archives:
```bash
tar -xzf casp15/predictions/ligands/H1171.tar.gz
tar -xzf casp16/pharma_ligands/pharma_ligands/L1000_exper_struct.tar.gz
```

List archive contents:
```bash
tar -tzf casp16/pharma_ligands/pharma_ligands/L1000.SMILES.tar.gz
```

## Data Analysis

View affinity data:
```bash
head casp16/pharma_ligands/pharma_ligands/L1000_exper_affinity.csv
```

Inspect results:
```bash
head casp15/ligand.csv
```
