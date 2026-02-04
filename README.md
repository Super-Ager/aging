This repository contains R scripts for analyzing protein expression variance decomposition using mixed-effects models based on longitudinal SomaScan 11K proteomics data and trajectory modeling.

## Scripts Overview

### 001-mix-effect-model.R

Mixed-effects model analysis for protein expression variance decomposition. Analyzes variance contributions from age, sex, individual differences, and residual components.

### 002-mix-effect-model-bysex.R

Sex-stratified analysis of protein expression variance decomposition. Performs separate analyses for male and female populations.

### 003-mix-effect-statics.R

Statistical analysis and visualization of variance decomposition results. Generates plots showing variance proportions for different components (age, sex, individual, residual).

### 004-crossBOA-validation.R

Cross-validation analysis comparing protein aging effects between USA and Iceland cohorts. Identifies proteins with consistent and significant age-related changes across populations.

### 005-crossBOA-validation-2.R

Analysis of top 100 genes from cross-validation results.

### 006-trajectoryanalysis.R

Protein trajectory modeling and survival analysis. Uses longitudinal data to predict visit 3 protein levels and analyzes the relationship between prediction gaps and mortality risk.

## Requirements

- R (>= 4.0)
- Required packages:
  - `parallel`, `dplyr`, `ggplot2`, `lme4`, `MuMIn`, `lmerTest`
  - `survival`, `glmnet`, `rms`
  - `tidyr`, `ggpubr`, `scales`, `ggrepel`, `ggalluvial`
  - `readxl`, `ggh4x`, `RColorBrewer`

## Usage

1. Update file paths in each script to match your data locations
2. Ensure input data files are available
3. Run scripts sequentially or independently based on your analysis needs

