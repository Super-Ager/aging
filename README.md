# Multi-Visit Longitudinal Proteomics Analysis of Aging and Longevity

This repository contains R scripts for analyzing protein expression variance decomposition using mixed-effects models based on longitudinal SomaScan 11K proteomics data and trajectory modeling. Pre-computed results can be interactively accessed at [super-ager.github.io](https://super-ager.github.io).

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

## Data Format

The input data should be a CSV file with the following columns:

- `ID`: Individual identifier
- `age`: Age at the time of measurement
- `Sex_F`: Sex indicator (0 = Male, 1 = Female)
- `Visit`: Visit number
- `seq.XXXX.XX`: Protein expression levels (one column per protein, columns starting with "seq.")

### Example Data Format

| ID | age | Sex_F | Visit | seq.3309.2 | seq.4374.45 | seq.1234.56 | ... |
|----|-----|-------|-------|------------|-------------|-------------|-----|
| 001 | 42.0 | 0 | 1 | 8.45 | 7.23 | 9.12 | ... |
| 001 | 47.0 | 0 | 2 | 8.52 | 7.31 | 9.08 | ... |
| 001 | 52.0 | 0 | 3 | 8.58 | 7.38 | 9.05 | ... |
| 001 | 62.0 | 0 | 4 | 8.65 | 7.45 | 9.01 | ... |
| 002 | 65.0 | 1 | 1 | 7.89 | 8.15 | 8.95 | ... |
| 002 | 70.0 | 1 | 2 | 7.95 | 8.22 | 8.87 | ... |
| 002 | 75.0 | 1 | 3 | 8.01 | 8.29 | 8.79 | ... |
| 002 | 85.0 | 1 | 4 | 8.08 | 8.36 | 8.71 | ... |

## Usage

1. Update file paths in each script to match your data locations
2. Ensure input data files are available
3. Run scripts sequentially or independently based on your analysis needs

## Contact

For questions regarding raw data and code, please contact ww20ya@163.com
