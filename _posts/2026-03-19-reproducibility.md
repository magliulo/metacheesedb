---
title: Reproducibility
author: Raffaele Magliulo
date: 2026-03-19
category: Jekyll
layout: post
---

This page provides a complete guide to the analytical code and computational resources associated with the MetaCheeseDB manuscript. All resources are publicly available to ensure full reproducibility of the reported findings.

---

## Bioinformatics pipelines

All pipelines for raw read preprocessing, taxonomic profiling, genome assembly, binning, quality assessment, and functional annotation are available at [SegataLab/MASTER-WP5-pipelines](https://github.com/SegataLab/MASTER-WP5-pipelines).

Please, refer to the following modules:

- `02-Preprocessing` for raw read quality control and trimming;
- `03-Assemblyfree_analysis` for taxonomic profiling;
- `04-Strainlevel_analysis` for strain-level phylogenomics
- `05-Assembly_pipeline` for assembly, indexing and mapping, and binning;

---

## Machine learning and explainable artificial intelligence analyses scripts

## Reproducibility notes

- Input data required to run the R script below is available [here](https://github.com/magliulo/metacheesedb/blob/main/data/).
- Script can be copied and pasted from below or downloadable [here](https://github.com/magliulo/metacheesedb/blob/main/code/R/ML_and_xAI.R) as R file.
- Update the `setwd()` and file paths in to match your local environment before running.

This [folder](https://github.com/magliulo/metacheesedb/blob/main/code/R/) contains the R script used to reproduce the revised machine learning analyses for [MetaCheeseDB](https://magliulo.github.io/metacheesedb/).

The analysis tests whether cheese microbiome profiles contain country-associated information after taking into account two major issues raised during peer review: possible data leakage between related samples, and confounding by cheesemaking metadata.

The script focuses on **core cheese samples** and should be interpreted as an exploratory sensitivity analysis, not as a ready-to-use traceability tool.

## Script

```text
code/R/ML_and_xAI.R
```

## Main analyses

The script runs four Random Forest models for the Italy versus Spain comparison:

1. **Optimistic sample-level model**  
   Each sample is treated as an independent unit. This is kept only as a comparator.

2. **Leakage-aware grouped model**  
   Samples from the same dataset are kept within the same cross-validation fold.

3. **Confounder-only model**  
   The model uses cheesemaking metadata only. This checks whether country can be predicted from production covariates.

4. **Confounder-adjusted microbiome model**  
   Microbial features are residualised against the available covariates within each training fold before prediction.

The same one-vs-one framework is also applied to the main represented countries: Italy, Spain, Ireland and Austria.

A SHAP analysis is included to identify the microbial taxa contributing most to the Italy versus Spain classification. This step can be slow and can be switched off.

## Required input files

Run the script from the repository root with these files available:

```text
data/MetaCheeseDB_metadata.tsv
data/MetaCheeseDB_taxonomy.tsv
```

The metadata table must include the columns used for filtering, grouping and covariate adjustment, including:

```text
Country
Part_of_cheese
Dataset
Milk_source
Animal_feeding
Milk_processing
Thermisation
Pasteurization
Skimming
Inoculation_of_moulds
Curd_cutting
Presence_of_rind
Backslopping
Starter_culture
Technology_of_production
Coagulation_method
Ripening_period
Rheological_properties
Temperature_of_curd_processing_Celsius_degree
```

The taxonomy table must contain species-level abundance profiles with sample IDs matching the metadata table.

## R packages

The script requires the following R packages:

```r
tidyverse
svglite
caret
pROC
vegan
ranger
writexl
ggsci
scico
kernelshap
shapviz
```

Install missing packages before running the script. For reproducibility, the script saves a `sessionInfo.txt` file in the output folder.

## How to run

Default command from the repository root:

```bash
Rscript code/R/ML_and_xAI.R
```

With explicit input and output paths:

```bash
Rscript code/R/ML_and_xAI.R \
  --metadata data/MetaCheeseDB_metadata.tsv \
  --taxonomy data/MetaCheeseDB_taxonomy.tsv \
  --outdir results/ML_reviewer_sensitivity
```

To skip SHAP:

```bash
Rscript code/R/ML_and_xAI.R --run-shap FALSE
```

## Main outputs

By default, results are written to:

```text
results/ML_reviewer_sensitivity/
```

Main output files include:

```text
06_model_performance_summary.tsv
performance_long_all_pairs.xlsx
multiroc.svg
performance_heatmap.svg
pred_prob.svg
p_varimp_lollipop.svg
p_feature_box.svg
9_3_shap_Country.svg
9_3_shap_Country_bee.svg
SHAP_Italy.svg
SHAP_Spain.svg
sessionInfo.txt
```

## Interpretation

The models are designed to test the robustness of country-associated microbiome patterns, not to prove that geography alone drives cheese microbiome composition.

High performance in the metadata-only model indicates that country labels are linked to cheesemaking variables. The confounder-adjusted microbiome model therefore provides a more cautious estimate of residual microbiome signal after accounting for the measured covariates.

All results should be interpreted together with the metadata covariation and PERMANOVA analyses reported in the manuscript and supplementary materials.

