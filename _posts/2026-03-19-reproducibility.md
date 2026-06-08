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

### Reproducibility notes

- Input data required to run the R script below is available [here](https://github.com/magliulo/metacheesedb/blob/main/data/s_Part_of_cheese.csv).
- Script can be copied and pasted from below or downloadable [here](https://github.com/magliulo/metacheesedb/blob/main/code/R/ML_and_xAI.R) as R file.
- Update the `setwd()` and file paths in to match your local environment before running.
