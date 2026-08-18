---
title: Reproducibility
author: Raffaele Magliulo
date: 2026-03-19
category: Jekyll
layout: post
---

This page provides a complete guide to the computational code resources associated with the MetaCheeseDB manuscript. All resources are publicly available to ensure full reproducibility, as reported below.

---

## Bioinformatics pipelines

All pipelines for raw read preprocessing, taxonomic profiling, genome assembly, binning, quality assessment, and functional annotation are available at the [`code/Bash/`](https://github.com/magliulo/metacheesedb/tree/main/code/Bash) directory, which contains minimal codes used for the MetaCheeseDB analyses. These scripts are intended to document the core commands, parameters, and database versions used in the study. Generic paths should therefore be replaced according to the local installation. Where relevant, the MetaCheeseDB scripts link directly to the corresponding upstream software. Reference genomes used for strain-level analyses correspond to those reported in Supplementary File 2.

---

## Data analysis and visualization

Downstream statistical analyses and figure generation are provided separately in the [`code/R/`](https://github.com/magliulo/metacheesedb/tree/main/code/R) directory. These scripts were from the original working analysis to the intermediate results while retaining the filtering rules, statistical tests, random seeds, and plotting logic used for the manuscript. Each figure script is designed to be run from the repository root, reads its required input files from `data/`, and writes outputs to the corresponding `output/FigX/` directory. The current repository contains standalone scripts for the main manuscript figures, allowing the reported analyses and visualisations to be reproduced without relying on the original local R workspace or absolute file paths.

---

## Data

The [`data/`](https://github.com/magliulo/metacheesedb/tree/main/data) directory contains the processed inputs required by these R scripts, including the harmonised MetaCheeseDB metadata, species-level taxonomic profiles, and alpha-diversity results. For analyses whose complete upstream reconstruction would require computationally intensive genome-resolved or functional annotation workflows, we additionally provide intermediate tables.

---
