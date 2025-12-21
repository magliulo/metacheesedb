---
layout: home
title: MetaCheeseDB
permalink: /
---

### Metadata-curated Cheese Database

**MetaCheeseDB** is an open source metadata-curated catalogue of cheese metagenomes for studying the cheese microbiome across production technologies, geography, and cheese style. It integrates **1,593** shotgun metagenomes spanning **156 cheese subtypes** from **19 countries** across **4 continents**, combining publicly available datasets with newly sequenced **PDO/PGI** cheeses.

Our goal is simple: make cheese metagenomics *comparable*. MetaCheeseDB pairs high-resolution microbiome data with **cheese-specific, production-aware metadata** so that researchers can test ecological and functional hypotheses at scale, and stakeholders can explore microbiome signatures relevant to quality, authenticity, and innovation.

---

## What makes MetaCheeseDB different?

### 1) Cheese-native metadata (not just sequencing descriptors)
Public repositories often lack standardised descriptors that matter for cheese ecology. For each sample we compiled **42 metadata fields**, including both technical descriptors and food-/process-related variables. We additionally provide a harmonised subset of **18 standardised metadata fields** designed for robust cross-study comparisons.

Examples of curated cheese-relevant fields include:
- milk source and animal management (where available)
- thermal treatment (e.g., pasteurisation)
- starter strategy (including backslopping / natural whey culture)
- curd treatment and ripening context
- rind presence and sampled compartment (core vs rind)
- texture class (soft / semi-hard / hard / very hard)
- geographical origin and certification status (PDO/PGI where applicable)

---

## What’s inside the database?

### Sample-level profiles (assembly-free)
- Species-level taxonomic profiles for all samples
- Diversity summaries and metadata-stratified comparisons
- Explainable machine learning outputs highlighting candidate microbial “biomarkers” for key metadata classes

### Genome-resolved catalogue (assembly-based)
- **4,170** high- and medium-quality MAGs (HQ/MQ)
- Dereplicated into **337** representative species-level genome bins (SGBs)
- Known vs unknown SGB designation (kSGB/uSGB), supporting novelty discovery
- Phylogenomics and functional annotation to connect taxa with metabolic potential

### Functional potential
Genome annotations include technologically relevant and health-/safety-adjacent traits, such as:
- fermentation/ripening-relevant functions (e.g., carbohydrate-active enzymes, proteases)
- vitamin biosynthesis pathways
- antimicrobial resistance (AMR) and virulence factor screening (as genetic potential)

---

## How to use MetaCheeseDB

### Explore
Use the [browsing interface](https://magliulo.github.io/raffaelemagliulo/files/MetaCheeseDB.html) to filter samples by:
- texture class, milk source, starter strategy, thermal treatment, rind/core, country, cheese subtype, and more.

### Download
MetaCheeseDB provides:
- curated metadata tables
- processed microbiome profiles
- genome-resolved outputs (MAG/SGB summaries and annotations)
- figures and intermediate results used in the analyses

### Reproduce
The full processing and analysis code is openly available (see the “Code” section). Raw reads remain hosted in their original public repositories; accession identifiers are provided in the metadata.

---

## Suggested applications

MetaCheeseDB is designed to support both fundamental and translational work, including:
1) **Defined-culture (“pitched”) cheese design**  
   Prioritise taxa/lineages associated with specific cheese styles and processing contexts, then validate candidates experimentally.

2) **PDO/PGI fingerprinting, traceability, and authenticity**  
   Link official production disciplinaries and cheese technologies to reproducible community signatures, supporting hypothesis generation for authenticity models.

3) **Health-oriented strain discovery (with appropriate caution)**  
   Identify strains encoding genetic potential for bioactive pathways as candidates for downstream validation, while recognising that gene presence does not imply expression or in-product activity.

---

## Important notes
- MetaCheeseDB reports **genetic potential**, not guaranteed activity. Functional claims require validation (e.g., culture work, transcriptomics/proteomics/metabolomics).
- As with all cross-study meta-analyses, residual batch effects may exist despite extensive harmonisation.

---

## Citation
A manuscript describing MetaCheeseDB, the harmonised metadata framework, and the main analyses is in preparation. In the meantime, please cite MetaCheeseDB as a resource and reference the associated manuscript when available.

---

## Contact / contributions
Feedback and contributions (new datasets, metadata corrections, feature requests) are welcome. Please use the repository issue tracker or the contact details provided on the site.

