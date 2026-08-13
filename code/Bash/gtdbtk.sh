#!/bin/bash

# 3. Assembly-based

## 3.5 GTDB-Tk taxonomic classification

gtdbtk classify_wf \
    --genome_dir /path/to/filtered_MAGs \
    --out_dir /path/to/gtdbtk_output \
    --extension fa \
    --cpus 20
