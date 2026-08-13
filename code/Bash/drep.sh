#!/bin/bash

# 3. Assembly-based

## 3.6 dRep dereplication

dRep dereplicate /path/to/drep_output \
    --processors 20 \
    --genomes /path/to/bins_path.txt \
    --completeness 50 \
    --contamination 10 \
    --checkM_method lineage_wf \
    --S_algorithm ANImf \
    --MASH_sketch 1000 \
    --P_ani 0.90 \
    --S_ani 0.95 \
    --cov_thresh 0.10
