#!/bin/bash

# 3. Assembly-based

## 3.7 PhyloPhlAn

phylophlan_write_default_configs.sh /path/to/configs

phylophlan \
    --input /path/to/filtered_MAGs \
    --output /path/to/phylophlan_output \
    -f /path/to/configs/supermatrix_aa.cfg \
    --diversity high \
    --fast \
    --min_num_markers 50 \
    -d phylophlan \
    --genome_ext fa \
    --nproc 20
