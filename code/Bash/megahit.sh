#!/bin/bash

# 3. Assembly-based

## 3.1 MEGAHIT

megahit \
    -1 SAMPLE_ID_R1.fastq.bz2 \
    -2 SAMPLE_ID_R2.fastq.bz2 \
    -m 0.85 \
    -t 8 \
    -o /path/to/output/SAMPLE_ID \
    --min-contig-len 1000 \
    --k-list 21,33,55,71,81,91
