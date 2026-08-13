#!/bin/bash

# 3. Assembly-based

## 3.3 MetaBAT2

jgi_summarize_bam_contig_depths \
    --outputDepth SAMPLE_ID_depth.txt \
    SAMPLE_ID.sorted.bam

metabat2 \
    -m 1000 \
    -i /path/to/SAMPLE_ID_contigs.fa \
    -a SAMPLE_ID_depth.txt \
    -o /path/to/output/SAMPLE_ID_metabat \
    --unbinned \
    --seed 0 \
    -t 20
