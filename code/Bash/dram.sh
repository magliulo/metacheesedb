#!/bin/bash

# 3. Assembly-based

## 3.11 DRAM

### 3.10.1 Functional annotation
DRAM.py annotate \
    -i /path/to/MAGs/SAMPLE_ID.fa \
    -o /path/to/dram_output/SAMPLE_ID_annotation \
    --threads 20 \
    --min_contig_size 300

### 3.10.2 Functional distillation
DRAM.py distill \
    -i /path/to/dram_output/SAMPLE_ID_annotation/annotations.tsv \
    -o /path/to/dram_output/SAMPLE_ID_distillation \
    --trna_path /path/to/dram_output/SAMPLE_ID_annotation/trnas.tsv \
    --rrna_path /path/to/dram_output/SAMPLE_ID_annotation/rrnas.tsv
    
