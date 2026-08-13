#!/bin/bash

# 2. Assembly-free

## 2.1 MetaPhlAn
metaphlan SAMPLE_ID_R1.fastq.bz2,SAMPLE_ID_R2.fastq.bz2 \
    --input_type fastq \
    --bowtie2db /path/to/db_folder \
    --index mpa_vOct22_CHOCOPhlAnSGB_202403 \
    --output_file SAMPLE_ID_profile.txt \
    --nproc 20 \
    --add_viruses \
    --unclassified_estimation \
    --bowtie2out SAMPLE_ID.bowtie2.bz2 \
    -s SAMPLE_ID.sam.bz2 \
    --force
