#!/bin/bash

# 1. Preprocessing of raw reads

# Link to 'preprocess.new.py' script: https://github.com/SegataLab/preprocessing/blob/master/preprocessing/preprocess.new.py
python3 preprocess.new.py \
    -i /path/to/fastq_folder \
    -s SAMPLE_ID \
    -e .fastq.gz \
    -f _R1 \
    -r _R2 \
    -n 20 \
    --rm_hsap
