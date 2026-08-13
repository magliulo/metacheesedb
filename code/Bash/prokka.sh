#!/bin/bash

# 3. Assembly-based

## 3.8 Prokka functional annotation

prokka \
    --addgenes \
    --compliant \
    --gcode 11 \
    --outdir /path/to/prokka_output/SAMPLE_ID \
    --prefix SAMPLE_ID \
    --force \
    --cpus 20 \
    /path/to/MAGs/SAMPLE_ID.fa
