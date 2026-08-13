#!/bin/bash

# 3. Assembly-based

## 3.9 ABRicate

### 3.9.1 Virulence factors
abricate \
    --db vfdb \
    /path/to/MAGs/SAMPLE_ID.fa \
    > /path/to/abricate_output/SAMPLE_ID_vfdb.txt

### 3.9.2 Antimicrobial resistance genes
abricate \
    --db card \
    /path/to/MAGs/SAMPLE_ID.fa \
    > /path/to/abricate_output/SAMPLE_ID_card.txt
