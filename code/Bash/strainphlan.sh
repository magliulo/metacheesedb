#!/bin/bash

# 2. Assembly-free

## 2.2 StrainPhlAn

### 2.2.1 sample2markers.py
sample2markers.py \
    -i sams/*.sam.bz2 \
    -o consensus_markers/ \
    -d /path/to/mpa_vOct22_CHOCOPhlAnSGB_202403.pkl \
    --nproc 16

### 2.2.2 extract_markers.py and tree generation

#### 2.2.2.1 Streptococcus thermophilus
extract_markers.py -c t__SGB8002 -o db_markers/ -d /path/to/mpa_vOct22_CHOCOPhlAnSGB_202403.pkl
mkdir -p output/Streptococcus_thermophilus
strainphlan \
    -s consensus_markers/*.bz2 \
    -m db_markers/t__SGB8002.fna \
    -r reference_genomes/s__Streptococcus_thermophilus/*.fna.gz \
    -o output/Streptococcus_thermophilus \
    -c t__SGB8002 \
    -d /path/to/mpa_vOct22_CHOCOPhlAnSGB_202403.pkl \
    --phylophlan_mode fast

#### 2.2.2.2 Lactococcus lactis
extract_markers.py -c t__SGB7985 -o db_markers/ -d /path/to/mpa_vOct22_CHOCOPhlAnSGB_202403.pkl
mkdir -p output/Lactococcus_lactis
strainphlan \
    -s consensus_markers/*.bz2 \
    -m db_markers/t__SGB7985.fna \
    -r reference_genomes/s__Lactococcus_lactis/*.fna.gz \
    -o output/Lactococcus_lactis \
    -c t__SGB7985 \
    -d /path/to/mpa_vOct22_CHOCOPhlAnSGB_202403.pkl \
    --phylophlan_mode fast

#### 2.2.2.3 Lacticaseibacillus paracasei
extract_markers.py -c t__SGB7142 -o db_markers/ -d /path/to/mpa_vOct22_CHOCOPhlAnSGB_202403.pkl
mkdir -p output/Lacticaseibacillus_paracasei
strainphlan \
    -s consensus_markers/*.bz2 \
    -m db_markers/t__SGB7142.fna \
    -r reference_genomes/s__Lacticaseibacillus_paracasei/*.fna.gz \
    -o output/Lacticaseibacillus_paracasei \
    -c t__SGB7142 \
    -d /path/to/mpa_vOct22_CHOCOPhlAnSGB_202403.pkl \
    --phylophlan_mode fast

#### 2.2.2.4 Lactobacillus delbrueckii
extract_markers.py -c t__SGB7020 -o db_markers/ -d /path/to/mpa_vOct22_CHOCOPhlAnSGB_202403.pkl
mkdir -p output/Lactobacillus_delbrueckii
strainphlan \
    -s consensus_markers/*.bz2 \
    -m db_markers/t__SGB7020.fna \
    -r reference_genomes/s__Lactobacillus_delbrueckii/*.fna.gz \
    -o output/Lactobacillus_delbrueckii \
    -c t__SGB7020 \
    -d /path/to/mpa_vOct22_CHOCOPhlAnSGB_202403.pkl \
    --phylophlan_mode fast

#### 2.2.2.5 Lactobacillus helveticus
extract_markers.py -c t__SGB7048 -o db_markers/ -d /path/to/mpa_vOct22_CHOCOPhlAnSGB_202403.pkl
mkdir -p output/Lactobacillus_helveticus
strainphlan \
    -s consensus_markers/*.bz2 \
    -m db_markers/t__SGB7048.fna \
    -r reference_genomes/s__Lactobacillus_helveticus/*.fna.gz \
    -o output/Lactobacillus_helveticus \
    -c t__SGB7048 \
    -d /path/to/mpa_vOct22_CHOCOPhlAnSGB_202403.pkl \
    --phylophlan_mode fast

#### 2.2.2.6 Lacticaseibacillus rhamnosus
extract_markers.py -c t__SGB7144 -o db_markers/ -d /path/to/mpa_vOct22_CHOCOPhlAnSGB_202403.pkl
mkdir -p output/Lacticaseibacillus_rhamnosus
strainphlan \
    -s consensus_markers/*.bz2 \
    -m db_markers/t__SGB7144.fna \
    -r reference_genomes/s__Lacticaseibacillus_rhamnosus/*.fna.gz \
    -o output/Lacticaseibacillus_rhamnosus \
    -c t__SGB7144 \
    -d /path/to/mpa_vOct22_CHOCOPhlAnSGB_202403.pkl \
    --phylophlan_mode fast
