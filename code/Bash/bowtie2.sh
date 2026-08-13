#!/bin/bash

# 3. Assembly-based

## 3.2 Bowtie2

### 3.2.1 Bowtie2 index
bowtie2-build \
    --threads 20 \
    /path/to/SAMPLE_ID_contigs.fa \
    /path/to/bowtie2_index/SAMPLE_ID

### 3.2.2 Bowtie2 mapping
bowtie2 \
    -x /path/to/bowtie2_index/SAMPLE_ID \
    -1 SAMPLE_ID_R1.fastq.bz2 \
    -2 SAMPLE_ID_R2.fastq.bz2 \
    -S SAMPLE_ID.sam \
    -p 20 \
    --very-sensitive-local \
    --no-unal

### 3.2.3 SAM to sorted BAM
samtools view \
    -b \
    -@ 20 \
    -o SAMPLE_ID.bam \
    SAMPLE_ID.sam
samtools sort \
    -@ 20 \
    -o SAMPLE_ID.sorted.bam \
    SAMPLE_ID.bam
samtools index \
    -@ 20 \
    SAMPLE_ID.sorted.bam
