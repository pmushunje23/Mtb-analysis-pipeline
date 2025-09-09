#!/bin/bash -ue
echo "Running Trimmomatic on ERR036221"

trimmomatic PE -threads 2 \
    ERR036221_1.fastq.gz ERR036221_2.fastq.gz \
    ERR036221_R1_paired.fastq.gz ERR036221_R1_unpaired.fastq.gz \
    ERR036221_R2_paired.fastq.gz ERR036221_R2_unpaired.fastq.gz \
    LEADING:3 TRAILING:3 \
    SLIDINGWINDOW:4:15 MINLEN:36
