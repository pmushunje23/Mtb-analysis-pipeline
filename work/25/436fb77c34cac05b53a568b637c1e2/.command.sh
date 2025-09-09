#!/bin/bash -ue
echo "Running Trimmomatic on ERR036223"

trimmomatic PE -threads 2 \
    ERR036223_1.fastq.gz ERR036223_2.fastq.gz \
    ERR036223_R1_paired.fastq.gz ERR036223_R1_unpaired.fastq.gz \
    ERR036223_R2_paired.fastq.gz ERR036223_R2_unpaired.fastq.gz \
    LEADING:3 TRAILING:3 \
    SLIDINGWINDOW:4:15 MINLEN:36
