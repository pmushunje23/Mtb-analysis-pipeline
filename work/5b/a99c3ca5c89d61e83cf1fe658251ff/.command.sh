#!/bin/bash -ue
echo "Running FastQC on raw reads: ERR036223"
fastqc ERR036223_1.fastq.gz ERR036223_2.fastq.gz
