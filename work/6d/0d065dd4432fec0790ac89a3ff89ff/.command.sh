#!/bin/bash -ue
echo "Running FastQC on raw reads: ERR036221"
fastqc ERR036221_1.fastq.gz ERR036221_2.fastq.gz
