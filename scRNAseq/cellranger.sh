#!/bin/bash

# SBATCH definitions ...

ref=/data/references/refdata-gex-GRCh38-2020-A/
cellranger=/data/tools/cellranger-7.1.0/cellranger
fastq_dir=/data/project/scRNAseq/raw_data

samples=("GB2-L1" "GB2-L2" "GB2-L3")

for sample in "${samples[@]}"; do
    echo "Running cell ranger count on ${sample}"
    ${cellranger} count --id="${sample}" --transcriptome="${ref}" --fastqs="${fastq_dir}/${sample}/" --sample="${sample}" --disable-ui
done
