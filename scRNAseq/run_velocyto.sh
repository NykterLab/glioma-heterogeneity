#!/bin/bash
# Script to produce loom files with velocyto for velocity analysis with scVelo
ref=/data/references/refdata-gex-GRCh38-2020-A/genes/genes.gtf

samples=("GB2-L1" "GB2-L2" "GB2-L3")

for sample in "${samples[@]}"; do
    echo "Running velocyto on ${sample}"
    velocyto run10x "/data/project/scRNAseq/cellranger/${sample}" "${ref}"
done
