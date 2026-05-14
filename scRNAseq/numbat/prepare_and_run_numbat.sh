#!/bin/bash
# Script for preparing Numbat inputs and running Numbat
# Author: Serafiina Jaatinen

if [ -z "$1" ]; then
    echo "Usage: $0 <sample_name>"
    exit 1
fi

SAMPLE="$1"

# Numbat preprocessing
# Preparing data
########################
# 1. Phasing and pileup
# 1.1. Perform phasing on the DNA-derived VCF. 
scRNAseq/numbat/run_phasing.sh ${SAMPLE}
# 1.2. Then run cellsnp-lite on scRNA-seq BAMs against the DNA-derived VCF to generate allele counts (only include heterozygous SNPs).
bcftools view -v snps \
    -i 'GT="0/1" || GT="0|1" || GT="1|0"' \
    -o /data/project/WGS/germline/${SAMPLE}.genotyped.annotated.filtered.funcotated.PASS.4.2.5.0.heterozygous.vcf \
    /data/project/WGS/germline/${SAMPLE}.genotyped.annotated.filtered.funcotated.PASS.4.2.5.0.vcf.gz 

bgzip /data/project/WGS/germline/${SAMPLE}.genotyped.annotated.filtered.funcotated.PASS.4.2.5.0.heterozygous.vcf
tabix /data/project/WGS/germline/${SAMPLE}.genotyped.annotated.filtered.funcotated.PASS.4.2.5.0.heterozygous.vcf.gz

for L in L1 L2 L3; do
    cellsnp-lite \
        -s /data/project/scRNAseq/cellranger/${SAMPLE}-${L}/outs/possorted_genome_bam.bam \
        -b /data/project/scRNAseq/cellranger/${SAMPLE}-${L}/outs/filtered_feature_bc_matrix/barcodes.tsv.gz \
        -O /data/project/scRNAseq/numbat/${SAMPLE}-${L}/pileup/${SAMPLE}-${L} \
        -R ${HET_VCF}.gz \
        -p 4 \
        --minMAF 0 \
        --minCOUNT 2 \
        --UMItag Auto \
        --cellTAG CB
done

# 1.3. Merge the phased GT fields (from phased DNA-derived VCF) with the obtained allele counts to produce an allele dataframe in the format of df_allele.
Rscript scRNAseq/numbat/merge_df_allele.R ${SAMPLE}-L1
Rscript scRNAseq/numbat/merge_df_allele.R ${SAMPLE}-L2
Rscript scRNAseq/numbat/merge_df_allele.R ${SAMPLE}-L3

# 2. Prepare the expression data &
# 3. Prepare the expression reference
Rscript scRNAseq/numbat/prepare_expression_data.R
# aggregate_counts in run_numbat.R

#######################
# Run Numbat with adequate resources
Rscript scRNAseq/numbat/run_numbat.R
