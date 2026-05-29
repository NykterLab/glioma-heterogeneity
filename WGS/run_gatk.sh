#!/bin/bash
# Script for GATK Mutect2 variant calling

DATA_DIR="/data"
REF="${DATA_DIR}/references/homo_sapiens/hg38.fna"

BAM_DIR="${DATA_DIR}/tumor_bams"
OUT_DIR="${DATA_DIR}/vcf"

NORMAL_BAM="${DATA_DIR}/germline/normal_sample.bam" # edit
NORMAL_NAME="normal_sample" # edit

PON="${DATA_DIR}/pon/panel_of_normals.vcf.gz"

GNOMAD="${DATA_DIR}/reference/gnomad.hg38.vcf.gz"

DOCKER_IMAGE="broadinstitute/gatk:latest"
CONTAINER="gatk_mutect2"

THREADS=4

mkdir -p "${OUT_DIR}"

sudo docker run \
    --name=${CONTAINER} \
    -v ${DATA_DIR}:/gatk/data \
    --cpus ${THREADS} \
    --log-driver=json-file \
    -td ${DOCKER_IMAGE}

for BAM in ${BAM_DIR}/*.bam; do

    BAM_BASENAME=$(basename "${BAM}")
    SAMPLE="${BAM_BASENAME%%.*}"
    echo "Processing ${SAMPLE}"

    NORMAL_BAM=$(find "${NORMAL_DIR}" -name "${SAMPLE%%-*}*.bam" | head -n 1)
    NORMAL_NAME=$(basename "${NORMAL_BAM}")
    echo "Matched normal: ${NORMAL_NAME}"

    OUTPUT_VCF="${OUT_DIR}/${SAMPLE}.mutect2.vcf.gz"
    FILTERED_VCF="${OUT_DIR}/${SAMPLE}.mutect2.filtered.vcf.gz"
    PASS_VCF="${OUT_DIR}/${SAMPLE}.mutect2.filtered.pass.vcf.gz"

    # Mutect2
    sudo docker exec ${CONTAINER} gatk Mutect2 \
        -R "/gatk${REF}" \
        -I "/gatk${BAM}" \
        -I "/gatk${NORMAL_BAM}" \
        -normal ${NORMAL_NAME} \
        -O "/gatk${OUTPUT_VCF}" \
        --pon "/gatk${PON}" \
        --germline-resource "/gatk${GNOMAD}" \
        --af-of-alleles-not-in-resource 0.0000625 \
        --disable-read-filter MateOnSameContigOrNoMappedMateReadFilter

    tabix -p vcf "${OUTPUT_VCF}"

    # Filter
    sudo docker exec ${CONTAINER} gatk FilterMutectCalls \
        -R "/gatk${REF}" \
        -V "/gatk${OUTPUT_VCF}" \
        -O "/gatk${FILTERED_VCF}"

    tabix -p vcf "${FILTERED_VCF}"

    bcftools view \
        -f PASS \
        -Oz \
        -o "${PASS_VCF}" \
        "${FILTERED_VCF}"

    tabix -p vcf "${PASS_VCF}"

    # Normalize
    bcftools norm \
        -m-both \
        -Oz \
        -o "${OUT_DIR}/${SAMPLE}.norm1.vcf.gz" \
        "${PASS_VCF}"

    bcftools norm \
        -f "${REF}" \
        -Oz \
        -o "${OUT_DIR}/${SAMPLE}.norm2.vcf.gz" \
        "${OUT_DIR}/${SAMPLE}.norm1.vcf.gz"

    tabix -p vcf "${OUT_DIR}/${SAMPLE}.norm2.vcf.gz"

    # Annotate
    perl /data/tools/annovar/table_annovar.pl \
        "${OUT_DIR}/${SAMPLE}.norm2.vcf.gz" \
        /data/tools/annovar/humandb/ \
        -out "${OUT_DIR}/${SAMPLE}" \
        -vcfinput \
        -buildver hg38 \
        -remove \
        -protocol refGene,avsnp150,exac03,dbnsfp42c,cosmic98_coding,cosmic98_noncoding,gnomad312_genome,clinvar_20221231 \
        -operation g,f,f,f,f,f,f,f \
        -nastring . \
        -polish

done
