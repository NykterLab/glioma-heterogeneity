#!/bin/bash

# Check input
if [ -z "$1" ]; then
    echo "Usage: $0 <sample_name>"
    exit 1
fi

SAMPLE="$1"
CONTAINER="gatk_${SAMPLE}_4.2.5.0"

# Start container
docker run --name=${CONTAINER} \
    -v /data:/gatk/data \
    --cpus 4 \
    --log-driver=json-file \
    -td broadinstitute/gatk:4.2.5.0

# GATK steps
docker exec ${CONTAINER} gatk HaplotypeCaller \
    --tmp-dir /gatk/data/project/WGS/germline/tmp_gatk \
    -R /gatk/data/homo_sapiens/hg38.fna \
    -I /gatk/data/project/WGS/germline/${SAMPLE}.hg38.rg.bqrecalib.compr.bam \
    -O /gatk/data/project/WGS/germline/${SAMPLE}.4.2.5.0.vcf.gz \
    -ERC GVCF

docker exec ${CONTAINER} gatk GenotypeGVCFs \
    -R /gatk/data/homo_sapiens/hg38.fna \
    -V /gatk/datasample/WGS/germline/${SAMPLE}.4.2.5.0.vcf.gz \
    -O /gatk/data/project/WGS/germline/${SAMPLE}.genotyped.4.2.5.0.vcf.gz

docker exec ${CONTAINER} gatk CNNScoreVariants \
    -V /gatk/data/project/WGS/germline/${SAMPLE}.genotyped.4.2.5.0.vcf.gz \
    -R /gatk/data/homo_sapiens/hg38.fna \
    -O /gatk/data/project/WGS/germline/${SAMPLE}.genotyped.annotated.4.2.5.0.vcf.gz

docker exec ${CONTAINER} gatk FilterVariantTranches \
    -V /gatk/data/project/WGS/germline/${SAMPLE}.genotyped.annotated.4.2.5.0.vcf.gz \
    --resource /gatk/data/broad/resources_broad_hg38_v0_1000G_phase1.snps.high_confidence.hg38.vcf.gz \
    --resource /gatk/data/broad/resources_broad_hg38_v0_Mills_and_1000G_gold_standard.indels.hg38.vcf.gz \
    --info-key CNN_1D \
    --snp-tranche 99.95 \
    --indel-tranche 99.4 \
    -O /gatk/data/project/WGS/germline/${SAMPLE}.genotyped.annotated.filtered.4.2.5.0.vcf.gz

docker exec ${CONTAINER} gatk Funcotator \
    --variant /gatk/data/project/WGS/germline/${SAMPLE}.genotyped.annotated.filtered.4.2.5.0.vcf.gz \
    --reference /gatk/data/homo_sapiens/hg38.fna \
    --ref-version hg38 \
    --data-sources-path /gatk/funcotator_dataSources.v1.6.20190124g \
    --output /gatk/data/project/WGS/germline/${SAMPLE}.genotyped.annotated.filtered.funcotated.vcf \
    --output-file-format VCF

# bcftools steps
bcftools view -i 'FILTER="PASS"' \
    ${SAMPLE}.genotyped.annotated.filtered.funcotated.4.2.5.0.vcf \
    > ${SAMPLE}.genotyped.annotated.filtered.funcotated.PASS.4.2.5.0.vcf

bcftools norm -m-both \
    -o ${SAMPLE}.genotyped.annotated.filtered.4.2.5.0.norm1.vcf \
    ${SAMPLE}.genotyped.annotated.filtered.4.2.5.0.vcf.gz

bcftools norm -f /data/homo_sapiens/hg38.fna \
    -o ${SAMPLE}.genotyped.annotated.filtered.4.2.5.0.norm2.vcf \
    ${SAMPLE}.genotyped.annotated.filtered.4.2.5.0.norm1.vcf

bcftools view -i 'FILTER="PASS"' \
    ${SAMPLE}.genotyped.annotated.filtered.4.2.5.0.norm2.vcf \
    > ${SAMPLE}.genotyped.annotated.filtered.4.2.5.0.norm2.PASS.vcf

# ANNOVAR
perl /data/tools/annovar/table_annovar.pl \
    ${SAMPLE}.genotyped.annotated.filtered.4.2.5.0.norm2.PASS.vcf \
    /data/tools/annovar/humandb/ \
    -buildver hg38 \
    -out ${SAMPLE}.genotyped.annotated.filtered.4.2.5.0 \
    -remove \
    -protocol refGene,avsnp150,gnomad312_genome,clinvar_20221231,dbnsfp42c \
    -operation g,f,f,f,f \
    -nastring . \
    -vcfinput \
    -polish
    