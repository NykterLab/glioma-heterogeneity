suppressPackageStartupMessages({
    library(glue)
    library(stringr)
    library(data.table)
    library(dplyr)
    library(vcfR)
    library(Matrix)
    library(numbat)
})

args = commandArgs(trailingOnly=TRUE)

# merge the phased GT fields (from phased DNA-derived VCF) with the obtained allele counts 
# to produce an allele dataframe in the format of df_allele (see section Preparing data).

gmap = "data/tools/Eagle_v2.4.1/tables/genetic_map_hg38_withX.txt.gz"
genome ="hg38"
sample = args[1]
outdir = paste0("/data/project/scRNAseq/numbat/", args[1])
label = sub("-L.*", "", args[1])

## Generate allele count dataframe
cat('Generating allele count dataframes\n')

if (genome == 'hg19') {
    gtf = gtf_hg19
} else {
    gtf = gtf_hg38
}

genetic_map = fread(gmap) %>% 
    setNames(c('CHROM', 'POS', 'rate', 'cM')) %>%
    group_by(CHROM) %>%
    mutate(
        start = POS,
        end = c(POS[2:length(POS)], POS[length(POS)])
    ) %>%
    ungroup()

# read in phased VCF
vcf_phased = lapply(1:22, function(chr) {
        vcf_file = glue('{outdir}/phasing/{label}_chr{chr}.phased.vcf.gz')
        if (file.exists(vcf_file)) {
            fread(vcf_file, skip = '#CHROM') %>%
                rename(CHROM = `#CHROM`) %>%   
                mutate(CHROM = str_remove(CHROM, 'chr'))
        } else {
            stop('Phased VCF not found')
        }
    }) %>%
    Reduce(rbind, .) %>%
    mutate(CHROM = factor(CHROM, unique(CHROM)))

# working with data.table
vcf_phased[, (label) := sub(":.*", "", .SD[[1]]), .SDcols=10]

pu_dir = glue('{outdir}/pileup/{sample}')

# pileup VCF
vcf_pu = fread(glue('{pu_dir}/cellSNP.base.vcf'), skip = '#CHROM') %>% 
    rename(CHROM = `#CHROM`) %>%
    mutate(CHROM = str_remove(CHROM, 'chr'))

# count matrices
AD = readMM(glue('{pu_dir}/cellSNP.tag.AD.mtx'))
DP = readMM(glue('{pu_dir}/cellSNP.tag.DP.mtx'))

cell_barcodes = fread(glue('{pu_dir}/cellSNP.samples.tsv'), header = F) %>% pull(V1)

df = numbat:::preprocess_allele(
    sample = label,
    vcf_pu = vcf_pu,
    vcf_phased = vcf_phased,
    AD = AD,
    DP = DP,
    barcodes = cell_barcodes,
    gtf = gtf,
    gmap = genetic_map
) %>%
filter(GT %in% c('1|0', '0|1'))
    
fwrite(df, glue('{outdir}/{sample}_allele_counts.tsv.gz'), sep = '\t')

cat('All done!\n')