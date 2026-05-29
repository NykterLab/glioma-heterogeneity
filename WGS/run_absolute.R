# Script for preparing ABSOLUTE inputs and running ABSOLUTE

library(ABSOLUTE)

samples = c("DA1-R1", "DA1-R2", "DA1-R3", "DA1-R4", 
	"GB1-R1", "GB1-R2", "GB1-R3", "GB1-R4", 
	"GB2-R1", "GB2-R2", "GB2-R3", "GB2-R4", 
	"GB2-L1", "GB2-L2", "GB2-L3")

# write loci files from vcf files
for (j in 1:length(samples)) {
	sample = samples[j]

	vcf = read.table(paste0(sample, "/", list.files(sample, ".pass.vcf$")), header=F, stringsAsFactors=F, sep="\t")

	loci = vcf[c(1,2,4,5)]
	colnames(loci) = c("chr", "pos")
	loci = loci[!(loci$chr %in% c("chrX", "chrY")),]

	loci_altref = loci

	loci = loci[1:2]

	write.table(loci, paste0(sample, "/", sample, ".single.noXY.loci"), sep="\t", quote=F, row.names=F, col.names=F)
	write.table(loci_altref, paste0(sample, "/", sample, ".single.noXY.altref.loci"), sep="\t", quote=F, row.names=F, col.names=F)

}


#####################################
# Run allelecounter here using separate loci file for each sample created above
# echo *.bam | alleleCounter -d -b $x -l ${x%%_*}.single.noXY.loci -o ${x/.bam/.allelecount}
#####################################


# Write MAF files from allelecounter files
for (j in 1:length(samples)) {
	sample = samples[j]

	allcount = read.table(paste0("allelecounter_singlesample/", list.files("allelecounter_singlesample", paste0(sample, ".*.allelecount"))), 
		sep="\t", header=F, stringsAsFactors=F)
	loci_altref = read.table(paste0(sample, "/", sample, ".single.noXY.altref.loci"), sep="\t", header=F, stringsAsFactors=F)
	colnames(allcount) = c("chr", "pos", "A", "C", "G", "T", "Good_depth")
	colnames(loci_altref) = c("chr", "pos", "ref", "alt")

	loci_altref$chrpos = paste0(loci_altref$chr,"_", loci_altref$pos)
	allcount$chrpos = paste0(allcount$chr, "_", allcount$pos)

	loci_altref = loci_altref[order(loci_altref$chrpos),]
	allcount = allcount[order(allcount$chrpos),]

	if (all(loci_altref$chrpos==allcount$chrpos)) {
		
		maf = allcount[1:2]
		colnames(maf) = c("Chromosome", "Start_position")
		maf$t_ref_count = sapply(1:nrow(loci_altref), function(i) allcount[i, loci_altref[i, "ref"]])
		maf$t_alt_count = sapply(1:nrow(loci_altref), function(i) allcount[i, loci_altref[i, "alt"]])
		maf$dbSNP_Val_Status = NA
		maf$Tumor_Sample_Barcode = NA
		maf$Hugo_Symbol = NA
	
		maf$Chromosome = sub("chr", "", maf$Chromosome)
		maf = maf[!(nchar(loci_altref$ref)!=1 | nchar(loci_altref$alt)!=1), ]
	
		maf$t_ref_count = unlist(maf$t_ref_count)
		maf$t_alt_count = unlist(maf$t_alt_count)
	
		write.table(maf, paste0(sample, "/", sample, ".maf"), sep="\t", quote=F, row.names=F)
	}

}

# Write seg files for ABSOLUTE input
for (j in 1:length(samples)) {
	sample = samples[j]
	
	logr = read.table(paste0(sample, "/", sample, ".logRsegmented.txt"), sep="\t", header=F, stringsAsFactors=F)	
	baf = read.table(paste0(sample, "/", sample, "_mutantBAF.tab"), sep="\t", header=T, stringsAsFactors=F)
	
	colnames(logr) = c("Chromosome", "Start", "Segment_Mean")
	
	segment_logr = unique(logr$Segment_Mean)
	
	df = data.frame(Chromosome=character(), Start=integer(), End=integer(), Num_Probes=integer(), Segment_Mean=numeric(),
		stringsAsFactors=F)
	
	for (i in 1:length(segment_logr)) {
	
		df[i, "Chromosome"] = logr[logr$Segment_Mean==segment_logr[i], "Chromosome"][1]
		df[i, "Start"] = logr[logr$Segment_Mean==segment_logr[i], "Start"][1]
		df[i, "End"] = logr[logr$Segment_Mean==segment_logr[i], "Start"][nrow(logr[logr$Segment_Mean==segment_logr[i],])]
		df[i, "Segment_Mean"] = segment_logr[i]
	
		if (df[i, "Start"] == df[i, "End"]) {df[i, "End"] = df[i, "End"]+1}
		
		df[i, "Num_Probes"] = sum(baf$Chromosome==df[i, "Chromosome"] & baf$Position>=df[i,"Start"] & baf$Position<df[i, "End"])
	
	}
	# write outputs because absolute wants file paths
	write.table(df, paste0(sample, "/", sample, ".seg"), sep="\t", quote=F, row.names=F)
}


# filter out SNVs with both alt and ref 0 from MAF
for (j in 1:length(samples)) {
	sample = samples[j]
	maf = read.table(paste0(sample, "/", sample, ".maf"), sep="\t", header=T, stringsAsFactors=F)
	zero.ix = which(maf$t_ref_count==0 & maf$t_alt_count==0)
	print(zero.ix)	
	if (length(zero.ix)>0) {maf = maf[-zero.ix,]}
	write.table(maf, paste0(sample, "/", sample, ".filtered.maf"), quote=F, row.names=F, sep="\t")

}

for (j in 1:length(samples)) {

	sample = samples[j]
	
	RunAbsolute(seg.dat.fn = paste0(sample, "/", sample, ".seg"),
		maf.fn = paste0(sample, "/", sample, ".maf"),
		min.mut.af=0.05,
		
		sigma.p=0,
		max.sigma.h=0.02,
		min.ploidy=0.8, # edit
		max.ploidy=10, # edit
		
		primary.disease="Glioblastoma",
		platform="Illumina_WES",
		sample.name=sample,
		results.dir=paste0(sample, "/" ),
	
		max.as.seg.count=1500, 
		max.neg.genome=0.2, # edit --> 0
		max.non.clonal=0.9, # edit --> 0
	
		copy_num_type="total",
		verbose=T)
}

