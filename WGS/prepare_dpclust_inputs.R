# Script for creating DPClust input files from AlleleCounter and Battenberg output

# Expects the following input files in the working directory: 
# - samplenames: sample names on rows
# - clustering_input_parameters: column called "filter_counts" for changing SNVs with that many counts to 0
# Expects the following folders in the working directory:
# - allelecounter: allelecounter output and loci files
# - rho_and_psi: Battenberg *rho_and_psi.txt files
# - subclones: Battenberg *subclones.txt files

source("/data/tools/multidimDPClust/R/interconvertMutationBurdens.R")
source("/data/tools/multidimDPClust/R/GetDirichletProcessInfo_multipleSamples_updated_NAP.R")

library(GenomicRanges)

error = function(message) {
  print(paste0("ERROR: ", message))
  quit(save="no")
}

read_cellularities <- function(sample_names) {

	cellularity_values = numeric()
	for (i in 1:length(sample_names)) {
	
	    rho_and_psi_filepath =  choose_file(samplenames[i], folder="rho_and_psi", file_filter=".*hg37_rho_and_psi.txt")
	    rho_and_psi_filepath = paste0("rho_and_psi/", rho_and_psi_filepath)
	    if ( is.null(rho_and_psi_filepath) ) { error(paste0("rho and psi file not found for sample: '", samplenames[i], "'.")) }	    
	    #cat("rho and psi filepath: ", rho_and_psi_filepath, "\n")
	    cellularity_values[i] = GetCellularity(rho_and_psi_filepath)
	}
	cat("cellularities: ", cellularity_values, "\n")
	return( cellularity_values)
}


nucleotides = c("A","C","G","T")

samplenames=read.table("samplenames")
samplenames=as.vector(samplenames[,1])

isMale = TRUE

clustering_input_parameters = read.table("clustering_input_parameters.txt", stringsAsFactors=FALSE, sep="\t", header=T)
filter_counts = clustering_input_parameters$filter_counts


GetWTandMutCount <- function(loci_file, allele_counts_file) {
  subs.data = read.table(loci_file, sep='\t', header=F, stringsAsFactors=F)
  subs.data = subs.data[order(subs.data[,1], subs.data[,2]),] # sort in alphapetical order (required by AlleleCounter)
  
  # Replace dinucleotides and longer with just the first base. Here we assume the depth of the second base is the same and the number of dinucleotides is so low that removing the second base is negligable
  subs.data[,3] = apply(as.data.frame(subs.data[,3]), 1, function(x) { substring(x, 1,1) })
  subs.data[,4] = apply(as.data.frame(subs.data[,4]), 1, function(x) { substring(x, 1,1) })
  
  subs.data.gr = GenomicRanges::GRanges(subs.data[,1], 
  										IRanges::IRanges(subs.data[,2], subs.data[,2]), 
  										rep('*', nrow(subs.data)))
  elementMetadata(subs.data.gr) = subs.data[,c(3,4)]
  
  alleleFrequencies = read.delim(allele_counts_file, sep='\t', header=T, quote=NULL, stringsAsFactors=F)
  alleleFrequencies = alleleFrequencies[order(alleleFrequencies[,1],alleleFrequencies[,2]),]


  print(head(alleleFrequencies))
  alleleFrequencies.gr = GenomicRanges::GRanges(alleleFrequencies[,1], 
  												IRanges::IRanges(alleleFrequencies[,2], alleleFrequencies[,2]), 
  												rep('*', nrow(alleleFrequencies)))
  elementMetadata(alleleFrequencies.gr) = alleleFrequencies[,3:7]
  
  # Subset the allele frequencies by the loci we would like to include
  overlap = findOverlaps(subs.data.gr, alleleFrequencies.gr)
  alleleFrequencies = alleleFrequencies[subjectHits(overlap),]
  
  nucleotides = c("A","C","G","T")
  ref.indices = match(subs.data[,3],nucleotides)
  alt.indices = match(subs.data[,4],nucleotides)
  WT.count = as.numeric(unlist(sapply(1:nrow(alleleFrequencies),function(v,a,i){v[i,a[i]+2]},v=alleleFrequencies,a=ref.indices)))
  mut.count = as.numeric(unlist(sapply(1:nrow(alleleFrequencies),function(v,a,i){v[i,a[i]+2]},v=alleleFrequencies,a=alt.indices)))
  
  combined = data.frame(chr=subs.data[,1],pos=subs.data[,2],WTCount=WT.count, mutCount=mut.count)
  colnames(combined) = c("chr","pos","WT.count","mut.count")
  
  combined.gr = GenomicRanges::GRanges(seqnames(subs.data.gr), ranges(subs.data.gr), rep('*', nrow(subs.data)))
  elementMetadata(combined.gr) = data.frame(WT.count=WT.count, mut.count=mut.count)
  
  combined.gr = sortSeqlevels(combined.gr)
  return(combined.gr)
}

GetCellularity <- function(rho_and_psi_file) {
  d = read.table(rho_and_psi_file, header=T, stringsAsFactors=F)
  return(d['FRAC_GENOME','rho'])
}

choose_file <- function( sample_name, folder=".", file_filter=".*", warn_if_multiple=TRUE ) {

     file_list = list.files( path = folder, pattern = file_filter,
							all.files = FALSE, full.names = FALSE, recursive = FALSE,
                            ignore.case = FALSE, include.dirs = FALSE, no.. = FALSE)
     retval = NULL
     n_found = 0

     for ( f in file_list ) {
         if ( grepl( sample_name, f)  ){
             n_found = n_found + 1
             if ( n_found == 1 ) {
                 retval = f
             }
         }
     }

     if ( n_found != 1 ) {
         if ( warn_if_multiple ) {
             cat("WARNING: ", n_found, " matching files found for sample'", sample_name, "'.\n")
         }
     }
     return( retval)
}


patient_id = sub("-.*", "", samplenames[1])
cat("patient: '", patient_id, "'\n")

loci_filepath = choose_file(patient_id ,"allelecounter", ".*\\.unique.noXY.loci")
if ( is.null(loci_filepath) ) {error(paste0("loci file not found for patient id ", patient_id))} 
loci_filepath = paste0("allelecounter/", loci_filepath)
cat("using loci file: '", loci_filepath, "'\n")

loci_info = read.table( loci_filepath, sep="\t", stringsAsFactors=FALSE, header=FALSE)
loci_info = loci_info[, 1:2]
colnames(loci_info) = c("chr", "pos")

mutCount = data.frame(matrix(NA, ncol=length(samplenames), nrow=nrow(loci_info)), stringsAsFactors=FALSE)
WTCount = data.frame(matrix(NA, ncol=length(samplenames), nrow=nrow(loci_info)), stringsAsFactors=FALSE)

cellularity_values = numeric()
subclone_filepaths = character()

for (i in 1:length(samplenames)) {
	allelecount_filepath = choose_file(samplenames[i], "allelecounter", file_filter=".*\\.allelecount")
	allelecount_filepath = paste0( "allelecounter/", allelecount_filepath)
	if ( is.null(allelecount_filepath) ) {error(paste0("allelecount file not found for sample ", samplenames[i]))} 
	cat("allelecount filepath: ", allelecount_filepath, "\n")

	count_file = GetWTandMutCount(loci_filepath, allelecount_filepath)
	mutCount[,i] = count_file$mut.count
	WTCount[,i] = count_file$WT.count

	rho_and_psi_filepath =  choose_file(samplenames[i], "rho_and_psi", file_filter=".*rho_and_psi.txt")
	rho_and_psi_filepath = paste0( "rho_and_psi/", rho_and_psi_filepath)
	if ( is.null(rho_and_psi_filepath) ) {error(paste0("rho and psi file not found for sample ", samplenames[i]))} 
	cat("rho and psi filepath: ", rho_and_psi_filepath, "\n")
	cellularity_values[i] =  GetCellularity(rho_and_psi_filepath)

	scfp = choose_file(samplenames[i], "subclones", file_filter=".*\\.txt")
	if ( is.null(scfp) ) {error(paste0("subclone file not found for sample ", samplenames[i]))} 
	subclone_filepaths[i] = paste0("subclones/", scfp)

}

cat("mutCount dimensions: ", dim(mutCount), "\n")
cat("WTCount dimensions: ", dim(WTCount), "\n")
cat("cellularities: ", cellularity_values, "\n")

# change low counts to zero
mutCount[mutCount<=filter_counts] = 0
WTCount[WTCount<=filter_counts] = 0


loci_info_right_order = data.frame(chr=as.data.frame(seqnames(count_file)), pos=start(count_file), stringsAsFactors=FALSE)
colnames(loci_info_right_order) = c("chr", "pos")


GetDirichletProcessInfo_multipleSamples(samplenames, cellularity_values, mutCount, WTCount, subclone_filepaths, 
    is.male = isMale, out.files = NULL, phase.dir = NULL, info = loci_info_right_order, 
    SNP.phase.file = NULL, mut.phase.file = NULL, keep.muts.not.explained.by.CN=T)

