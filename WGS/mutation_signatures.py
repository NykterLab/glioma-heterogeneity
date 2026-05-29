# Script for running SigProfilerAssignment

from SigProfilerAssignment import Analyzer as Analyze
Analyze.cosmic_fit(samples="path/to/vcf", 
                   output="output_vcf",
                   input_type="vcf",
                   context_type="96",
                   genome_build="GRCh38")
