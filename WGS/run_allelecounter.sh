#!/bin/bash
# Script for running allelecounter, output to be used in DPClust
# Run for each patient separately, in this example for DA1

subject=DA1
zgrep -h -v '^#\|^chrX\|^chrY' $subject*.vcf | /data/scripts/sort_vcf.py -s -I > $subject.hg38.somatic.union.noXY.vcf;
/data/scripts/comp_vcf.py -m $subject.hg38.somatic.union.noXY.vcf $subject.hg38.somatic.union.noXY.vcf > $subject.tmp;
mv -f $subject.tmp $subject.hg38.somatic.union.noXY.vcf;
cut -f 1,2,4,5 $subject.hg38.somatic.union.noXY.vcf > $subject.hg38.somatic.union.noXY.altref.loci;
cut -f 1,2 $subject.hg38.somatic.union.noXY.vcf > $subject.hg38.somatic.union.noXY.loci;

echo *.bam | alleleCounter -d -b $x -l ${x%%-*}.hg38.somatic.union.noXY.loci -o ${x/.bam/.allelecount}
