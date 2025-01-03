#! /bin/bash
path_re="/home/renshida/jieduan1/result1"
path_da="/home/renshida/jieduan1/data"
#提取gene
grep ENSG00000188976.6 $path_da/gencode.v19.annotation.gtf_withproteinids > $path_re/gene4.gtf
awk '$3=="transcript" {OFS="\t";print $1,$4-1,$5,$8,0,$7}' $path_re/gene4.gtf > $path_re/gene4.bed 
grep ENSG00000188976.6 $path_re/v19.intron.bed | awk '{OFS="\t";$1 = "chr" $1;print $1,$2,$3,".",0,$4}' > $path_re/4intron.bed
bedtools subtract -a $path_re/gene4.bed -b $path_re/4intron.bed -s >$path_re/trans4.bed
bedtools getfasta -fi $path_da/GRCh37.p13.genome.fa -bed $path_re/trans4.bed -s -name > $path_re/trans4.fa



awk '$3=="CDS" {OFS="\t";print $1,$4-1,$5,$6,0,$7}' $path_re/gene4.gtf > $path_re/CDS4.bed
bedtools getfasta -fi $path_da/GRCh37.p13.genome.fa -bed $path_re/CDS4.bed -s -name > $path_re/CDS4.fa
