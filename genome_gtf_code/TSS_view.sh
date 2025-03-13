#! /bin/bash
path_re="/home/renshida/jieduan1/result1"
path_da="/home/renshida/jieduan1/data"
cat $path_re/intron.sort.geneco.bed | awk '$6=="+" {OFS="\t";print $1,$2-10,$2+11,$4,$5,$6}' > $path_re/5ss.pos.bed
cat $path_re/intron.sort.geneco.bed | awk '$6=="-" {OFS="\t";print $1,$3-10,$3+11,$4,$5,$6}' > $path_re/5ss.neg.bed
cat $path_re/intron.sort.geneco.bed | awk '$6=="+" {OFS="\t";print $1,$3-10,$3+11,$4,$5,$6}' > $path_re/3ss.pos.bed
cat $path_re/intron.sort.geneco.bed | awk '$6=="-" {OFS="\t";print $1,$2-10,$2+11,$4,$5,$6}' > $path_re/3ss.neg.bed

bedtools getfasta -fi $path_da/GRCh37.p13.genome.fa -bed $path_re/5ss.pos.bed -s -name > $path_re/5ss.pos.fa
bedtools getfasta -fi $path_da/GRCh37.p13.genome.fa -bed $path_re/5ss.neg.bed -s -name > $path_re/5ss.neg.fa
bedtools getfasta -fi $path_da/GRCh37.p13.genome.fa -bed $path_re/3ss.pos.bed -s -name > $path_re/3ss.pos.fa
bedtools getfasta -fi $path_da/GRCh37.p13.genome.fa -bed $path_re/3ss.neg.bed -s -name > $path_re/3ss.neg.fa
 
grep -v ">" $path_re/5ss.pos.fa > $path_re/5ss.pos.txt
grep -v ">" $path_re/5ss.neg.fa > $path_re/5ss.neg.txt
grep -v ">" $path_re/3ss.pos.fa > $path_re/3ss.pos.txt
grep -v ">" $path_re/3ss.neg.fa > $path_re/3ss.neg.txt



cat $path_re/v19_trans.gtf | awk '$7=="+" {OFS="\t";print $1,$5-50,$5+51,$6,0,$7}' > $path_re/trans.5.pos.bed
cat $path_re/v19_trans.gtf | awk '$7=="-" {OFS="\t";print $1,$4-50,$4+51,$6,0,$7}' > $path_re/trans.5.neg.bed
bedtools getfasta -fi $path_da/GRCh37.p13.genome.fa -bed $path_re/trans.5.pos.bed -s -name > $path_re/trans.5.pos.fa
bedtools getfasta -fi $path_da/GRCh37.p13.genome.fa -bed $path_re/trans.5.neg.bed -s -name > $path_re/trans.5.neg.fa
grep -v ">" $path_re/trans.5.pos.fa > $path_re/trans.5.pos.txt
grep -v ">" $path_re/trans.5.neg.fa > $path_re/trans.5.neg.txt
