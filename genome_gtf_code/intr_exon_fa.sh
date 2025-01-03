#! /bin/bash
#cat /home/renshida/jieduan1/result/v19.intron.bed | awk '{OFS="\t";print $1,$2,$3,$5,$4}' > /home/renshida/jieduan1/result/intron.sort.bed
#cat /home/renshida/jieduan1/result/intron.sort.bed | awk '{OFS="\t"; $1 = "chr" $1; print }' > /home/renshida/jieduan1/result/intron.sort.geneco.bed
bedtools getfasta -fi /home/renshida/jieduan1/data/GRCh37.p13.genome.fa -bed /home/renshida/jieduan1/result1/intron.sort.geneco.bed -s -name >/home/renshida/jieduan1/result1/intron.fa
#rm /home/renshida/jieduan1/result/intron.sort.merge.ga.bed


#cat /home/renshida/jieduan1/result/exon.bed | awk '{OFS="\t";print $1,$2,$3,$5,$4}' > /home/renshida/jieduan1/result/exon.sort1.bed
bedtools getfasta -fi /home/renshida/jieduan1/data/GRCh37.p13.genome.fa -bed /home/renshida/jieduan1/result1/exon.bed -s -name >/home/renshida/jieduan1/result1/exon.fa

