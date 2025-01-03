#! /bin/bash
#cat /home/renshida/jieduan1/result/utr3.bed | awk '{OFS="\t"} {print $4,$5,$6,$10,$8}' > /home/renshida/jieduan1/result/utr3_fa.bed
#awk '{OFS="\t"} {gsub(/"/, "", $1); gsub(/"/, "", $5); print}' /home/renshida/jieduan1/result/utr3_fa.bed > /home/renshida/jieduan1/result/temp.bed && mv /home/renshida/jieduan1/result/temp.bed /home/renshida/jieduan1/result/utr3_fa.bed
bedtools getfasta -fi /home/renshida/jieduan1/data/GRCh37.p13.genome.fa -bed /home/renshida/jieduan1/result1/gencode_3prime_UTRs.bed -s -name >/home/renshida/jieduan1/result1/utr3.fa

bedtools getfasta -fi /home/renshida/jieduan1/data/GRCh37.p13.genome.fa -bed /home/renshida/jieduan1/result1/gencode_5prime_UTRs.bed -s -name >/home/renshida/jieduan1/result1/utr5.fa
