#! /bin/bash
cat /home/renshida/rna-seq/data/SRR_Acc_List.txt|while read id;
do
        samtools sort -o /home/renshida/rna-seq/data/hisat2/${id}.sort.bam /home/renshida/rna-seq/data/hisat2/${id}.sam
done
