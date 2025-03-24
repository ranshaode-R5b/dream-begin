#! /bin/bash
cat /home/renshida/rna-seq/data/SRR_Acc_List.txt|while read id;
do
	featureCounts -a /home/renshida/rna-seq/data/ref/rn6.ncbiRefSeq.gtf -o /home/renshida/rna-seq/data/subread/${id}_final_counts.txt -T 8 -p /home/renshida/rna-seq/data/hisat2/trimmomatic/sort/${id}.sort.bam
done
