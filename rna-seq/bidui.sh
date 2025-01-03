#! /bin/bash
#for filename in /home/renshida/rna-seq/data/trim
#cat hisat2_id | while read id;
#do
#       name=($id)
#        pos=${name[0]}
#       rev=${neme[1]}
#	hisat2 -t -p 8 -x /home/renshida/rna-seq/data/ref/rn6/genome \
#	-1 /home/renshida/rna-seq/data/trim/$pos \
#	-2 /home/renshida/rna-seq/data/trim/$rev -S ${id}.sam	
#done

cd /home/renshida/rna-seq/data
cat SRR_Acc_List.txt|while read id;
do
hisat2 -t -p 8 -x /home/renshida/rna-seq/data/ref/rn6/genome \
        -1 /home/renshida/rna-seq/data/trim/${id}_1_val_1.fq.gz\
	-2 /home/renshida/rna-seq/data/trim/${id}_2_val_2.fq.gz -S /home/renshida/rna-seq/data/hisat2/${id}.sam
done
