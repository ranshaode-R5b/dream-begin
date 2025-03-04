#! /bin/bash
cat trim_id | while read id;
do
	arr=($id)
	fq1=${arr[0]}
	fq2=${arr[1]}
	java -jar /home/renshida/apps/Trimmomatic-0.39/trimmomatic-0.39.jar PE \
	-threads 8 \
	-phred33 \
	/home/renshida/rna-seq/data/fa/$fq1 /home/renshida/rna-seq/data/fa/$fq2 \
	-baseout /home/renshida/rna-seq/data/trim/trimmomatic/$fq1 \
	ILLUMINACLIP:/home/renshida/apps/Trimmomatic-0.39/adapters/TruSeq3-PE.fa:2:30:10 \
	LEADING:3 \
	TRAILING:3 \
	SLIDINGWINDOW:4:15 \
	MINLEN:36;
done      		

	

