#! /bin/bash
java -jar /home/renshida/apps/Trimmomatic-0.39/trimmomatic-0.39.jar PE \
	 -basein SRR27175483_1.fastq.gz \
	 -baseout SRR27175483_1.paired.fastq.gz \
	 ILLUMINACLIP:/home/renshida/apps/Trimmomatic-0.39/adapters/TruSeq3-PE.fa:2:30:10 \
	 LEADING:3 \
	 TRAILING:3 \
	 SLIDINGWINDOW:4:15 \
	 MINLEN:36				


