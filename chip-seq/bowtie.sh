#！ /bin/bash
###index
mkdir -p index/bowtie2/hg38_res
path="/home/renshida/chipseq"
# --threads设置线程数
bowtie2-build ${path}/reference/hg38/hg38.fa  ${path}/index/bowtie2/hg38_res/hg38 --threads 8

###usage
path="/home/renshida/chipseq"
bowtie2_index="/home/lm/Z_Projects/chipseq/index/bowtie2/hg38_res/hg38"
mkdir -p ${path}/align

nohup sh -c 'cat SRR_Acc_List.txt | while read id; do
  bowtie2 -p 8 -x ${bowtie2_index} -U ${path}/clean/${id}.fq.gz | samtools sort -O bam -@ 8 -o - > ${path}/align/${id}.bam
  samtools flagstat ${path}/align/${id}.bam > ${path}/align/${id}_flagstat.txt
done' > ${path}/align/bowtie2_run.log 2>&1 &  		#生成比对统计信息