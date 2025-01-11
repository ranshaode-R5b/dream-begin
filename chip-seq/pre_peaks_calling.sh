#! /bin/bash
#########在peaks calling分析之前可以先重命名样本，这样做的目的是1、用生物学特征去替换SRA数据库ID号，有助于我们了解每个样本是什么并做过什么处理。2、也可以发现有多少样本存在生物学重复，为后续合并合并做准备。


cd /home/renshida/chipseq/rename
vim rename.tsv
# 如下内容，中间的空格是tab键
SRR30273118     control
SRR30273119     p300
SRR30273120     BRD4

path="/home/renshida/chipseq"

for i in {1..3}
do
  id=$(cat rename.tsv | awk '{print $1}'| head -n $i | tail -1)
  filename=$(cat rename.tsv | awk '{print $2}'| head -n $i | tail -1)
  mv ${path}/rmdup/${id}_rmdup.bam ${path}/rename/${filename}.bam
  samtools sort ${path}/rename/${filename}.bam -o ${path}/rename/${filename}_sort.bam
done


#or
#awk '{print $1, $2}' rename.tsv | while read id filename; do
  mv ${path}/rmdup/${id}_rmdup.bam ${path}/rename/${filename}.bam
  samtools sort ${path}/rename/${filename}.bam -o ${path}/rename/${filename}_sort.bam
done
#########如果一个样品分成了多个lane进行测序，那么进行peaks calling的时候，需要把bam进行合并，使用samtools merge, 知道即可，这个数据集不需要,以下代码在正式使用时需要修改！！！merge之后的文件需要重新构建索引。


####merge
mkdir mergeBam
cd align

# 不用循环
# samtools merge control.merge.bam control_1_trimmed.bam control_2_trimmed.bam

# 通常的循环处理方式
# 这里使用“_”进行分割，但是针对不同的数据需要采用不同的分割方式！ 
ls *.bam | cut -d "_" -f 1 |sort -u |while read id; do samtools merge 
../mergeBam/$id.merge.bam $id.bam; done


#######无需合并
mkdir intersect
cd intersect
cat id_norep.tsv 
cp -r ./rename/*_sort.bam ./intersect/

path="/home/renshida/chipseq"

cat ${path}/intersect/id_norep.tsv | while read id;
do
  index_file="${path}/intersect/${id}_sort.bam"  # Target index file path

  # Index the BAM file and specify the index file's location
  samtools index $index_file
done > ${path}/indexing_output.log 2>&1 &

#######需合并
cat ${path}/intersect/id_norep.tsv | while read id;
do
  bedtools intersect -a ${path}/rename/${id}_rep1_sort.bam -b ${path}/rename/${id}_rep2_sort.bam > ${path}/intersect/${id}_intersect.bam       # -a 需要合并的文件1  -b 需要合并的文件
  samtools index ${path}/intersect/${id}_intersect.bam
done