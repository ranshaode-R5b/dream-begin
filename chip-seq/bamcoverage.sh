mkdir bamCoverage
path="/home/renshida/chipseq"
for i in {1..3}
do
  id_p300=$(cat ${path}/intersect/id_norep.tsv | grep "p300" | head -n ${i} | tail -1)
  id_BRD4=$(cat ${path}/intersect/id_norep.tsv | grep "BRD4" | head -n ${i} | tail -1)
  
  bamCoverage --binSize 10 -p 15 --normalizeUsing RPKM -b ${path}/intersect/${id_p300}_sort.bam -o ${path}/bamCoverage/${id_p300}.bw
  bamCoverage --binSize 10 -p 15 --normalizeUsing RPKM -b ${path}/intersect/${id_BRD4}_sort.bam -o ${path}/bamCoverage/${id_BRD4}.bw
done > ${path}/bamCoverage.log 2>&1 &



#--binSize 10：设置 bin 的大小为 10 bp，这是计算覆盖度时使用的数据区块大小。
#-p 15：线程数为 15
#--normalizeUsing RPKM：一共有RPKM,CPM,BPM,RPGC,None这几种方式，标准化数据
#-b：指定输入的 BAM 文件，根据提取的 ID 动态指定。
#-o：指定输出的 BigWig 文件路径和名称。
