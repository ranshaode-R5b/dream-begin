path="/home/renshida/chipseq"

cd ${path}
mkdir macs2

for i in {1..3}
do
  id_p300=$(cat ${path}/intersect/id_norep.tsv | grep "p300" | head -n ${i} | tail -1)
  id_BRD4=$(cat ${path}/intersect/id_norep.tsv | grep "BRD4" | head -n ${i} | tail -1)
  
  macs2 callpeak -t ${path}/intersect/${id_p300}_sort.bam \
    --control ${path}/intersect/control_sort.bam \
    -g hs --outdir ${path}/macs2 \
    -n ${id_p300}_intersect -q 0.01
  
  macs2 callpeak -t ${path}/intersect/${id_BRD4}_sort.bam \
    --control ${path}/intersect/control_sort.bam \
    -g hs --outdir ${path}/macs2 \
    -n ${id_BRD4}_intersect -q 0.01

done > ${path}/macs2_output.log 2>&1 &

# macs2 callpeak参数：
# -t 输入文件
# --control 对照组文件。
# -g 基因组大小。输入具体的数字（如1.0e+9或1000000000）。
# 对于人、小鼠、线虫或果蝇，可以分别用hs、mm、ce或dm。
# --outdir 输出文件的路径
# -n 文件名前缀
# -q Q值