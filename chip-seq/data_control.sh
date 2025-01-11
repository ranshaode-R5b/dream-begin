#！ /bin/bash
#### fastp 质控
path="/home/renshida/chipseq"
cd ${path}
mkdir clean

cat SRR_Acc_List.txt | while read id;
do
 nohup fastp \
    -i ${path}/fastq/${id}.fastq.gz \
    -o ${path}/clean/${id}.fq.gz \
    -j ${path}/clean/${id}.fastp.json \
    -h ${path}/clean/${id}.fastp.html &
done


