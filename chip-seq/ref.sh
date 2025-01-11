#！ /bin/bash
###参考基因组下载
cd
mkdir reference
mkdir reference/hg19 reference/hg38 reference/mm10 reference/mm39

cd reference/hg19
## hg19的基因组文件
wget https://hgdownload.soe.ucsc.edu/goldenPath/hg19/bigZips/hg19.fa.gz
gzip -d hg19.fa.gz
## hg19的基因组注释文件
wget https://ftp.ebi.ac.uk/pub/databases/gencode/Gencode_human/release_42/GRCh37_mapping/gencode.v42lift37.annotation.gtf.gz
gzip -d gencode.v42lift37.annotation.gtf.gz

# 下载hg38（GRCh38）
cd
cd reference/hg38
## hg38的基因组文件
wget http://hgdownload.cse.ucsc.edu/goldenPath/hg38/bigZips/hg38.fa.gz
gzip -d hg38.fa.gz
## hg38的基因组注释文件
wget https://ftp.ebi.ac.uk/pub/databases/gencode/Gencode_human/release_42/gencode.v42.annotation.gtf.gz
gzip -d gencode.v42.annotation.gtf.gz

# 下载mm10（GRCm38）
cd
cd reference/mm10
## mm10的基因组文件
wget ftp://hgdownload.soe.ucsc.edu/goldenPath/mm10/bigZips/mm10.fa.gz
gzip -d mm10.fa.gz
## mm10的基因组注释文件
wget https://ftp.ebi.ac.uk/pub/databases/gencode/Gencode_mouse/release_M25/gencode.vM25.annotation.gtf.gz
gzip -d gencode.vM25.annotation.gtf.gz

# 下载mm39（GRCm39）
cd
cd reference/mm39
## mm39的基因组文件
wget ftp://hgdownload.soe.ucsc.edu/goldenPath/mm39/bigZips/mm39.fa.gz
gzip -d mm39.fa.gz
## mm39的基因组注释文件
wget https://ftp.ebi.ac.uk/pub/databases/gencode/Gencode_mouse/release_M32/gencode.vM32.annotation.gtf.gz
gzip -d gencode.vM32.annotation.gtf.gz



###########基因组数据下载
path="/home/renshida/chipseq"
mkdir -p ${path}/reference/hg38  
#nohub wget -O ${path}/reference/hg38/hg38.fa.gz http://hgdownload.cse.ucsc.edu/goldenPath/hg38/bigZips/hg38.fa.gz &
gzip -d hg38.fa.gz



path="/home/lm/Z_Projects/chipseq"
mkdir -p ${path}/reference/hg38  

#nohup wget -O ${path}/reference/hg38/gencode.v47.annotation.gtf.gz https://ftp.ebi.ac.uk/pub/databases/gencode/Gencode_human/release_47/gencode.v47.annotation.gtf.gz
gzip -d gencode.v47.annotation.gtf.gz

# gff3文件
# nohup wget -O ${path}/reference/hg38/gencode.v47.annotation.gff3.gz https://ftp.ebi.ac.uk/pub/databases/gencode/Gencode_human/release_47/gencode.v47.annotation.gff3.gz

# 解压缩
gzip -d gencode.v47.annotation.gtf.gz