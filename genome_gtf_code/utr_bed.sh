#! /bin/bash
gtf_file="/home/renshida/jieduan1/data/gencode.v19.annotation.gtf_withproteinids"
awk -F '\t' '$3 != "gene"' ${gtf_file} |grep -v "##" > tmp.txt
awk -F '\t' '{print $1,$4,$5,$3}' OFS='\t' 

#提取stop_codon
awk -F '\t' '$3 == "stop_codon"' ${gtf_file} > stop_codon.txt
awk -F '"' '{print $8,$4,$6}' OFS='|' stop_codon.txt >tmp.txt
paste stop_codon.txt tmp.txt | awk -F '\t' '{print $1,$4,$5,$3,$NF,$7}' OFS="\t" > stop_codon.bed
rm tmp.txt stop_codon.txt
#提取start_codon
awk -F '\t' '$3 == "start_codon"' ${gtf_file} > start_codon.txt
awk -F '"' '{print $8,$4,$6}' OFS='|' start_codon.txt >tmp.txt
paste start_codon.txt tmp.txt | awk -F '\t' '{print $1,$4,$5,$3,$NF,$7}' OFS="\t" > start_codon.bed
rm tmp.txt start_codon.txt
#提取transcript
awk -F '\t' '$3 == "transcript"' ${gtf_file} > transcript.txt
awk -F '"' '{print $8,$4,$6}' OFS='|' transcript.txt >tmp.txt
paste transcript.txt tmp.txt | awk -F '\t' '{print $1,$4,$5,$3,$NF,$7}' OFS="\t" > transcript.bed
rm tmp.txt transcript.txt
#提取CDS
awk -F '\t' '$3 == "CDS"' ${gtf_file} > CDS.txt
awk -F '"' '{print $8,$4,$6}' OFS='|' CDS.txt >tmp.txt
paste CDS.txt tmp.txt | awk -F '\t' '{print $1,$4,$5,$3,$NF,$7}' OFS="\t" > CDS.bed
rm tmp.txt CDS.txt
#提取exon
awk -F '\t' '$3 == "exon"' ${gtf_file} > exon.txt
awk -F '"' '{print $8,$4,$6}' OFS='|' exon.txt >tmp.txt
paste exon.txt tmp.txt | awk -F '\t' '{print $1,$4,$5,$3,$NF,$7}' OFS="\t" > exon.bed
rm tmp.txt exon.txt

#提取3pUTR
awk -F '\t' 'NR == FNR {a[$5,$6]=$0; next} ($5,$6) in a {print a[$5,$6], $0}' OFS='\t' stop_codon.bed transcript.bed > stop_trans.bed 
awk -F '\t' '$6 == "+" {print $1,$2,$9, "3'\''UTR",$5,$6}' OFS='\t' stop_trans.bed > zheng.bed
awk -F '\t' '$6 == "-" {print $1,$8,$3, "3'\''UTR",$5,$6}' OFS='\t' stop_trans.bed > fu.bed
cat zheng.bed fu.bed | sort -k1,1 -k2,2n > 3pUTR.bed
rm zheng.bed fu.bed
#提取5pUTR
awk -F '\t' 'NR == FNR {a[$5,$6]=$0; next} ($5,$6) in a {print a[$5,$6], $0}' OFS='\t' start_codon.bed transcript.bed > start_trans.bed 
awk -F '\t' '$6 == "+" {print $1,$2,$9, "5'\''UTR",$5,$6}' OFS='\t' start_trans.bed > zheng.bed
awk -F '\t' '$6 == "-" {print $1,$8,$3, "5'\''UTR",$5,$6}' OFS='\t' start_trans.bed > fu.bed
cat zheng.bed fu.bed | sort -k1,1 -k2,2n > 5pUTR.bed
rm zheng.bed fu.bed

