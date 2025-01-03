#! /bin/bash
protein_coding=$("protein_coding" "IG_C_gene" "IG_D_gene" "IG_J_gene" "IG_LV_gene" "IG_V_gene" "TR_C_gene" "TR_J_gene" "TR_V_gene" "TR_D_gene" "antisense")
protein_coding_str=$(IFS=' ';echo "${protein_coding[*]}")


awk -v arr="$protein_coding_str" 'BEGIN{for(i in split(arr, a, " ")) if (a[i] != "") b[a[i]]=1} $14 in b' /home/renshida/jieduan1/data/gencode.v19.annotation.gtf_withproteinids | sort | uniq -c
