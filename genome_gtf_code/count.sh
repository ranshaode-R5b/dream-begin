#! /bin/bash
protein_coding=("protein_coding" "IG_C_gene" "IG_D_gene" "IG_J_gene" "IG_LV_gene" "IG_V_gene" "TR_C_gene" "TR_J_gene" "TR_V_gene" "TR_D_gene" "antisense")

protein_coding_str=$(IFS=' '; echo "${protein_coding[*]}")

# 使用 awk 过滤 GTF 文件
 awk -v arr="$protein_coding_str" '
 BEGIN {
     n = split(arr, a, " ");
         for (i = 1; i <= n; i++) if (a[i] != "") b[a[i]] = 1;
	 }
$3 == "gene" {
    for (i = 1; i <= NF; i++) {
	        if ($i ~ /gene_type/) {
		split($(i+1), gene_type_arr, /"/);
		if (gene_type_arr[2] in b) {
		gene_type = gene_type_arr[2];
		count[gene_type]++;
		               }
		        }
		}
	}
	END {
	for (type in count) {
	print type, count[type];
	}
}' /home/renshida/jieduan1/data/gencode.v19.annotation.gtf_withproteinids



protein_coding=("protein_coding" "IG_C_gene" "IG_D_gene" "IG_J_gene" "IG_LV_gene" "IG_V_gene" "TR_C_gene" "TR_J_gene" "TR_V_gene" "TR_D_gene" "antisense")

protein_coding_str=$(IFS=' '; echo "${protein_coding[*]}")

# 使用 awk 过滤 GTF 文件
 awk -v arr="$protein_coding_str" '
  BEGIN {
       n = split(arr, a, " ");
                for (i = 1; i <= n; i++) if (a[i] != "") b[a[i]] = 1;
			         }
				 $3 == "transcript" {
				     for (i = 1; i <= NF; i++) {
					                     if ($i ~ /transcript_type/) {
								                     split($(i+1), transcript_type_arr, /"/);
										                     if (transcript_type_arr[2] in b) {
													                     transcript_type = transcript_type_arr[2];
															                     count[transcript_type]++;
																	                                    }
																					                            }
																								                    }
																										            }
																											            END {
																												            for (type in count) {
																														            print type, count[type];
																															            }
																															    }' /home/renshida/jieduan1/data/gencode.v19.annotation.gtf_withproteinids
