#! /bin/bash
protein_coding=("Mt_rRNA" "Mt_tRNA" "miRNA" "misc_RNA" "rRNA" "scRNA" "snRNA" "snoRNA" "ribozyme" "sRNA" "3prime_overlapping_ncRNA" "scaRNA")

protein_coding_str=$(IFS=' '; echo "${protein_coding[*]}")

# 使用 awk 过滤 GTF 文件
 awk -v arr="$protein_coding_str" '
  BEGIN {
       n = split(arr, a, " ");
                for (i = 1; i <= n; i++) if (a[i] != "") b[a[i]] = 1;
			         }
				 $3 == "transcript" {
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
