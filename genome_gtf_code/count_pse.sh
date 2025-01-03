#! /bin/bash
pse=("IG_pseudogene" "IG_C_pseudogene" "IG_J_pseudogene" "IG_V_pseudogene" "TR_V_pseudogene" "TR_J_pseudogene" "Mt_tRNA_pseudogene" "tRNA_pseudogene" "snoRNA_pseudogene" "snRNA_pseudogene" "scRNA_pseudogene" "rRNA_pseudogene" "misc_RNA_pseudogene" "miRNA_pseudogene" "pseudogene" "processed_pseudogene" "polymorphic_pseudogene" "transcribed_processed_pseudogene" "transcribed_unprocessed_pseudogene" "transcribed_unitary_pseudogene" "translated_processed_pseudogene" "translated_unprocessed_pseudogene" "unitary_pseudogene" "unprocessed_pseudogene")

pse_str=$(IFS=' '; echo "${pse[*]}")

# 使用 awk 过滤 GTF 文件
 awk -v arr="$pse_str" '
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





																		awk -v arr="$pse_str" '
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
