#! /bin/bash
ln=("lncRNA" "sense_intronic" "sense_overlapping" "lincRNA" "bidirectional_promoter_lncRNA")

ln_str=$(IFS=' '; echo "${ln[*]}")

# 使用 awk 过滤 GTF 文件
 awk -v arr="$ln_str" '
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













																		awk -v arr="$ln_str" '
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
