#######fasta  提取
write_fasta <- function(data, file) {
  # 打开文件连接
  con <- file(file, "w")
  
  # 循环遍历数据框的每一行
  for (i in 1:nrow(test_merge)) {
    # 写入序列ID行
    cat(paste0(">", test_merge$allele[i]), "\n", file = con)
    # 写入序列行
    cat(test_merge$peptide[i], "\n", file = con)
  }
  
  # 关闭文件连接
  close(con)
}

# 调用函数并指定输出文件
write_fasta(test_merge, "output.fasta")



########fasta 匹配
unique_ids <- unique(test_mhcii$hlatype)

# 遍历每个唯一 ID
for (id in unique_ids) {
  # 获取当前 ID 对应的序列]
  sequences <- test_mhcii$antigen_peptide[test_mhcii$hlatype == id]
  id1 <- gsub("[\\*/:]","",id)
  # 构造文件名
  fasta_file <- paste0("~/mhcii-data/IEDB_data/", id1, ".fa")
  
  # 打开一个文件进行写入
  fileConn <- file(fasta_file,open = "a")
  #print(paste("Writing ID:", id))
  # 写入 FASTA 格式
  for (seq in sequences) {
    writeLines(paste0(">", id), fileConn)
    writeLines(seq, fileConn)               # 写入对应的序列
  }
  close(fileConn)
}
