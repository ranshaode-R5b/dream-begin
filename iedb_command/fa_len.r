#raw_data to iedb_data format
unique_ids <- unique(test_mhcii$hlatype)

# 遍历每个唯一 ID
for (id in unique_ids) {
  # 获取当前 ID 对应的序列和长度
  sequences <- test_mhcii$antigen_peptide[test_mhcii$hlatype == id]
  lengths <- test_mhcii$antigen_peptide_length[test_mhcii$hlatype == id]
  
  id1 <- gsub("[\\*/:]", "", id)  # 去除特定字符
  
  # 按长度分组序列
  split_sequences <- split(sequences, lengths)
  
  # 为每个长度创建对应的文件
  for (length in names(split_sequences)) {
    # 构造文件名
    fasta_file <- paste0("~/mhcii-data/IEDB_data_length/", id1, "_length_", length, ".fa")
    
    # 打开一个文件进行写入（覆盖模式）
    fileConn <- file(fasta_file, open = "a")
    
    # 写入 FASTA 格式，先写入ID，然后写入所有序列
    for (seq in split_sequences[[length]]) {
      writeLines(paste0(">", id), fileConn)
      writeLines(seq, fileConn)  # 写入对应的序列
    }
    
    # 关闭文件连接
    close(fileConn)
  }
}
