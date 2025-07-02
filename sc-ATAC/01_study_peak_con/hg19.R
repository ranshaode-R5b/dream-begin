library(ArchR)
library(rtracklayer)
library(dplyr)
library(readr)
library(GenomicRanges)

addArchRThreads(threads = 20)
addArchRGenome("hg19")

datasets <- list(
  GSE162798 = list(
    path = "/data1/renshida/scATAC/data/01_study_Archr/GSE162798/ArchR",
    projRds = "peak_res_0.3/GSE162798_peak_celltype.rds"
  ),
  GSE214085_Con = list(
    path = "/data1/renshida/scATAC/data/01_study_Archr/GSE214085_Con/ArchR",
    projRds = "peak_res_0.8/GSE214085_Con_peak_celltype.rds"
  ),
  GSE214085_IPA = list(
    path = "/data1/renshida/scATAC/data/01_study_Archr/GSE214085_IPA/ArchR",
    projRds = "peak_res_0.8/GSE214085_IPA_peak_celltype.rds"
  ),
  GSE272504 = list(
    path = "/data1/renshida/scATAC/data/01_study_Archr/GSE272504/ArchR",
    projRds = "peak_res_0.3/GSE272504_peak_celltype.rds"
  )
)
options(scipen = 999)
# 循环处理每个数据集
for (ds_name in names(datasets)) {
  message("Processing ", ds_name)

  ds_info <- datasets[[ds_name]]
  rds_path <- file.path(ds_info$path, ds_info$projRds)
  
  proj <- readRDS(rds_path)
  peakSet <- getPeakSet(proj)
  
  # 提取细胞类型
  celltypes <- unique(names(peakSet))
  
  # 输出路径
  output_dir <- file.path(ds_info$path, "peak_bed_files")
  dir.create(output_dir, showWarnings = FALSE, recursive = TRUE)
  
  # 按细胞类型输出 BED 文件
  for (ct in celltypes) {
    peaks_ct <- peakSet[names(peakSet) == ct]
    df_ct <- data.frame(
      chr = as.character(seqnames(peaks_ct)),
      start = start(peaks_ct) - 1,  # BED is 0-based
      end = end(peaks_ct),
      score = mcols(peaks_ct)$score
    )
    
    # 文件名安全处理
    safe_ct <- gsub("[^A-Za-z0-9]", "_", ct)
    bed_path <- file.path(output_dir, paste0(safe_ct, ".bed"))
    
    write.table(df_ct, file = bed_path, sep = "\t", quote = FALSE, row.names = FALSE, col.names = FALSE)
  }
}