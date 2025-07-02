library(ArchR)
library(rtracklayer)
library(dplyr)
library(readr)
library(GenomicRanges)
library(BSgenome.Hsapiens.UCSC.hg38)
addArchRThreads(threads = 20)
addArchRGenome("hg38")

datasets <- list(
  GSE158013 = list(
    path = "/data1/renshida/scATAC/data/01_study_Archr/GSE158013/ArchR",
    markerSubdir = "BED_markerPeaks_by_CellType",
    projRds = "peak_res_0.1/GSE158013_peak_celltype.rds"
  ),
  `38956385_Naive` = list(
    path = "/data1/renshida/scATAC/data/01_study_Archr/38956385_Naive/ArchR",
    markerSubdir = "BED_markerPeaks_by_CellType",
    projRds = "peak_res_0.3/38956385_Naive_peak_celltype.rds"
  ),
  `38956385_Prime` = list(
    path = "/data1/renshida/scATAC/data/01_study_Archr/38956385_Prime/ArchR",
    markerSubdir = "BED_markerPeaks_by_CellType",
    projRds = "peak_res_0.3/38956385_Prime_peak_celltype.rds"
  ),
  GSE231708_hg_Lung = list(
    path = "/data1/renshida/scATAC/data/01_study_Archr/GSE231708_hg_Lung/ArchR",
    markerSubdir = "BED_markerPeaks_by_CellType",
    projRds = "peak_res_0.5/GSE231708_hg_Lung_peak_celltype.rds"
  ),
  GSE184386_retina_PF_treat = list(
    path = "/data1/renshida/scATAC/data/01_study_Archr/GSE184386_retina_PF_treat/ArchR",
    markerSubdir = "BED_markerPeaks_by_CellType",
    projRds = "peak_res_0.1/GSE184386_retina_PF_treat_peak_celltype.rds"
  ),
  GSE184386_retina_con = list(
    path = "/data1/renshida/scATAC/data/01_study_Archr/GSE184386_retina_con/ArchR",
    markerSubdir = "BED_markerPeaks_by_CellType",
    projRds = "peak_res_0.1/GSE184386_retina_con_peak_celltype.rds"
  ),
  GSE184386_retina_4D_treat = list(
    path = "/data1/renshida/scATAC/data/01_study_Archr/GSE184386_retina_4D_treat/ArchR",
    markerSubdir = "BED_markerPeaks_by_CellType",
    projRds = "peak_res_0.1/GSE184386_retina_4D_treat_peak_celltype.rds"
  ),
  GSE184386_org_H7_ctrl = list(
    path = "/data1/renshida/scATAC/data/01_study_Archr/GSE184386_org_H7_ctrl/ArchR",
    markerSubdir = "BED_markerPeaks_by_CellType",
    projRds = "peak_res_0.1/GSE184386_org_H7_ctrl_peak_celltype.rds"
  ),
  GSE184386_org_H7_24_treat = list(
    path = "/data1/renshida/scATAC/data/01_study_Archr/GSE184386_org_H7_24_treat/ArchR",
    markerSubdir = "BED_markerPeaks_by_CellType",
    projRds = "peak_res_0.1/GSE184386_org_H7_24_treat_peak_celltype.rds"
  ),
  GSE168669_RESA = list(
    path = "/data1/renshida/scATAC/data/01_study_Archr/GSE168669_RESA/ArchR",
    markerSubdir = "BED_markerPeaks_by_CellType",
    projRds = "peak_res_1/GSE168669_RESA_peak_celltype.rds"
  ),
  GSE183771_P = list(
    path = "/data1/renshida/scATAC/data/01_study_Archr/GSE183771_P/ArchR",
    markerSubdir = "BED_markerPeaks_by_CellType",
    projRds = "peak_res_0.3/GSE183771_P_peak_celltype.rds"
  ),
  GSE183771_S = list(
    path = "/data1/renshida/scATAC/data/01_study_Archr/GSE183771_S/ArchR",
    markerSubdir = "BED_markerPeaks_by_CellType",
    projRds = "peak_res_0.3/GSE183771_S_peak_celltype.rds"
  ),
  OMIX001616 = list(
    path = "/data1/renshida/scATAC/data/01_study_Archr/OMIX001616/ArchR",
    markerSubdir = "BED_markerPeaks_by_CellType",
    projRds = "peak_res_0.5/OMIX001616_peak_celltype.rds"
  ),
  GSE227265 = list(
    path = "/data1/renshida/scATAC/data/01_study_Archr/GSE227265/ArchR",
    markerSubdir = "BED_markerPeaks_by_CellType",
    projRds = "peak_res_0.1/GSE227265_peak_celltype.rds"
  )
)


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