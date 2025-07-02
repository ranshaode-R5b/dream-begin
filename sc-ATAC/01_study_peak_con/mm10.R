library(ArchR)
library(rtracklayer)
library(dplyr)
library(readr)
library(GenomicRanges)
library(BSgenome.Mmusculus.UCSC.mm10)
addArchRThreads(threads = 20)
addArchRGenome("mm10")

datasets <- list(
  GSE162798 = list(
    path = "/data1/renshida/scATAC/data/01_study_Archr/GSE162798/ArchR_mm",
    markerSubdir = "BED_markerPeaks_by_CellType",
    projRds = "peak_res_0.3/GSE162798_peak_celltype.rds"
  ),
  GSE231708_mm_Liver = list(
    path = "/data1/renshida/scATAC/data/01_study_Archr/GSE231708_mm_Liver/ArchR_mm",
    markerSubdir = "BED_markerPeaks_by_CellType",
    projRds = "peak_res_0.1/GSE231708_mm_Liver_peak_celltype.rds"
  ),
  GSE231708_mm_Lung = list(
    path = "/data1/renshida/scATAC/data/01_study_Archr/GSE231708_mm_Lung/ArchR_mm",
    markerSubdir = "BED_markerPeaks_by_CellType",
    projRds = "peak_res_0.1/GSE231708_mm_Lung_peak_celltype.rds"
  ),
  GSE242466_mm_COL = list(
    path = "/data1/renshida/scATAC/data/01_study_Archr/GSE242466_mm_COL/ArchR_mm",
    markerSubdir = "BED_markerPeaks_by_CellType",
    projRds = "peak_res_0.1/GSE242466_mm_COL_peak_celltype.rds"
  ),
  GSE242466_mm_SKI = list(
    path = "/data1/renshida/scATAC/data/01_study_Archr/GSE242466_mm_SKI/ArchR_mm",
    markerSubdir = "BED_markerPeaks_by_CellType",
    projRds = "peak_res_0.1/GSE242466_mm_SKI_peak_celltype.rds"
  ),
  GSE242466_mm_SPL = list(
    path = "/data1/renshida/scATAC/data/01_study_Archr/GSE242466_mm_SPL/ArchR_mm",
    markerSubdir = "BED_markerPeaks_by_CellType",
    projRds = "peak_res_0.1/GSE242466_mm_SPL_peak_celltype.rds"
  ),
  GSE242466_mm_VAT = list(
    path = "/data1/renshida/scATAC/data/01_study_Archr/GSE242466_mm_VAT/ArchR_mm",
    markerSubdir = "BED_markerPeaks_by_CellType",
    projRds = "peak_res_0.1/GSE242466_mm_VAT_peak_celltype.rds"
  ),
  `E-MTAB-11264` = list(
    path = "/data1/renshida/scATAC/data/01_study_Archr/E-MTAB-11264/ArchR_mm",
    markerSubdir = "BED_markerPeaks_by_CellType",
    projRds = "peak_res_0.1/E-MTAB-11264_peak_celltype.rds"
  ),
  GSE164439 = list(
    path = "/data1/renshida/scATAC/data/01_study_Archr/GSE164439/ArchR_mm",
    markerSubdir = "BED_markerPeaks_by_CellType",
    projRds = "peak_res_0.5/GSE164439_peak_celltype.rds"
  ),
  GSE149622_IL23 = list(
    path = "/data1/renshida/scATAC/data/01_study_Archr/GSE149622_IL23/ArchR_mm",
    markerSubdir = "BED_markerPeaks_by_CellType",
    projRds = "peak_res_0.3/GSE149622_IL23_peak_celltype.rds"
  ),
  GSE149622_naive = list(
    path = "/data1/renshida/scATAC/data/01_study_Archr/GSE149622_naive/ArchR_mm",
    markerSubdir = "BED_markerPeaks_by_CellType",
    projRds = "peak_res_0.3/GSE149622_naive_peak_celltype.rds"
  ),
  GSE211077_Ascl = list(
    path = "/data1/renshida/scATAC/data/01_study_Archr/GSE211077_Ascl/ArchR_mm",
    markerSubdir = "BED_markerPeaks_by_CellType",
    projRds = "peak_res_0.3/GSE211077_Ascl_peak_celltype.rds"
  ),
  GSE211077_IPA = list(
    path = "/data1/renshida/scATAC/data/01_study_Archr/GSE211077_IPA/ArchR_mm",
    markerSubdir = "BED_markerPeaks_by_CellType",
    projRds = "peak_res_0.3/GSE211077_IPA_peak_celltype.rds"
  ),
  GSE211077_Retinal = list(
    path = "/data1/renshida/scATAC/data/01_study_Archr/GSE211077_Retinal/ArchR_mm",
    markerSubdir = "BED_markerPeaks_by_CellType",
    projRds = "peak_res_0.3/GSE211077_Retinal_peak_celltype.rds"
  ),
  GSE139950 = list(
    path = "/data1/renshida/scATAC/data/01_study_Archr/GSE139950/ArchR_mm",
    markerSubdir = "BED_markerPeaks_by_CellType",
    projRds = "peak_res_0.1/GSE139950_peak_celltype.rds"
  ),
  GSE214082_CD = list(
    path = "/data1/renshida/scATAC/data/01_study_Archr/GSE214082_CD/ArchR_mm",
    markerSubdir = "BED_markerPeaks_by_CellType",
    projRds = "peak_res_0.1/GSE214082_CD_peak_celltype.rds"
  ),
  GSE214082_HFD = list(
    path = "/data1/renshida/scATAC/data/01_study_Archr/GSE214082_HFD/ArchR_mm",
    markerSubdir = "BED_markerPeaks_by_CellType",
    projRds = "peak_res_0.1/GSE214082_HFD_peak_celltype.rds"
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