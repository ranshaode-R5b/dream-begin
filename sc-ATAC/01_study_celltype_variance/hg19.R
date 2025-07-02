library(ArchR)
library(rtracklayer)
library(ggplot2)
library(dplyr)
library(BSgenome.Hsapiens.UCSC.hg19)

addArchRThreads(threads = 20)
addArchRGenome("hg19")

# 数据设置：每组数据集的目录、分辨率、UMAP/Cluster命名、以及cluster -> celltype 映射
datasets <- list(
#  GSE162798 = list(
#    path = "/data1/renshida/scATAC/data/01_study_Archr/GSE162798/ArchR",
#    resolution = "0.3",
#    cluster_col = "Clusters_peak_res0.3",
#    umap_col = "UMAP_peak_res0.3",
#    cluster_to_celltype = c(
#      "C1" = "Myoepithelial cell",
#      "C2" = "Ductal cell",
#      "C3" = "Neutrophil",
#      "C4" = "Fibroblast",
#      "C5" = "Stem cell",
#      "C6" = "Adipocyte",
#      "C7" = "Endothelial cell",
#      "C8" = "Luminal epithelial cell",
#      "C9" = "Natural killer cell",
#      "C10" = "Macrophage"
#    )
#  ),
  GSE214085_Con = list(
    path = "/data1/renshida/scATAC/data/01_study_Archr/GSE214085_Con/ArchR",
    resolution = "0.8",
    cluster_col = "Clusters_peak_res0.8",
    umap_col = "UMAP_peak_res0.8",
    cluster_to_celltype = c(
"C1" = "Unknown1",
"C2" = "Basal epithelial cell",
"C3" = "Alveolar type II cell",
"C4" = "Pulmonary neuroendocrine cell",
"C5" = "Endothelial cell",
"C6" = "Fibroblast",
"C7" = "Club cell",
"C8" = "B cell",
"C9" = "Unknown2",
"C10" = "Unknown3",
"C11" = "T cell",
"C12" = "Natural killer cell",
"C13" = "Unknown4",
"C14" = "Unknown4",
"C15" = "Dendritic cell",
"C16" = "Macrophage",
"C17" = "Dendritic cell"
    )
  )
#    GSE214085_IPA = list(
#    path = "/data1/renshida/scATAC/data/01_study_Archr/GSE214085_IPA/ArchR",
#    resolution = "0.8",
#    cluster_col = "Clusters_peak_res0.8",
#    umap_col = "UMAP_peak_res0.8",
#    cluster_to_celltype = c(
#"C1" = "Alveolar epithelial cell",
#"C2" = "Endothelial cell",
#"C3" = "Endothelial cell",
#"C4" = "Ciliated cell",
#"C5" = "Alveolar type II cell",
#"C6" = "Basal cell",
#"C7" = "Fibroblast",
#"C8" = "Alveolar type II cell",
#"C9" = "Dendritic cell",
#"C10" = "Dendritic cell",
#"C11" = "Unknown",
#"C12" = "Monocyte",
#"C13" = "B cell",
#"C14" = "B cell",
#"C15" = "Mast cell",
#"C16" = "T cell",
#"C17" = "Natural killer cell"
#    )
#  )
#  GSE272504 =  list(
#    path = "/data1/renshida/scATAC/data/01_study_Archr/GSE272504/ArchR",
#    resolution = "0.3",
#    cluster_col = "Clusters_peak_res0.3",
#    umap_col = "UMAP_peak_res0.3",
#    cluster_to_celltype = c(
#"C1" = "Endothelial cell",
#"C2" = "T cell",
#"C3" = "Myoepithelial cell",
#"C4" = "Fibroblast",
#"C5" = "Adipocyte",
#"C6" = "B cell",
#"C7" = "Luminal epithelial cell",
#"C8" = "Basal cell",
#"C9" = "Basal cell",
#"C10" = "Ductal cell"
#    )
#  )
)

macs2_path <- "/home/renshida/miniconda3/envs/scatac_r_410/bin/macs2"

# 遍历每个数据集
for (dataset_name in names(datasets)) {
  info <- datasets[[dataset_name]]
  rds_path <- file.path(info$path, paste0("peak_res_", info$resolution),
                        paste0(dataset_name, "_proj_peak_umap_res_", info$resolution, ".rds"))
  
  # 读取 RDS
  message("Reading: ", rds_path)
  proj <- readRDS(rds_path)
  setwd(info$path)
  # 映射 cluster 到 celltype
  cluster_vec <- getCellColData(proj, select = info$cluster_col)[,1]
  proj$CellType <- mapLabels(cluster_vec, newLabels = info$cluster_to_celltype)

  # UMAP 图 - Sample
  p1 <- plotEmbedding(
    ArchRProj = proj,
    colorBy   = "cellColData",
    name      = "Sample",
    embedding = info$umap_col,
    order     = TRUE,
    randomize = FALSE
  ) + ggtitle(paste0("UMAP Peak (res=", info$resolution, ") - Sample")) +
    theme(text = element_text(size = 16))

  # UMAP 图 - CellType
  p2 <- plotEmbedding(
    ArchRProj = proj,
    colorBy   = "cellColData",
    name      = "CellType",
    embedding = info$umap_col,
    order     = TRUE,
    randomize = FALSE
  ) + ggtitle(paste0("UMAP Peak (res=", info$resolution, ") - CellType")) +
    theme(text = element_text(size = 16))

  # 保存 PDF 图像
  pdf_file <- file.path(info$path, paste0("peak_res_", info$resolution), "UMAP_CellType.pdf")
  pdf(pdf_file, width = 6, height = 6)
  print(p1)
  print(p2)
  dev.off()

  # 生成 group coverage
  proj <- addGroupCoverages(
    ArchRProj = proj,
    groupBy   = "CellType",
    force     = TRUE
  )

  # 生成 reproducible peak set
  proj <- addReproduciblePeakSet(
    ArchRProj   = proj,
    groupBy     = "CellType",
    pathToMacs2 = macs2_path,
    force       = TRUE
  )
  proj <- addPeakMatrix(proj)

# 差异可及性分析（每个 CellType vs 其他）
  markers_peaks <- getMarkerFeatures(
    ArchRProj = proj,
    useMatrix = "PeakMatrix",
    groupBy = "CellType",
    bias = c("TSSEnrichment", "log10(nFrags)"),
    testMethod = "wilcoxon"
  )

  markerGRList <- getMarkers(
  seMarker  = markers_peaks,
  cutOff    = "FDR <= 0.05 & Log2FC >= 1",
  returnGR  = TRUE
)

# 输出目录
bed_dir <- file.path(proj@projectMetadata$outputDirectory, "BED_markerPeaks_by_CellType")
dir.create(bed_dir, showWarnings = FALSE)

# 循环导出每个细胞类型的峰集
for (ct in names(markerGRList)) {
  gr <- markerGRList[[ct]]
  
  # 将 Log₂FC 值写入 score 字段
  mcols(gr)$score <- mcols(gr)$Log2FC
  # 清理不必要的元信息
  mcols(gr) <- mcols(gr)[, "score", drop = FALSE]
  ct_safe <- gsub("[ \\./]", "_", ct)
  
  # 导出为标准 BED 格式
  bed_file <- file.path(bed_dir, paste0(ct_safe, "_markerPeaks.bed"))
  export(gr, con = bed_file, format = "bed")
  
  message("导出 ", ct, " - ", bed_file)
}
  
  bw_dir <- file.path(info$path, paste0("peak_res_", info$resolution), "BigWig_CellType")
  dir.create(bw_dir, showWarnings = FALSE)

  # 生成并导出 bigWig 使用 CPM 归一化，bin 大小为 50bp
  getGroupBW(
    ArchRProj = proj,
    groupBy = "CellType",
    normMethod = "ReadsInTSS",
    tileSize = 50,
    threads = 4
  )

  # 保存更新后的 proj 对象
  saveRDS(proj, file = file.path(info$path, paste0("peak_res_", info$resolution),
                                 paste0(dataset_name, "_peak_celltype.rds")))
}

readRDS()
