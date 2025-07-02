library(ArchR)
library(rtracklayer)
library(ggplot2)
library(dplyr)
library(BSgenome.Mmusculus.UCSC.mm10)

addArchRThreads(threads = 20)
addArchRGenome("mm10")

# 数据设置：每组数据集的目录、分辨率、UMAP/Cluster命名、以及cluster -> celltype 映射
datasets <- list(
   GSE162798 = list(
    path = "/data1/renshida/scATAC/data/01_study_Archr/GSE162798/ArchR_mm",
    resolution = "0.3",
    cluster_col = "Clusters_peak_res0.3",
    umap_col = "UMAP_peak_res0.3",
    cluster_to_celltype = c(
  "C1" = "Endothelial cell",
  "C2" = "Neuroendocrine Cell",
  "C3" = "Alveolar type II cell",
  "C4" = "Basal cell",
  "C5" = "Fibroblast",
  "C6" = "Dendritic cell",
  "C7" = "Natural killer cell",
  "C8" = "Monocyte",
  "C9" = "Macrophage"
    )
  ),
  GSE231708_mm_Liver = list(
    path = "/data1/renshida/scATAC/data/01_study_Archr/GSE231708_mm_Liver/ArchR_mm",
    resolution = "0.1",
    cluster_col = "Clusters_peak_res0.1",
    umap_col = "UMAP_peak_res0.1",
    cluster_to_celltype = c(
"C1" = "Hepatic progenitor cell",
"C2" = "Hepatocyte",
"C3" = "Venous endothelial cell",
"C4" = "Kupffer cell",
"C5" = "T cell",
"C6" = "Stromal cell",
"C7" = "Endothelial cell"
    )
  ),
  GSE231708_mm_Lung = list(
    path = "/data1/renshida/scATAC/data/01_study_Archr/GSE231708_mm_Lung/ArchR_mm",
    resolution = "0.1",
    cluster_col = "Clusters_peak_res0.1",
    umap_col = "UMAP_peak_res0.1",
    cluster_to_celltype = c(
"C1" = "Alveolar Type II Cell",
"C2" = "Club cell",
"C3" = "Ependymal cell",
"C4" = "Club cell",
"C5" = "Alveolar Type II Cell",
"C6" = "Macrophage",
"C7" = "Neutrophil",
"C8" = "Unknown",
"C9" = "Smooth Muscle Cell",
"C10" = "Fibroblast"
    )
  ),
  GSE242466_mm_COL = list(
    path = "/data1/renshida/scATAC/data/01_study_Archr/GSE242466_mm_COL/ArchR_mm",
    resolution = "0.1",
    cluster_col = "Clusters_peak_res0.1",
    umap_col = "UMAP_peak_res0.1",
    cluster_to_celltype = c(
"C1" = "Innate lymphoid cell",
"C2" = "Neutrophil",
"C3" = "Dendritic cell",
"C4" = "Mast cell",
"C5" = "Macrophage",
"C6" = "Eosinophil",
"C7" = "CD8 T cell",
"C8" = "Th17 cell",
"C9" = "Eosinophil",
"C10" = "B cell",
"C11" = "Myeloid",
"C12" = "Dendritic cell",
"C13" = "Unknown",
"C14" = "Unknown"
    )
  ),
  GSE242466_mm_SKI = list(
    path = "/data1/renshida/scATAC/data/01_study_Archr/GSE242466_mm_SKI/ArchR_mm",
    resolution = "0.1",
    cluster_col = "Clusters_peak_res0.1",
    umap_col = "UMAP_peak_res0.1",
    cluster_to_celltype = c(
"C1" = "Mast cell",
"C2" = "B cell",
"C3" = "Natural killer cell",
"C4" = "T cell",
"C5" = "Macrophage",
"C6" = "Dendritic cell",
"C7" = "Macrophage",
"C8" = "Monocyte",
"C9" = "Basophil"
    )
  ),
  GSE242466_mm_SPL = list(
    path = "/data1/renshida/scATAC/data/01_study_Archr/GSE242466_mm_SPL/ArchR_mm",
    resolution = "0.1",
    cluster_col = "Clusters_peak_res0.1",
    umap_col = "UMAP_peak_res0.1",
    cluster_to_celltype = c(
"C1" = "Monocyte",
"C2" = "Neutrophil",
"C3" = "Dendritic cell",
"C4" = "Natural killer cell",
"C5" = "Regulatory T Cell",
"C6" = "T memory cell",
"C7" = "B cell",
"C8" = "Macrophage",
"C9" = "CD4 T cell",
"C10" = "Unknown"
    )
  ),
  GSE242466_mm_VAT = list(
    path = "/data1/renshida/scATAC/data/01_study_Archr/GSE242466_mm_VAT/ArchR_mm",
    resolution = "0.1",
    cluster_col = "Clusters_peak_res0.1",
    umap_col = "UMAP_peak_res0.1",
    cluster_to_celltype = c(
"C1" = "Neutrophil",
"C2" = "Neutrophil",
"C3" = "Dendritic cell",
"C4" = "Dendritic cell",
"C5" = "Macrophage",
"C6" = "B cell",
"C7" = "Monocyte",
"C8" = "Macrophage",
"C9" = "B cell",
"C10" = "T cell",
"C11" = "Natural Killer Cell",
"C12" = "Innate lymphoid cell",
"C13" = "T memory cell",
"C14" = "Macrophage"
    )
  ),
  `E-MTAB-11264` = list(
    path = "/data1/renshida/scATAC/data/01_study_Archr/E-MTAB-11264/ArchR_mm",
    resolution = "0.1",
    cluster_col = "Clusters_peak_res0.1",
    umap_col = "UMAP_peak_res0.1",
    cluster_to_celltype = c(
"C1" = "Astrocyte",
"C2" = "Oligodendrocyte",
"C3" = "Oligodendrocyte",
"C4" = "Interneuron",
"C5" = "Pyramidal neuron",
"C6" = "Ependymal Cell",
"C7" = "unknown",
"C8" = "Microglia",
"C9" = "Endothelial cell"
    )
  ),
  GSE164439 = list(
    path = "/data1/renshida/scATAC/data/01_study_Archr/GSE164439/ArchR_mm",
    resolution = "0.5",
    cluster_col = "Clusters_peak_res0.5",
    umap_col = "UMAP_peak_res0.5",
    cluster_to_celltype = c(
"C1" = "Germ cell",
"C2" = "Blood vessel",
"C3" = "Endothelial cell",
"C4" = "Peritubular myoid cell",
"C5" = "Stromal cell",
"C6" = "Peritubular myoid cell",
"C7" = "Sertoli cell",
"C8" = "Leydig cell"
    )
  ),
  GSE149622_IL23 = list(
    path = "/data1/renshida/scATAC/data/01_study_Archr/GSE149622_IL23/ArchR_mm",
    resolution = "0.3",
    cluster_col = "Clusters_peak_res0.3",
    umap_col = "UMAP_peak_res0.3",
    cluster_to_celltype = c(
"C1" = "Th17 cell",
"C2" = "T cell",
"C3" = "Dendritic cell",
"C4" = "Keratinocyte",
"C5" = "Regulatory T cell",
"C6" = "T memory cell"
    )
  ),
  GSE149622_naive = list(
    path = "/data1/renshida/scATAC/data/01_study_Archr/GSE149622_naive/ArchR_mm",
    resolution = "0.3",
    cluster_col = "Clusters_peak_res0.3",
    umap_col = "UMAP_peak_res0.3",
    cluster_to_celltype = c(
"C1" = "B cell",
"C2" = "Dendritic cell",
"C3" = "T cell",
"C4" = "Natural killer cell",
"C5" = "Fibroblast",
"C6" = "Lymphatic endothelial cell"
    )
  ),
  GSE211077_Ascl = list(
    path = "/data1/renshida/scATAC/data/01_study_Archr/GSE211077_Ascl/ArchR_mm",
    resolution = "0.3",
    cluster_col = "Clusters_peak_res0.3",
    umap_col = "UMAP_peak_res0.3",
    cluster_to_celltype = c(
"C1" = "Microglia",
"C2" = "Retinal ganglion cell",
"C3" = "Müller glia",
"C4" = "Müller glia",
"C5" = "Neuron",
"C6" = "Bipolar cell",
"C7" = "Photoreceptor Cell"
    )
  ),
  GSE211077_IPA = list(
    path = "/data1/renshida/scATAC/data/01_study_Archr/GSE211077_IPA/ArchR_mm",
    resolution = "0.3",
    cluster_col = "Clusters_peak_res0.3",
    umap_col = "UMAP_peak_res0.3",
    cluster_to_celltype = c(
"C1" = "Microglia",
"C2" = "Müller glia",
"C3" = "Photoreceptor",
"C4" = "Bipolar Cell",
"C5" = "Neuron",
"C6" = "Müller glia",
"C7" = "Retinal ganglion cell"
    )
  ),
  GSE211077_Retinal = list(
    path = "/data1/renshida/scATAC/data/01_study_Archr/GSE211077_Retinal/ArchR_mm",
    resolution = "0.3",
    cluster_col = "Clusters_peak_res0.3",
    umap_col = "UMAP_peak_res0.3",
    cluster_to_celltype = c(
"C1" = "Müller glia",
"C2" = "Bipolar cell",
"C3" = "Neuron",
"C4" = "Photoreceptor cell",
"C5" = "Retinal ganglion cell",
"C6" = "Retinal ganglion cell"
    )
  ),
  GSE139950 = list(
    path = "/data1/renshida/scATAC/data/01_study_Archr/GSE139950/ArchR_mm",
    resolution = "0.1",
    cluster_col = "Clusters_peak_res0.1",
    umap_col = "UMAP_peak_res0.1",
    cluster_to_celltype = c(
"C1" = "Podocyte",
"C2" = "Proximal tubule cell",
"C3" = "Proximal tubule cell-Inj",
"C4" = "Intercalated cell",
"C5" = "Distal tubule cell",
"C6" = "Loop of henle cell",
"C7" = "Intercalated cell",
"C8" = "Macrophage",
"C9" = "T cell",
"C10" = "Endothelial cell",
"C11" = "Myofibroblast"
    )
  ),
  GSE214082_CD = list(
    path = "/data1/renshida/scATAC/data/01_study_Archr/GSE214082_CD/ArchR_mm",
    resolution = "0.1",
    cluster_col = "Clusters_peak_res0.1",
    umap_col = "UMAP_peak_res0.1",
    cluster_to_celltype = c(
"C1" = "Endothelial cell",
"C2" = "Smooth muscle cell",
"C3" = "Fibroblast",
"C4" = "Macrophage",
"C5" = "T cell"
    )
  ),
  GSE214082_HFD = list(
    path = "/data1/renshida/scATAC/data/01_study_Archr/GSE214082_HFD/ArchR_mm",
    resolution = "0.1",
    cluster_col = "Clusters_peak_res0.1",
    umap_col = "UMAP_peak_res0.1",
    cluster_to_celltype = c(
"C1" = "Fibroblast",
"C2" = "Fibroblast",
"C3" = "Smooth muscle cell",
"C4" = "Endothelial cell",
"C5" = "T cell",
"C6" = "B cell",
"C7" = "Macrophage",
"C8" = "Dendritic cell"
    )
  )
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

  # 可选：保存更新后的 proj 对象
  saveRDS(proj, file = file.path(info$path, paste0("peak_res_", info$resolution),
                                 paste0(dataset_name, "_peak_celltype.rds")))
  message(">>>  ", dataset_name, " 已完成...")
}

