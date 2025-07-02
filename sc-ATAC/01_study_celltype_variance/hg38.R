library(ArchR)
library(rtracklayer)
library(ggplot2)
library(dplyr)
library(BSgenome.Hsapiens.UCSC.hg38)

addArchRThreads(threads = 20)
addArchRGenome("hg38")

# 数据设置：每组数据集的目录、分辨率、UMAP/Cluster命名、以及cluster -> celltype 映射
datasets <- list(
  GSE158013 = list(
    path = "/data1/renshida/scATAC/data/01_study_Archr/GSE158013/ArchR",
    resolution = "0.1",
    cluster_col = "Clusters_peak_res0.1",
    umap_col = "UMAP_peak_res0.1",
    cluster_to_celltype = c(
"C1" = "Monocyte",
"C2" = "B cell",
"C3" = "Natural killer cell",
"C4" = "T cell",
"C5" = "Neutrophil",
"C6" = "Dendritic cell"
    )
  ),
  `38956385_Naive` = list(
    path = "/data1/renshida/scATAC/data/01_study_Archr/38956385_Naive/ArchR" ,
    resolution = "0.3",
    cluster_col = "Clusters_peak_res0.3",
    umap_col = "UMAP_peak_res0.3",
    cluster_to_celltype = c(
"C1" = "B cell",
"C2" = "Hepatocyte",
"C3" = "Germ cell",
"C4" = "Neutrophil",
"C5" = "Extraembryonic cell",
"C6" = "T cell",
"C7" = "Enteroendocrine cell",
"C8" = "Fibroblast",
"C9" = "Unknown",
"C10" = "Epithelial cell",
"C11" = "Neuron",
"C12" = "Placental cell"
    )
  ),
  `38956385_Prime` = list(
    path = "/data1/renshida/scATAC/data/01_study_Archr/38956385_Prime/ArchR" ,
    resolution = "0.3",
    cluster_col = "Clusters_peak_res0.3",
    umap_col = "UMAP_peak_res0.3",
    cluster_to_celltype = c(
"C1" = "B cell",
"C2" = "Hepatocyte",
"C3" = "Mesenchymal stem cell",
"C4" = "Germ cell",
"C5" = "Dendritic cell",
"C6" = "Smooth muscle cell",
"C7" = "Monocyte",
"C8" = "Neuron",
"C9" = "Syncytiotrophoblast",
"C10" = "Goblet cell",
"C11" = "Fibroblast",
"C12" = "Hepatocyte",
"C13" = "Endothelial cell"
    )
  ),
  GSE231708_hg_Lung = list(
    path = "/data1/renshida/scATAC/data/01_study_Archr/GSE231708_hg_Lung/ArchR" ,
    resolution = "0.5",
    cluster_col = "Clusters_peak_res0.5",
    umap_col = "UMAP_peak_res0.5",
    cluster_to_celltype = c(
"C1" = "Alveolar epithelial cell",
"C2" = "Ciliated cell",
"C3" = "Basal epithelial cell",
"C4" = "Dendritic cell",
"C5" = "Neutrophil",
"C6" = "Macrophage",
"C7" = "Natural killer cell",
"C8" = "Endothelial cell",
"C9" = "Fibroblast"
    )
  ),
  GSE184386_retina_PF_treat = list(
    path = "/data1/renshida/scATAC/data/01_study_Archr/GSE184386_retina_PF_treat/ArchR" ,
    resolution = "0.1",
    cluster_col = "Clusters_peak_res0.1",
    umap_col = "UMAP_peak_res0.1",
    cluster_to_celltype = c(
"C1" = "Rod cell",
"C2" = "Cone cell",
"C3" = "Horizontal cell",
"C4" = "Microglia",
"C5" = "Amacrine cell",
"C6" = "Bipolar cell",
"C7" = "Müller glia"
    )
  ),
  GSE184386_retina_con = list(
    path = "/data1/renshida/scATAC/data/01_study_Archr/GSE184386_retina_con/ArchR" ,
    resolution = "0.1",
    cluster_col = "Clusters_peak_res0.1",
    umap_col = "UMAP_peak_res0.1",
    cluster_to_celltype = c(
"C1" = "Rod cell",
"C2" = "Horizontal cell",
"C3" = "Amacrine cell",
"C4" = "Müller glia",
"C5" = "Microglia",
"C6" = "Bipolar cell",
"C7" = "Astrocyte",
"C8" = "Cone cell",
"C9" = "Retinal ganglion cell"
    )
  ),
  GSE184386_retina_4D_treat = list(
    path = "/data1/renshida/scATAC/data/01_study_Archr/GSE184386_retina_4D_treat/ArchR" ,
    resolution = "0.1",
    cluster_col = "Clusters_peak_res0.1",
    umap_col = "UMAP_peak_res0.1",
    cluster_to_celltype = c(
"C1" = "Horizontal cell",
"C2" = "Lens epithelial cell",
"C3" = "Retinal ganglion cell",
"C4" = "Microglia",
"C5" = "Astrocyte",
"C6" = "Bipolar cell",
"C7" = "Müller glia",
"C8" = "Cone cell",
"C9" = "Amacrine cell"
    )
  ),
  GSE184386_org_H7_ctrl = list(
    path = "/data1/renshida/scATAC/data/01_study_Archr/GSE184386_org_H7_ctrl/ArchR" ,
    resolution = "0.1",
    cluster_col = "Clusters_peak_res0.1",
    umap_col = "UMAP_peak_res0.1",
    cluster_to_celltype = c(
"C1" = "Bipolar cell",
"C2" = "Rod photoreceptor",
"C3" = "Retinal ganglion cell",
"C4" = "Amacrine cell",
"C5" = "Microglia",
"C6" = "Horizontal cell"
    )
  ),
  GSE184386_org_H7_24_treat = list(
    path = "/data1/renshida/scATAC/data/01_study_Archr/GSE184386_org_H7_24_treat/ArchR" ,
    resolution = "0.1",
    cluster_col = "Clusters_peak_res0.1",
    umap_col = "UMAP_peak_res0.1",
    cluster_to_celltype = c(
"C1" = "Horizontal cell",
"C2" = "Microglia",
"C3" = "Bipolar cell",
"C4" = "Müller glia",
"C5" = "Cone cell",
"C6" = "Amacrine cell",
"C7" = "Astrocyte",
"C8" = "Retinal ganglion cell"
    )
  ),
  GSE168669_RESA = list(
    path = "/data1/renshida/scATAC/data/01_study_Archr/GSE168669_RESA/ArchR" ,
    resolution = "1",
    cluster_col = "Clusters_peak_res1",
    umap_col = "UMAP_peak_res1",
    cluster_to_celltype = c(
"C1" = "Luminal epithelial cell",
"C2" = "Basal epithelial cell",
"C3" = "T cell",
"C4" = "Stromal cell",
"C5" = "Luminal epithelial cell",
"C6" = "Neuroendocrine cell"
    )
  ),
  GSE183771_P = list(
    path = "/data1/renshida/scATAC/data/01_study_Archr/GSE183771_P/ArchR" ,
    resolution = "0.3",
    cluster_col = "Clusters_peak_res0.3",
    umap_col = "UMAP_peak_res0.3",
    cluster_to_celltype = c(
"C1" = "Dendritic cell",
"C2" = "Natural killer cell",
"C3" = "Luminal epithelial cell",
"C4" = "Stromal cell",
"C5" = "Mast cell",
"C6" = "Secretory epithelial cell",
"C7" = "Ciliated epithelial cell",
"C8" = "Unknown",
"C9" = "Secretory epithelial cell"
    )
  ),
  GSE183771_S = list(
    path = "/data1/renshida/scATAC/data/01_study_Archr/GSE183771_S/ArchR" ,
    resolution = "0.3",
    cluster_col = "Clusters_peak_res0.3",
    umap_col = "UMAP_peak_res0.3",
    cluster_to_celltype = c(
"C1" = "Ductal cell",
"C2" = "Fibroblast",
"C3" = "Secretory epithelial cell",
"C4" = "Ciliated epithelial cell",
"C5" = "Natural killer cell",
"C6" = "Dendritic cell",
"C7" = "Endothelial cell",
"C8" = "T cell",
"C9" = "Stromal cell",
"C10" = "Glandular epithelial cell"
    )
  ),
  OMIX001616 = list(
    path = "/data1/renshida/scATAC/data/01_study_Archr/OMIX001616/ArchR" ,
    resolution = "0.5",
    cluster_col = "Clusters_peak_res0.5",
    umap_col = "UMAP_peak_res0.5",
    cluster_to_celltype = c(
"C1" = "Beta cell",
"C2" = "Pancreatic stellate cell",
"C3" = "Mast cell",
"C4" = "Endocrine cell",
"C5" = "Ductal epithelial cell",
"C6" = "Acinar cell",
"C7" = "Acinar cell",
"C8" = "Dendritic cell"
    )
  ),
  GSE227265 = list(
    path = "/data1/renshida/scATAC/data/01_study_Archr/GSE227265/ArchR" ,
    resolution = "0.1",
    cluster_col = "Clusters_peak_res0.1",
    umap_col = "UMAP_peak_res0.1",
    cluster_to_celltype = c(
"C1" = "Endothelial cell",
"C2" = "Hepatic stellate cell",
"C3" = "T cell",
"C4" = "B cell",
"C5" = "Malignant epithelial cell",
"C6" = "Monocyte",
"C7" = "Mast cell",
"C8" = "Liver cancer stem cell",
"C9" = "Epithelial cell",
"C10" = "Hepatocyte",
"C11" = "Malignant cell",
"C12" = "Erythrocyte",
"C13" = "Hepatocyte"
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
}