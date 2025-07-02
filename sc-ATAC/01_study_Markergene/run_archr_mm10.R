#!/usr/bin/env Rscript
suppressPackageStartupMessages({
    library(optparse)
    library(ArchR)
    library(ggplot2)
    library(dplyr)
    library(BSgenome.Mmusculus.UCSC.mm10)
})

# -----------------------
# 1. 解析命令行参数
# -----------------------
option_list <- list(
    make_option(c("--study_id"), type="character", default=NULL,
                help="Study ID ", metavar="character"),
    make_option(c("--base_dir"), type="character", default=NULL,
                help="Base directory path", metavar="character"),
    make_option(c("--fragments"), type="character", default=NULL,
                help="Path to file listing all fragment.tsv.gz paths", metavar="character")
)
opt_parser <- OptionParser(option_list = option_list)
opt <- parse_args(opt_parser)

if (is.null(opt$study_id) || is.null(opt$base_dir) || is.null(opt$fragments)) {
    print_help(opt_parser)
    stop("请提供 --study_id, --base_dir, --fragments 三个参数。", call.=FALSE)
}

study_id   <- opt$study_id
base_dir   <- normalizePath(opt$base_dir)
frag_list  <- normalizePath(opt$fragments)

# 研究目录
study_dir <- file.path(base_dir, study_id)
if (!dir.exists(study_dir)) {
    message(">>> 研究目录 ", study_dir, " 不存在，正在创建...")
    dir.create(study_dir, recursive = TRUE)
}

# ArchR 输出目录
archr_outdir <- file.path(study_dir, "ArchR_mm")
if (!dir.exists(archr_outdir)) {
    message(">>> 输出目录 ", archr_outdir, " 不存在，正在创建...")
    dir.create(archr_outdir, recursive = TRUE)
}
message("==== 开始处理研究：", study_id, " ====")

# -----------------------
# 2. 读取 fragment 列表并生成 sample_names 和 fragment_paths
# -----------------------
all_frags <- readLines(frag_list)
if (length(all_frags) == 0) {
    stop("在 fragments 文件中未读取到任何路径，请检查。")
}

# 提取每个 fragment 文件的文件名，再去掉后缀 ".fragment.tsv.gz" 或 ".fragments.tsv.gz"
extractSampleName <- function(filepath) {
    fname <- basename(filepath)
    sub(pattern = "(\\.?fragments?\\.tsv\\.gz)$", replacement = "", x = fname)
}
sample_names   <- sapply(all_frags, extractSampleName)
fragment_paths <- all_frags

# 如果有重复的样本名，添加后缀
if (any(duplicated(sample_names))) {
    dup_idx <- which(duplicated(sample_names) | duplicated(sample_names, fromLast = TRUE))
    sample_names[dup_idx] <- paste0(sample_names[dup_idx], "_", seq_along(sample_names)[dup_idx])
}

message("共检测到 ", length(fragment_paths), " 个 fragment 文件，样本名示例：", 
        paste0(sample_names[1:min(5, length(sample_names))], collapse=","), " ...")

# -----------------------
# 3. ArchR 初始化
# -----------------------
addArchRThreads(threads = 20)
addArchRGenome("mm10")

# -----------------------
# 4. 创建 ArrowFiles
# -----------------------
dir.create(archr_outdir, recursive = TRUE, showWarnings = FALSE)
setwd(archr_outdir)

arrow_files <- paste0(sample_names, ".arrow")
files_exist <- all(file.exists(arrow_files))

if (!files_exist) {
    message(">>> 正在创建 ArrowFiles ...")
    ArrowFiles <- createArrowFiles(
        inputFiles       = fragment_paths,
        sampleNames      = sample_names,
        minTSS           = 4,
        minFrags         = 1000,
        addTileMat       = TRUE,
        addGeneScoreMat  = TRUE,
        geneAnnotation   = geneAnnoMm10,
        genomeAnnotation = genomeAnnoMm10
    )
} else {
    message(">>> ArrowFiles 已存在，跳过创建。")
    ArrowFiles <- file.path(archr_outdir, arrow_files)
}

# -----------------------
# 5. 创建或加载 ArchRProject
# -----------------------
proj_file <- file.path(archr_outdir, paste0(study_id, "_ArchRProject.rds"))
if (file.exists(proj_file)) {
    message(">>> ArchRProject 文件已存在，正在加载：", proj_file)
    proj <- readRDS(proj_file)
} else {
    message(">>> 创建 ArchRProject ...")
    proj <- ArchRProject(
        ArrowFiles      = ArrowFiles,
        copyArrows      = FALSE,
        outputDirectory = archr_outdir
    )
}
# -----------------------
# 6. 去除双细胞 (Doublets)
# -----------------------
doub_scores_file <- file.path(archr_outdir, paste0("doublet_scores_", study_id, ".rds"))
if (!file.exists(doub_scores_file)) {
    message(">>> 添加双细胞评分 (DoubletScores) ...")
    proj <- addDoubletScores(
        input      = proj,
        k          = 10,
        knnMethod  = "UMAP",
        LSIMethod  = 1
    )
} else {
    message(">>> DoubletScores 已存在，跳过计算。")
}

proj_filt_file <- file.path(archr_outdir, paste0("proj_filtered_", study_id, ".rds"))
if (!file.exists(proj_filt_file)) {
    message(">>> 进行双细胞过滤 (filterDoublets) ...")
    proj <- filterDoublets(
        ArchRProj   = proj,
        filterRatio = 1.5
    )
} else {
    message(">>> 已经完成双细胞过滤，直接读取项目。")
    proj <- readRDS(proj_filt_file)
}

# -----------------------
# 7. 绘制 FragmentSizes 和 QC 图
# -----------------------
qc_pdf <- file.path(archr_outdir, "QC_Fragment_TSS_nFrags.pdf")
if (!file.exists(qc_pdf)) {
    message(">>> 绘制 FragmentSizes 和 TSS/nFrags 分布图 ...")
    p1 <- plotFragmentSizes(ArchRProj = proj) +
    theme(text = element_text(size = 14))  # 调整字体大小
    p2 <- plotGroups(
        ArchRProj = proj,
        groupBy   = "Sample",
        colorBy   = "cellColData",
        name      = "TSSEnrichment",
        plotAs    = "ridges"
    ) + theme(text = element_text(size = 14))
    p3 <- plotGroups(
        ArchRProj = proj,
        groupBy    = "Sample",
        colorBy    = "cellColData",
        name       = "TSSEnrichment",
        plotAs     = "violin",
        alpha      = 0.4,
        addBoxPlot = TRUE
    ) + theme(text = element_text(size = 14))
    p4 <- plotGroups(
        ArchRProj = proj,
        groupBy   = "Sample",
        colorBy   = "cellColData",
        name      = "log10(nFrags)",
        plotAs    = "ridges"
    ) + theme(text = element_text(size = 14))
    p5 <- plotGroups(
        ArchRProj = proj,
        groupBy    = "Sample",
        colorBy    = "cellColData",
        name       = "log10(nFrags)",
        plotAs     = "violin",
        alpha      = 0.4,
        addBoxPlot = TRUE
    ) + theme(text = element_text(size = 14))
    pdf(qc_pdf, width = 10, height = 12)
    print(p1)
    print(p2+p3)
    print(p4+p5)
    print(p5)
    dev.off()
} else {
    message(">>> QC 图已存在，跳过绘制。")
}

# -----------------------
# 8. TileMatrix 上的 IterativeLSI
# -----------------------
res_tile_lsi <- 0.5
lsi_tile_name <- paste0("IterativeLSI_tile_res", "c")

message(">>> 进行 IterativeLSI（TileMatrix） ...")
proj <- addIterativeLSI(
    ArchRProj     = proj,
    useMatrix     = "TileMatrix",
    name          = lsi_tile_name,
    iterations    = 3,
    clusterParams = list(
        resolution  = c(0.2, 0.4),
        sampleCells = 10000,
        n.start     = 10
    ),
    varFeatures   = 25000,
    dimsToUse     = 1:30,
    force         = TRUE
)

# -----------------------
# 9. TileMatrix 聚类和 UMAP
# -----------------------
tile_cluster_name <- paste0("Clusters_tile_res", res_tile_lsi)
tile_umap_name    <- paste0("UMAP_tile_res", res_tile_lsi)

# 聚类
if (!(tile_cluster_name %in% names(proj@cellColData))) {
    message(">>> 进行聚类（TileMatrix, resolution=0.5） ...")
    proj <- addClusters(
        input        = proj,
        reducedDims  = lsi_tile_name,
        method       = "Seurat",
        name         = tile_cluster_name,
        resolution   = res_tile_lsi,
        force        = TRUE,
        nOutlier     = 25,
        maxClusters  = 60
    )
} else {
    message(">>> 聚类（TileMatrix, resolution=0.5）已存在，跳过。")
}

# UMAP
if (!(tile_umap_name %in% names(proj@embeddings))) {
    message(">>> 进行 UMAP（TileMatrix, resolution=0.5） ...")
    proj <- addUMAP(
        ArchRProj    = proj,
        reducedDims  = lsi_tile_name,
        name         = tile_umap_name,
        nNeighbors   = 15,
        minDist      = 0.8,
        metric       = "cosine",
        force        = TRUE
    )
} else {
    message(">>> UMAP（TileMatrix, resolution=0.5）已存在，跳过。")
}

resolutions <- c(0.1, 0.3, 0.5, 0.8, 1.0, 1.3, 1.5)
# macs2 可执行文件路径
macs2_path <- "/home/renshida/miniconda3/envs/scatac_r_410/bin/macs2"
for (res in resolutions) {
    subdir <- file.path(archr_outdir, paste0("peak_res_", res))
    final_proj_file <- file.path(
        subdir, paste0(study_id, "_proj_peak_umap_res_", res, ".rds")
    )
  if (file.exists(final_proj_file)) {
    message("[res=", res, "] 最终 .rds 已存在，跳过该分辨率。")
    next
  }
  # 11.1 为当前分辨率创建子目录
  subdir <- file.path(archr_outdir, paste0("peak_res_", res))
  if (!dir.exists(subdir)) {
    dir.create(subdir, recursive = TRUE)
  }
  message("\n\n>>> [res=", res,
          "] 开始 Peak calling 与后续分析，输出到：", subdir)

  # 11.2 克隆主 proj 为 proj_res
  cloned_proj_file <- file.path(
    subdir, paste0(study_id, "_archr_project_res_", res, ".rds")
  )
  if (!file.exists(cloned_proj_file)) {
    saveRDS(proj, file = cloned_proj_file)
  }
  proj_res <- readRDS(cloned_proj_file)

  # 11.3 计算 GroupCoverages（基于 clusters_tile_res{res}）
  group_cov_name <- paste0("GroupCoverages_tile_res", res)
  cluster_col    <- paste0("Clusters_tile_res", res)
  if (!(group_cov_name %in% names(proj_res@projectMetadata))) {
    message(
      "[res=", res,
      "] 计算 GroupCoverages（",
      cluster_col,
      "） …"
    )
    proj_res <- addGroupCoverages(
      ArchRProj = proj_res,
      groupBy   = "Clusters_tile_res0.5",  # 必须与 cellColData 中的列名完全一致
      force     = TRUE
    )
  } else {
    message("[res=", res, "] GroupCoverages 已存在，跳过。")
  }

  # 11.4 调用 MACS2 生成可复现 PeakSet
  peak_set_name <- paste0("PeakSet_tile_res", res)
  if (!(peak_set_name %in% names(proj_res@projectMetadata))) {
    message("[res=", res, "] 调用 MACS2 进行 Peak calling …")
    proj_res <- addReproduciblePeakSet(
      ArchRProj   = proj_res,
      groupBy     = "Clusters_tile_res0.5",
      pathToMacs2 = macs2_path,
      force       = TRUE
    )
  } else {
    message("[res=", res, "] PeakSet 已存在，跳过。")
  }

    # 11.5 构建 PeakMatrix
    message("[res=", res, "] 构建 PeakMatrix …")
    proj_res <- addPeakMatrix(proj_res)

  # 11.6 保存含有 PeakMatrix 的中间项目
  peak_proj_file <- file.path(
    subdir, paste0(study_id, "_proj_peak_res_", res, ".rds")
  )
  if (!file.exists(peak_proj_file)) {
    saveRDS(proj_res, file = peak_proj_file)
  }

  # 11.7 在 PeakMatrix 上做 IterativeLSI（注意首字母大写，完全匹配 ArchR 内部的命名）
  lsi_peak_name <- paste0("IterativeLSI_peak_res", res)
  message("[res=", res, "] 对 PeakMatrix 进行 IterativeLSI …")
  proj_res <- addIterativeLSI(
        ArchRProj     = proj_res,
        useMatrix     = "PeakMatrix",
        name          = lsi_peak_name,
        iterations    = 3,
        clusterParams = list(
        resolution  = res,
        sampleCells = NULL,
        n.start     = 10
        ),
        varFeatures   = 25000,
        dimsToUse     = 1:30,
        force         = TRUE
    )


  # 11.8 PeakMatrix 上做聚类
  clusters_peak_name <- paste0("Clusters_peak_res", res)
  if (!(clusters_peak_name %in% colnames(proj_res@cellColData))) {
    message("[res=", res, "] 对 PeakMatrix 做聚类（resolution=", res, "）…")
    proj_res <- addClusters(
      input       = proj_res,
      reducedDims = lsi_peak_name,
      method      = "Seurat",
      name        = clusters_peak_name,
      resolution  = res,
      force       = TRUE,
      nOutlier    = 25,
      maxClusters = 50
    )
  } else {
    message("[res=", res, "] ", clusters_peak_name, " 已存在，跳过。")
  }

  # 11.9 PeakMatrix 上做 UMAP
  umap_peak_name <- paste0("UMAP_peak_res", res)
  if (!(umap_peak_name %in% names(proj_res@embeddings))) {
    message("[res=", res, "] 对 PeakMatrix 做 UMAP（res=", res, "）…")
    proj_res <- addUMAP(
      ArchRProj   = proj_res,
      reducedDims = lsi_peak_name,
      name        = umap_peak_name,
      nNeighbors  = 15,
      minDist     = 0.5,
      metric      = "cosine",
      force       = TRUE
    )
  } else {
    message("[res=", res, "] ", umap_peak_name, " 已存在，跳过。")
  }

  # 11.10 保存最终包含 PeakMatrix LSI/聚类/UMAP 的项目
  final_proj_file <- file.path(
    subdir, paste0(study_id, "_proj_peak_umap_res_", res, ".rds")
  )
  if (!file.exists(final_proj_file)) {
    saveRDS(proj_res, file = final_proj_file)
  }

  # 11.11 可选：将每个分辨率的 PeakMatrix UMAP 图输出至 PDF（放在子目录里）
  umap_pdf_res <- file.path(subdir, paste0("umap_peak_res_", res, ".pdf"))
  if (!file.exists(umap_pdf_res)) {
    message("[res=", res, "] 绘制 PeakMatrix UMAP 图，输出到：", umap_pdf_res)
    pdf(umap_pdf_res, width = 8, height = 6)

    p1 <- plotEmbedding(
      ArchRProj = proj_res,
      colorBy   = "cellColData",
      name      = "Sample",
      embedding = umap_peak_name
    ) + ggtitle(paste0("UMAP Peak (res=", res, ") - Sample")) + theme(text = element_text(size = 14))

    p2 <- plotEmbedding(
      ArchRProj = proj_res,
      colorBy   = "cellColData",
      name      = clusters_peak_name,
      embedding = umap_peak_name
    ) + ggtitle(paste0("UMAP Peak (res=", res, ") - Cluster")) + theme(text = element_text(size = 14))

    print(p1)
    print(p2)
    dev.off()
  }
}

# -----------------------
# 12. 多分辨率 UMAP 总览（可选）
# -----------------------
umap_all_pdf <- file.path(
  archr_outdir, "umap_peak_multi_res_summary.pdf"
)
if (!file.exists(umap_all_pdf)) {
  message(">>> 绘制横跨所有分辨率的 Peak UMAP 总览，输出到：", umap_all_pdf)
  pdf(umap_all_pdf, width = 12, height = 4 * length(resolutions))

  for (res in resolutions) {
    proj_res <- readRDS(
      file.path(
        archr_outdir, paste0("peak_res_", res),
        paste0(study_id, "_proj_peak_umap_res_", res, ".rds")
      )
    )
    umap_peak_name    <- paste0("UMAP_peak_res", res)
    clusters_peak_name <- paste0("Clusters_peak_res", res)

    p_s <- plotEmbedding(
      ArchRProj = proj_res,
      colorBy   = "cellColData",
      name      = "Sample",
      embedding = umap_peak_name
    ) + ggtitle(paste0("UMAP PeakRes", res, " - Sample")) + theme(text = element_text(size = 14))

    p_c <- plotEmbedding(
      ArchRProj = proj_res,
      colorBy   = "cellColData",
      name      = clusters_peak_name,
      embedding = umap_peak_name
    ) + ggtitle(paste0("UMAP PeakRes", res, " - Cluster")) + theme(text = element_text(size = 14))

    print(p_s)
    print(p_c)
  }

  dev.off()
} else {
  message(">>> UMAP 多分辨率总览已存在，跳过。")
}
# -----------------------
# 13. 提取 Marker 基因并输出 Top10（针对每个分辨率）
# -----------------------
for (res in resolutions) {
    # 读取对应分辨率的项目
    proj_res <- readRDS(file.path(archr_outdir, paste0("peak_res_", res),
                                   paste0(study_id, "_proj_peak_umap_res_", res, ".rds")))
    clusters_peak <- paste0("Clusters_peak_res", res)
    message(paste0(">>> [res=", res, "] 提取 Marker 基因（PeakMatrix, resolution=", res, "） ..."))

    markersGS <- getMarkerFeatures(
        ArchRProj = proj_res,
        useMatrix = "GeneScoreMatrix",
        groupBy   = clusters_peak,
        bias      = c("TSSEnrichment", "log10(nFrags)"),
        testMethod= "wilcoxon"
    )
    markerList <- getMarkers(markersGS, cutOff = "FDR <= 0.05 & Log2FC >= 0.25")
    marker <- markerList@listData
    names(marker) <- as.character(seq_along(markerList))

    use <- data.frame(
        cluster = integer(),
        gene    = character(),
        Log2FC  = numeric(),
        stringsAsFactors = FALSE
    )
    colnames(use) <- c("cluster", "gene", "Log2FC")

    for (i in seq_along(markerList)) {
        if (nrow(marker[[i]]) > 0) {
            tmp <- data.frame(
                cluster = i,
                gene    = marker[[i]]$name,
                Log2FC  = marker[[i]]$Log2FC,
                stringsAsFactors = FALSE
            )
            tmp <- tmp[order(-tmp$Log2FC), ]
        } else {
            tmp <- data.frame(
                cluster = i,
                gene    = "null",
                Log2FC  = NA,
                stringsAsFactors = FALSE
            )
        }
        use <- rbind(use, tmp)
    }

    # 写出所有簇的 Marker 列表
    write.csv(use, file = file.path(archr_outdir, paste0("marker_res", res, ".csv")), row.names = FALSE)

    # 打印并保存 Top10
    top10_per_cluster <- use %>%
        group_by(cluster) %>%
        arrange(desc(Log2FC)) %>%
        slice_head(n = 10) %>%
        ungroup()
    print(top10_per_cluster)

    out <- top10_per_cluster %>%
        group_by(cluster) %>%
        summarise(
            genes = paste(gene, collapse = ","),
            .groups = "drop"
        )
    for (i in seq_len(nrow(out))) {
        cat("[res=", res, "] Cluster", out$cluster[i], ":\n")
        cat(out$genes[i], "\n\n")
    }
}

message("==== 研究 ", study_id, " 处理完成 ====")
