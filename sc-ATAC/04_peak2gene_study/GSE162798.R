setwd("/data1/renshida/scATAC/data/01_study_Archr/GSE162798/ArchR")

library(ArchR)
library(rtracklayer)
library(dplyr)
library(readr)
library(ggplot2)
library(GenomicRanges)
library(tibble)
library(BSgenome.Hsapiens.UCSC.hg19)
addArchRThreads(threads = 20)
addArchRGenome("hg19")

# 读取含 CellType 注释的 ArchRProject
proj <- readRDS("peak_res_0.3/GSE162798_peak_celltype.rds")

# 获取 marker peaks BED 路径列表
marker_dir1 <- file.path(proj@projectMetadata$outputDirectory, "BED_markerPeaks_by_CellType")
marker_dir  <- file.path("/data1/renshida/scATAC/data/04_co_2gene", "GSE162798")

# 确保输出目录存在
if (!dir.exists(marker_dir)) dir.create(marker_dir, recursive = TRUE)

# 列出所有 cell type 的 marker BED 文件
bed_files <- list.files(marker_dir1, pattern = "_markerPeaks.bed$", full.names = TRUE)

# 只取第一个细胞类型进行测试
bed_file <- bed_files[1]
ct <- gsub("_markerPeaks.bed$", "", basename(bed_file))
message("Processing TEST cell type: ", ct)

# 导入 marker peaks
gr_marker <- import(bed_file)

# 获取该 cell type 的细胞条码
cells_ct <- proj$cellNames[proj$CellType == ct]
if (length(cells_ct) == 0) {
  stop("No cells found for ", ct, "; stopping test.")
}

# 指定子项目目录
subdir <- file.path(marker_dir, paste0("subproj_", ct))
if (dir.exists(subdir)) unlink(subdir, recursive = TRUE)
dir.create(subdir, recursive = TRUE)

# 创建子项目
subProj <- subsetArchRProject(
  ArchRProj = proj,
  cells = cells_ct,
  outputDirectory = subdir,
  dropCells = TRUE,
  force = TRUE
)

# 添加特异 peak 并生成 PeakMatrix
subProj <- addPeakSet(subProj, peakSet = gr_marker,
                      genomeAnnotation = getGenomeAnnotation(subProj),
                      force = TRUE)
subProj <- addPeakMatrix(subProj)



# 2. Peak-to-Gene Linkage 分析
subProj <- addPeak2GeneLinks(
  ArchRProj = subProj,
  reducedDims = "IterativeLSI_peak_res0.3",
  useMatrix = "GeneScoreMatrix",
  corCutOff = 0.45,
  maxDist = 250000,
  addPermutedPval = FALSE
)
p2g_df <- getPeak2GeneLinks(subProj, returnLoops = FALSE)

# 提取 gene/peak 名称
pset <- metadata(p2g_df)$peakSet
gset <- metadata(p2g_df)$geneSet
p2g_df$geneName <- mcols(gset)$name[p2g_df$idxRNA]
p2g_df$peakName <- paste(
  seqnames(pset)[p2g_df$idxATAC],
  start(pset)[p2g_df$idxATAC],
  end(pset)[p2g_df$idxATAC],
  sep = "_"
)

# 创建一个只包含 Chr, Start, End, geneName, Correlation, FDR, VarQATAC, VarQRNA 的新表
# 转换为标准 data.frame
df <- as.data.frame(p2g_df)

# 添加 Chr、Start、End，并调整 Start 为 0-based
pset <- metadata(p2g_df)$peakSet
df2 <- df %>%
  mutate(
    Chr   = as.character(seqnames(pset)[idxATAC]),
    Start = start(pset)[idxATAC] - 1,
    End   = end(pset)[idxATAC]
  ) %>%
  select(
    Chr, Start, End,
    geneName,
    Correlation,
    FDR,
    VarQATAC,
    VarQRNA
  )

# 查看并保存
head(df2)
write.csv(df2, file = file.path(subdir, paste0("p2g_", ct, "_formatted.csv")), row.names = FALSE)
message("Processing completed for cell type: ", ct)

