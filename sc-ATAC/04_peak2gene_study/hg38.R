library(ArchR)
library(rtracklayer)
library(dplyr)
library(readr)
library(GenomicRanges)
library(BSgenome.Hsapiens.UCSC.hg38)
library(GenomeInfoDb)
library(stringr)
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

for (ds_name in names(datasets)) {
  ds <- datasets[[ds_name]]
  message("Processing dataset: ", ds_name)
  setwd(ds$path)
  proj <- readRDS(ds$projRds)

  # 提取分辨率，构造 reducedDims 名称
  res_match <- str_match(ds$projRds, "peak_res[_]?([0-9]+(?:\\.[0-9]+)?)")
  if (is.na(res_match[1,2])) {
    stop("无法从 projRds 路径中提取峰分辨率: ", ds$projRds)
  }
  resolution <- res_match[1,2]
  reduced_name <- paste0("IterativeLSI_peak_res", resolution)
  message("  → Using reducedDims: ", reduced_name)

  marker_dir <- file.path(proj@projectMetadata$outputDirectory, ds$markerSubdir)
  out_base <- file.path("/data1/renshida/scATAC/data/04_p2g", ds_name, "ArchR_mm")
  dir.create(out_base, recursive = TRUE, showWarnings = FALSE)

  bed_files <- list.files(marker_dir, pattern = "_markerPeaks.bed$", full.names = TRUE)
  for (bed in bed_files) {
    ct <- sub("_markerPeaks\\.bed$", "", basename(bed))
    ct_clean <- gsub("_", " ", ct)
    message("  ➤ CellType: '", ct_clean, "' from file: '", ct, "'")

    cells_ct <- proj$cellNames[proj$CellType == ct_clean]
    message("    → cells_ct length: ", length(cells_ct))
    if (length(cells_ct) == 0) {
      warning("   ⚠️ no cells for ", ct, "; skipped")
      next
    }

    gr_marker <- import(bed)

    # 跳过空 BED
    if (length(gr_marker) == 0) {
      message("    ⚠️ bed empty for ", ct, "; skipped")
      next
    }

    # 修复染色体信息
    seqlevelsStyle(gr_marker) <- "UCSC"
    gr_marker <- keepStandardChromosomes(gr_marker, pruning.mode = "coarse")
    gr_marker <- sortSeqlevels(gr_marker)

    subproj_dir <- file.path(out_base, paste0("subproj_", ct))
    unlink(subproj_dir, recursive = TRUE, force = TRUE)
    dir.create(subproj_dir, recursive = TRUE)

    subProj <- subsetArchRProject(
      ArchRProj = proj,
      cells = cells_ct,
      outputDirectory = subproj_dir,
      dropCells = TRUE,
      force = TRUE
    )

    subProj <- addPeakSet(
      subProj,
      peakSet = gr_marker,
      genomeAnnotation = getGenomeAnnotation(subProj),
      force = TRUE
    ) %>% addPeakMatrix()

    n_cells <- length(cells_ct)
    k_use <- min(100, floor(length(cells_ct) * 0.1))

  if (n_cells < 30) {
  message("⚠️ Skipping ", ct, ": too few cells (", n_cells, ").")
  } else {
  subProj <- addPeak2GeneLinks(
    ArchRProj = subProj,
    reducedDims = reduced_name,
    useMatrix   = "GeneScoreMatrix",
    corCutOff   = 0.45,
    maxDist     = 250000,
    k           = k_use,
    addPermutedPval = FALSE
  )
  p2g <- getPeak2GeneLinks(subProj, returnLoops = FALSE)
  # 输出逻辑省略
  }

    if (nrow(p2g) == 0) {
      message("     no peak–gene links for ", ct, "; skipped CSV")
      unlink(subproj_dir, recursive = TRUE, force = TRUE)
      next
    }

    pset <- metadata(p2g)$peakSet
    gset <- metadata(p2g)$geneSet
    p2g$geneName <- mcols(gset)$name[p2g$idxRNA]

    df <- as.data.frame(p2g) %>%
      mutate(
        Chr = as.character(seqnames(pset)[idxATAC]),
        Start = start(pset)[idxATAC] - 1,
        End = end(pset)[idxATAC]
      ) %>%
      select(Chr, Start, End, geneName, Correlation, FDR, VarQATAC, VarQRNA)

    out_csv <- file.path(out_base, paste0("p2g_", ct, "_links.csv"))
    write_csv(df, out_csv)
    message("    -----------------------------------wrote ", nrow(df), " links → ", out_csv)

    unlink(subproj_dir, recursive = TRUE, force = TRUE)
    message("    ----------------------------------removed temp folder")
  }
}