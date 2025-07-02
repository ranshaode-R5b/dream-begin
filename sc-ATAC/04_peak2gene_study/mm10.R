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