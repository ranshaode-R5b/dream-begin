#!/usr/bin/env bash
set -euo pipefail

# 🔧 配置路径
BED_DIR="/data1/renshida/scATAC/data/01_study_Archr/GSE162798/ArchR/BED_markerPeaks_by_CellType_hg_38"
MOTIF_ROOT="/data1/renshida/scATAC/data/03_motif/GSE162798/Motif_hommer"
OUT_DIR="/data1/renshida/scATAC/data/03_motif/GSE162798/Annotated_Peaks_full"
GENOME="hg38"

mkdir -p "${OUT_DIR}"

for BED in "${BED_DIR}"/*_markerPeaks.bed; do
  CT=$(basename "${BED}" "_markerPeaks.bed")
  MOTIF_DIR="${MOTIF_ROOT}/${CT}_homer"
  echo "🔍 正在处理：${CT}"

  ANNO_STATS="${OUT_DIR}/${CT}_annotation_stats.txt"
  ANNO_TSV="${OUT_DIR}/${CT}_peaks_annotated.tsv"
  MATRIX_PREFIX="${OUT_DIR}/${CT}_motif_matrix"
  MOTIF_FILES=("${MOTIF_DIR}/knownResults"/*.motif)

  annotatePeaks.pl \
    "${BED}" "${GENOME}" \
    -annStats "${ANNO_STATS}" \
    -m "${MOTIF_FILES[@]}" \
    -mscore \
    -nmotifs \
    -mbed "${OUT_DIR}/${CT}_motif_positions.bed" \
    -matrix "${MATRIX_PREFIX}" \
    > "${ANNO_TSV}"

  echo "✅ 完成 ${CT} 注释"
done

echo "🎯 全部 cell type 注释完成，结果保存在：${OUT_DIR}"
