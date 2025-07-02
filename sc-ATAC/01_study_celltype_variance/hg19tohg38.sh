#!/usr/bin/env bash
set -euo pipefail

# ====== 配置区域，根据环境修改 ======
BED_DIR="/data1/renshida/scATAC/data/01_study_Archr/GSE214085_Con/ArchR/BED_markerPeaks_by_CellType"
CHAIN="/data1/renshida/scATAC/data/01_study_Archr/hg19ToHg38.over.chain.gz"   
LIFTOVER="/home/renshida/miniconda3/envs/scatac_r_410/bin/liftOver"                # UCSC liftOver 
OUT_DIR="/data1/renshida/scATAC/data/02_study_peaks/GSE214085_Con/ArchR/BED_markerPeaks_by_CellType/"
OUT_DIR_unmap="/data1/renshida/scATAC/data/02_study_peaks/GSE214085_Con/ArchR/unmap"
mkdir -p "$OUT_DIR"
mkdir -p "$OUT_DIR_unmap"
echo " Processing all .bed files in $BED_DIR"

for bed in "$BED_DIR"/*.bed; do
  name=$(basename "$bed" .bed)
  out_bed="${OUT_DIR}/${name}.bed"
  unmap="${OUT_DIR_unmap}/${name}_unmapped.bed"

  echo "name.bed → $name".bed

  # 使用 UCSC liftOver 保留原始6列以上信息
  "$LIFTOVER" -bedPlus=6 "$bed" "$CHAIN" "$out_bed" "$unmap"

done

echo "All done! Results in: $OUT_DIR"
