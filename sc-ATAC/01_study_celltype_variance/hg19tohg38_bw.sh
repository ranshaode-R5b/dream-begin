#!/usr/bin/env bash
set -euo pipefail

IN_DIR="/data1/renshida/scATAC/data/01_study_Archr/GSE214085_IPA/ArchR/GroupBigWigs/CellType"
CHAIN="/data1/renshida/scATAC/data/01_study_Archr/hg19ToHg38.over.chain.gz"
OUT_DIR="/data1/renshida/scATAC/data/02_study_peaks/GSE214085_IPA/ArchR/bigwig"
mkdir -p "$OUT_DIR"

for bw in "$IN_DIR"/*.bw; do
  name=$(basename "$bw" .bw)
  out="${OUT_DIR}/${name}.bw"
  echo " Converting $name.bw → $name.bw"
  CrossMap bigwig "$CHAIN" "$bw" "$out"
  echo " Done: $out"
done

