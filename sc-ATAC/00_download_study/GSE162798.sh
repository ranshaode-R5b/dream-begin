#!/usr/bin/env bash

BASE_URL="ftp://ftp.ncbi.nlm.nih.gov/geo/series/GSE162nnn/GSE162798/"

# 本地保存目录
DEST_DIR="/data1/renshida/scATAC/data/00_rawdata_study"
dirs=(
  matrix
  miniml
  soft
  suppl
)

for d in "${dirs[@]}"; do
  echo "==> 下载目录 ${d} ..."
  wget -r -c -np -nH --cut-dirs=3 -R "index.html*" --no-parent \
     -P "$DEST_DIR" "${BASE_URL}/${d}/"
done