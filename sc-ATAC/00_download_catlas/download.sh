#!/usr/bin/env bash

BASE_URL="https://decoder-genetics.wustl.edu/data/catlas/"

# 本地保存目录
DEST_DIR="/data1/renshida/scATAC/data/00_rawdata"
dirs=(
  Hic
)

for d in "${dirs[@]}"; do
  echo "==> 下载目录 ${d} （跳过已存在文件）…"
  wget -r -c -np -nH --cut-dirs=2 -R "index.html*" \
     -P "$DEST_DIR" "${BASE_URL}/${d}/"
done