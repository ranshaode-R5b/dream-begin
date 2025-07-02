#!/usr/bin/env bash

DATASETS=("humanneo" "macaqueneo" "humanbraindev" "humanfetal" "marmosetneo" "mouseneo" "mousetissues" "pairedTag" "humanbrain")  # 根据需要添加多个子目录
BASE="/data1/renshida/scATAC/data/00_rawdata"

for ds in "${DATASETS[@]}"; do
  dir="$BASE/$ds/cCREs"
  for file in "$dir"/*; do
    filename="$(basename "$file")"
    prefix="${filename%.*}"
    out="${BASE}/$ds/cCREs_awk"
    mkdir -p "$out"
    out_file="${out}/${prefix}.bed"
    # 使用 tail 跳过首行，再 cut 提取指定列
    tail -n +2 "$file" | cut -f1,2,3,5 > "$out_file"

    echo "提取完成（跳过首行）：$file -> $out_file"
  done
done
