#!/bin/bash

# 用户参数
IN_DIR="/data1/renshida/scATAC/data/00_rawdata/humanheart/cCREs"
OUT_DIR="/data1/renshida/scATAC/data/00_rawdata/humanheart/cCREs_merge"

# 检查目录
if [ ! -d "$IN_DIR" ]; then
  echo "Error: 输入目录 $IN_DIR 不存在" >&2
  exit 1
fi

mkdir -p "$OUT_DIR"

# 主循环
for file in "$IN_DIR"/*.bed; do
  base=$(basename "$file" .bed)
  out="$OUT_DIR/${base}.bed"
  echo "Processing $file → $out"

  sort -k1,1 -k2,2n "$file" \
    | bedtools merge -i - \
    > "$out"
done

echo "All BED files in $IN_DIR have been merged into $OUT_DIR"
