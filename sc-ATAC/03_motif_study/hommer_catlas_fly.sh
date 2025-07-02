#!/usr/bin/env bash
set -euo pipefail

# 配置区

DATASETS=("flybrain") 
ROOT_DIR="/data1/renshida/scATAC/data/00_rawdata"
OUT_ROOT="/data1/renshida/scATAC/data/03_motif"
GENOME="dm6"
THREADS=16

for DS in "${DATASETS[@]}"; do
  echo -e "\n 开始处理数据集：$DS"

  INPUT_DIR="${ROOT_DIR}/${DS}/cCREs"
  OUTPUT_BASE="${OUT_ROOT}/${DS}/Motif_homer"
  BACKGROUND_DIR="${OUTPUT_BASE}/background"
  BACKGROUND="${BACKGROUND_DIR}/all_peaks.bed"

  if [[ ! -d "${INPUT_DIR}" ]]; then
    echo " 输入目录不存在：${INPUT_DIR}, 跳过 $DS"
    continue
  fi

  mkdir -p "${BACKGROUND_DIR}"

  echo " 生成背景 peaks 文件：${BACKGROUND}"
  cat "${INPUT_DIR}"/*.bed \
    | sort -k1,1 -k2,2n \
    | bedtools merge -i - \
    > "${BACKGROUND}"

  for BED in "${INPUT_DIR}"/*.bed; do
    CT="$(basename "${BED}" ".bed")"
    OUT_DIR="${OUTPUT_BASE}/${CT}_homer"
    mkdir -p "${OUT_DIR}"

    echo " [${DS}] 正在处理 CellType：${CT}"
    findMotifsGenome.pl \
      "${BED}" \
      "${GENOME}" \
      "${OUT_DIR}" \
      -bg "${BACKGROUND}" \
      -size given \
      -p "${THREADS}"

    echo " [${DS}/${CT}] 完成，结果在：${OUT_DIR}"
  done

done

echo -e "\n 所有指定数据集 motif 富集任务完成！"


