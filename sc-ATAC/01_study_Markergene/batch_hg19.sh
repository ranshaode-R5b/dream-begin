#!/bin/bash

RAW_BASE_DIR="/data1/renshida/scATAC/data/00_rawdata_study"
ARCHR_BASE_DIR="/data1/renshida/scATAC/data/01_study_Archr"
SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
R_SCRIPT="${SCRIPT_DIR}/run_archr_hg19.R"

# 要处理的 StudyID 列表
#round test
#round1
#STUDIES=("GSE272504" "GSE162798" "GSE214085_IPA" "GSE214085_Con")
#round2
STUDIES=("GSE162798")
for STUDY in "${STUDIES[@]}"; do
    RAW_STUDY_DIR="${RAW_BASE_DIR}/${STUDY}"
    DATA_CYCLE_DIR="${RAW_STUDY_DIR}/data_cycle"

    if [[ ! -d "${DATA_CYCLE_DIR}" ]]; then
        echo "跳过：${STUDY} 没有 data_cycle 子目录"
        continue
    fi

    echo "=== 处理研究 ${STUDY} ==="
    TMP_LIST="/tmp/${STUDY}_fragments.list"

    # 在 data_cycle 下搜 *.fragments.tsv.gz，写入临时文件
    find "${DATA_CYCLE_DIR}" -type f -name "*fragments.tsv.gz" > "${TMP_LIST}"

    if [[ ! -s "${TMP_LIST}" ]]; then
        echo "未找到任何 fragments.tsv.gz，跳过 ${STUDY}"
        rm -f "${TMP_LIST}"
        continue
    fi

    # 调用 R 脚本：--base_dir 指向 ARCHR_BASE_DIR
    Rscript "${R_SCRIPT}" \
        --study_id "${STUDY}" \
        --base_dir "${ARCHR_BASE_DIR}" \
        --fragments "${TMP_LIST}"

    rm -f "${TMP_LIST}"
done

echo "全部研究处理完成"
