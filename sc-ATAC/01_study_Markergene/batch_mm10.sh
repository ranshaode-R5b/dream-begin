RAW_BASE_DIR="/data1/renshida/scATAC/data/00_rawdata_study"
ARCHR_BASE_DIR="/data1/renshida/scATAC/data/01_study_Archr"
SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
R_SCRIPT="${SCRIPT_DIR}/run_archr_mm10.R"

"STUDIES=(
        "GSE162798" 
        "GSE174226" 
        "GSE231708_mm_Liver" 
        "GSE231708_mm_Lung" 
        "GSE245957" 
        "GSE242466_mm_COL" 
        "GSE242466_mm_SKI" 
        "GSE242466_mm_SPL" 
        "GSE242466_mm_VAT" 
        "E-MTAB-11264" 
        "GSE164439" 
        "GSE149622_IL23" 
        "GSE149622_naive" 
        "GSE211077_Ascl" 
        "GSE211077_IPA" 
        "GSE211077_Retinal"
        "GSE139950"
        "GSE214082_CD"
        "GSE214082_HFD"
        )"
#STUDIES=("GSE245957" "GSE231708_mm_Liver" "GSE231708_mm_Lung")
STUDIES=("GSE164439" "GSE231708_mm_Lung")
for STUDY in "${STUDIES[@]}"; do
    RAW_STUDY_DIR="${RAW_BASE_DIR}/${STUDY}"
    DATA_CYCLE_DIR="${RAW_STUDY_DIR}/data_cycle_mm"

    if [[ ! -d "${DATA_CYCLE_DIR}" ]]; then
        echo "跳过：${STUDY} 没有 data_cycle_mm 子目录"
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


