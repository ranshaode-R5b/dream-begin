#!/usr/bin/env bash
set -euo pipefail

DATASETS=("mouseneo" "mousetissues" "pairedTag")
BASE="/data1/renshida/scATAC/data/00_rawdata"

for ds in "${DATASETS[@]}"; do
  echo "=== 处理 dataset: $ds ==="
  BIGWIG_DIR="$BASE/$ds/bigwig"
  CCRE_DIR="$BASE/$ds/cCREs_auto"

  [ -d "$BIGWIG_DIR" ] || { echo "⚠️ 目录不存在：$BIGWIG_DIR，跳过"; continue; }
  mkdir -p "$CCRE_DIR"
  cd "$BIGWIG_DIR"

  # 获取前三个 bigWig 文件
  mapfile -t files < <(ls *.bw 2>/dev/null | head -n3)
  if (( ${#files[@]} == 0 )); then
    echo "❌ 未找到 .bw 文件，跳过"
    continue
  fi

  cutoffs=()
  for bw in "${files[@]}"; do
    base="${bw%.bw}"
    echo "➤ 计算 $bw 的 cutoff"

    # 转 bedGraph 并排序
    bigWigToBedGraph "$bw" temp.bedGraph
    sort --parallel=2 -k1,1 -k2,2n temp.bedGraph > temp.sorted.bedGraph

    # 生成 cutoff 统计
    stats="$CCRE_DIR/${ds}_${base}_cutoff_stats.txt"
    macs3 bdgpeakcall \
      -i temp.sorted.bedGraph \
      --min-length 200 --max-gap 30 \
      --cutoff-analysis \
      --cutoff-analysis-max 50 \
      --cutoff-analysis-steps 25 \
      -o "$stats"

    # Python 脚本：尾部剪裁 + 拐点检测

    cutoff=$(python3 <<EOF
import numpy as np, os
import matplotlib.pyplot as plt
from kneed import KneeLocator

stats = "$stats"
assert os.path.exists(stats), f"文件不存在: {stats}"

data = np.loadtxt(stats, skiprows=1)
x_all, y_all = data[:,0], data[:,1]

# 检测 abrupt "cliff" 片段
sec = y_all[2:] + y_all[:-2] - 2 * y_all[1:-1]
cliff_mask = np.concatenate(([True], sec < np.percentile(sec, 90), [True]))
peak_mask = y_all > 1000
mask = cliff_mask & peak_mask
x, y = x_all[mask], y_all[mask]

if len(x) < 3:
    x, y = x_all, y_all

# 拐点检测
kl = KneeLocator(
    x, y,
    curve="convex", direction="increasing",
    S=5.0, interp_method="polynomial", online=True
)

cutoff = kl.knee if kl.knee is not None else x[np.argmax(np.abs(y[2:] + y[:-2] - 2 * y[1:-1]))]
print(cutoff)

basename = os.path.basename(stats).rsplit("_cutoff_stats.txt", 1)[0]
# 带上后缀 .png，确保唯一
filename = f"{basename}_cutoff_plot.png"
save_path = os.path.join(os.path.dirname(stats), filename)
# 绘图
plt.figure(figsize=(8,6))
kl.plot_knee(figsize=(8,6), title="Cutoff Knee Plot", xlabel="Score", ylabel="Peak Count")
plt.savefig(save_path, dpi=150)
plt.close()
EOF
)

    echo "    => 得到 cutoff = $cutoff"
    cutoffs+=("$cutoff")

    rm temp.bedGraph temp.sorted.bedGraph
  done

  # 计算平均 cutoff
  avg_cutoff=$(printf '%s\n' "${cutoffs[@]}" | awk '{sum+=$1; cnt+=1} END{printf "%.2f", sum/cnt}')
  echo "✳️ $ds 平均 cutoff = $avg_cutoff"

  # 按平均 cutoff 统一调用所有 bigWig 文件
  for bw in *.bw; do
    base="${bw%.bw}"
    echo "♦️ 处理 $bw"
    bigWigToBedGraph "$bw" "${base}.bedGraph"
    sort -k1,1 -k2,2n "${base}.bedGraph" > "${base}.sorted.bedGraph"

    macs3 bdgpeakcall \
      -i "${base}.sorted.bedGraph" \
      --cutoff "$avg_cutoff" \
      --min-length 200 --max-gap 30 \
      -o "$CCRE_DIR/${base}.bed"

    rm "${base}.bedGraph" "${base}.sorted.bedGraph"
  done

  echo "✅ 完成 dataset：$ds 的所有样本调用峰"
done

echo "🎉 所有 dataset 完成！"
