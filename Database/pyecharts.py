import pandas as pd
import numpy as np
from pyecharts.charts import Scatter
from pyecharts import options as opts
from pyecharts.commons.utils import JsCode

# 染色体长度定义
chrom_len = {
    1: 249250621, 2: 243199373, 3: 198022430, 4: 191154276,
    5: 180915260, 6: 171115067, 7: 159138663, 8: 146364022,
    9: 141213431, 10: 135534747, 11: 135006516, 12: 133851895,
    13: 115169878, 14: 107349540, 15: 102531392, 16: 90354753,
    17: 81195210, 18: 78077248, 19: 59128983, 20: 63025520,
    21: 48129895, 22: 51304566, 23: 156040895, 24: 57227415
}

# 计算染色体起始位置
chrom_start = {}
cumulative = 0
for chr_num in sorted(chrom_len.keys()):
    chrom_start[chr_num] = cumulative
    cumulative += chrom_len[chr_num]

# 读取数据
data = pd.read_csv("/home/renshida/eqtm/data/manhattan/data/31791327-cis-LLS.txt", sep="\s+")

# 处理染色体编号
def parse_chromosome(chr_str):
    chr_type = chr_str[3:]
    if chr_type == "X": return 23
    if chr_type == "Y": return 24
    return int(chr_type)

data["chr_num"] = data["C_Chr"].apply(parse_chromosome)

# 计算基因组位置和-log10(p)
data["x"] = data.apply(lambda row: chrom_start[row["chr_num"]] + row["C_Pos"], axis=1)
data["y"] = -np.log10(data["P.value"])

# 颜色配置
colors = [
    "#40E0D0", "#5ED5B6", "#7CCB9C", "#9AC082", "#B8B668", "#D6AB4E",
    "#F4A134", "#F69A2A", "#F79320", "#F98C16", "#FA850C", "#FB7E02",
    "#FA6E1D", "#F95E38", "#F84E53", "#F73E6E", "#F62E89", "#F51EA4",
    "#F40EBF", "#E40FC1", "#D410C3", "#C411C5", "#B412C7", "#A413C9"
]
data["color"] = data["chr_num"].apply(lambda x: colors[(x-1)%24])

# 计算染色体刻度和标签
chrom_ticks = []
chrom_labels = []
cumulative = 0
for chr_num in sorted(chrom_len.keys()):
    if chr_num > 22: continue  # 可选是否显示XY染色体
    end = cumulative + chrom_len[chr_num]
    chrom_ticks.append(end)
    chrom_labels.append(str(chr_num))
    cumulative += chrom_len[chr_num]



# 创建散点图
scatter = Scatter(init_opts=opts.InitOpts(width="100%", height="500px"))

# 自定义x轴标签格式
xaxis_formatter = JsCode(f"""
function (value) {{
    var ticks = {chrom_ticks};
    var labels = {chrom_labels};
    for (var i = 0; i < ticks.length; i++) {{
        if (Math.abs(value - ticks[i]) < 1500000) {{
            return labels[i];
        }}
    }}
    return '';
}}
""")

# 配置图表选项
scatter.set_global_opts(
    title_opts=opts.TitleOpts(title="Manhattan Plot"),
    xaxis_opts=opts.AxisOpts(
        name="Chromosome",
        type_="value",
        splitline_opts=opts.SplitLineOpts(
            is_show=False
        ),
       axistick_opts=opts.AxisTickOpts(is_show=False),
        axislabel_opts=opts.LabelOpts(
            formatter=xaxis_formatter,
            rotate=0,
            interval = 0,
        ),
        min_=0,
        max_=cumulative,
        split_number=1000,
    ),
    yaxis_opts=opts.AxisOpts(
        name="-log10(p value)",
        type_="value",
        splitline_opts=opts.SplitLineOpts(is_show=False)
    ),
    tooltip_opts=opts.TooltipOpts(
        formatter=JsCode(
            "function (params) {"
            "var data = params.data;"
            "return 'Chromosome: ' + params.value[2] + '<br/>'"
            "+ 'Position: ' + params.value[3] + '<br/>'"
            "+ 'CpG: ' + params.value[5] + '<br/>'"
            "+ 'Gene: ' + params.value[6] + '<br/>'"
            "+ '-log10(p): ' + params.value[1].toFixed(2);}"
        )
    )
)

significant_data = data[data["y"] > -np.log10(0.05)]
non_significant_data = data[data["y"] <= -np.log10(0.05)]

# 添加显著数据系列
scatter.add_xaxis(significant_data["x"].tolist())
scatter.add_yaxis(
    series_name="Significant",
    y_axis=[list(z) for z in zip(significant_data["y"], significant_data["chr_num"], significant_data["C_Pos"], significant_data["color"], significant_data["CpG"], significant_data["Gene"])],
    symbol_size=5,
    symbol="circle",
    itemstyle_opts=opts.ItemStyleOpts(
        color=JsCode("function(params) { return params.data[4]; }")
    ),
    label_opts=opts.LabelOpts(is_show=False),
)

# 添加不显著数据系列
scatter.add_xaxis(non_significant_data["x"].tolist())
scatter.add_yaxis(
    series_name="Non-Significant",
    y_axis=[list(z) for z in zip(non_significant_data["y"], non_significant_data["chr_num"], non_significant_data["C_Pos"], non_significant_data["color"], non_significant_data["CpG"], non_significant_data["Gene"])],
    symbol_size=5,
    symbol="triangle",
    itemstyle_opts=opts.ItemStyleOpts(
        color=JsCode("function(params) { return params.data[4]; }")
    ),
    label_opts=opts.LabelOpts(is_show=False),
)
# 添加显著性线
scatter.set_series_opts(
    markline_opts=opts.MarkLineOpts(
        data=[opts.MarkLineItem(y=-np.log10(0.05))],
        label_opts=opts.LabelOpts(position="end", formatter="0.05"),
        linestyle_opts=opts.LineStyleOpts(color="red", type_="dashed")
    )
)

# 生成HTML文件
scatter.render("/home/renshida/eqtm/data/manhattan/data/data-manhattan/31791327-LLS.html")