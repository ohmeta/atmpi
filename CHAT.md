# ATMPI 下游分析完整报告

**分析日期**: 2026 年 3 月 20 日  
**项目**: ASV-to-MAG 整合分析 (ATMPI)

---

## 📋 目录

1. [流程概述](#流程概述)
2. [运行结果汇总](#运行结果汇总)
3. [下游分析发现](#下游分析发现)
4. [Python 下游分析脚本](#python-下游分析脚本)
5. [R 下游分析脚本](#r-下游分析脚本)
6. [输出文件说明](#输出文件说明)

---

## 流程概述

### ATMPI 工作流程

ATMPI (ASV-to-MAG Integration Pipeline) 是一个用于将 16S rRNA 扩增子序列变异 (ASV) 与宏基因组组装基因组 (MAG) 进行匹配的生物信息学流程。

**核心步骤：**

| 步骤 | 规则 | 工具 | 功能 |
|------|------|------|------|
| 0 | `validate_inputs` | Python | 验证输入文件格式 |
| 1 | `extract_16s_from_mags` | Barrnap | 从 MAG 中提取 16S rRNA 基因 |
| 2 | `build_mag_16s_db` | makeblastdb | 构建 BLAST 数据库 |
| 3 | `blast_asv_to_mag_16s` | blastn | ASV 与预测的 16S 基因比对 |
| 4 | `blast_asv_to_mag_contigs` | blastn | ASV 与 MAG 全序列比对 |
| 5 | `process_blast_results` | Python | 处理 BLAST 结果并分配 |
| 6 | `filter_target_species` | Python | 过滤目标物种 |
| 7 | `generate_summary` | Python | 生成统计摘要和 HTML 报告 |

### 运行命令

```bash
# 激活环境
mamba activate env-atmpi

# 运行流程
snakemake --cores 16 --configfile config.yaml
```

---

## 运行结果汇总

### 整体统计

| 指标 | 数值 |
|------|------|
| 总 ASV 数 | 1,923 |
| 总 MAG 数 | 1,316 |
| ASV-MAG 配对数 | 537 |
| 有匹配的 ASV | 537 (27.9%) |
| 有匹配的 MAG | 61 (4.6%) |
| 平均匹配身份 | 98.86% |
| 匹配身份范围 | 97.0% - 100.0% |
| 具有菌株多样性的 MAG | 50 个 (82.0%) |

### 菌门水平分布

| 菌门 | ASV 数量 | MAG 数量 | 平均身份 |
|------|---------|---------|---------|
| p__Bacillota | 305 | 33 | 98.98% |
| p__Pseudomonadota | 96 | 6 | 98.28% |
| p__Actinomycetota | 66 | 13 | 99.00% |
| p__Bacteroidota | 58 | 6 | 98.93% |
| p__Verrucomicrobiota | 11 | 2 | 99.29% |
| p__un_k__Bacteria | 1 | 1 | 100.00% |

---

## 下游分析发现

### 1. 高置信度配对 (Identity ≥ 99%)

- **总数**: 304 条配对
- **占比**: 56.6% 的所有配对

**Top 10 高置信度 MAG:**

| MAG ID | ASV 数量 | 平均身份 | 物种 |
|--------|---------|---------|------|
| HOLA.SRR23604321.bin.34.fa.gz | 24 | 99.28% | Staphylococcus petrasii |
| HOLA.SRR23604292.bin.11.fa.gz | 21 | 99.30% | Veillonella parvula |
| HOLA.SRR23604272.bin.18.fa.gz | 18 | 99.66% | Streptococcus oralis |
| HOLA.SRR23604277.bin.14.fa.gz | 18 | 99.74% | Escherichia/Shigella coli |
| HOLA.SRR23604287.bin.63.fa.gz | 15 | 99.30% | Enterococcus avium |
| HOLA.SRR23604316.bin.15.fa.gz | 13 | 99.66% | Blautia wexlerae |
| HOLA.SRR23604311.bin.27.fa.gz | 13 | 99.70% | un_f__Lachnospiraceae |
| HOLA.SRR23604272.bin.12.fa.gz | 10 | 99.75% | Clostridium_sensu_stricto |
| HOLA.SRR23604288.bin.24.fa.gz | 9 | 99.64% | un_f__Lachnospiraceae |
| HOLA.SRR23604294.bin.18.fa.gz | 9 | 99.67% | Hungatella effluvii |

### 2. 菌株多样性分析

**关键发现：**
- 50 个 MAG (82.0%) 显示多 ASV 匹配
- 提示这些 MAG 代表包含多个菌株的种群（微多样性）

**菌株多样性 Top 5 MAG:**

| MAG ID | ASV 数量 | 平均身份 | 属 |
|--------|---------|---------|-----|
| HOLA.SRR23604299.bin.24.fa.gz | 59 | 97.60% | Enterobacter |
| HOLA.SRR23604321.bin.34.fa.gz | 35 | 99.08% | Staphylococcus |
| HOLA.SRR23604328.bin.21.fa.gz | 29 | 98.64% | Parabacteroides |
| HOLA.SRR23604272.bin.18.fa.gz | 29 | 99.11% | Streptococcus |
| HOLA.SRR23604292.bin.11.fa.gz | 28 | 98.99% | Veillonella |

### 3. 生物学洞见

1. **匹配质量优秀**: 平均 98.86% 的身份表明 ASV 和 MAG 之间有良好的对应关系
2. **丰富的菌株多样性**: 82% 的 MAG 显示多 ASV 匹配，提示样本中存在显著的微多样性
3. **优势菌门**: Bacillota 占主导（305 ASVs, 33 MAGs），符合肠道微生物组的典型组成
4. **高置信度配对**: 56.6% 的配对达到 99% 以上身份，可用于精确的菌株追踪

### 4. 后续研究方向

1. 对高置信度 MAG 进行功能注释，链接分类与功能
2. 分析具有多 ASV 的 MAG 的微多样性模式
3. 调查分类不一致的配对，可能发现新的分类单元
4. 整合样本元数据进行关联分析

---

## Python 下游分析脚本

**文件**: `downstream/creative_downstream_analysis.py`

**依赖安装**:
```bash
# 使用 pip
pip install pandas numpy matplotlib seaborn scipy

# 或使用 conda
conda install -c conda-forge pandas numpy matplotlib seaborn scipy
```

**运行**:
```bash
# 方法 1: 使用 mamba run
mamba run -n env-atmpi python downstream/creative_downstream_analysis.py

# 方法 2: 激活环境后运行
mamba activate env-atmpi
python downstream/creative_downstream_analysis.py
```

**功能特点**:
- ✅ 整合 ASV 和 MAG 分类学信息
- ✅ 计算匹配质量统计
- ✅ 识别高置信度配对 (≥99% identity)
- ✅ 分析菌株多样性
- ✅ 评估分类学一致性
- ✅ 生成 5 张可视化图表
- ✅ 生成 Markdown 综合报告

**输出文件**:
- `phylum_statistics.tsv` - 菌门水平统计
- `high_confidence_pairs.tsv` - 高置信度配对
- `strain_diversity.tsv` - 菌株多样性分析
- `taxonomy_consistency.tsv` - 分类一致性分析
- `identity_distribution.png` - 身份分布图
- `phylum_identity.png` - 菌门身份图
- `strain_diversity.png` - 菌株多样性图
- `taxonomy_consistency.png` - 分类一致性热图
- `high_confidence_genus.png` - 高置信度属分布图
- `analysis_report.md` - 综合分析报告

---

## R 下游分析脚本

**文件**: `downstream/creative_downstream_analysis.R`

**依赖安装**:
```r
# 在 R 中运行
install.packages(c("tidyverse", "ggplot2", "pheatmap", "data.table", 
                   "scales", "gridExtra", "patchwork"))

# 或使用 conda
conda install -c conda-forge r-tidyverse r-ggplot2 r-pheatmap r-datarable
```

**运行**:
```bash
# 方法 1: 直接运行
Rscript downstream/creative_downstream_analysis.R

# 方法 2: 在 RStudio 中打开并运行
# 或在 R 会话中 source
source("downstream/creative_downstream_analysis.R")
```

**功能特点**:
- ✅ 数据加载与整合 (使用 data.table 加速)
- ✅ 统计分析与可视化
- ✅ 使用 ggplot2 生成高质量图表
- ✅ 生成组合图 (all_plots_combined.png)
- ✅ 计算多样性指数 (Shannon, Evenness)
- ✅ 生成 Markdown 综合报告 (R 版本)

**输出文件**:
- `phylum_statistics.tsv` - 菌门水平统计
- `high_confidence_pairs.tsv` - 高置信度配对
- `strain_diversity.tsv` - 菌株多样性分析
- `taxonomy_consistency.tsv` - 分类一致性分析
- `mag_diversity_indices.tsv` - MAG 多样性指数
- `identity_distribution.png` - 身份分布图
- `phylum_identity.png` - 菌门身份箱线图
- `strain_diversity.png` - 菌株多样性柱状图
- `taxonomy_consistency.png` - 分类一致性热图
- `high_confidence_genus.png` - 高置信度属分布图
- `diversity_indices.png` - 多样性指数分布图
- `all_plots_combined.png` - 所有图表组合
- `analysis_report_R.md` - 综合分析报告 (R 版本)

---

## Python vs R 版本对比

| 特性 | Python 版本 | R 版本 |
|------|-----------|--------|
| 数据加载 | pandas | data.table |
| 可视化 | matplotlib + seaborn | ggplot2 |
| 热图 | seaborn.heatmap | pheatmap |
| 组合图 | 手动排列 | patchwork |
| 多样性指数 | ❌ | ✅ Shannon, Evenness |
| 报告生成 | Markdown | Markdown |
| 运行时间 | ~30 秒 | ~45 秒 |
| 输出图表数 | 5 | 7 |

---

## 输出文件说明

### 主要输出文件

| 文件路径 | 描述 |
|---------|------|
| `assignments/asv_to_mag_assignments.tsv` | ASV-MAG 配对结果 |
| `assignments/asv_strain_resolution.tsv` | 菌株分辨率分析 |
| `assignments/target_species_assignments.tsv` | 目标物种过滤结果 |
| `summary/asv_mag_summary.tsv` | 统计摘要 |
| `summary/visualization_report.html` | HTML 可视化报告 |
| `summary/per_mag_summary.tsv` | 每 MAG 统计 |

### 下游分析输出

**目录**: `downstream/analysis_results/`

| 文件 | 描述 |
|------|------|
| `analysis_report.md` | 综合分析报告 |
| `phylum_statistics.tsv` | 菌门水平统计 |
| `high_confidence_pairs.tsv` | 高置信度配对 (≥99%) |
| `strain_diversity.tsv` | 菌株多样性分析 |
| `taxonomy_consistency.tsv` | 分类一致性分析 |
| `identity_distribution.png` | 身份分布直方图 |
| `phylum_identity.png` | 菌门身份箱线图 |
| `strain_diversity.png` | 菌株多样性柱状图 |
| `taxonomy_consistency.png` | 分类一致性热图 |
| `high_confidence_genus.png` | 高置信度属分布图 |

---

## 配置文件说明

**config.yaml** 关键参数：

```yaml
# 输入文件
asv_table: "asv_table.tsv"
mag_table: "mag_table.tsv"

# BLAST 参数
blast_evalue: 0.0001
blast_max_target_seqs: 10

# 阈值
identity_threshold: 97    # 物种水平
coverage_threshold: 80    # 覆盖度
strain_threshold: 98.5    # 菌株水平

# 目标物种
target_species:
  - "Bifidobacterium longum"
  - "Gemella taiwanensis"
  - "Streptococcus oralis"
  - "Klebsiella pneumoniae"

# 计算资源
threads: 16
memory: "64G"
```

---

## 数据可用性

**输入数据**:
- ASV 表：`asv_table.tsv` (1,923 ASVs)
- MAG 表：`mag_table.tsv` (1,316 MAGs)
- ASV 分类学：`downstream/taxa_rdp_details.tsv`
- MAG 分类学：`downstream/MAGs_hmq.megahit.semibin_single.gtdbtk.all.tsv`

**分析环境**:
- Conda 环境：`env-atmpi`
- Snakemake 版本：最新
- Python 依赖：pandas, numpy, matplotlib, seaborn, scipy
- R 依赖：tidyverse, ggplot2, pheatmap, vegan

---

*报告生成时间：2026-03-20*  
*分析完成状态：✅ 成功*
