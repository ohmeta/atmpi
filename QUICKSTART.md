# ATMPI 快速启动指南

## 🚀 快速开始

### 1. 激活环境并运行流程

```bash
# 激活 conda 环境
mamba activate env-atmpi

# 运行 Snakemake 流程
snakemake --cores 16 --configfile config.yaml

# 或先进行 dry-run 检查
snakemake -n --cores 16 --configfile config.yaml
```

### 2. 运行下游分析

#### Python 版本
```bash
mamba run -n env-atmpi python downstream/creative_downstream_analysis.py
```

#### R 版本
```bash
Rscript downstream/creative_downstream_analysis.R
```

---

## 📁 项目结构

```
atmpi/
├── Snakefile                    # Snakemake 流程定义
├── config.yaml                  # 配置文件
├── asv_table.tsv               # ASV 输入表
├── mag_table.tsv               # MAG 输入表
├── scripts/                    # Python 脚本
│   ├── validate_inputs.py
│   ├── extract_16s.py
│   ├── process_blast_results.py
│   └── ...
├── downstream/                 # 下游分析
│   ├── creative_downstream_analysis.py   # Python 分析脚本
│   ├── creative_downstream_analysis.R    # R 分析脚本
│   ├── taxa_rdp_details.tsv             # ASV 分类学
│   └── MAGs_hmq.megahit.semibin_single.gtdbtk.all.tsv  # MAG 分类学
├── assignments/                # 输出：ASV-MAG 分配
├── summary/                    # 输出：统计摘要
└── downstream/analysis_results/ # 下游分析结果
```

---

## 📊 主要输出文件

### 流程输出

| 文件 | 描述 |
|------|------|
| `assignments/asv_to_mag_assignments.tsv` | ASV-MAG 配对结果 |
| `assignments/asv_strain_resolution.tsv` | 菌株分辨率分析 |
| `summary/asv_mag_summary.tsv` | 统计摘要 |
| `summary/visualization_report.html` | HTML 报告 |

### 下游分析输出

**目录**: `downstream/analysis_results/`

| 文件 | 描述 |
|------|------|
| `analysis_report.md` | Python 版综合报告 |
| `analysis_report_R.md` | R 版综合报告 |
| `identity_distribution.png` | 身份分布图 |
| `phylum_identity.png` | 菌门身份图 |
| `strain_diversity.png` | 菌株多样性图 |
| `taxonomy_consistency.png` | 分类一致性热图 |
| `high_confidence_genus.png` | 高置信度属分布 |
| `all_plots_combined.png` | 所有图表组合 |
| `diversity_indices.png` | 多样性指数图 |

---

## 🔧 配置说明

编辑 `config.yaml` 调整参数：

```yaml
# 输入文件
asv_table: "asv_table.tsv"
mag_table: "mag_table.tsv"

# 阈值设置
identity_threshold: 97      # 物种水平 (97%)
coverage_threshold: 80      # 覆盖度 (80%)
strain_threshold: 98.5      # 菌株水平 (98.5%)

# 目标物种
target_species:
  - "Bifidobacterium longum"
  - "Streptococcus oralis"

# 计算资源
threads: 16
```

---

## 📈 结果解读

### ASV-MAG 分配文件

**asv_to_mag_assignments.tsv** 列说明：

| 列名 | 描述 |
|------|------|
| `asv_id` | ASV 标识符 |
| `mag_id` | MAG 标识符 |
| `identity` | 匹配身份 (%) |
| `alignment_length` | 比对长度 |
| `evalue` | BLAST E 值 |
| `source` | 来源 (contig/16s) |

### 关键指标

- **Identity ≥ 97%**: 物种水平匹配
- **Identity ≥ 99%**: 菌株水平匹配
- **高置信度**: Identity ≥ 99% 的配对
- **菌株多样性**: 一个 MAG 匹配多个 ASV

---

## 🐛 常见问题

### 1. 流程运行失败

```bash
# 检查输入文件格式
head asv_table.tsv
head mag_table.tsv

# 清理中间文件重新运行
snakemake --cores 16 --configfile config.yaml --rerun-incomplete
```

### 2. 依赖缺失

```bash
# Python 依赖
pip install pandas numpy matplotlib seaborn scipy

# R 依赖
R -e "install.packages(c('tidyverse', 'ggplot2', 'pheatmap', 'data.table'))"
```

### 3. 内存不足

```bash
# 减少线程数
snakemake --cores 8 --configfile config.yaml
```

---

## 📚 参考资源

- [Snakemake 文档](https://snakemake.readthedocs.io/)
- [GTDB-Tk 分类](https://ecogenomics.github.io/GTDBTk/)
- [16S rRNA 基因数据库](https://www.arb-silva.de/)

---

## 📞 联系与支持

如有问题，请查看：
- `README.md` - 项目详细说明
- `CHAT.md` - 完整分析报告
- `DEBUG.md` - 调试指南

---

*最后更新：2026-03-20*
