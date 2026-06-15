# bulkseq_auto

Bulk RNA-seq 一键自动化分析流程，基于 Snakemake 构建。从 FASTQ 到 PDF 报告全自动。

## 前置依赖

- [conda](https://docs.conda.io/en/latest/) / [mamba](https://mamba.readthedocs.io/)
- Snakemake >= 9

## 快速开始

```bash
# 1. 创建 conda 环境
conda env create -f bioenv/environment.yml
# 或使用便携版（不锁定子依赖版本）：
conda env create -f bioenv/environment_portable.yml

# 2. 激活环境
conda activate bioinfo_env

# 3. 构建参考基因组（仅 fastq 模式需要，首次一次性操作）
bash scripts/build_reference.sh human   # 可选 human / mouse / rat

# 4. 创建项目目录，准备输入数据，放入 config.yaml 和 sample_sheet.csv

# 5. 运行
bash run_pipeline.sh /path/to/project/ --mode counts --cores 16
```

## 构建参考基因组

fastq 模式需要 HISAT2 索引和基因注释 GTF。`scripts/build_reference.sh` 一键完成下载→校验→构建。

```bash
conda activate bioinfo_env
bash scripts/build_reference.sh <species>
```

| 参数 | 基因组 | 下载大小 | 构建耗时 | 内存需求 |
|------|--------|----------|----------|----------|
| `human` | GRCh38.p14 | ~3 GB | 2-3 小时 | 64 GB+ |
| `mouse` | GRCm39 | ~3 GB | 1-2 小时 | 32 GB+ |
| `rat` | GRCr8 | ~2 GB | 1-2 小时 | 32 GB+ |

> 下载源和输出路径在 `config/reference_sources.yaml` 中配置。

产物目录结构（`{reference_root}/{genome_name}/`）：

```
reference/
├── GRCh38.p14/
│   ├── GRCh38.p14_index.1.ht2 ... .8.ht2
│   └── GRCh38.p14.gtf
├── GRCm39/
│   ├── GRCm39_index.1.ht2 ... .8.ht2
│   └── GRCm39.gtf
└── GRCr8/
    ├── GRCr8_index.1.ht2 ... .8.ht2
    └── GRCr8.gtf
```

项目 config 中只需指定 `species` 和 `reference.root`，索引/GTF 路径自动推导：

```yaml
species: "human"
reference:
  root: "/home/zhuzp/reference"
  # index 和 gtf 留空，自动推导为:
  #   {root}/GRCh38.p14/GRCh38.p14_index
  #   {root}/GRCh38.p14/GRCh38.p14.gtf
```

手动指定 `index`/`gtf` 的优先级高于自动推导。

## 两种分析模式

| 模式 | 起点 | 输入要求 | 运行内容 |
|------|------|----------|----------|
| `counts` | 定量矩阵 | `04_Counts/counts.txt` | 定量转换 → 样本过滤 → DEG → GO/KEGG → 报告 |
| `fastq` | 原始测序数据 | `RawData/*.fq.gz` | fastp → FastQC → HISAT2 → featureCounts → ... → 报告 |

## 项目目录结构

```
project/
├── config.yaml          # 项目配置（必填）
├── sample_sheet.csv     # 样本信息表（必填）
├── RawData/             # FASTQ 文件（fastq 模式输入）
├── 01_CleanData/        # fastp 质控过滤后 FASTQ（fastq 模式输出）
├── 02_QC_Reports/       # FastQC + fastp + MultiQC 质控报告（fastq 模式输出）
├── 03_Alignment/        # HISAT2 比对 BAM 文件（fastq 模式输出）
├── 04_Counts/           # 定量矩阵及标准化表达量
│   └── filtered/        # 过滤后定量矩阵（下游分析输入）
├── 05_DEG/              # 差异表达分析结果
├── 06_GO_KEGG/          # GO/KEGG 富集分析结果
├── 07_GSEA/             # GSEA 分析结果
├── 00_Config/           # 配置及附属文件（报告后移入）
├── Report.pdf           # 最终 PDF 报告
```

> 流程完成后自动清理原始 `counts.txt`、`counts.txt.summary` 和 `logs/`，这些文件包含内部路径，不适合保留在输出目录中。定量信息已保留在 `04_Counts/filtered/` 中。

## 配置说明

### sample_sheet.csv

```csv
sample,group
sample1,Control
sample2,Control
sample3,Treatment1
sample4,Treatment1
```

- `sample`: 样本名，须与 counts 列名一致
- `group`: 分组名，第一组默认为对照（可通过 `comparison_mode` 调整）
- `batch`（可选）: 批次信息，启用 `batch.correct` 时需要

### config.yaml 关键参数

| 参数 | 说明 | 默认值 |
|------|------|--------|
| `project_dir` | 项目根目录绝对路径 | - |
| `project_id` | 项目标识符 | - |
| `mode` | 分析模式: `fastq` \| `counts` | `fastq` |
| `species` | 物种: `human` \| `mouse` \| `rat` | `human` |
| `reference.root` | 参考基因组根目录（自动推导 index/gtf 时使用） | `""` |
| `reference.index` | HISAT2 索引前缀（手动指定，优先于自动推导） | `""` |
| `reference.gtf` | GTF 注释文件路径（手动指定，优先于自动推导） | `""` |
| `deg.logFC_cutoff` | 差异倍数阈值 | `0.5` |
| `deg.p_cutoff` | 显著性阈值 | `0.05` |
| `comparison_mode` | 比较策略: `auto` \| `pairwise` \| `control_vs_rest` \| `all_vs_all` | `auto` |
| `comparisons` | 手动指定比较列表（优先级高于自动模式） | `[]` |
| `exclude_samples` | 剔除离群样本列表 | `[]` |
| `enrichment.go_kegg` | 是否运行 GO/KEGG 富集 | `true` |
| `enrichment.gsea` | 是否运行 GSEA | `false` |

## 输出说明

**上游（fastq 模式）：**
- `01_CleanData/` — fastp 质控过滤后的 FASTQ 文件
- `02_QC_Reports/` — fastp HTML + FastQC HTML + MultiQC 汇总报告
- `03_Alignment/` — HISAT2 比对 BAM + BAI 索引
- `total_mapping_summary.txt` — 各样本比对率汇总

**下游（counts/fastq 模式）：**
- `04_Counts/` — counts/fpkm/tpm 标准化表达矩阵
- `04_Counts/filtered/` — 样本过滤后的矩阵及样本表（下游分析输入，安全保留）
- `05_DEG/` — 差异基因 CSV + 火山图 + 热图（PDF/PNG）
- `06_GO_KEGG/` — GO/KEGG 富集结果 CSV + 富集图（PDF/PNG）
- `07_GSEA/` — GSEA 分析结果
- `Report.pdf` — 整合 PDF 报告（含表号和图号，带目录）

> 报告完成后，`config.yaml`、`total_mapping_summary.txt` 和 `Report.tex` 移入 `00_Config/`；原始 `counts.txt`、`counts.txt.summary` 和 `logs/` 自动删除（包含内部路径）。

## License

GPL-3.0
