# Bulkseq 自动化分析报告使用说明

## 运行完整流程

```bash
cd /path/to/bulkseq_auto

# 从 counts 矩阵开始（下游分析）
bash run_pipeline.sh /path/to/project

# 从 FASTQ 开始（完整流程）
bash run_pipeline.sh /path/to/project --mode fastq

# 预览运行计划
bash run_pipeline.sh /path/to/project --dry

# 指定 CPU 核心数
bash run_pipeline.sh /path/to/project --cores 8
```

## 仅重新生成报告

修改 `scripts/generate_report.Rmd` 后：

```bash
cd /path/to/bulkseq_auto
conda activate bioinfo_env
snakemake --configfile /path/to/project/config.yaml --cores 16 -R report
```

## 重新生成 GO/KEGG 富集图

修改 `scripts/go_kegg.R` 后：

```bash
cd /path/to/bulkseq_auto
conda activate bioinfo_env
snakemake --configfile /path/to/project/config.yaml --cores 16 -R go_kegg report
```

## 项目目录结构

```
<项目目录>/
  config.yaml              # 项目配置（物种、模式、阈值、参考路径）
  sample_sheet.csv          # 样本分组信息（sample, group）
  05_Counts/                # 定量结果（counts, FPKM, TPM）
  06_DEG/                   # 差异分析结果 + 火山图/热图
  07_GO_KEGG/               # GO/KEGG 富集结果 + 组合图
  08_GSEA/                  # GSEA 结果（可选）
  Report.pdf                # 最终 PDF 报告
```

## config.yaml 关键参数

```yaml
project_dir: "/path/to/project"
project_id: "Project-XXX"
mode: "counts"               # counts（下游）或 fastq（完整流程）
species: "human"             # human 或 mouse
deg:
  logFC_cutoff: 0.5
  p_cutoff: 0.05
enrichment:
  go_kegg: true
  gsea: false
```
