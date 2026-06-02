#!/bin/bash
# ============================================================
# Bulkseq 自动化分析流程 — 入口脚本
# ============================================================
# 用法:
#   Mode 2 (从定量矩阵开始):
#     bash run_pipeline.sh /path/to/project
#
#   Mode 1 (从原始 FASTQ 开始):
#     bash run_pipeline.sh /path/to/project --mode fastq
#
#   预览 DAG:
#     bash run_pipeline.sh /path/to/project --dry
# ============================================================

set -e

PIPELINE_DIR="$(cd "$(dirname "$0")" && pwd)"
CONDA_ENV="bioinfo_env"

# ---- 解析参数 ----
if [ $# -lt 1 ]; then
    echo "用法: bash run_pipeline.sh <项目目录> [选项]"
    echo "选项:"
    echo "  --mode fastq|counts   分析模式（默认从 config.yaml 读取）"
    echo "  --dry                 预览运行计划（不执行）"
    echo "  --cores N             CPU 核心数（默认 16）"
    echo "  --unlock              解锁工作目录（异常退出后使用）"
    exit 1
fi

PROJECT_DIR=$(realpath "$1")
shift

MODE=""
DRY=""
CORES=16
UNLOCK=""

while [ $# -gt 0 ]; do
    case "$1" in
        --mode) MODE="$2"; shift 2 ;;
        --dry) DRY="-n"; shift ;;
        --cores) CORES="$2"; shift 2 ;;
        --unlock) UNLOCK="--unlock"; shift ;;
        *) echo "未知参数: $1"; exit 1 ;;
    esac
done

# ---- 验证项目目录 ----
if [ ! -d "$PROJECT_DIR" ]; then
    echo "错误: 项目目录不存在: $PROJECT_DIR"
    exit 1
fi

# ---- 查找/生成 config.yaml ----
CONFIG_FILE="$PROJECT_DIR/config.yaml"
if [ ! -f "$CONFIG_FILE" ]; then
    echo "未找到 config.yaml，正在从模板生成..."
    if [ -f "$PIPELINE_DIR/config/config.yaml" ]; then
        cp "$PIPELINE_DIR/config/config.yaml" "$CONFIG_FILE"
        echo "请编辑 $CONFIG_FILE 并填写项目参数，然后重新运行。"
        exit 1
    else
        echo "错误: 模板 config.yaml 不存在"
        exit 1
    fi
fi

# ---- 如果指定了 --mode，覆盖 config ----
EXTRA_CONFIG=""
if [ -n "$MODE" ]; then
    EXTRA_CONFIG="mode=${MODE}"
fi

# ---- 验证 sample_sheet.csv ----
if [ ! -f "$PROJECT_DIR/sample_sheet.csv" ]; then
    if [ -f "$PIPELINE_DIR/config/sample_sheet.csv" ]; then
        cp "$PIPELINE_DIR/config/sample_sheet.csv" "$PROJECT_DIR/sample_sheet.csv"
        echo "已创建 sample_sheet.csv 模板，请编辑后重新运行。"
        exit 1
    fi
fi

# ---- 激活环境 ----
source ~/.bashrc
echo "激活 conda 环境: $CONDA_ENV"
conda activate "$CONDA_ENV"

# ---- 运行 Snakemake ----
echo "========================================"
echo "项目: $(basename "$PROJECT_DIR")"
echo "目录: $PROJECT_DIR"
echo "模式: $MODE"
echo "核心: $CORES"
echo "========================================"

cd "$PIPELINE_DIR"

CMD="snakemake --configfile '$CONFIG_FILE' --cores $CORES $DRY $UNLOCK"
if [ -n "$EXTRA_CONFIG" ]; then
    CMD="$CMD --config $EXTRA_CONFIG"
fi

echo "运行: $CMD"
eval "$CMD"
EXIT_CODE=$?

if [ $EXIT_CODE -eq 1 ]; then
    echo ""
    echo "========================================"
    echo "流程在检查点暂停，请按提示检查中间结果，"
    echo "修改 config.yaml 中的 checkpoint 标志后重新运行。"
    echo "========================================"
elif [ $EXIT_CODE -eq 0 ]; then
    echo "========================================"
    echo "流程完成！报告: $PROJECT_DIR/Report.pdf"
    echo "========================================"
else
    echo "错误: Snakemake 退出码 $EXIT_CODE"
    exit $EXIT_CODE
fi
