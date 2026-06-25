#!/bin/bash
#SBATCH --job-name=bulkseq_upstream
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=1
#SBATCH --cpus-per-task=26
#SBATCH --mem=64G
#SBATCH --time=24:00:00
#SBATCH --output=bulkseq_upstream_%j.log

# ============================================================
# bulkseq_auto — 独立上游分析脚本 (fastp → HISAT2 → featureCounts)
# ============================================================
# 用法:
#   1. 编辑下方 === 配置区 === 的参数
#   2. 直接运行: bash scripts/run_upstream.sh
#   3. 或提交集群: sbatch scripts/run_upstream.sh
#
# 依赖: conda 环境 bioinfo_env
# ============================================================

set -e

# =================================================================
# === 配置区 ===
# =================================================================

# --- 项目路径 ---
PROJECT_DIR="/home/zhuzp/bulkseq/FQ260512001"

# --- 参考基因组 ---
# 方式一: 自动推导（推荐，需先运行 build_reference.sh 构建索引）
SPECIES="human"                             # human / mouse / rat
REF_ROOT="/home/zhuzp/reference"            # 参考基因组根目录

# 方式二: 手动指定（优先级高于自动推导，留空则使用方式一）
INDEX_PREFIX=""                             # HISAT2 索引前缀（不含 .1.ht2）
GTF_FILE=""                                 # 基因注释 GTF 路径

# --- 流程控制 ---
SKIP_QC="NO"                                # YES=跳过质控，直接比对
QC_ONLY="NO"                                # YES=仅运行质控，完成后退出
SKIP_ALIGN="NO"                             # YES=跳过比对，直接定量（需已有 BAM）

# --- 参数 ---
TRIM_FRONT=15                               # fastp 5' 端切除碱基数
MULTIMAPPING="NO"                           # YES=计入多比对 reads（-M --fraction）

# --- 环境 ---
CONDA_ENV="bioinfo_env"

# =================================================================
# === 初始化 ===
# =================================================================

source ~/.bashrc
conda activate "$CONDA_ENV"

# 自动推导参考路径
if [ -z "$INDEX_PREFIX" ] || [ -z "$GTF_FILE" ]; then
    case "$SPECIES" in
        human)  GENOME_NAME="GRCh38.p14" ;;
        mouse)  GENOME_NAME="GRCm39" ;;
        rat)    GENOME_NAME="GRCr8" ;;
        *)      echo "错误: 未知物种 '$SPECIES'"; exit 1 ;;
    esac
    REF_DIR="${REF_ROOT}/${GENOME_NAME}"
    IDX_DIR="${REF_ROOT}/${GENOME_NAME}_His2_refrence"
    [ -z "$INDEX_PREFIX" ] && INDEX_PREFIX="${IDX_DIR}/${GENOME_NAME}_index"
    [ -z "$GTF_FILE" ]     && GTF_FILE=$(ls "${REF_DIR}"/*.gtf 2>/dev/null | head -1)
    echo "自动推导参考路径:"
    echo "  索引: $INDEX_PREFIX"
    echo "  GTF:  $GTF_FILE"
fi

# 验证
if [ ! -f "${INDEX_PREFIX}.1.ht2" ]; then
    echo "错误: HISAT2 索引不存在: ${INDEX_PREFIX}.1.ht2"
    echo "请先构建: bash scripts/build_reference.sh $SPECIES"
    exit 1
fi
if [ ! -f "$GTF_FILE" ]; then
    echo "错误: GTF 文件不存在: $GTF_FILE"
    exit 1
fi

# 目录
CLEAN_DIR="${PROJECT_DIR}/01_CleanData"
QC_DIR="${PROJECT_DIR}/02_QC_Reports"
ALIGN_DIR="${PROJECT_DIR}/03_Alignment"
COUNT_DIR="${PROJECT_DIR}/05_Counts"
LOG_DIR="${PROJECT_DIR}/logs"
SUMMARY_FILE="${PROJECT_DIR}/total_mapping_summary.txt"

mkdir -p "$CLEAN_DIR" "$QC_DIR" "$ALIGN_DIR" "$COUNT_DIR" "$LOG_DIR"

BOLD="\033[1m"
GREEN="\033[32m"
YELLOW="\033[33m"
RED="\033[31m"
RESET="\033[0m"

banner() { echo -e "\n${BOLD}${GREEN}==== $1 ====${RESET}\n"; }
warn()  { echo -e "${YELLOW}[WARN]${RESET} $*"; }
err()   { echo -e "${RED}[ERROR]${RESET} $*"; }

# =================================================================
# === 样本发现 ===
# =================================================================

banner "样本发现"

# 优先子目录布局: RawData/<样本名>/*_1.fq.gz
SAMPLES=()
RAW_DIR="${PROJECT_DIR}/RawData"
if [ ! -d "$RAW_DIR" ]; then
    err "RawData 目录不存在: $RAW_DIR"
    exit 1
fi

for dir in "$RAW_DIR"/*/; do
    dir_name=$(basename "$dir")
    r1=$(ls "$dir"/*_1.fq.gz 2>/dev/null | head -1)
    if [ -n "$r1" ]; then
        SAMPLES+=("$dir_name")
    fi
done

# 回退: 扁平布局 RawData/<样本名>_1.fq.gz
if [ ${#SAMPLES[@]} -eq 0 ]; then
    for r1 in "$RAW_DIR"/*_1.fq.gz; do
        [ -f "$r1" ] || continue
        base=$(basename "$r1" _1.fq.gz)
        SAMPLES+=("$base")
    done
fi

if [ ${#SAMPLES[@]} -eq 0 ]; then
    err "未找到任何 *_1.fq.gz 文件（子目录或扁平布局）"
    exit 1
fi

echo "发现 ${#SAMPLES[@]} 个样本: ${SAMPLES[*]}"

# 辅助函数：获取样本的 R1/R2 路径
get_r1() {
    local s="$1"
    local sub="${RAW_DIR}/${s}/"*_1.fq.gz
    local flat="${RAW_DIR}/${s}_1.fq.gz"
    ls $sub 2>/dev/null | head -1 || echo "$flat"
}
get_r2() {
    local r1; r1=$(get_r1 "$1")
    echo "${r1//_1.fq.gz/_2.fq.gz}"
}

# =================================================================
# === Step 1: 质控 ===
# =================================================================

if [ "$SKIP_QC" == "YES" ]; then
    warn "SKIP_QC=YES, 跳过质控"
else
    banner "Step 1/3: fastp + FastQC"

    echo "fastp 参数: --trim_front1 ${TRIM_FRONT} --trim_front2 ${TRIM_FRONT} -w 4"

    for sample in "${SAMPLES[@]}"; do
        R1=$(get_r1 "$sample")
        R2=$(get_r2 "$sample")

        if [ ! -f "$R2" ]; then
            err "缺少配对 R2: $sample"
            exit 1
        fi

        echo "  [$sample] fastp..."
        fastp -i "$R1" -I "$R2" \
              -o "${CLEAN_DIR}/${sample}_1.clean.fq.gz" \
              -O "${CLEAN_DIR}/${sample}_2.clean.fq.gz" \
              --trim_front1 "$TRIM_FRONT" --trim_front2 "$TRIM_FRONT" \
              -w 4 -h "${QC_DIR}/${sample}_fastp.html" -j "${QC_DIR}/${sample}_fastp.json" \
              2>&1 | tail -1

        echo "  [$sample] FastQC..."
        fastqc "${CLEAN_DIR}/${sample}_1.clean.fq.gz" \
               "${CLEAN_DIR}/${sample}_2.clean.fq.gz" \
               -t 2 -o "$QC_DIR" --quiet
    done

    echo ""
    echo "MultiQC..."
    multiqc "$QC_DIR" -o "$QC_DIR" -n "multiqc_final_report" --force --quiet

    echo "质控完成"
fi

if [ "$QC_ONLY" == "YES" ]; then
    banner "QC_ONLY=YES, 脚本退出"
    exit 0
fi

# =================================================================
# === Step 2: 比对 ===
# =================================================================

if [ "$SKIP_ALIGN" == "YES" ]; then
    warn "SKIP_ALIGN=YES, 跳过 HISAT2 比对"
else
    banner "Step 2/3: HISAT2 比对"

    for sample in "${SAMPLES[@]}"; do
        R1="${CLEAN_DIR}/${sample}_1.clean.fq.gz"
        R2="${CLEAN_DIR}/${sample}_2.clean.fq.gz"
        BAM="${ALIGN_DIR}/${sample}.sorted.bam"
        HISAT2_LOG="${LOG_DIR}/${sample}_hisat2.log"

        echo "  [$sample] 比对中..."
        hisat2 -p 10 -x "$INDEX_PREFIX" \
               -1 "$R1" -2 "$R2" \
               2> "$HISAT2_LOG" | \
        samtools view -@ 2 -bS - | \
        samtools sort -@ 4 -m 4G -o "$BAM" -

        samtools index "$BAM"

        rate=$(grep "overall alignment rate" "$HISAT2_LOG" | awk '{print $1}')
        echo "    -> $rate"
    done

    echo "比对完成"
fi

# =================================================================
# === Step 3: 定量 & 汇总 ===
# =================================================================

banner "Step 3/3: featureCounts 定量"

MULTI_FLAG=""
[ "$MULTIMAPPING" == "YES" ] && MULTI_FLAG="-M --fraction"

BAMS=("${ALIGN_DIR}"/*.sorted.bam)
if [ ${#BAMS[@]} -eq 0 ]; then
    err "未找到 BAM 文件: ${ALIGN_DIR}/*.sorted.bam"
    exit 1
fi

echo "定量 ${#BAMS[@]} 个 BAM 文件..."
featureCounts -T 16 -p -B -C -s 0 $MULTI_FLAG \
              -a "$GTF_FILE" \
              -o "${COUNT_DIR}/counts.txt" \
              "${BAMS[@]}" 2> "${LOG_DIR}/fc.log"

# 比对率汇总
echo -e "SampleID\tHISAT2_Mapping_Rate" > "$SUMMARY_FILE"
for log in "${LOG_DIR}"/*_hisat2.log; do
    [ -f "$log" ] || continue
    s_name=$(basename "$log" _hisat2.log)
    h_rate=$(grep "overall alignment rate" "$log" | awk '{print $1}')
    echo -e "${s_name}\t${h_rate}" >> "$SUMMARY_FILE"
done

echo ""
echo "比对率汇总:"
cat "$SUMMARY_FILE"

banner "上游分析完成"
echo "counts.txt:  ${COUNT_DIR}/counts.txt"
echo "BAM 目录:   ${ALIGN_DIR}/"
echo "质控报告:   ${QC_DIR}/multiqc_final_report.html"
echo "日志目录:   ${LOG_DIR}/"
