#!/bin/bash
# ============================================================
# 参考基因组自动构建脚本
# ============================================================
# 用法:
#   conda activate bioinfo_env
#   bash scripts/build_reference.sh <species_key>
#
# 支持的 species_key: human, mouse, rat
#
# 流程:
#   1. 从 config/reference_sources.yaml 读取配置
#   2. wget -c 断点续传下载 FASTA + GTF
#   3. wget -c 下载 md5checksums.txt 并校验
#   4. gunzip 解压
#   5. hisat2_extract_exons.py / hisat2_extract_splice_sites.py 提取注释
#   6. hisat2-build 构建索引
# ============================================================

set -e

SCRIPT_DIR="$(cd "$(dirname "$0")" && pwd)"
PIPELINE_DIR="$(dirname "$SCRIPT_DIR")"
CONFIG_FILE="$PIPELINE_DIR/config/reference_sources.yaml"

if [ $# -lt 1 ]; then
    echo "用法: bash scripts/build_reference.sh <species_key>"
    echo "支持的 species_key: human, mouse, rat"
    echo "示例: bash scripts/build_reference.sh human"
    exit 1
fi

SPECIES_KEY="$1"

# ---- 检查 conda 环境 ----
if ! command -v hisat2-build &>/dev/null; then
    echo "错误: hisat2-build 未找到，请先激活 bioinfo_env 环境"
    echo "  conda activate bioinfo_env"
    exit 1
fi

# ---- 从 YAML 解析配置 ----
parse_yaml_val() {
    local key="$1"
    grep -A10 "^${SPECIES_KEY}:" "$CONFIG_FILE" | grep "^  ${key}:" | head -1 | sed "s/.*${key}: *\"\?\([^\"]*\)\"\?/\1/"
}

GENOME_NAME=$(parse_yaml_val "genome_name")
FASTA_URL=$(parse_yaml_val "fasta_url")
GTF_URL=$(parse_yaml_val "gtf_url")
OUTPUT_DIR=$(parse_yaml_val "output_dir")

if [ -z "$GENOME_NAME" ]; then
    echo "错误: 未找到物种 '$SPECIES_KEY' 的配置，请检查 $CONFIG_FILE"
    exit 1
fi

echo "============================================"
echo " 参考基因组构建"
echo " 物种: $SPECIES_KEY"
echo " 版本: $GENOME_NAME"
echo " 输出: $OUTPUT_DIR"
echo "============================================"

# ---- 创建工作目录 ----
WORK_DIR="$OUTPUT_DIR/build_work"
mkdir -p "$WORK_DIR"
mkdir -p "$OUTPUT_DIR"
cd "$WORK_DIR"

# ---- Step 1: 下载 FASTA ----
FASTA_GZ=$(basename "$FASTA_URL")
FASTA_FILE="${FASTA_GZ%.gz}"

echo ""
echo "[1/6] 下载基因组 FASTA..."
wget -c "$FASTA_URL" -O "$FASTA_GZ"

# ---- Step 2: 下载 GTF ----
GTF_GZ=$(basename "$GTF_URL")
GTF_FILE="${GTF_GZ%.gz}"

echo ""
echo "[2/6] 下载基因注释 GTF..."
wget -c "$GTF_URL" -O "$GTF_GZ"

# ---- Step 3: 下载并校验 md5 ----
echo ""
echo "[3/6] 校验文件完整性..."
MD5_URL="${FASTA_URL%/*}/md5checksums.txt"
wget -c "$MD5_URL" -O md5checksums.txt 2>/dev/null || {
    echo "  警告: md5checksums.txt 下载失败，跳过校验"
}
if [ -f md5checksums.txt ]; then
    echo "  ✓ 校验 FASTA..."
    grep "$(basename "$FASTA_URL")" md5checksums.txt | md5sum -c 2>/dev/null && echo "    FASTA 校验通过" || echo "    FASTA 校验失败或警告"
    echo "  ✓ 校验 GTF..."
    grep "$(basename "$GTF_URL")" md5checksums.txt | md5sum -c 2>/dev/null && echo "    GTF 校验通过" || echo "    GTF 校验失败或警告"
fi

# ---- Step 4: 解压 ----
echo ""
echo "[4/6] 解压文件..."
if [ ! -f "$FASTA_FILE" ]; then
    echo "  解压 $FASTA_GZ ..."
    gunzip -k "$FASTA_GZ"  # -k 保留 .gz 文件
fi
if [ ! -f "$GTF_FILE" ]; then
    echo "  解压 $GTF_GZ ..."
    gunzip -k "$GTF_GZ"
fi

# ---- Step 5: 提取注释 ----
echo ""
echo "[5/6] 提取外显子和剪接位点..."
EXONS_FILE="$WORK_DIR/exons.txt"
SPLICE_FILE="$WORK_DIR/splice_sites.txt"

hisat2_extract_exons.py "$GTF_FILE" > "$EXONS_FILE"
echo "  ✓ exons.txt ($(wc -l < "$EXONS_FILE") lines)"

hisat2_extract_splice_sites.py "$GTF_FILE" > "$SPLICE_FILE"
echo "  ✓ splice_sites.txt ($(wc -l < "$SPLICE_FILE") lines)"

# ---- Step 6: 构建索引 ----
echo ""
echo "[6/6] 构建 HISAT2 索引（预计 1-3 小时）..."

INDEX_PREFIX="$OUTPUT_DIR/${GENOME_NAME}_index"

hisat2-build -p 8 \
    --ss "$SPLICE_FILE" \
    --exon "$EXONS_FILE" \
    "$FASTA_FILE" \
    "$INDEX_PREFIX" 2>&1 | tee build.log

echo ""
echo "============================================"
echo " 构建完成！"
echo " 索引路径: $INDEX_PREFIX"
echo ""
echo " 在项目 config.yaml 中配置:"
echo "   reference:"
echo "     index: \"$INDEX_PREFIX\""
echo "     gtf: \"$WORK_DIR/$GTF_FILE\""
echo "============================================"
