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
# 产物（直接放在 output_dir 下）:
#   {genome_name}_index.*.ht2  — HISAT2 索引（8 个文件）
#   {genome_name}.gtf          — 基因注释
# ============================================================

set -e

SCRIPT_DIR="$(cd "$(dirname "$0")" && pwd)"
PIPELINE_DIR="$(dirname "$SCRIPT_DIR")"
CONFIG_FILE="$PIPELINE_DIR/config/reference_sources.yaml"

if [ $# -lt 1 ]; then
    echo "用法: bash scripts/build_reference.sh <species_key>"
    echo "支持的 species_key: human, mouse, rat"
    exit 1
fi

SPECIES_KEY="$1"

# ---- 检查工具 ----
if ! command -v hisat2-build &>/dev/null; then
    echo "错误: hisat2-build 未找到，请先激活 bioinfo_env 环境"
    exit 1
fi

# ---- 解析 YAML ----
parse_yaml_val() {
    local key="$1"
    grep -A10 "^${SPECIES_KEY}:" "$CONFIG_FILE" | grep "^  ${key}:" | head -1 | sed "s/.*${key}: *\"\?\([^\"]*\)\"\?/\1/"
}

GENOME_NAME=$(parse_yaml_val "genome_name")
FASTA_URL=$(parse_yaml_val "fasta_url")
GTF_URL=$(parse_yaml_val "gtf_url")
OUTPUT_DIR=$(parse_yaml_val "output_dir")

if [ -z "$GENOME_NAME" ]; then
    echo "错误: 未找到物种 '$SPECIES_KEY' 的配置"
    exit 1
fi

echo "============================================"
echo " 参考基因组构建"
echo " 物种: $SPECIES_KEY ($GENOME_NAME)"
echo " 输出: $OUTPUT_DIR"
echo "============================================"

mkdir -p "$OUTPUT_DIR"

# ---- 下载 & 构建工作目录（临时，完成后清理） ----
WORK_DIR="$OUTPUT_DIR/.build_tmp"
rm -rf "$WORK_DIR"
mkdir -p "$WORK_DIR"
cd "$WORK_DIR"

cleanup() {
    echo ""
    echo "清理临时文件..."
    rm -rf "$WORK_DIR"
}
trap cleanup EXIT

FASTA_GZ=$(basename "$FASTA_URL")
FASTA_FILE="${FASTA_GZ%.gz}"
GTF_GZ=$(basename "$GTF_URL")
GTF_FILE="${GTF_GZ%.gz}"

# ---- Step 1-2: 下载 ----
echo ""
echo "[1/6] 下载 FASTA..."
wget -c "$FASTA_URL" -O "$FASTA_GZ"

echo ""
echo "[2/6] 下载 GTF..."
wget -c "$GTF_URL" -O "$GTF_GZ"

# ---- Step 3: MD5 校验 ----
echo ""
echo "[3/6] 校验完整性..."
MD5_URL="${FASTA_URL%/*}/md5checksums.txt"
wget -c "$MD5_URL" -O md5checksums.txt 2>/dev/null || echo "  警告: 无法下载 md5checksums.txt"
if [ -f md5checksums.txt ]; then
    grep "$(basename "$FASTA_URL")" md5checksums.txt | md5sum -c 2>/dev/null && echo "  FASTA 校验通过 ✓" || echo "  FASTA 校验失败 ✗"
    grep "$(basename "$GTF_URL")" md5checksums.txt | md5sum -c 2>/dev/null && echo "  GTF 校验通过 ✓" || echo "  GTF 校验失败 ✗"
fi

# ---- Step 4: 解压 ----
echo ""
echo "[4/6] 解压..."
gunzip -k "$FASTA_GZ"
gunzip -k "$GTF_GZ"

# ---- Step 5: 提取注释 ----
echo ""
echo "[5/6] 提取外显子和剪接位点..."
EXONS_FILE="$WORK_DIR/exons.txt"
SPLICE_FILE="$WORK_DIR/splice_sites.txt"

hisat2_extract_exons.py "$GTF_FILE" > "$EXONS_FILE"
echo "  exons.txt: $(wc -l < "$EXONS_FILE") lines"

hisat2_extract_splice_sites.py "$GTF_FILE" > "$SPLICE_FILE"
echo "  splice_sites.txt: $(wc -l < "$SPLICE_FILE") lines"

# ---- Step 6: 构建索引 ----
echo ""
echo "[6/6] 构建 HISAT2 索引（预计 1-3 小时）..."

INDEX_PREFIX="$WORK_DIR/${GENOME_NAME}_index"

hisat2-build -p 8 \
    --ss "$SPLICE_FILE" \
    --exon "$EXONS_FILE" \
    "$FASTA_FILE" \
    "$INDEX_PREFIX" 2>&1 | tee build.log

# ---- 产物搬运到 output_dir ----
echo ""
echo "搬运产物到 $OUTPUT_DIR ..."
mv "${GENOME_NAME}_index"*.ht2 "$OUTPUT_DIR/"
cp "$GTF_FILE" "$OUTPUT_DIR/${GENOME_NAME}.gtf"

echo ""
echo "============================================"
echo " 构建完成！"
echo ""
echo " 索引: $OUTPUT_DIR/${GENOME_NAME}_index.*.ht2"
echo " GTF:  $OUTPUT_DIR/${GENOME_NAME}.gtf"
echo ""
echo " 项目 config.yaml 配置:"
echo "   species: \"$SPECIES_KEY\""
echo "   reference:"
echo "     root: \"$(dirname "$OUTPUT_DIR")\""
echo "============================================"
