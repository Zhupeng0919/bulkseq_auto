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
# 目录结构（以 human 为例）:
#   ref/GRCh38.p14/                      — 原始参考文件 (FASTA, GTF)
#   ref/GRCh38.p14_His2_refrence/        — HISAT2 索引及注释提取产物
#       ├── GRCh38.p14_index.*.ht2       (HISAT2 索引)
#       ├── exons.txt
#       └── splice_sites.txt
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

# 原始参考文件目录 / HISAT2 索引目录
REF_DIR="$OUTPUT_DIR"
IDX_DIR="${OUTPUT_DIR}_His2_refrence"

echo "============================================"
echo " 参考基因组构建"
echo " 物种: $SPECIES_KEY ($GENOME_NAME)"
echo " 参考文件: $REF_DIR"
echo " 索引目录: $IDX_DIR"
echo "============================================"

mkdir -p "$REF_DIR" "$IDX_DIR"

FASTA_GZ=$(basename "$FASTA_URL")
FASTA_FILE="${FASTA_GZ%.gz}"
GTF_GZ=$(basename "$GTF_URL")
GTF_FILE="${GTF_GZ%.gz}"

# ---- Step 1-2: 下载（跳过已存在且非空的文件） ----
echo ""
echo "[1/6] 下载 FASTA..."
if [ -s "$REF_DIR/$FASTA_GZ" ]; then
    echo "  已存在，跳过: $REF_DIR/$FASTA_GZ"
else
    wget -c "$FASTA_URL" -O "$REF_DIR/$FASTA_GZ"
fi

echo ""
echo "[2/6] 下载 GTF..."
if [ -s "$REF_DIR/$GTF_GZ" ]; then
    echo "  已存在，跳过: $REF_DIR/$GTF_GZ"
else
    wget -c "$GTF_URL" -O "$REF_DIR/$GTF_GZ"
fi

# ---- Step 3: MD5 校验（如果未校验过） ----
echo ""
echo "[3/6] 校验完整性..."
MD5_URL="${FASTA_URL%/*}/md5checksums.txt"
MD5_FILE="$REF_DIR/.md5_checked"
if [ -f "$MD5_FILE" ]; then
    echo "  已校验，跳过"
else
    wget -c "$MD5_URL" -O "$REF_DIR/md5checksums.txt" 2>/dev/null || echo "  警告: 无法下载 md5checksums.txt"
    if [ -f "$REF_DIR/md5checksums.txt" ]; then
        (cd "$REF_DIR" && grep "$FASTA_GZ" md5checksums.txt | md5sum -c 2>/dev/null) && echo "  FASTA 校验通过 ✓" || echo "  FASTA 校验失败 ✗"
        (cd "$REF_DIR" && grep "$GTF_GZ" md5checksums.txt | md5sum -c 2>/dev/null) && echo "  GTF 校验通过 ✓" || echo "  GTF 校验失败 ✗"
    fi
    touch "$MD5_FILE"
fi

# ---- Step 4: 解压 ----
echo ""
echo "[4/6] 解压..."
if [ -s "$REF_DIR/$FASTA_FILE" ]; then
    echo "  FASTA 已解压，跳过: $REF_DIR/$FASTA_FILE"
else
    gunzip -k "$REF_DIR/$FASTA_GZ"
fi
if [ -s "$REF_DIR/$GTF_FILE" ]; then
    echo "  GTF 已解压，跳过: $REF_DIR/$GTF_FILE"
else
    gunzip -k "$REF_DIR/$GTF_GZ"
fi

# ---- Step 5: 提取注释 ----
echo ""
echo "[5/6] 提取外显子和剪接位点..."
EXONS_FILE="$IDX_DIR/exons.txt"
SPLICE_FILE="$IDX_DIR/splice_sites.txt"

hisat2_extract_exons.py "$REF_DIR/$GTF_FILE" > "$EXONS_FILE"
echo "  exons.txt: $(wc -l < "$EXONS_FILE") lines"

hisat2_extract_splice_sites.py "$REF_DIR/$GTF_FILE" > "$SPLICE_FILE"
echo "  splice_sites.txt: $(wc -l < "$SPLICE_FILE") lines"

# ---- Step 6: 构建索引 ----
echo ""
echo "[6/6] 构建 HISAT2 索引（预计 1-3 小时）..."

INDEX_PREFIX="$IDX_DIR/${GENOME_NAME}_index"

# 如果索引已完整（至少 6 个 ht2 文件），跳过构建
EXISTING_HT2=$(find "$IDX_DIR" -maxdepth 1 -name "${GENOME_NAME}_index*.ht2" 2>/dev/null | wc -l)
if [ "$EXISTING_HT2" -ge 6 ]; then
    echo "  索引已存在 ($EXISTING_HT2 个 .ht2 文件)，跳过构建"
else
    hisat2-build -p 8 \
        --ss "$SPLICE_FILE" \
        --exon "$EXONS_FILE" \
        "$REF_DIR/$FASTA_FILE" \
        "$INDEX_PREFIX" 2>&1 | tee "$IDX_DIR/build.log"
fi

echo ""
echo "============================================"
echo " 构建完成！"
echo ""
echo " 参考文件: $REF_DIR/"
echo " 索引目录: $IDX_DIR/"
echo " GTF:      $REF_DIR/$GTF_FILE"
echo ""
echo " 项目 config.yaml 配置:"
echo "   species: \"$SPECIES_KEY\""
echo "   reference:"
echo "     root: \"$(dirname "$OUTPUT_DIR")\""
echo "============================================"
