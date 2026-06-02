#!/usr/bin/env Rscript
# ============================================================
# 批次效应校正：使用 limma::removeBatchEffect 对 TPM 进行校正
# 校正后的 TPM 用于 PCA 和相关性热图等可视化
# DEG 分析中 batch 通过 DESeq2 设计公式 (~ batch + group) 处理
# ============================================================
# 用法: Rscript batch_correct.R \
#   --tpm    <tpm_filtered.xls> \
#   --sample <sample_sheet_filtered.csv> \
#   --outdir <输出目录> \
#   --project <项目ID>
# ============================================================

suppressPackageStartupMessages(library(limma))

# ---- 解析参数 ----
args <- commandArgs(trailingOnly = TRUE)
if (length(args) < 8) {
  stop("Usage: Rscript batch_correct.R --tpm <file> --sample <file> --outdir <dir> --project <id>")
}

tpm_file    <- args[which(args == "--tpm") + 1]
sample_file <- args[which(args == "--sample") + 1]
outdir      <- args[which(args == "--outdir") + 1]
proj_id     <- args[which(args == "--project") + 1]

# ---- 读取数据 ----
tpm <- read.table(tpm_file, header = TRUE, sep = "\t", row.names = 1, check.names = FALSE)
has_length <- "Length" %in% colnames(tpm)
gene_lengths <- NULL
if (has_length) {
  gene_lengths <- tpm$Length
  tpm <- tpm[, !colnames(tpm) %in% "Length", drop = FALSE]
}

samples <- read.csv(sample_file, stringsAsFactors = FALSE)

# ---- 校验 batch 列 ----
if (!"batch" %in% colnames(samples)) {
  stop("sample_sheet 中未找到 'batch' 列，无法进行批次校正")
}

# ---- 样本对齐 ----
common <- intersect(colnames(tpm), samples$sample)
if (length(common) == 0) {
  stop("样本名在 TPM 和 sample_sheet 中不匹配！")
}
tpm <- tpm[, common, drop = FALSE]
samples <- samples[match(common, samples$sample), , drop = FALSE]

batch <- factor(samples$batch)
group <- factor(samples$group)

cat(sprintf("样本数: %d, 基因数: %d\n", ncol(tpm), nrow(tpm)))
cat(sprintf("批次: %s\n", paste(levels(batch), collapse = ", ")))
cat(sprintf("分组: %s\n", paste(levels(group), collapse = ", ")))
cat("批次分布:\n")
print(table(batch, group))

# ---- 批次校正 ----
if (length(levels(batch)) < 2) {
  cat("警告: 仅检测到 1 个批次，批次校正将无效果，直接复制原始 TPM\n")
  tpm_corrected <- tpm
} else {
  log_tpm <- log2(as.matrix(tpm) + 1)
  design <- model.matrix(~ group, data = samples)
  log_tpm_corrected <- removeBatchEffect(log_tpm, batch = batch, design = design)
  tpm_corrected <- 2^log_tpm_corrected - 1
  tpm_corrected[tpm_corrected < 0] <- 0
  cat("批次校正完成\n")
}

# ---- 写出结果 ----
tpm_out <- as.data.frame(tpm_corrected)
tpm_out$Geneid <- rownames(tpm_out)
if (has_length) {
  tpm_out$Length <- gene_lengths
  tpm_out <- tpm_out[, c("Geneid", "Length", setdiff(colnames(tpm_out), c("Geneid", "Length")))]
} else {
  tpm_out <- tpm_out[, c("Geneid", setdiff(colnames(tpm_out), "Geneid"))]
}

outfile <- file.path(outdir, paste0(proj_id, "_tpm_batch_corrected.xls"))
write.table(tpm_out, file = outfile, sep = "\t", quote = FALSE, row.names = FALSE)
cat(sprintf("✓ 批次校正完成 → %s\n", outfile))
