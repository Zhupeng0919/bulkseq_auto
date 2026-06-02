#!/usr/bin/env Rscript
# ============================================================
# 将 featureCounts 输出（counts.txt）转换为标准化的
# counts.xls / fpkm.xls / tpm.xls
# ============================================================
# 用法: Rscript counts_to_tpm_fpkm.R \
#   --input   <counts.txt> \
#   --outdir  <输出目录> \
#   --project <项目ID>
# ============================================================

suppressPackageStartupMessages(library(readr))
suppressPackageStartupMessages(library(dplyr))
suppressPackageStartupMessages(library(tidyr))

# ---- 解析命令行参数 ----
args <- commandArgs(trailingOnly = TRUE)
if (length(args) < 6) {
  stop("Usage: Rscript counts_to_tpm_fpkm.R --input <counts.txt> --outdir <dir> --project <id>")
}

input_file <- args[which(args == "--input") + 1]
outdir     <- args[which(args == "--outdir") + 1]
proj_id    <- args[which(args == "--project") + 1]

dir.create(outdir, showWarnings = FALSE, recursive = TRUE)

# ---- 读取数据 ----
counts <- read_tsv(input_file, skip = 1, show_col_types = FALSE) %>% as.data.frame()
rownames(counts) <- counts$Geneid

# 自动检测计数列起始位置（Chr,Start,End,Strand,Length 之后）
# featureCounts 标准格式前 6 列为：Geneid, Chr, Start, End, Strand, Length
# 之后的所有列为样本计数列
annot_cols <- 6  # Geneid + Chr + Start + End + Strand + Length
count_matrix <- counts[, (annot_cols + 1):ncol(counts), drop = FALSE]

# ---- 清洗列名（去除路径、.sorted.bam 后缀、Illumina 测序前缀） ----
clean_colnames <- gsub("^.*/|\\.sorted\\.bam$", "", colnames(count_matrix))
clean_colnames <- gsub("^[A-Z]\\d+_L\\d+_", "", clean_colnames)
colnames(count_matrix) <- clean_colnames

# 按字母顺序排列列
count_matrix <- count_matrix[, gtools::mixedorder(colnames(count_matrix)), drop = FALSE]

# ---- 计算 FPKM / TPM ----
kb <- counts$Length / 1000
count_mat <- as.matrix(count_matrix)

# FPKM
fpkm <- t(t(count_mat) / colSums(count_mat) * 10^6) / kb
# TPM
rpk <- count_mat / kb
tpm <- t(t(rpk) / colSums(rpk) * 10^6)

# ---- 写出带 Geneid 和 Length 的表 ----
write_output <- function(mat, label) {
  df <- as.data.frame(mat)
  df$Geneid <- rownames(df)
  df$Length <- counts$Length
  df <- df[, c("Geneid", "Length", setdiff(colnames(df), c("Geneid", "Length")))]
  outfile <- file.path(outdir, paste0(proj_id, "_", label, ".xls"))
  write.table(df, file = outfile, sep = "\t", quote = FALSE, row.names = FALSE)
  cat(sprintf("  ✓ %s → %s\n", label, outfile))
}

write_output(count_matrix, "counts")
write_output(fpkm, "fpkm")
write_output(tpm, "tpm")

cat(sprintf(">>> 转换完成，文件已保存至: %s\n", outdir))
