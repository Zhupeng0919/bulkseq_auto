#!/usr/bin/env Rscript
# ============================================================
# 样本剔除：根据 exclude_samples 列表过滤 counts/fpkm/tpm 和 sample_sheet
# ============================================================
# 用法: Rscript filter_samples.R \
#   --counts    <counts.xls> \
#   --sample    <sample_sheet.csv> \
#   --outdir    <输出目录> \
#   --project   <项目ID> \
#   --exclude   "sample1,sample2"
# ============================================================

# ---- 解析参数 ----
args <- commandArgs(trailingOnly = TRUE)
if (length(args) < 10) {
  stop("Usage: Rscript filter_samples.R --counts <file> --sample <file> --outdir <dir> --project <id> --exclude <list>")
}

counts_file <- args[which(args == "--counts") + 1]
sample_file <- args[which(args == "--sample") + 1]
outdir      <- args[which(args == "--outdir") + 1]
proj_id     <- args[which(args == "--project") + 1]
exclude_str <- args[which(args == "--exclude") + 1]

exclude_samples <- if (nchar(exclude_str) > 0) unlist(strsplit(exclude_str, ",")) else character(0)

# ---- 过滤函数 ----
filter_matrix <- function(file_path, label) {
  if (!file.exists(file_path)) return()
  mat <- read.table(file_path, header = TRUE, sep = "\t", row.names = 1, check.names = FALSE)
  has_length <- "Length" %in% colnames(mat)
  if (has_length) {
    len_col <- mat$Length
    mat <- mat[, !colnames(mat) %in% "Length", drop = FALSE]
  }
  keep_cols <- setdiff(colnames(mat), exclude_samples)
  mat <- mat[, keep_cols, drop = FALSE]
  if (has_length) {
    mat$Geneid <- rownames(mat)
    mat$Length <- len_col
    mat <- mat[, c("Geneid", "Length", setdiff(colnames(mat), c("Geneid", "Length")))]
  }
  outfile <- file.path(outdir, paste0(proj_id, "_", label, "_filtered.xls"))
  write.table(mat, file = outfile, sep = "\t", quote = FALSE, row.names = FALSE)
  cat(sprintf("  ✓ %s_filtered → %s\n", label, outfile))
}

# ---- 主流程 ----
if (length(exclude_samples) > 0) {
  cat(sprintf("剔除样本: %s\n", paste(exclude_samples, collapse = ", ")))
} else {
  cat("无样本需要剔除，直接复制原始文件\n")
}

# 过滤 counts
filter_matrix(counts_file, "counts")

# 过滤 fpkm
fpkm_file <- file.path(outdir, paste0(proj_id, "_fpkm.xls"))
filter_matrix(fpkm_file, "fpkm")

# 过滤 tpm
tpm_file <- file.path(outdir, paste0(proj_id, "_tpm.xls"))
filter_matrix(tpm_file, "tpm")

# 过滤 sample_sheet
samples <- read.csv(sample_file, stringsAsFactors = FALSE)
if (length(exclude_samples) > 0) {
  samples <- samples[!samples$sample %in% exclude_samples, , drop = FALSE]
  cat(sprintf("sample_sheet: %d → %d 样本\n", nrow(samples) + length(exclude_samples), nrow(samples)))
}
write.csv(samples, file.path(outdir, "sample_sheet_filtered.csv"), row.names = FALSE, quote = FALSE)
cat("  ✓ sample_sheet_filtered.csv 已保存\n")

# 同时复制原始的 counts.txt.summary（若有）
summary_file <- file.path(outdir, "counts.txt.summary")
if (file.exists(summary_file)) {
  file.copy(summary_file, file.path(outdir, "counts_filtered.txt.summary"), overwrite = TRUE)
}

cat(">>> 样本过滤完成\n")
