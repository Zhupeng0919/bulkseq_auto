#!/usr/bin/env Rscript
# ============================================================
# GSEA 分析 — 基于 msigdbr C2 CP 通路
# ============================================================
# 用法: Rscript gsea.R \
#   --counts    <counts.xls> \
#   --deg       <All_Diff.csv> \
#   --outdir    <输出目录> \
#   --project   <项目ID> \
#   --species   human|mouse|rat \
#   [--genes    GENE1,GENE2,...]
# ============================================================

suppressPackageStartupMessages(library(clusterProfiler))
suppressPackageStartupMessages(library(dplyr))
suppressPackageStartupMessages(library(ggplot2))
suppressPackageStartupMessages(library(msigdbr))
suppressPackageStartupMessages(library(enrichplot))

# ---- 解析参数 ----
args <- commandArgs(trailingOnly = TRUE)
if (length(args) < 10) {
  stop("Usage: Rscript gsea.R --counts <file> --deg <file> --outdir <dir> --project <id> --species human|mouse|rat [--genes GENE1,GENE2,...]")
}

counts_file <- args[which(args == "--counts") + 1]
deg_file    <- args[which(args == "--deg") + 1]
outdir      <- args[which(args == "--outdir") + 1]
proj_id     <- args[which(args == "--project") + 1]
species     <- args[which(args == "--species") + 1]

target_genes <- NULL
if ("--genes" %in% args) {
  genes_str <- args[which(args == "--genes") + 1]
  target_genes <- trimws(unlist(strsplit(genes_str, ",")))
}

dir.create(outdir, showWarnings = FALSE, recursive = TRUE)

# ---- 物种配置 ----
msigdbr_species <- if (species == "mouse") "Mus musculus" else if (species == "rat") "Rattus norvegicus" else "Homo sapiens"

# ---- 读取表达矩阵（FPKM 或 TPM，用于相关性排序） ----
expr <- read.table(counts_file, header = TRUE, sep = "\t", row.names = 1, check.names = FALSE)
if ("Length" %in% colnames(expr)) {
  expr <- expr[, !colnames(expr) %in% "Length", drop = FALSE]
}

# ---- 确定目标基因 ----
if (is.null(target_genes) || length(target_genes) == 0) {
  # 从 DEG 结果中取 top 基因
  all_diff <- read.csv(deg_file, stringsAsFactors = FALSE)
  degs <- all_diff %>% filter(change != "NOT") %>% arrange(padj)
  if (nrow(degs) == 0) {
    cat("!!! 无差异基因，跳过 GSEA\n")
    file.create(file.path(outdir, "SKIPPED"))
    quit(save = "no")
  }
  target_genes <- head(degs$Symbol, 3)  # 最多取 TOP 3
  cat(sprintf("自动选取 TOP %d 差异基因进行 GSEA\n", length(target_genes)))
}

# 只保留在表达矩阵中的基因
target_genes <- intersect(target_genes, rownames(expr))
if (length(target_genes) == 0) {
  cat("!!! 目标基因不在表达矩阵中，跳过 GSEA\n")
  file.create(file.path(outdir, "SKIPPED"))
  quit(save = "no")
}

cat(sprintf("GSEA 目标基因: %s\n", paste(target_genes, collapse = ", ")))

# ---- 基因集 ----
m_df <- msigdbr(species = msigdbr_species, category = "C2", subcategory = "CP")
geneset <- m_df %>% dplyr::select(gs_name, gene_symbol) %>% distinct()

# ---- 逐基因 GSEA ----
idx <- 1
for (gene in target_genes) {
  cat(sprintf("\n>>> [%02d] GSEA for: %s\n", idx, gene))

  gene_expr <- as.numeric(expr[gene, ])
  cor_results <- apply(expr, 1, function(x) cor(gene_expr, x, method = "spearman"))
  genelist <- sort(cor_results, decreasing = TRUE)

  gsea_res <- tryCatch(
    GSEA(genelist, TERM2GENE = geneset,
         pvalueCutoff = 0.05, pAdjustMethod = "BH", eps = 0, verbose = FALSE),
    error = function(e) NULL
  )

  if (!is.null(gsea_res) && nrow(gsea_res) > 0) {
    prefix <- sprintf("%02d_%s", idx, gene)

    write.csv(as.data.frame(gsea_res), file.path(outdir, paste0(prefix, "_GSEA_Result.csv")), row.names = FALSE)

    # Top 5 通路
    top_pathways <- gsea_res@result %>%
      arrange(desc(abs(NES))) %>%
      head(5) %>%
      pull(ID)

    if (length(top_pathways) > 0) {
      p <- gseaplot2(gsea_res, geneSetID = top_pathways,
                     title = paste("GSEA for", gene),
                     rel_heights = c(1.5, 0.5, 1),
                     pvalue_table = FALSE,
                     color = RColorBrewer::brewer.pal(min(5, length(top_pathways)), "Set1"))
      plot_file <- file.path(outdir, paste0(prefix, "_GSEA_Plot"))
      ggsave(paste0(plot_file, ".pdf"), p, width = 10, height = 8)
      ggsave(paste0(plot_file, ".png"), p, width = 10, height = 8, dpi = 300)
      cat(sprintf("  ✓ %s: %d 通路\n", gene, nrow(gsea_res)))
    }
    idx <- idx + 1
  } else {
    cat(sprintf("  ! %s: 无显著富集\n", gene))
  }
}

cat(sprintf("\n>>> GSEA 分析完成，结果保存至: %s\n", outdir))
