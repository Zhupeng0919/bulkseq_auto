#!/usr/bin/env Rscript
# ============================================================
# GSEA 富集分析 — gseGO + gseKEGG + GSEA(MSigDB Hallmark)
# ============================================================
# 用法: Rscript gsea.R \
#   --deg        <per-comparison Diff.csv> \
#   --outdir     <输出目录> \
#   --project    <项目ID> \
#   --species    human|mouse|rat \
#   --comparison <比较标签>
# ============================================================

suppressPackageStartupMessages({
  library(clusterProfiler)
  library(dplyr)
  library(ggplot2)
  library(enrichplot)
  library(msigdbr)
})

# ---- 解析参数 ----
args <- commandArgs(trailingOnly = TRUE)
deg_file   <- args[which(args == "--deg") + 1]
outdir     <- args[which(args == "--outdir") + 1]
proj_id    <- args[which(args == "--project") + 1]
species    <- args[which(args == "--species") + 1]
comparison <- args[which(args == "--comparison") + 1]

dir.create(outdir, showWarnings = FALSE, recursive = TRUE)
set.seed(12345)

# ---- 物种配置 ----
species_config <- list(
  human = list(OrgDb = "org.Hs.eg.db", kegg = "hsa", msigdbr = "Homo sapiens"),
  mouse = list(OrgDb = "org.Mm.eg.db", kegg = "mmu", msigdbr = "Mus musculus"),
  rat   = list(OrgDb = "org.Rn.eg.db", kegg = "rno", msigdbr = "Rattus norvegicus")
)
cfg <- species_config[[species]]
if (is.null(cfg)) stop("Unknown species: ", species)

library(cfg$OrgDb, character.only = TRUE)
org_db <- get(cfg$OrgDb)

# ---- 读取 DEG，构建排序基因列表 ----
deg <- read.csv(deg_file, stringsAsFactors = FALSE)
# 仅保留有 log2FoldChange 的行
deg <- deg[!is.na(deg$log2FoldChange), ]
# 按 log2FoldChange 降序排列
deg <- deg[order(deg$log2FoldChange, decreasing = TRUE), ]

gene_list <- deg$log2FoldChange
names(gene_list) <- deg$Symbol

cat(sprintf("基因数: %d\n", length(gene_list)))

# ---- 公共 GSEA 参数 ----
gsea_params <- list(
  minGSSize     = 10,
  maxGSSize     = 500,
  pvalueCutoff  = 0.05,
  pAdjustMethod = "BH",
  eps           = 1e-10,
  by            = "fgsea",
  seed          = TRUE
)

prefix <- paste0(proj_id, "_", comparison)

# ============================================================
# 1. GO 富集 (gseGO, SYMBOL)
# ============================================================
cat("\n>>> GO BP GSEA\n")
go_res <- tryCatch(
  suppressMessages(gseGO(
    geneList     = gene_list,
    OrgDb        = org_db,
    ont          = "BP",
    keyType      = "SYMBOL",
    minGSSize    = gsea_params$minGSSize,
    maxGSSize    = gsea_params$maxGSSize,
    pvalueCutoff = gsea_params$pvalueCutoff,
    pAdjustMethod = gsea_params$pAdjustMethod,
    eps          = gsea_params$eps,
    by           = gsea_params$by,
    seed         = gsea_params$seed
  )),
  error = function(e) { cat(sprintf("  ! GO Error: %s\n", e$message)); NULL }
)

if (!is.null(go_res) && nrow(go_res) > 0) {
  cat(sprintf("  ✓ GO: %d 通路\n", nrow(go_res)))
  go_res <- setReadable(go_res, OrgDb = org_db, keyType = "SYMBOL")
  write.csv(go_res@result, file.path(outdir, paste0(prefix, "_GSEA_GO.csv")), row.names = FALSE)

  # Up Top5: NES > 0, p.adjust < 0.05
  top5_up <- go_res@result %>% filter(p.adjust < 0.05, NES > 0) %>% arrange(desc(NES)) %>% head(5) %>% pull(ID)
  if (length(top5_up) > 0) {
    p_gsea <- gseaplot2(go_res, geneSetID = top5_up, title = paste(prefix, "- GO BP Up Top5"), pvalue_table = TRUE)
    ggsave(file.path(outdir, paste0(prefix, "_GSEA_GO_Top5_Up.pdf")), p_gsea, width = 14, height = 10)
    ggsave(file.path(outdir, paste0(prefix, "_GSEA_GO_Top5_Up.png")), p_gsea, width = 14, height = 10, dpi = 300)
  }

  # Down Top5: NES < 0, p.adjust < 0.05
  top5_down <- go_res@result %>% filter(p.adjust < 0.05, NES < 0) %>% arrange(NES) %>% head(5) %>% pull(ID)
  if (length(top5_down) > 0) {
    p_gsea <- gseaplot2(go_res, geneSetID = top5_down, title = paste(prefix, "- GO BP Down Top5"), pvalue_table = TRUE)
    ggsave(file.path(outdir, paste0(prefix, "_GSEA_GO_Top5_Down.pdf")), p_gsea, width = 14, height = 10)
    ggsave(file.path(outdir, paste0(prefix, "_GSEA_GO_Top5_Down.png")), p_gsea, width = 14, height = 10, dpi = 300)
  }
} else {
  cat("  ! GO: 无显著富集\n")
}

# ============================================================
# 2. KEGG 富集 (gseKEGG, ENTREZID)
# ============================================================
cat("\n>>> KEGG GSEA\n")

# SYMBOL → ENTREZID
id_map <- tryCatch(
  bitr(names(gene_list), fromType = "SYMBOL", toType = "ENTREZID", OrgDb = org_db),
  error = function(e) { cat(sprintf("  ! ID转换失败: %s\n", e$message)); NULL }
)

if (!is.null(id_map) && nrow(id_map) > 0) {
  # 合并并去重
  deg_kegg <- data.frame(SYMBOL = names(gene_list), logFC = gene_list, stringsAsFactors = FALSE)
  deg_kegg <- merge(deg_kegg, id_map, by = "SYMBOL")
  deg_kegg <- deg_kegg[order(deg_kegg$logFC, decreasing = TRUE), ]
  deg_kegg <- deg_kegg[!duplicated(deg_kegg$ENTREZID), ]

  gene_list_kegg <- deg_kegg$logFC
  names(gene_list_kegg) <- deg_kegg$ENTREZID

  cat(sprintf("  转换后基因数: %d (丢失 %.1f%%)\n",
              length(gene_list_kegg),
              (1 - length(gene_list_kegg)/length(gene_list)) * 100))

  kegg_res <- tryCatch(
    suppressMessages(gseKEGG(
      geneList      = gene_list_kegg,
      organism      = cfg$kegg,
      keyType       = "kegg",
      minGSSize     = gsea_params$minGSSize,
      maxGSSize     = gsea_params$maxGSSize,
      pvalueCutoff  = gsea_params$pvalueCutoff,
      pAdjustMethod = gsea_params$pAdjustMethod,
      eps           = gsea_params$eps,
      by            = gsea_params$by,
      seed          = gsea_params$seed
    )),
    error = function(e) { cat(sprintf("  ! KEGG Error: %s\n", e$message)); NULL }
  )

  if (!is.null(kegg_res) && nrow(kegg_res) > 0) {
    cat(sprintf("  ✓ KEGG: %d 通路\n", nrow(kegg_res)))
    kegg_res <- setReadable(kegg_res, OrgDb = org_db, keyType = "ENTREZID")
    write.csv(kegg_res@result, file.path(outdir, paste0(prefix, "_GSEA_KEGG.csv")), row.names = FALSE)

    # Up Top5: NES > 0, p.adjust < 0.05
    top5_up <- kegg_res@result %>% filter(p.adjust < 0.05, NES > 0) %>% arrange(desc(NES)) %>% head(5) %>% pull(ID)
    if (length(top5_up) > 0) {
      p_gsea <- gseaplot2(kegg_res, geneSetID = top5_up, title = paste(prefix, "- KEGG Up Top5"), pvalue_table = TRUE)
      ggsave(file.path(outdir, paste0(prefix, "_GSEA_KEGG_Top5_Up.pdf")), p_gsea, width = 14, height = 10)
      ggsave(file.path(outdir, paste0(prefix, "_GSEA_KEGG_Top5_Up.png")), p_gsea, width = 14, height = 10, dpi = 300)
    }

    # Down Top5: NES < 0, p.adjust < 0.05
    top5_down <- kegg_res@result %>% filter(p.adjust < 0.05, NES < 0) %>% arrange(NES) %>% head(5) %>% pull(ID)
    if (length(top5_down) > 0) {
      p_gsea <- gseaplot2(kegg_res, geneSetID = top5_down, title = paste(prefix, "- KEGG Down Top5"), pvalue_table = TRUE)
      ggsave(file.path(outdir, paste0(prefix, "_GSEA_KEGG_Top5_Down.pdf")), p_gsea, width = 14, height = 10)
      ggsave(file.path(outdir, paste0(prefix, "_GSEA_KEGG_Top5_Down.png")), p_gsea, width = 14, height = 10, dpi = 300)
    }
  } else {
    cat("  ! KEGG: 无显著富集\n")
  }
} else {
  cat("  ! ID转换失败，跳过KEGG\n")
}

# ============================================================
# 3. MSigDB Hallmark 富集 (GSEA, SYMBOL)
# ============================================================
cat("\n>>> Hallmark GSEA\n")

hallmark_gs <- tryCatch(
  msigdbr(species = cfg$msigdbr, collection = "H") %>%
    dplyr::select(gs_name, gene_symbol) %>%
    distinct(),
  error = function(e) { cat(sprintf("  ! msigdbr Error: %s\n", e$message)); NULL }
)

if (!is.null(hallmark_gs) && nrow(hallmark_gs) > 0) {
  h_res <- tryCatch(
    suppressMessages(GSEA(
      geneList      = gene_list,
      TERM2GENE     = hallmark_gs,
      minGSSize     = gsea_params$minGSSize,
      maxGSSize     = gsea_params$maxGSSize,
      pvalueCutoff  = gsea_params$pvalueCutoff,
      pAdjustMethod = gsea_params$pAdjustMethod,
      eps           = gsea_params$eps,
      by            = gsea_params$by,
      seed          = gsea_params$seed
    )),
    error = function(e) { cat(sprintf("  ! Hallmark Error: %s\n", e$message)); NULL }
  )

  if (!is.null(h_res) && nrow(h_res) > 0) {
    cat(sprintf("  ✓ Hallmark: %d 通路\n", nrow(h_res)))
    write.csv(h_res@result, file.path(outdir, paste0(prefix, "_GSEA_Hallmark.csv")), row.names = FALSE)

    # Up Top5: NES > 0, p.adjust < 0.05
    top5_up <- h_res@result %>% filter(p.adjust < 0.05, NES > 0) %>% arrange(desc(NES)) %>% head(5) %>% pull(ID)
    if (length(top5_up) > 0) {
      p_gsea <- gseaplot2(h_res, geneSetID = top5_up, title = paste(prefix, "- Hallmark Up Top5"), pvalue_table = TRUE)
      ggsave(file.path(outdir, paste0(prefix, "_GSEA_Hallmark_Top5_Up.pdf")), p_gsea, width = 14, height = 10)
      ggsave(file.path(outdir, paste0(prefix, "_GSEA_Hallmark_Top5_Up.png")), p_gsea, width = 14, height = 10, dpi = 300)
    }

    # Down Top5: NES < 0, p.adjust < 0.05
    top5_down <- h_res@result %>% filter(p.adjust < 0.05, NES < 0) %>% arrange(NES) %>% head(5) %>% pull(ID)
    if (length(top5_down) > 0) {
      p_gsea <- gseaplot2(h_res, geneSetID = top5_down, title = paste(prefix, "- Hallmark Down Top5"), pvalue_table = TRUE)
      ggsave(file.path(outdir, paste0(prefix, "_GSEA_Hallmark_Top5_Down.pdf")), p_gsea, width = 14, height = 10)
      ggsave(file.path(outdir, paste0(prefix, "_GSEA_Hallmark_Top5_Down.png")), p_gsea, width = 14, height = 10, dpi = 300)
    }
  } else {
    cat("  ! Hallmark: 无显著富集\n")
  }
} else {
  cat("  ! 无法加载 Hallmark 基因集\n")
}

cat(sprintf("\n>>> GSEA 完成: %s\n", outdir))
