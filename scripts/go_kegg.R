#!/usr/bin/env Rscript
# ============================================================
# GO/KEGG 富集分析（clusterProfiler）— v2 上下调分别富集
# ============================================================
# 用法: Rscript go_kegg.R \
#   --deg         <Diff.csv> \
#   --outdir      <输出目录> \
#   --project     <项目ID> \
#   --species     human|mouse|rat \
#   [--comparison AB3_vs_SC] \
#   [--pval       0.05]
# ============================================================

suppressPackageStartupMessages(library(clusterProfiler))
suppressPackageStartupMessages(library(dplyr))
suppressPackageStartupMessages(library(ggplot2))

# ---- 解析参数 ----
args <- commandArgs(trailingOnly = TRUE)
if (length(args) < 8) {
  stop("Usage: Rscript go_kegg.R --deg <file> --outdir <dir> --project <id> --species human|mouse|rat [--comparison <label>] [--pval 0.05]")
}

deg_file    <- args[which(args == "--deg") + 1]
outdir      <- args[which(args == "--outdir") + 1]
proj_id     <- args[which(args == "--project") + 1]
species     <- args[which(args == "--species") + 1]
p_cut       <- ifelse("--pval" %in% args, as.numeric(args[which(args == "--pval") + 1]), 0.05)
comp_label  <- ifelse("--comparison" %in% args, args[which(args == "--comparison") + 1], "")

dir.create(outdir, showWarnings = FALSE, recursive = TRUE)

# 解析比较标签用于标题
if (nchar(comp_label) > 0) {
  comp_parts <- strsplit(comp_label, "_vs_")[[1]]
  if (length(comp_parts) == 2) {
    title_kd <- comp_parts[1]
    title_ctrl <- comp_parts[2]
  } else {
    title_kd <- comp_label
    title_ctrl <- ""
  }
} else {
  title_kd <- proj_id
  title_ctrl <- ""
}

# ---- 物种配置 ----
if (species == "human") {
  OrgDb <- "org.Hs.eg.db"
  kegg_org <- "hsa"
} else if (species == "mouse") {
  OrgDb <- "org.Mm.eg.db"
  kegg_org <- "mmu"
} else if (species == "rat") {
  OrgDb <- "org.Rn.eg.db"
  kegg_org <- "rno"
} else {
  stop("species 必须为 human、mouse 或 rat")
}

suppressPackageStartupMessages(library(OrgDb, character.only = TRUE))

# ---- 读取 DEG 结果 ----
all_diff <- read.csv(deg_file, stringsAsFactors = FALSE)

# ---- 对上下调分别富集 ----
run_enrichment <- function(gene_list, direction_label) {
  prefix <- paste0(comp_label, "_", direction_label)

  genes <- unique(gene_list)
  cat(sprintf("\n>>> %s 差异基因数: %d\n", direction_label, length(genes)))

  if (length(genes) < 3) {
    cat(sprintf("!!! %s 差异基因太少，跳过富集分析\n", direction_label))
    return(invisible())
  }

  # ID 转换
  gene_convert <- tryCatch(
    bitr(genes, fromType = "SYMBOL", toType = "ENTREZID", OrgDb = OrgDb),
    error = function(e) NULL
  )

  if (is.null(gene_convert) || nrow(gene_convert) < 3) {
    cat(sprintf("!!! %s ID 转换后基因太少，跳过富集分析\n", direction_label))
    return(invisible())
  }

  entrz <- gene_convert$ENTREZID

  # GO 富集
  go_results <- enrichGO(
    gene          = entrz,
    OrgDb         = OrgDb,
    ont           = "ALL",
    pAdjustMethod = "BH",
    pvalueCutoff  = p_cut,
    qvalueCutoff  = p_cut,
    readable      = TRUE
  )

  go_file <- file.path(outdir, paste0(prefix, "_GO_enrich_results.csv"))
  if (!is.null(go_results) && nrow(go_results) > 0) {
    write.csv(as.data.frame(go_results), go_file, row.names = FALSE)
    cat(sprintf("  ✓ GO 富集: %d 条通路\n", nrow(go_results)))
  } else {
    cat("  ! GO 无显著富集结果\n")
  }

  # KEGG 富集
  kegg_results <- tryCatch(
    enrichKEGG(gene = entrz, organism = kegg_org,
               pvalueCutoff = p_cut, qvalueCutoff = p_cut),
    error = function(e) NULL
  )

  kegg_file <- file.path(outdir, paste0(prefix, "_KEGG_enrich_results.csv"))
  if (!is.null(kegg_results) && nrow(kegg_results) > 0) {
    write.csv(as.data.frame(kegg_results), kegg_file, row.names = FALSE)
    cat(sprintf("  ✓ KEGG 富集: %d 条通路\n", nrow(kegg_results)))
  } else {
    cat("  ! KEGG 无显著富集结果\n")
  }

  # ---- 绘图：组合柱状图（Category 标签 + Count 气泡 + -log10(p.adjust) 柱）----
  type_color <- c(
    "BP"   = "#e87d8d",
    "CC"   = "#f8c387",
    "MF"   = "#a3d393",
    "KEGG" = "#7cc8e9"
  )

  plot_data_list <- list()

  if (!is.null(go_results) && nrow(go_results) > 0) {
    go_df <- as.data.frame(go_results) %>% mutate(Category = ONTOLOGY)
    plot_data_list[["GO"]] <- go_df
  }
  if (!is.null(kegg_results) && nrow(kegg_results) > 0) {
    kegg_df <- as.data.frame(kegg_results) %>% mutate(Category = "KEGG")
    plot_data_list[["KEGG"]] <- kegg_df
  }

  if (length(plot_data_list) > 0) {
    plot_data <- bind_rows(plot_data_list) %>%
      filter(p.adjust < p_cut) %>%
      group_by(Category) %>%
      slice_max(order_by = Count, n = 5, with_ties = FALSE) %>%
      mutate(logP = -log10(p.adjust)) %>%
      ungroup() %>%
      mutate(Description = ifelse(nchar(as.character(Description)) > 32,
             paste0(substr(Description, 1, 29), "..."), as.character(Description))) %>%
      mutate(Category = factor(Category, levels = c("KEGG", "MF", "CC", "BP"))) %>%
      arrange(Category, logP) %>%
      mutate(Description = factor(Description, levels = unique(Description)))

    if (nrow(plot_data) > 0) {
      title_text <- if (nchar(title_ctrl) > 0) {
        paste0(title_kd, " VS ", title_ctrl, " ", direction_label)
      } else {
        paste0(proj_id, " ", direction_label)
      }

      xaxis_max <- max(plot_data$logP, na.rm = TRUE) * 1.2
      left_w <- xaxis_max * 0.5

      rect_data <- plot_data %>%
        group_by(Category) %>%
        summarise(n = n(), .groups = "drop") %>%
        mutate(
          ymax = cumsum(n) + 0.4,
          ymin = lag(ymax, default = 0) + 0.6,
          xmin = -left_w, xmax = -left_w * 0.55
        )

      p <- ggplot(plot_data, aes(x = logP, y = Description, fill = Category)) +
        geom_col(width = 0.6, alpha = 0.85) +
        geom_text(aes(x = xaxis_max * 0.85, label = Description),
                  hjust = 0, size = 3.5, fontface = "bold") +
        geom_point(aes(x = -left_w * 0.3, size = Count, fill = Category),
                   shape = 21, color = "black", stroke = 1) +
        geom_text(aes(x = -left_w * 0.3, label = Count), size = 3, fontface = "bold") +
        scale_size_continuous(range = c(6, 12)) +
        geom_rect(data = rect_data,
                  aes(xmin = xmin, xmax = xmax, ymin = ymin, ymax = ymax, fill = Category),
                  inherit.aes = FALSE) +
        geom_text(data = rect_data,
                  aes(x = (xmin + xmax) / 2, y = (ymin + ymax) / 2, label = Category),
                  inherit.aes = FALSE, angle = 90, color = "white", fontface = "bold", size = 5.5) +
        annotate(geom = "segment", x = 0, y = 0.4, xend = xaxis_max, yend = 0.4,
                 linewidth = 1, color = "black") +
        scale_fill_manual(values = type_color) +
        scale_x_continuous(expand = expansion(mult = c(0, 0.05))) +
        coord_cartesian(xlim = c(-left_w * 1.05, xaxis_max * 1.35), clip = "off") +
        theme_minimal() +
        theme(
          axis.text.y = element_blank(), axis.ticks.y = element_blank(),
          axis.line.x = element_line(linewidth = 0.8, color = "black"),
          panel.grid = element_blank(), legend.position = "right",
          plot.title = element_text(hjust = 0.5, face = "bold", size = 16),
          plot.margin = margin(l = 10, r = 200, t = 10, b = 10)
        ) +
        labs(x = "-log10(p.adjust)", y = NULL,
             title = paste0(title_text, " GO & KEGG Enrichment"))

      fig_base <- file.path(outdir, paste0("Fig_Enrichment_", direction_label, "_", comp_label))
      ggsave(paste0(fig_base, ".pdf"), p, width = 14, height = 10)
      ggsave(paste0(fig_base, ".png"), p, width = 14, height = 10, dpi = 300)
      cat(sprintf("  ✓ 富集图已保存 (PDF + PNG): %s\n", fig_base))
    }
  }
}

# 上调基因富集
up_genes <- all_diff %>% filter(change == "UP") %>% pull(Symbol)
run_enrichment(up_genes, "UP")

# 下调基因富集
down_genes <- all_diff %>% filter(change == "DOWN") %>% pull(Symbol)
run_enrichment(down_genes, "DOWN")

cat(sprintf("\n>>> GO/KEGG 分析完成，结果保存至: %s\n", outdir))
