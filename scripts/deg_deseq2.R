#!/usr/bin/env Rscript
# ============================================================
# DESeq2 差异分析 + 火山图 + 热图
# 支持三种比较模式：pairwise / control_vs_rest / all_vs_all
# 每个比较独立输出 CSV、火山图(PDF+PNG)、热图(PDF+PNG)
# ============================================================
# 用法: Rscript deg_deseq2.R \
#   --counts    <counts.xls> \
#   --sample    <sample_sheet.csv> \
#   --outdir    <输出目录> \
#   --project   <项目ID> \
#   [--lfc      0.5] \
#   [--pval     0.05] \
#   [--comparison_mode auto] \
#   [--comparisons "AB3_vs_SC,AB7_vs_SC,All_vs_SC"]
# ============================================================

suppressPackageStartupMessages(library(DESeq2))
suppressPackageStartupMessages(library(dplyr))
suppressPackageStartupMessages(library(tibble))
suppressPackageStartupMessages(library(ggplot2))
suppressPackageStartupMessages(library(pheatmap))
suppressPackageStartupMessages(library(ggrepel))

# ---- 解析参数 ----
args <- commandArgs(trailingOnly = TRUE)
if (length(args) < 8) {
  stop("Usage: Rscript deg_deseq2.R --counts <file> --sample <file> --outdir <dir> --project <id> [--lfc 0.5] [--pval 0.05] [--comparison_mode auto]")
}

counts_file <- args[which(args == "--counts") + 1]
sample_file <- args[which(args == "--sample") + 1]
outdir      <- args[which(args == "--outdir") + 1]
proj_id     <- args[which(args == "--project") + 1]

lfc_cut     <- ifelse("--lfc" %in% args, as.numeric(args[which(args == "--lfc") + 1]), 0.5)
p_cut       <- ifelse("--pval" %in% args, as.numeric(args[which(args == "--pval") + 1]), 0.05)
comp_mode   <- ifelse("--comparison_mode" %in% args, args[which(args == "--comparison_mode") + 1], "auto")
manual_comps_raw <- ifelse("--comparisons" %in% args, args[which(args == "--comparisons") + 1], "")
use_manual_comps  <- nchar(manual_comps_raw) > 0

dir.create(outdir, showWarnings = FALSE, recursive = TRUE)

# ---- 读取数据 ----
counts <- read.table(counts_file, header = TRUE, sep = "\t", row.names = 1, check.names = FALSE)
if ("Length" %in% colnames(counts)) {
  counts <- counts[, !colnames(counts) %in% "Length", drop = FALSE]
}

samples <- read.csv(sample_file, stringsAsFactors = FALSE)
original_group_levels <- unique(samples$group)

# 对齐样本
common <- intersect(colnames(counts), samples$sample)
if (length(common) == 0) {
  stop("样本名在 counts 和 sample_sheet 中不匹配！")
}
counts <- counts[, common, drop = FALSE]
samples <- samples[match(common, samples$sample), , drop = FALSE]

samples$group <- factor(samples$group, levels = original_group_levels)
group_levels <- levels(samples$group)

cat(sprintf("样本数: %d\n", ncol(counts)))
cat(sprintf("基因数: %d\n", nrow(counts)))
cat(sprintf("分组: %s\n", paste(group_levels, collapse = ", ")))
cat("分组信息:\n"); print(table(samples$group))

# ---- 批次检测 ----
has_batch <- "batch" %in% colnames(samples)
if (has_batch) {
  samples$batch <- factor(samples$batch)
  cat(sprintf("检测到 batch 列，批次: %s\n", paste(levels(samples$batch), collapse = ", ")))
  cat("批次分布:\n"); print(table(samples$batch, samples$group))
  design_formula <- as.formula("~ batch + group")
} else {
  design_formula <- as.formula("~ group")
}
cat(sprintf("DESeq2 设计公式: %s\n", deparse(design_formula)))

# ---- 构建比较列表 ----
comparisons <- list()

if (use_manual_comps) {
  # ---- 手动模式: 解析逗号分隔的比较标签 ----
  comp_labels <- trimws(unlist(strsplit(manual_comps_raw, ",")))
  comp_labels <- comp_labels[comp_labels != ""]
  cat(sprintf("手动指定的比较 (%d 个): %s\n", length(comp_labels),
              paste(comp_labels, collapse = ", ")))

  for (label in comp_labels) {
    parts <- strsplit(label, "_vs_")[[1]]
    if (length(parts) != 2) {
      stop(sprintf("无法解析比较标签: '%s'，格式应为 {分子}_{vs}_{分母}", label))
    }
    kd   <- parts[1]
    ctrl <- parts[2]
    combined <- (kd == "All" || ctrl == "All")

    # 校验非 All 组名是否存在于 sample_sheet
    check_groups <- c(kd, ctrl)
    check_groups <- check_groups[check_groups != "All"]
    for (g in check_groups) {
      if (!(g %in% group_levels)) {
        stop(sprintf("手动比较 '%s' 中的组 '%s' 不在 sample_sheet 中！已知分组: %s",
                     label, g, paste(group_levels, collapse = ", ")))
      }
    }

    if (combined && !("All" %in% c(kd, ctrl))) {
      stop(sprintf("组合比较 '%s' 必须包含 'All' 标签在 _vs_ 的某一侧", label))
    }

    comparisons[[label]] <- list(kd = kd, ctrl = ctrl, combined = combined)
  }

} else {
  # ---- 自动模式 ----
  if (comp_mode == "auto") {
    comp_mode <- if (length(group_levels) == 2) "pairwise" else "control_vs_rest"
  }
  cat(sprintf("比较模式: %s\n", comp_mode))

  if (comp_mode == "pairwise") {
    ctrl <- group_levels[1]
    kd   <- group_levels[2]
    label <- paste0(kd, "_vs_", ctrl)
    comparisons[[label]] <- list(kd = kd, ctrl = ctrl, combined = FALSE)

  } else if (comp_mode == "control_vs_rest") {
    ctrl <- group_levels[1]
    for (kd in group_levels[-1]) {
      label <- paste0(kd, "_vs_", ctrl)
      comparisons[[label]] <- list(kd = kd, ctrl = ctrl, combined = FALSE)
    }
    if (length(group_levels) > 2) {
      label <- paste0("All_vs_", ctrl)
      comparisons[[label]] <- list(kd = "All", ctrl = ctrl, combined = TRUE)
    }

  } else if (comp_mode == "all_vs_all") {
    for (i in 1:(length(group_levels) - 1)) {
      for (j in (i + 1):length(group_levels)) {
        label <- paste0(group_levels[j], "_vs_", group_levels[i])
        comparisons[[label]] <- list(kd = group_levels[j], ctrl = group_levels[i], combined = FALSE)
      }
    }
  }
}

# ---- DESeq2 ----
dds <- DESeqDataSetFromMatrix(
  countData = round(counts),
  colData   = samples,
  design    = design_formula
)

keep <- rowSums(counts(dds) >= 10) >= 1
dds <- dds[keep, ]
cat(sprintf("过滤后基因数: %d\n", nrow(dds)))

dds <- DESeq(dds)

# 为 combined 比较准备重标分组后的 DESeq2 对象（支持 All 在分子或分母）
dds_combined_map <- list()
has_combined <- any(sapply(comparisons, function(x) x$combined))

if (has_combined) {
  for (cn in names(comparisons)) {
    comp <- comparisons[[cn]]
    if (!comp$combined) next

    # 非 All 的组名（作为 factor 的 reference level）
    ref_group <- if (comp$kd == "All") comp$ctrl else comp$kd

    if (!is.null(dds_combined_map[[ref_group]])) next  # 已构建

    if (!(ref_group %in% group_levels)) {
      stop(sprintf("组合比较 '%s' 中的参考组 '%s' 不在样本分组中！", cn, ref_group))
    }

    samples_combined <- samples
    samples_combined$group_combined <- ifelse(
      as.character(samples_combined$group) == ref_group, ref_group, "All"
    )
    samples_combined$group_combined <- factor(samples_combined$group_combined,
                                              levels = c(ref_group, "All"))

    if (has_batch) {
      samples_combined$batch <- factor(samples_combined$batch)
      dds_c <- DESeqDataSetFromMatrix(
        countData = round(counts),
        colData   = samples_combined,
        design    = ~ batch + group_combined
      )
    } else {
      dds_c <- DESeqDataSetFromMatrix(
        countData = round(counts),
        colData   = samples_combined,
        design    = ~ group_combined
      )
    }
    keep_c <- rowSums(counts(dds_c) >= 10) >= 1
    dds_c <- dds_c[keep_c, ]
    dds_c <- suppressMessages(DESeq(dds_c))
    dds_combined_map[[ref_group]] <- dds_c
    cat(sprintf("已为组合比较构建合并 DESeq2 对象 (参考组: %s)\n", ref_group))
  }
}

# ---- 逐比较分析 ----
all_res_list <- list()
all_degs <- c()

for (comp_name in names(comparisons)) {
  comp <- comparisons[[comp_name]]
  cat(sprintf("\n=== 比较: %s ===\n", comp_name))

  if (comp$combined) {
    ref_group <- if (comp$kd == "All") comp$ctrl else comp$kd
    dds_use <- dds_combined_map[[ref_group]]
    # All 在分子侧: contrast = (All, ref), All 在分母侧: contrast = (ref, All)
    if (comp$kd == "All") {
      res <- results(dds_use, contrast = c("group_combined", "All", comp$ctrl))
    } else {
      res <- results(dds_use, contrast = c("group_combined", comp$kd, "All"))
    }
  } else {
    res <- results(dds, contrast = c("group", comp$kd, comp$ctrl))
  }

  res_df <- as.data.frame(res) %>%
    rownames_to_column("Symbol") %>%
    mutate(
      comparison = comp_name,
      change = case_when(
        padj < p_cut & log2FoldChange > lfc_cut  ~ "UP",
        padj < p_cut & log2FoldChange < -lfc_cut ~ "DOWN",
        TRUE ~ "NOT"
      )
    ) %>%
    arrange(padj)

  n_up <- sum(res_df$change == "UP")
  n_down <- sum(res_df$change == "DOWN")
  cat(sprintf("上调基因: %d, 下调基因: %d\n", n_up, n_down))

  all_res_list[[comp_name]] <- res_df
  all_degs <- union(all_degs, res_df$Symbol[res_df$change != "NOT"])

  # 保存 per-comparison CSV
  csv_file <- file.path(outdir, paste0(proj_id, "_", comp_name, "_Diff.csv"))
  write.csv(res_df, csv_file, row.names = FALSE)
  cat(sprintf("  ✓ CSV → %s\n", csv_file))

  # ---- 火山图 ----
  kd_label <- comp$kd
  ctrl_label <- comp$ctrl

  top10 <- res_df %>%
    filter(change != "NOT") %>%
    slice_min(padj, n = 10)

  p_volcano <- ggplot(res_df, aes(x = log2FoldChange, y = -log10(padj), color = change)) +
    geom_point(alpha = 0.4, size = 1.5) +
    scale_color_manual(values = c("UP" = "#d73027", "DOWN" = "#4575b4", "NOT" = "grey80")) +
    geom_vline(xintercept = c(-lfc_cut, lfc_cut), lty = 4, col = "black", lwd = 0.8) +
    geom_hline(yintercept = -log10(p_cut), lty = 4, col = "black", lwd = 0.8) +
    geom_text_repel(data = top10, aes(label = Symbol), size = 3, max.overlaps = 50, box.padding = 0.5) +
    theme_bw() +
    labs(x = "log2(Fold Change)", y = "-log10(Adjusted P-value)")

  volcano_base <- file.path(outdir, paste0("Fig_Volcano_", comp_name))
  ggsave(paste0(volcano_base, ".pdf"), p_volcano, width = 8, height = 6)
  ggsave(paste0(volcano_base, ".png"), p_volcano, width = 8, height = 6, dpi = 300)
  cat("  ✓ 火山图已保存 (PDF + PNG)\n")

  # ---- 热图 ----
  # 上下调各取 top 10，仅展示参与比较的分组，不做样本聚类
  up_genes <- res_df %>% filter(change == "UP") %>% slice_min(padj, n = 10) %>% pull(Symbol)
  down_genes <- res_df %>% filter(change == "DOWN") %>% slice_min(padj, n = 10) %>% pull(Symbol)
  top_genes <- unique(c(up_genes, down_genes))

  if (length(top_genes) > 1) {
    vsd <- vst(dds, blind = FALSE)
    vsd_mat <- assay(vsd)

    # 确定参与比较的分组
    if (comp$combined) {
      comp_groups <- group_levels  # 组合比较展示所有组
    } else {
      comp_groups <- c(comp$kd, comp$ctrl)
    }

    # 仅保留参与比较的样本
    keep_samples <- samples$group %in% comp_groups
    vsd_mat_sub <- vsd_mat[, keep_samples, drop = FALSE]
    samples_sub <- samples[keep_samples, , drop = FALSE]
    samples_sub$group <- factor(samples_sub$group, levels = comp_groups)

    common_genes <- intersect(top_genes, rownames(vsd_mat_sub))

    if (length(common_genes) > 1) {
      heat_data <- vsd_mat_sub[common_genes, , drop = FALSE]
      heat_data_scaled <- t(scale(t(heat_data)))

      annotation_col <- data.frame(Group = samples_sub$group)
      rownames(annotation_col) <- colnames(vsd_mat_sub)

      p_heatmap <- pheatmap(heat_data_scaled,
                            show_colnames = TRUE,
                            annotation_col = annotation_col,
                            fontsize = 6,
                            fontsize_col = 5,
                            color = colorRampPalette(c("navy", "white", "firebrick3"))(100),
                            cluster_cols = FALSE,
                            cluster_rows = TRUE,
                            main = "")

      heatmap_width <- max(6, ncol(vsd_mat_sub) * 0.6)
      heatmap_height <- max(6, length(common_genes) * 0.25)

      heatmap_base <- file.path(outdir, paste0("Fig_Heatmap_", comp_name))
      pdf(paste0(heatmap_base, ".pdf"), width = heatmap_width, height = heatmap_height)
      print(p_heatmap)
      dev.off()
      png(paste0(heatmap_base, ".png"), width = heatmap_width, height = heatmap_height, units = "in", res = 300)
      print(p_heatmap)
      dev.off()
      cat("  ✓ 热图已保存 (PDF + PNG)\n")
    }
  }
}

# ---- 合并保存（向后兼容） ----
all_res <- bind_rows(all_res_list)
write.csv(all_res, file.path(outdir, paste0(proj_id, "_All_Diff.csv")), row.names = FALSE)
cat(sprintf("\n合并 DEG 结果已保存，总差异基因数: %d\n", length(all_degs)))

cat(sprintf(">>> DEG 分析完成，结果保存至: %s\n", outdir))
