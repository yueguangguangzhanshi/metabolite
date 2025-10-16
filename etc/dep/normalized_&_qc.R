run_qc_pre_by_user_style <- function(
    metabo.data.fill,
    metabo.data.rsd,
    sample_df,
    ion_type = "POS",
    out_root = "3.QC",
    file_suffix = ""   # ← 新增：文件名后缀，如 "pre" / "norm" / "qcrlsc"
){
  suppressPackageStartupMessages({
    library(dplyr); library(tidyr); library(stringr); library(tibble)
    library(ggplot2); library(plotly); library(ggrepel); library(viridis); library(scales)
  })
  have_psych       <- requireNamespace("psych", quietly = TRUE)
  have_corrplot    <- requireNamespace("corrplot", quietly = TRUE)
  have_highcharter <- requireNamespace("highcharter", quietly = TRUE)
  have_htmlwidgets <- requireNamespace("htmlwidgets", quietly = TRUE)
  have_webshot     <- requireNamespace("webshot", quietly = TRUE)
  have_ggord       <- requireNamespace("ggord", quietly = TRUE)
  
  # 小工具：把后缀拼在文件名末尾（带下划线）
  suf <- if (!is.null(file_suffix) && nzchar(file_suffix)) paste0("_", file_suffix) else ""
  
  # ---- 标准化样本信息 ----
  stopifnot(all(c("File name","Condition","injection.order") %in% colnames(sample_df)))
  sample_df <- sample_df %>%
    mutate(
      `File name` = as.character(`File name`),
      Condition   = Condition,
      injection.order = as.numeric(injection.order)
    )
  keep_samp <- intersect(colnames(metabo.data.fill), sample_df$`File name`)
  if (length(keep_samp) < 3) stop("样本交集少于3，无法生成QC图。")
  metabo.data.fill <- metabo.data.fill[, keep_samp, drop = FALSE]
  sample_df <- sample_df[match(keep_samp, sample_df$`File name`), , drop = FALSE]
  
  # -------- 2. QC相关性 --------
  dir.create(file.path(out_root,"2.correlation"), recursive = TRUE, showWarnings = FALSE)
  qc_names <- sample_df$`File name`[sample_df$Condition == "QC"]
  if (length(qc_names) >= 2) {
    X_qc <- as.data.frame(metabo.data.fill[, qc_names, drop = FALSE])
    X_qc[] <- lapply(X_qc, function(v) suppressWarnings(as.numeric(v)))
    
    if (have_psych) {
      result  <- psych::corr.test(X_qc, method = "pearson", adjust = "none", alpha = .05)
      rmatrix <- result$r; pmatrix <- result$p
    } else {
      rmatrix <- cor(X_qc, use = "pairwise.complete.obs", method = "pearson")
      pmatrix <- matrix(NA_real_, nrow = ncol(X_qc), ncol = ncol(X_qc),
                        dimnames = list(colnames(X_qc), colnames(X_qc)))
    }
    
    colpal <- colorRampPalette(c("darkblue", "white", "red"))(200)
    if (have_corrplot) {
      pdf(file.path(out_root,"2.correlation", sprintf("correlation_%s%s.pdf", ion_type, suf)))
      corrplot::corrplot(rmatrix, method="circle", col=colpal, tl.col="black", type="upper")
      dev.off()
      
      png(file.path(out_root,"2.correlation", sprintf("correlation_%s%s.png", ion_type, suf)))
      corrplot::corrplot(rmatrix, method="circle", col=colpal, tl.col="black", type="upper")
      dev.off()
    } else {
      rlong <- as.data.frame(rmatrix) |>
        tibble::rownames_to_column("row") |>
        tidyr::pivot_longer(-row, names_to="col", values_to="corr")
      p <- ggplot(rlong, aes(col, row, fill=corr)) +
        geom_tile() + scale_fill_gradient2(limits=c(-1,1), midpoint=0) +
        theme_bw(base_size=10) + theme(axis.text.x=element_text(angle=90,vjust=.5,hjust=1))
      ggsave(file.path(out_root,"2.correlation", sprintf("correlation_%s%s.png", ion_type, suf)),
             p, width=6, height=5, dpi=150)
    }
    
    write.csv(rmatrix, file.path(out_root,"2.correlation", sprintf("correlation_%s%s.csv", ion_type, suf)), row.names = TRUE)
    write.csv(pmatrix, file.path(out_root,"2.correlation", sprintf("correlation_pvalue_%s%s.csv", ion_type, suf)), row.names = TRUE)
  } else {
    message("QC样本少于2，跳过相关性图。")
  }
  
  # =========================================================================================
  # 3) RSD —— 仅“累加阶梯图”（X=%RSD，Y=#Peaks）
  # =========================================================================================
  dir.create(file.path(out_root,"3.RSD"), recursive = TRUE, showWarnings = FALSE)
  rsd_list_csv <- file.path(out_root,"3.RSD", sprintf("RSD_%s%s.csv", ion_type, suf))
  rsd_cum_csv  <- file.path(out_root,"3.RSD", sprintf("RSD_cumulative_%s%s.csv", ion_type, suf))
  rsd_cum_png  <- file.path(out_root,"3.RSD", sprintf("RSD_cumulative_%s%s.png", ion_type, suf))
  rsd_cum_pdf  <- file.path(out_root,"3.RSD", sprintf("RSD_cumulative_%s%s.pdf", ion_type, suf))
  
  qc_cols <- sample_df$`File name`[sample_df$Condition=="QC"]
  rsd_perc <- function(v){ mu <- mean(v, na.rm=TRUE); if(!is.finite(mu)||mu<=0) return(NA_real_); 100*stats::sd(v,na.rm=TRUE)/mu }
  
  # 读取 RSD：优先用 metabo.data.rsd 的列（兼容 qc_rsd / QC.nor.rsd）
  qc_rsd_vec <- NULL
  if (!is.null(metabo.data.rsd) && is.data.frame(metabo.data.rsd)) {
    cand <- c("qc_rsd", "QC.nor.rsd", "QC.RSD", "qc.RSD")
    hit  <- intersect(cand, colnames(metabo.data.rsd))
    if (length(hit) >= 1) {
      qc_rsd_vec <- metabo.data.rsd[[ hit[1] ]]
    }
  }
  # 若没有，就从表达矩阵现算
  if (is.null(qc_rsd_vec)) {
    if (length(qc_cols) >= 2) {
      qc_rsd_vec <- apply(metabo.data.fill[, qc_cols, drop=FALSE], 1, rsd_perc)
    } else {
      warning("QC样本少于2，无法计算 RSD。RSD 累加图将为空。")
      qc_rsd_vec <- rep(NA_real_, nrow(metabo.data.fill))
    }
  }
  # 标准化为百分比
  if (is.numeric(qc_rsd_vec) && max(qc_rsd_vec, na.rm = TRUE) <= 1.5) qc_rsd_vec <- qc_rsd_vec * 100
  
  # 保存 RSD 列表
  write.csv(qc_rsd_vec, rsd_list_csv, row.names = FALSE)
  
  # 累加数据
  rsd_valid <- qc_rsd_vec[is.finite(qc_rsd_vec)]
  rsd_cum_df <- tibble(
    rsd_pct   = sort(rsd_valid),
    cum_peaks = seq_along(rsd_valid)
  )
  
  # ---- 关键：补“水平尾巴”到 100 ----
  if (nrow(rsd_cum_df) > 0 && max(rsd_cum_df$rsd_pct) < 100) {
    rsd_cum_df <- dplyr::bind_rows(
      rsd_cum_df,
      tibble::tibble(rsd_pct = 100, cum_peaks = max(rsd_cum_df$cum_peaks))
    )
  }
  
  write.csv(rsd_cum_df, rsd_cum_csv, row.names = FALSE)
  
  # 画累加阶梯图
  p_rsd_cum <- ggplot(rsd_cum_df, aes(x = rsd_pct, y = cum_peaks)) +
    geom_step(direction = "hv") +
    geom_vline(xintercept = c(15, 20, 30), linetype = "dashed", color = "gray40") +
    labs(x = "% RSD", y = "# of Peaks", title = sprintf("QC RSD Cumulative (%s)", ion_type)) +
    theme_bw(base_size = 12)
  ggsave(rsd_cum_pdf, p_rsd_cum, width=7.5, height=5)
  ggsave(rsd_cum_png, p_rsd_cum, width=7.5, height=5, dpi=150)
  
  # -------- 4. PCA --------
  dir.create(file.path(out_root,"4.PCA"), recursive = TRUE, showWarnings = FALSE)
  
  # 1) 转数值 + 统一非有限值为 NA
  pca_qc <- as.data.frame(metabo.data.fill, check.names = FALSE)
  pca_qc[] <- lapply(pca_qc, function(x) suppressWarnings(as.numeric(x)))
  # 把 Inf/-Inf/"NaN" 等全部置为 NA
  pca_qc_mat <- as.matrix(pca_qc)
  pca_qc_mat[!is.finite(pca_qc_mat)] <- NA
  pca_qc <- as.data.frame(pca_qc_mat, check.names = FALSE)
  
  # 2) 按特征（行）用中位数填补 NA（仅用于 PCA，不改原数据）
  feat_median <- apply(pca_qc, 1, function(v) median(v, na.rm = TRUE))
  na_rows <- which(rowSums(is.na(pca_qc)) > 0)
  for (i in na_rows) {
    pca_qc[i, is.na(pca_qc[i, ])] <- feat_median[i]
  }
  
  # 3) 去掉“跨样本方差=0”的特征（否则 scale.=TRUE 会除以 0 → NaN）
  var_feat <- apply(pca_qc, 1, function(v) stats::var(v, na.rm = TRUE))
  keep_feat <- is.finite(var_feat) & var_feat > 0
  removed_feat_n <- sum(!keep_feat)
  pca_qc <- pca_qc[keep_feat, , drop = FALSE]
  
  # 4) （可选）去掉在剩余特征上方差=0的样本
  sd_samp <- apply(pca_qc, 2, function(v) stats::sd(v, na.rm = TRUE))
  keep_samp <- is.finite(sd_samp) & sd_samp > 0
  removed_samp_n <- sum(!keep_samp)
  if (removed_samp_n > 0) {
    pca_qc <- pca_qc[, keep_samp, drop = FALSE]
    # 同步 sample_df（保持分组向量对齐）
    sample_df <- sample_df[match(colnames(pca_qc), sample_df$`File name`), , drop = FALSE]
  }
  
  # 5) 真正做 PCA
  pc.cr <- prcomp(t(pca_qc), scale. = TRUE)
  
  # 分组
  pca_group <- sample_df$Condition
  
  if (have_ggord && have_htmlwidgets) {
    g <- ggord::ggord(pc.cr, pca_group, obslab=FALSE, arrow=NULL, txt=NULL, ellipse=TRUE, alpha=0.5,
                      poly=FALSE, coord_fix=FALSE, facet=FALSE, size=4, repel=TRUE) +
      geom_text(aes(label=lab), vjust=-0.5, hjust="inward", check_overlap=TRUE, size=3, color="black", show.legend=TRUE) +
      theme(plot.title = element_text(hjust = 0.5, size=12), text = element_text(size = 8), panel.grid = element_blank()) +
      labs(title="PCA Plot", fill=NULL, colour=NULL, shape=NULL) +
      scale_color_brewer(palette="Set2") + scale_fill_brewer(palette="Set2")
    
    g_pca <- plotly::ggplotly(g)
    if (length(g_pca$x$data) >= 5) {
      g_pca$x$data[[3]]$showlegend <- FALSE
      g_pca$x$data[[4]]$showlegend <- FALSE
      g_pca$x$data[[5]]$showlegend <- TRUE
      g_pca$x$data[[5]]$name <- "label"
    }
    htmlwidgets::saveWidget(g_pca, file = file.path(out_root,"4.PCA", sprintf("PCA_%s%s.html", ion_type, suf)), selfcontained = TRUE)
    
    g2 <- ggord::ggord(pc.cr, pca_group, obslab=FALSE, arrow=NULL, txt=NULL, ellipse=TRUE, alpha=0.5,
                       poly=FALSE, coord_fix=FALSE, facet=FALSE, size=4, repel=TRUE) +
      ggrepel::geom_text_repel(aes(label=lab)) +
      theme(plot.title = element_text(hjust = 0.5, size=12), text = element_text(size = 8), panel.grid = element_blank()) +
      labs(title="PCA Plot", fill=NULL, colour=NULL, shape=NULL) +
      scale_color_brewer(palette="Set2") + scale_fill_brewer(palette="Set2")
    
    ggsave(file.path(out_root,"4.PCA", sprintf("PCA_%s%s.png", ion_type, suf)), g2, width=7, height=5.5, dpi=150)
    ggsave(file.path(out_root,"4.PCA", sprintf("PCA_%s%s.pdf", ion_type, suf)), g2, width=7, height=5.5)
  } else {
    pca_df <- as.data.frame(pc.cr$x[,1:2]) |>
      mutate(Sample = rownames(pc.cr$x), Group = pca_group)
    gp <- ggplot(pca_df, aes(PC1, PC2, color=Group)) + geom_point(size=2) +
      theme_bw(base_size=12) + labs(title="PCA Plot")
    ggsave(file.path(out_root,"4.PCA", sprintf("PCA_%s%s.png", ion_type, suf)), gp, width=7, height=5.5, dpi=150)
  }
  
  # =========================================================================================
  # 5) CV 累积分布：单图（0–1 全覆盖，突出 0–0.1 / 0.1–0.3 / 0.3–1）
  # =========================================================================================
  dir.create(file.path(out_root,"5.CV"), recursive = TRUE, showWarnings = FALSE)
  metabo.data.cv <- metabo.data.fill
  cv_group <- data.frame(Sample = sample_df$`File name`, Group = sample_df$Condition, stringsAsFactors = FALSE)
  
  df_long <- metabo.data.cv %>%
    tibble::rownames_to_column("Metabolite") %>%
    tidyr::pivot_longer(-Metabolite, names_to = "Sample", values_to = "Intensity") %>%
    dplyr::left_join(cv_group, by = "Sample")
  
  cv_data <- df_long %>%
    dplyr::group_by(Metabolite, Group) %>%
    dplyr::summarise(CV = sd(Intensity, na.rm=TRUE) / mean(Intensity, na.rm=TRUE), .groups = "drop") %>%
    dplyr::filter(is.finite(CV))
  
  cumulative_plot <- cv_data %>%
    dplyr::group_by(Group) %>%
    dplyr::arrange(CV, .by_group = TRUE) %>%
    dplyr::mutate(Cumulative = seq_along(CV)/n()) %>%
    dplyr::ungroup()
  
  cumulative_plot_ext <- cumulative_plot %>%
    dplyr::group_split(Group) %>%
    lapply(function(df) {
      if (nrow(df) > 0 && max(df$CV, na.rm = TRUE) < 1) {
        dplyr::bind_rows(
          df,
          tibble::tibble(Metabolite = NA_character_, Group = df$Group[1], CV = 1.0, Cumulative = 1.0)
        )
      } else df
    }) %>% dplyr::bind_rows()
  
  # 背景区间带：0–0.1 / 0.1–0.3 / 0.3–1
  seg_bands <- tibble::tibble(
    xmin = c(0.00, 0.10, 0.30),
    xmax = c(0.10, 0.30, 1.00),
    seg  = factor(c("0–0.1","0.1–0.3","0.3–1"), levels = c("0–0.1","0.1–0.3","0.3–1"))
  )
  
  p_cv <- ggplot(cumulative_plot_ext, aes(x = CV, y = Cumulative, color = Group)) +
    geom_rect(data = seg_bands, inherit.aes = FALSE,
              aes(xmin = xmin, xmax = xmax, ymin = -Inf, ymax = Inf, fill = seg),
              alpha = 0.08, color = NA) +
    scale_fill_manual(values = c("0–0.1" = "#ffd166", "0.1–0.3" = "#06d6a0", "0.3–1" = "#118ab2"), guide = "none") +
    geom_step(direction = "hv", linewidth = 0.9) +
    geom_vline(xintercept = c(0.1, 0.3), linetype = "dashed", color = "gray40") +
    scale_x_continuous(
      limits = c(0, 1), expand = c(0, 0),
      breaks = c(0, 0.05, 0.10, 0.20, 0.30, 0.50, 1.00),
      minor_breaks = c(seq(0, 0.30, by = 0.01), seq(0.30, 1.00, by = 0.05))
    ) +
    scale_y_continuous(labels = scales::percent_format(), breaks = c(0, 0.25, 0.5, 0.75, 1)) +
    labs(x = "Coefficient of Variation (CV)", y = "Cumulative Proportion",
         color = "Sample Group", title = sprintf("CV Cumulative (%s)", ion_type)) +
    theme_minimal(base_size = 12) +
    theme(legend.position = "top", panel.grid.minor = element_line(linewidth = 0.2, colour = "grey90"))
  
  ggsave(file.path(out_root,"5.CV", sprintf("cv_%s%s.pdf", ion_type, suf)), p_cv, width=8, height=5.5)
  ggsave(file.path(out_root,"5.CV", sprintf("cv_%s%s.png", ion_type, suf)), p_cv, width=8, height=5.5, dpi=150)
  
  invisible(list(
    qc_n = length(qc_names),
    paths = list(
      corr_dir = file.path(out_root,"2.correlation"),
      rsd_dir  = file.path(out_root,"3.RSD"),
      pca_dir  = file.path(out_root,"4.PCA"),
      cv_dir   = file.path(out_root,"5.CV")
    )
  ))
}

# 归一化前
run_qc_pre_by_user_style(
  metabo.data.fill = metabo.data.fill.pre,
  metabo.data.rsd  = metabo.data.rsd,                  # 若含 qc_rsd 列则用你的高图流程
  sample_df        = sample.data[[ion_type]],
  ion_type         = ion_type,
  out_root         = "3.QC",
  file_suffix      = "pre"                             # ← 关键
)

run_qc_pre_by_user_style(
  metabo.data.fill = metabo.data.fill,                 # 归一化后的矩阵
  metabo.data.rsd  = metabo.data.rsd,            # 对应 RSD 表（可同结构）
  sample_df        = sample.data[[ion_type]],
  ion_type         = ion_type,
  out_root         = "3.QC",
  file_suffix      = "norm_svr"                            # ← 关键
)

# 归一化后（比如 qcrlsc 或 svr）
run_qc_pre_by_user_style(
  metabo.data.fill = metabo.data.fill,                 # 归一化后的矩阵
  metabo.data.rsd  = metabo.data.rsd,            # 对应 RSD 表（可同结构）
  sample_df        = sample.data[[ion_type]],
  ion_type         = ion_type,
  out_root         = "3.QC",
  file_suffix      = "norm_rlsc"                            # ← 关键
)

run_qc_pre_by_user_style(
  metabo.data.fill = metabo.data.fill,                 # 归一化后的矩阵
  metabo.data.rsd  = metabo.data.rsd,            # 对应 RSD 表（可同结构）
  sample_df        = sample.data[[ion_type]],
  ion_type         = ion_type,
  out_root         = "3.QC",
  file_suffix      = "norm_normae"                            # ← 关键
)
