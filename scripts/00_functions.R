# ============================================================
# 00_functions.R — Shared utility functions for PST analysis
# ============================================================

suppressPackageStartupMessages({
  library(readxl)
  library(dplyr)
  library(tidyr)
  library(tibble)
  library(stringr)
  library(ggplot2)
  library(ggrepel)
  library(patchwork)
  library(pls)
})

# Color palette (colorblind-friendly)
group_colors <- c(
  "0-mo"  = "#E64B35",
  "6-mo"  = "#4DBBD5",
  "QC"    = "#00A087",
  "Blank" = "#B09C85"
)

# Group shapes
group_shapes <- c(
  "0-mo"  = 16,   # filled circle
  "6-mo"  = 16,
  "QC"    = 17,   # filled triangle
  "Blank" = 15    # filled square
)


# ----------------------------------------------------------
# Parse MS-DIAL Area export xlsx
# Returns: list(metadata, intensity, class_info)
# ----------------------------------------------------------
parse_msdial_export <- function(filepath) {
  # Read all data as text to preserve structure
  # Row 1: Class, Row 2: File type, Row 3: Injection order, Row 4: Batch ID
  # Row 5: Column names, Row 6+: Data
  #
  # IMPORTANT: read_excel with n_max drops columns that are entirely NA,
  # so we must read the full file at once to preserve column alignment.

  # Read first 5 rows to get headers (include row 5 = column names)
  header_block <- read_excel(filepath, sheet = 1, col_names = FALSE,
                             n_max = 5, .name_repair = "minimal")
  header_block <- as.data.frame(header_block)
  total_cols <- ncol(header_block)

  # Extract header rows
  class_row <- as.character(header_block[1, ])   # Row 1: Class
  batch_row <- as.character(header_block[4, ])   # Row 4: Batch ID
  col_names <- as.character(header_block[5, ])   # Row 5: Column names

  # Find the boundary: "Class" in row 1 marks the start of sample info
  class_label_idx <- which(class_row == "Class")
  if (length(class_label_idx) == 0) {
    stop("Cannot find 'Class' label in row 1")
  }
  meta_end <- class_label_idx[1]  # "Class" column is the last metadata col
  sample_start <- meta_end + 1

  # Identify summary columns (last columns with "Average"/"Stdev" in batch row)
  summary_idx <- which(batch_row %in% c("Average", "Stdev"))
  if (length(summary_idx) > 0) {
    data_end <- min(summary_idx) - 1
  } else {
    data_end <- total_cols
  }

  # Sample column indices
  sample_idx <- sample_start:data_end
  n_samples <- length(sample_idx)

  cat(sprintf("  Metadata columns: 1-%d (%d cols)\n", meta_end, meta_end))
  cat(sprintf("  Sample columns: %d-%d (%d cols)\n", sample_start, data_end, n_samples))
  cat(sprintf("  Summary columns excluded: %d cols\n", length(summary_idx)))

  # Read the actual data (skip 5 header rows)
  data_raw <- read_excel(filepath, sheet = 1, col_names = FALSE,
                         skip = 5, .name_repair = "minimal")
  data_raw <- as.data.frame(data_raw)

  # Separate metadata and intensity
  metadata <- data_raw[, 1:meta_end, drop = FALSE]
  colnames(metadata) <- col_names[1:meta_end]

  intensity <- data_raw[, sample_idx, drop = FALSE]
  # Convert to numeric
  intensity <- as.data.frame(lapply(intensity, function(x) as.numeric(as.character(x))))
  colnames(intensity) <- col_names[sample_idx]

  # Class info for sample columns
  sample_classes <- class_row[sample_idx]
  names(sample_classes) <- col_names[sample_idx]

  list(
    metadata = metadata,
    intensity = intensity,
    class_info = sample_classes
  )
}


# ----------------------------------------------------------
# Load sample grouping file
# ----------------------------------------------------------
load_grouping <- function(filepath) {
  grp <- read_excel(filepath, col_names = TRUE)
  grp <- as.data.frame(grp)
  colnames(grp) <- c("Prefix", "Sample_number", "Sample_name", "Group")

  # Reassign 12-mo to 6-mo
  grp$Group[grp$Group == "12-mo"] <- "6-mo"

  grp
}


# ----------------------------------------------------------
# Match MS-DIAL columns to grouping file
# Returns a data.frame with: col_name, sample_number, sample_name, group
# ----------------------------------------------------------
match_samples <- function(col_names, class_info, grouping_df) {
  # Extract sample number from MS-DIAL column names
  # Sample numbers: B001-B168, QC01-QC10, BK##_## patterns
  # Strategy: try to match each grouping sample number within the column name

  result <- data.frame(
    col_name = col_names,
    class = class_info,
    sample_number = NA_character_,
    sample_name = NA_character_,
    group = NA_character_,
    stringsAsFactors = FALSE
  )

  # For each column, find matching sample number from grouping file
  for (i in seq_along(col_names)) {
    cn <- col_names[i]
    for (j in seq_len(nrow(grouping_df))) {
      sn <- grouping_df$Sample_number[j]
      # Check if column name ends with the sample number (possibly after underscore/prefix)
      pattern <- paste0("(^|_)", sn, "$")
      if (grepl(pattern, cn)) {
        result$sample_number[i] <- sn
        result$sample_name[i] <- grouping_df$Sample_name[j]
        result$group[i] <- grouping_df$Group[j]
        break
      }
    }
  }

  # For unmatched columns, use class_info from Row 1 as fallback
  unmatched <- is.na(result$group)
  if (any(unmatched)) {
    valid_classes <- c("0-mo", "6-mo", "12-mo", "QC", "Blank")
    for (i in which(unmatched)) {
      cls <- result$class[i]
      if (!is.na(cls) && cls %in% valid_classes) {
        result$group[i] <- ifelse(cls == "12-mo", "6-mo", cls)
        result$sample_name[i] <- result$col_name[i]
        cat(sprintf("  Fallback: %s assigned to group '%s' via Class info\n",
                    result$col_name[i], result$group[i]))
      }
    }
    still_unmatched <- is.na(result$group)
    if (any(still_unmatched)) {
      cat(sprintf("  Unmatched columns: %d (Test/NA, will be excluded)\n",
                  sum(still_unmatched)))
    }
  }

  result
}


# ----------------------------------------------------------
# Blank filtering
# Remove features where mean(Blank)/mean(QC) > 0.3
# Also remove features where mean(Blank) > mean(all Samples)
# ----------------------------------------------------------
blank_filter <- function(intensity, sample_info) {
  blank_cols <- sample_info$col_name[sample_info$group == "Blank"]
  qc_cols <- sample_info$col_name[sample_info$group == "QC"]
  sample_cols <- sample_info$col_name[sample_info$group %in% c("0-mo", "6-mo")]

  blank_mat <- intensity[, blank_cols, drop = FALSE]
  qc_mat <- intensity[, qc_cols, drop = FALSE]
  sample_mat <- intensity[, sample_cols, drop = FALSE]

  mean_blank <- rowMeans(blank_mat, na.rm = TRUE)
  mean_qc <- rowMeans(qc_mat, na.rm = TRUE)
  mean_sample <- rowMeans(sample_mat, na.rm = TRUE)

  # Avoid division by zero
  ratio <- ifelse(mean_qc == 0, Inf, mean_blank / mean_qc)

  keep <- (ratio <= 0.3) & (mean_blank <= mean_sample)

  cat(sprintf("  Blank filter: %d -> %d features (removed %d)\n",
              nrow(intensity), sum(keep), sum(!keep)))

  keep
}


# ----------------------------------------------------------
# Missing value assessment
# ----------------------------------------------------------
missing_value_filter <- function(intensity, sample_info) {
  qc_cols <- sample_info$col_name[sample_info$group == "QC"]
  group_0mo_cols <- sample_info$col_name[sample_info$group == "0-mo"]
  group_6mo_cols <- sample_info$col_name[sample_info$group == "6-mo"]

  n_features <- nrow(intensity)

  # Detection rate in QC (non-zero proportion)
  qc_mat <- intensity[, qc_cols, drop = FALSE]
  qc_detect <- rowMeans(qc_mat > 0 & !is.na(qc_mat), na.rm = TRUE)
  keep_qc <- qc_detect >= 0.8

  cat(sprintf("  QC detection rate filter (>=80%%): %d -> %d features (removed %d)\n",
              n_features, sum(keep_qc), sum(!keep_qc)))

  # Detection rate per group (only check 0-mo and 6-mo)
  mat_0mo <- intensity[, group_0mo_cols, drop = FALSE]
  mat_6mo <- intensity[, group_6mo_cols, drop = FALSE]

  detect_0mo <- rowMeans(mat_0mo > 0 & !is.na(mat_0mo), na.rm = TRUE)
  detect_6mo <- rowMeans(mat_6mo > 0 & !is.na(mat_6mo), na.rm = TRUE)

  # Remove features with detection rate < 50% in ALL sample groups
  keep_group <- (detect_0mo >= 0.5) | (detect_6mo >= 0.5)

  cat(sprintf("  Group detection rate filter (>=50%% in any group): %d -> %d features (removed %d)\n",
              sum(keep_qc), sum(keep_qc & keep_group), sum(keep_qc & !keep_group)))

  keep_qc & keep_group
}


# ----------------------------------------------------------
# QC RSD filtering (on RAW intensity, no normalization)
# ----------------------------------------------------------
qc_rsd_filter <- function(intensity, sample_info, rsd_threshold = 30) {
  qc_cols <- sample_info$col_name[sample_info$group == "QC"]
  qc_mat <- intensity[, qc_cols, drop = FALSE]

  qc_mean <- rowMeans(qc_mat, na.rm = TRUE)
  qc_sd <- apply(qc_mat, 1, sd, na.rm = TRUE)

  # RSD = CV% = SD/mean * 100
  rsd <- ifelse(qc_mean == 0, Inf, qc_sd / qc_mean * 100)

  keep <- rsd <= rsd_threshold

  cat(sprintf("  QC RSD filter (<=%.0f%%): %d -> %d features (removed %d)\n",
              rsd_threshold, nrow(intensity), sum(keep), sum(!keep)))

  list(keep = keep, rsd = rsd)
}


# ----------------------------------------------------------
# Missing value imputation (half-minimum per feature)
# ----------------------------------------------------------
impute_half_min <- function(intensity) {
  for (j in seq_len(ncol(intensity))) {
    col <- intensity[, j]
    missing <- is.na(col) | col == 0
    if (any(missing)) {
      positive_vals <- col[!missing & col > 0]
      if (length(positive_vals) > 0) {
        half_min <- min(positive_vals) / 2
      } else {
        half_min <- 1  # fallback if all zeros
      }
      intensity[missing, j] <- half_min
    }
  }
  intensity
}


# ----------------------------------------------------------
# PCA score plot for QC assessment
# ----------------------------------------------------------
plot_pca_qc <- function(intensity, sample_info, output_path) {
  # Use samples + QC only (exclude Blank for cleaner plot)
  keep_groups <- c("0-mo", "6-mo", "QC")
  keep_cols <- sample_info$col_name[sample_info$group %in% keep_groups]
  keep_groups_vec <- sample_info$group[sample_info$group %in% keep_groups]

  mat <- t(as.matrix(intensity[, keep_cols, drop = FALSE]))

  # Log2 transform for PCA visualization (add 1 to avoid log(0))
  mat_log <- log2(mat + 1)

  # Center and scale
  mat_scaled <- scale(mat_log, center = TRUE, scale = TRUE)
  # Remove constant columns (zero variance)
  zero_var <- apply(mat_scaled, 2, function(x) all(is.na(x)))
  mat_scaled <- mat_scaled[, !zero_var, drop = FALSE]
  # Replace any remaining NAs with 0
  mat_scaled[is.na(mat_scaled)] <- 0

  pca_res <- prcomp(mat_scaled, center = FALSE, scale. = FALSE)

  # Variance explained
  var_explained <- summary(pca_res)$importance[2, ] * 100

  pca_df <- data.frame(
    PC1 = pca_res$x[, 1],
    PC2 = pca_res$x[, 2],
    Group = factor(keep_groups_vec, levels = c("0-mo", "6-mo", "QC"))
  )

  p <- ggplot(pca_df, aes(x = PC1, y = PC2, color = Group, shape = Group)) +
    geom_point(size = 3, alpha = 0.7) +
    scale_color_manual(values = group_colors) +
    scale_shape_manual(values = group_shapes) +
    labs(
      x = sprintf("PC1 (%.1f%%)", var_explained[1]),
      y = sprintf("PC2 (%.1f%%)", var_explained[2]),
      title = "PCA Score Plot — QC Assessment"
    ) +
    theme_bw() +
    theme(
      plot.title = element_text(hjust = 0.5, face = "bold"),
      legend.position = "right"
    )

  pdf(output_path, width = 8, height = 6)
  print(p)
  dev.off()
  cat(sprintf("  Saved: %s\n", output_path))

  invisible(p)
}


# ----------------------------------------------------------
# QC RSD histogram
# ----------------------------------------------------------
plot_rsd_histogram <- function(rsd_values, rsd_threshold, output_path) {
  rsd_df <- data.frame(RSD = rsd_values[is.finite(rsd_values)])

  p <- ggplot(rsd_df, aes(x = RSD)) +
    geom_histogram(binwidth = 2, fill = "#4DBBD5", color = "white", alpha = 0.8) +
    geom_vline(xintercept = rsd_threshold, linetype = "dashed", color = "#E64B35",
               linewidth = 1) +
    annotate("text", x = rsd_threshold + 2, y = Inf, vjust = 2, hjust = 0,
             label = sprintf("Threshold = %d%%", rsd_threshold),
             color = "#E64B35", fontface = "bold") +
    labs(
      x = "QC RSD (%)",
      y = "Number of Features",
      title = "Distribution of QC RSD Values (Raw Intensity)"
    ) +
    theme_bw() +
    theme(plot.title = element_text(hjust = 0.5, face = "bold"))

  pdf(output_path, width = 8, height = 6)
  print(p)
  dev.off()
  cat(sprintf("  Saved: %s\n", output_path))

  invisible(p)
}


# ============================================================
# Phase 2/3 Functions
# ============================================================

# ----------------------------------------------------------
# Data normalization: sum-norm -> log2 -> Pareto scaling
# Returns list(normalized, log2_data, scaled)
#   normalized: sum-normalized (for FC calculation reference)
#   log2_data:  sum-normalized + log2 (for FC and univariate)
#   scaled:     sum-normalized + log2 + Pareto (for PLS-DA)
# ----------------------------------------------------------
normalize_data <- function(intensity, sample_info) {
  sample_groups <- c("0-mo", "6-mo")
  sample_cols <- sample_info$col_name[sample_info$group %in% sample_groups]
  mat <- as.matrix(intensity[, sample_cols, drop = FALSE])

  # 1. Sum normalization (per sample = per column)
  col_sums <- colSums(mat, na.rm = TRUE)
  median_sum <- median(col_sums)
  mat_norm <- sweep(mat, 2, col_sums, "/") * median_sum

  # 2. Log2 transformation
  mat_log2 <- log2(mat_norm)

  # 3. Pareto scaling (per feature = per row): mean-center, divide by sqrt(SD)
  row_means <- rowMeans(mat_log2, na.rm = TRUE)
  row_sds <- apply(mat_log2, 1, sd, na.rm = TRUE)
  row_sds[row_sds == 0] <- 1  # avoid division by zero
  mat_scaled <- sweep(mat_log2, 1, row_means, "-")
  mat_scaled <- sweep(mat_scaled, 1, sqrt(row_sds), "/")

  groups <- sample_info$group[sample_info$group %in% sample_groups]

  list(
    normalized = mat_norm,
    log2_data  = mat_log2,
    scaled     = mat_scaled,
    groups     = factor(groups, levels = c("0-mo", "6-mo")),
    col_names  = sample_cols
  )
}


# ----------------------------------------------------------
# PLS-DA using pls::plsr with dummy Y
# Returns list(model, scores, vip, r2, q2)
# ----------------------------------------------------------
run_plsda <- function(scaled_mat, groups, ncomp = 2) {
  X <- t(scaled_mat)  # samples x features
  Y <- ifelse(groups == "6-mo", 1, 0)

  df <- data.frame(Y = I(matrix(Y, ncol = 1)), X = I(X))
  model <- plsr(Y ~ X, data = df, ncomp = ncomp, validation = "LOO",
                method = "oscorespls")

  scores <- model$scores[, 1:ncomp, drop = FALSE]
  colnames(scores) <- paste0("Comp", 1:ncomp)

  # VIP calculation
  vip <- compute_vip(model)

  # R2 and Q2 from cross-validation
  r2_vals <- R2(model)$val[1, 1, ]  # intercept + ncomp values
  q2_vals <- 1 - MSEP(model)$val[1, 1, ] / var(Y)

  list(
    model  = model,
    scores = scores,
    vip    = vip,
    R2     = r2_vals,
    Q2     = q2_vals,
    groups = groups
  )
}


# ----------------------------------------------------------
# VIP (Variable Importance in Projection) calculation
# ----------------------------------------------------------
compute_vip <- function(model) {
  W <- model$loading.weights  # p x ncomp
  T_mat <- model$scores       # n x ncomp
  Q <- model$Yloadings        # 1 x ncomp

  ncomp <- ncol(W)
  p <- nrow(W)

  # SS for each component
  ss <- numeric(ncomp)
  for (a in seq_len(ncomp)) {
    ss[a] <- (t(T_mat[, a]) %*% T_mat[, a]) * (Q[1, a]^2)
  }
  ss_total <- sum(ss)

  vip <- numeric(p)
  for (j in seq_len(p)) {
    weighted_sum <- sum(W[j, ]^2 * ss)
    vip[j] <- sqrt(p * weighted_sum / ss_total)
  }

  vip
}


# ----------------------------------------------------------
# PLS-DA permutation test
# ----------------------------------------------------------
run_permutation_test <- function(scaled_mat, groups, ncomp = 2,
                                 n_perm = 1000) {
  X <- t(scaled_mat)
  Y <- ifelse(groups == "6-mo", 1, 0)

  # Original model Q2
  df <- data.frame(Y = I(matrix(Y, ncol = 1)), X = I(X))
  orig_model <- plsr(Y ~ X, data = df, ncomp = ncomp, validation = "LOO",
                     method = "oscorespls")
  orig_q2 <- 1 - MSEP(orig_model)$val[1, 1, ncomp + 1] / var(Y)
  orig_r2 <- R2(orig_model)$val[1, 1, ncomp + 1]

  # Permutation
  perm_q2 <- numeric(n_perm)
  perm_r2 <- numeric(n_perm)
  cor_vals <- numeric(n_perm)

  for (i in seq_len(n_perm)) {
    Y_perm <- sample(Y)
    cor_vals[i] <- abs(cor(Y, Y_perm))
    df_perm <- data.frame(Y = I(matrix(Y_perm, ncol = 1)), X = I(X))
    perm_model <- tryCatch({
      plsr(Y ~ X, data = df_perm, ncomp = ncomp, validation = "LOO",
           method = "oscorespls")
    }, error = function(e) NULL)

    if (!is.null(perm_model)) {
      perm_q2[i] <- 1 - MSEP(perm_model)$val[1, 1, ncomp + 1] / var(Y_perm)
      perm_r2[i] <- R2(perm_model)$val[1, 1, ncomp + 1]
    } else {
      perm_q2[i] <- NA
      perm_r2[i] <- NA
    }
  }

  p_value <- mean(perm_q2 >= orig_q2, na.rm = TRUE)

  list(
    orig_r2  = orig_r2,
    orig_q2  = orig_q2,
    perm_r2  = perm_r2,
    perm_q2  = perm_q2,
    cor_vals = cor_vals,
    p_value  = p_value
  )
}


# ----------------------------------------------------------
# PLS-DA score plot
# ----------------------------------------------------------
plot_plsda_scores <- function(plsda_result, title_suffix, output_path,
                              paired_lines = FALSE, subject_ids = NULL) {
  scores_df <- data.frame(
    Comp1 = plsda_result$scores[, 1],
    Comp2 = plsda_result$scores[, 2],
    Group = plsda_result$groups
  )

  model <- plsda_result$model
  xvar <- explvar(model)

  p <- ggplot(scores_df, aes(x = Comp1, y = Comp2, color = Group)) +
    stat_ellipse(level = 0.95, linetype = "dashed", linewidth = 0.5) +
    geom_point(size = 3, alpha = 0.7) +
    scale_color_manual(values = group_colors) +
    labs(
      x = sprintf("Component 1 (%.1f%%)", xvar[1]),
      y = sprintf("Component 2 (%.1f%%)", xvar[2]),
      title = paste("PLS-DA Score Plot —", title_suffix)
    ) +
    theme_bw() +
    theme(plot.title = element_text(hjust = 0.5, face = "bold"))

  # Add paired lines if requested
  if (paired_lines && !is.null(subject_ids)) {
    scores_df$Subject <- subject_ids
    paired_df <- scores_df |>
      select(Subject, Group, Comp1, Comp2) |>
      pivot_wider(names_from = Group, values_from = c(Comp1, Comp2),
                  names_sep = "_")
    if (all(c("Comp1_0-mo", "Comp1_6-mo", "Comp2_0-mo", "Comp2_6-mo") %in%
            colnames(paired_df))) {
      p <- p + geom_segment(
        data = paired_df,
        aes(x = `Comp1_0-mo`, y = `Comp2_0-mo`,
            xend = `Comp1_6-mo`, yend = `Comp2_6-mo`),
        color = "grey60", linewidth = 0.3, alpha = 0.5,
        inherit.aes = FALSE
      )
    }
  }

  pdf(output_path, width = 8, height = 6)
  print(p)
  dev.off()
  cat(sprintf("  Saved: %s\n", output_path))
  invisible(p)
}


# ----------------------------------------------------------
# Permutation test plot
# ----------------------------------------------------------
plot_permutation <- function(perm_result, output_path) {
  df <- data.frame(
    Correlation = perm_result$cor_vals,
    R2 = perm_result$perm_r2,
    Q2 = perm_result$perm_q2
  ) |>
    pivot_longer(cols = c(R2, Q2), names_to = "Metric", values_to = "Value")

  # Add original values
  orig_df <- data.frame(
    Correlation = c(1, 1),
    Metric = c("R2", "Q2"),
    Value = c(perm_result$orig_r2, perm_result$orig_q2)
  )

  p <- ggplot(df, aes(x = Correlation, y = Value, color = Metric)) +
    geom_point(alpha = 0.3, size = 1) +
    geom_point(data = orig_df, size = 4, shape = 18) +
    geom_smooth(method = "lm", se = FALSE, linetype = "dashed", linewidth = 0.5) +
    scale_color_manual(values = c("R2" = "#E64B35", "Q2" = "#4DBBD5")) +
    labs(
      x = "Correlation with Original Y",
      y = "Value",
      title = sprintf("Permutation Test (n=1000, p=%.4f)", perm_result$p_value)
    ) +
    theme_bw() +
    theme(plot.title = element_text(hjust = 0.5, face = "bold"))

  pdf(output_path, width = 8, height = 6)
  print(p)
  dev.off()
  cat(sprintf("  Saved: %s\n", output_path))
  invisible(p)
}


# ----------------------------------------------------------
# Volcano plot
# ----------------------------------------------------------
plot_volcano <- function(stats_df, output_path, title_suffix = "All Samples") {
  stats_df <- stats_df |>
    mutate(
      Significance = case_when(
        FDR < 0.05 & log2FC > 1  ~ "Up in 6-mo",
        FDR < 0.05 & log2FC < -1 ~ "Down in 6-mo",
        TRUE ~ "Not significant"
      )
    )

  volcano_colors <- c(
    "Up in 6-mo"       = "#E64B35",
    "Down in 6-mo"     = "#4DBBD5",
    "Not significant"  = "grey70"
  )

  n_up <- sum(stats_df$Significance == "Up in 6-mo")
  n_down <- sum(stats_df$Significance == "Down in 6-mo")

  p <- ggplot(stats_df, aes(x = log2FC, y = -log10(FDR), color = Significance)) +
    geom_point(size = 1.5, alpha = 0.7) +
    geom_vline(xintercept = c(-1, 1), linetype = "dashed", color = "grey40") +
    geom_hline(yintercept = -log10(0.05), linetype = "dashed", color = "grey40") +
    scale_color_manual(values = volcano_colors) +
    labs(
      x = "log2(Fold Change) [6-mo / 0-mo]",
      y = "-log10(FDR)",
      title = paste("Volcano Plot —", title_suffix),
      subtitle = sprintf("Up: %d | Down: %d", n_up, n_down)
    ) +
    theme_bw() +
    theme(
      plot.title = element_text(hjust = 0.5, face = "bold"),
      plot.subtitle = element_text(hjust = 0.5)
    )

  # Label top significant features
  top_features <- stats_df |>
    filter(Significance != "Not significant") |>
    arrange(FDR) |>
    head(15)

  if (nrow(top_features) > 0) {
    p <- p + geom_text_repel(
      data = top_features,
      aes(label = Metabolite_name),
      size = 2.5, max.overlaps = 15, show.legend = FALSE
    )
  }

  pdf(output_path, width = 8, height = 6)
  print(p)
  dev.off()
  cat(sprintf("  Saved: %s\n", output_path))
  invisible(p)
}


# ----------------------------------------------------------
# Extract subject ID from sample name (number before underscore)
# ----------------------------------------------------------
extract_subject_id <- function(sample_name) {
  as.integer(str_extract(sample_name, "^\\d+"))
}
