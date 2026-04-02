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

  # Read header rows (rows 1-4) with no col_names
  header_raw <- read_excel(filepath, sheet = 1, col_names = FALSE,
                           n_max = 4, .name_repair = "minimal")
  header_raw <- as.data.frame(header_raw)

  # Row 1 = Class info
  class_row <- as.character(header_raw[1, ])
  # Row 4 = Batch ID (used to detect summary columns: "Average"/"Stdev")
  batch_row <- as.character(header_raw[4, ])

  # Read column names from row 5
  col_names_raw <- read_excel(filepath, sheet = 1, col_names = FALSE,
                              skip = 4, n_max = 1, .name_repair = "minimal")
  col_names <- as.character(col_names_raw[1, ])

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
    data_end <- length(col_names)
  }

  # Sample column indices
  sample_idx <- sample_start:data_end
  n_samples <- length(sample_idx)

  cat(sprintf("  Metadata columns: 1-%d (%d cols)\n", meta_end, meta_end))
  cat(sprintf("  Sample columns: %d-%d (%d cols)\n", sample_start, data_end, n_samples))
  cat(sprintf("  Summary columns excluded: %d cols\n", length(summary_idx)))

  # Read the actual data (skip 5 header rows)
  # Specify column types: metadata as text, samples as numeric
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

  # For unmatched columns, use class_info if available
  unmatched <- is.na(result$group)
  if (any(unmatched)) {
    # Test and NA columns stay unmatched
    cat(sprintf("  Unmatched columns: %d (will be excluded if Test/NA)\n",
                sum(unmatched)))
    cat(sprintf("  Unmatched classes: %s\n",
                paste(unique(result$class[unmatched]), collapse = ", ")))
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
