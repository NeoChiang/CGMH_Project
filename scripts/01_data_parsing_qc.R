# ============================================================
# 01_data_parsing_qc.R — Phase 1: Data Parsing & QC Assessment
# PST Untargeted Metabolomics & Lipidomics Analysis
# ============================================================

# ============================================================
# CONFIGURATION — Change this block to switch datasets
# ============================================================
dataset_name   <- "Meta_Pos"   # Options: "Meta_Pos", "Meta_Neg", "Lip_Pos", "Lip_Neg"

input_files <- list(
  Meta_Pos = "Meta_Pos_Area_2_2026_03_17_15_52_18.xlsx",
  Meta_Neg = "Meta_Neg_Area_1_2026_03_17_14_48_29.xlsx",
  Lip_Pos  = "Lip_Pos_Area_1_2026_03_17_14_42_20.xlsx",
  Lip_Neg  = "Lip_Neg_Area_1_2025_10_08_21_17_38_raw.xlsx"
)

input_file     <- input_files[[dataset_name]]
grouping_file  <- "Sample grouping.xlsx"
output_dir     <- file.path("output", dataset_name)
dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)

set.seed(42)

# ============================================================
# Load shared functions
# ============================================================
source("scripts/00_functions.R")

cat("============================================================\n")
cat(sprintf("Processing dataset: %s\n", dataset_name))
cat(sprintf("Input file: %s\n", input_file))
cat("============================================================\n\n")

# ============================================================
# Step 1: Parse MS-DIAL Export
# ============================================================
cat("Step 1: Parsing MS-DIAL export...\n")

parsed <- parse_msdial_export(input_file)
metadata  <- parsed$metadata
intensity <- parsed$intensity
class_info <- parsed$class_info

cat(sprintf("  Total features: %d\n", nrow(intensity)))
cat(sprintf("  Total sample columns: %d\n", ncol(intensity)))

# Load grouping
grouping <- load_grouping(grouping_file)

# Match samples
sample_info <- match_samples(
  col_names = colnames(intensity),
  class_info = class_info,
  grouping_df = grouping
)

# Keep only valid groups: 0-mo, 6-mo, QC, Blank
valid_groups <- c("0-mo", "6-mo", "QC", "Blank")
keep_cols <- sample_info$group %in% valid_groups
excluded_count <- sum(!keep_cols)

cat(sprintf("  Excluding %d columns (Test/NA/unmatched)\n", excluded_count))

# Filter columns
intensity <- intensity[, keep_cols, drop = FALSE]
sample_info <- sample_info[keep_cols, , drop = FALSE]

# Group summary
cat("\n  Group summary:\n")
grp_summary <- table(sample_info$group)
for (g in names(grp_summary)) {
  cat(sprintf("    %s: %d samples\n", g, grp_summary[g]))
}

n_total <- nrow(intensity)
filtering_log <- data.frame(
  Step = "Raw (after parsing)",
  Features = n_total,
  stringsAsFactors = FALSE
)

# ============================================================
# Step 2: Blank Filtering
# ============================================================
cat("\nStep 2: Blank filtering...\n")

keep_blank <- blank_filter(intensity, sample_info)
intensity <- intensity[keep_blank, , drop = FALSE]
metadata  <- metadata[keep_blank, , drop = FALSE]

filtering_log <- rbind(filtering_log, data.frame(
  Step = "After blank filtering",
  Features = nrow(intensity)
))

# ============================================================
# Step 3: Missing Value Assessment
# ============================================================
cat("\nStep 3: Missing value assessment...\n")

keep_mv <- missing_value_filter(intensity, sample_info)
intensity <- intensity[keep_mv, , drop = FALSE]
metadata  <- metadata[keep_mv, , drop = FALSE]

filtering_log <- rbind(filtering_log, data.frame(
  Step = "After missing value filtering",
  Features = nrow(intensity)
))

# ============================================================
# Step 4: QC RSD Filtering (RAW intensity)
# ============================================================
cat("\nStep 4: QC RSD filtering (raw intensity)...\n")

rsd_result <- qc_rsd_filter(intensity, sample_info, rsd_threshold = 60)
rsd_values <- rsd_result$rsd
keep_rsd <- rsd_result$keep

# Plot RSD histogram BEFORE filtering (to show distribution)
plot_rsd_histogram(
  rsd_values = rsd_values,
  rsd_threshold = 60,
  output_path = file.path(output_dir, "02_qc_rsd_distribution.pdf")
)

intensity <- intensity[keep_rsd, , drop = FALSE]
metadata  <- metadata[keep_rsd, , drop = FALSE]

filtering_log <- rbind(filtering_log, data.frame(
  Step = "After QC RSD filtering",
  Features = nrow(intensity)
))

cat(sprintf("\n  Final feature count: %d\n", nrow(intensity)))

# ============================================================
# Step 5: Missing Value Imputation
# ============================================================
cat("\nStep 5: Missing value imputation (half-minimum)...\n")

intensity <- impute_half_min(intensity)

n_zeros <- sum(intensity == 0, na.rm = TRUE)
n_na <- sum(is.na(intensity))
cat(sprintf("  Remaining zeros: %d, NAs: %d\n", n_zeros, n_na))

# ============================================================
# Save filtering summary
# ============================================================
cat("\nFiltering summary:\n")
print(filtering_log, row.names = FALSE)

write.csv(filtering_log,
          file.path(output_dir, "01_filtering_summary.csv"),
          row.names = FALSE, fileEncoding = "UTF-8")
cat(sprintf("  Saved: %s\n", file.path(output_dir, "01_filtering_summary.csv")))

# ============================================================
# Step 6: QC PCA (raw intensity after filtering)
# ============================================================
cat("\nStep 6: QC PCA assessment...\n")

plot_pca_qc(
  intensity = intensity,
  sample_info = sample_info,
  output_path = file.path(output_dir, "03_pca_qc_scoreplot.pdf")
)

# ============================================================
# Save filtered data matrix
# ============================================================
cat("\nSaving filtered data matrix...\n")

# Combine metadata + intensity for export
output_matrix <- cbind(metadata, intensity)
write.csv(output_matrix,
          file.path(output_dir, "04_filtered_data_matrix.csv"),
          row.names = FALSE, fileEncoding = "UTF-8")
cat(sprintf("  Saved: %s\n", file.path(output_dir, "04_filtered_data_matrix.csv")))

# Also save sample_info for use by Phase 2/3 scripts
saveRDS(list(
  metadata = metadata,
  intensity = intensity,
  sample_info = sample_info
), file = file.path(output_dir, "phase1_data.rds"))
cat(sprintf("  Saved: %s\n", file.path(output_dir, "phase1_data.rds")))

cat("\n============================================================\n")
cat(sprintf("Phase 1 complete for %s!\n", dataset_name))
cat("============================================================\n")
