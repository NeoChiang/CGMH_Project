# ============================================================
# 03_statistics_paired.R — Phase 3: Paired Analysis
# PST Untargeted Metabolomics & Lipidomics Analysis
# ============================================================

# ============================================================
# CONFIGURATION — Change this block to switch datasets
# ============================================================
dataset_name   <- "Meta_Pos"   # Options: "Meta_Pos", "Meta_Neg", "Lip_Pos", "Lip_Neg"

output_dir     <- file.path("output", dataset_name)

# Paired subjects (n=40)
paired_subjects <- c(63, 64, 65, 66, 67, 71, 76, 77, 81, 88,
                     91, 92, 94, 95, 96, 97, 99, 101, 103, 105,
                     107, 108, 111, 113, 114, 116, 119, 120, 121, 122,
                     123, 124, 125, 126, 127, 128, 131, 133, 135, 138)

set.seed(42)

# ============================================================
# Load shared functions and Phase 2 data
# ============================================================
source("scripts/00_functions.R")

cat("============================================================\n")
cat(sprintf("Phase 3: Paired Analysis — %s\n", dataset_name))
cat("============================================================\n\n")

phase1 <- readRDS(file.path(output_dir, "phase1_data.rds"))
norm_data <- readRDS(file.path(output_dir, "phase2_norm_data.rds"))
metadata <- phase1$metadata
sample_info <- phase1$sample_info

# ============================================================
# Step 10: Extract Paired Subset
# ============================================================
cat("Step 10: Extracting paired subset...\n")

# Get sample info for 0-mo and 6-mo only
sample_cols <- norm_data$col_names
sample_groups <- as.character(norm_data$groups)

# Build mapping: col_name -> sample_name -> subject_id -> group
sample_mapping <- sample_info |>
  filter(col_name %in% sample_cols) |>
  mutate(subject_id = extract_subject_id(sample_name))

# Filter to paired subjects only
paired_mapping <- sample_mapping |>
  filter(subject_id %in% paired_subjects)

cat(sprintf("  Paired subjects found: %d\n", length(unique(paired_mapping$subject_id))))

# Verify each subject has both 0-mo and 6-mo
paired_check <- paired_mapping |>
  group_by(subject_id) |>
  summarize(
    has_0mo = any(group == "0-mo"),
    has_6mo = any(group == "6-mo"),
    .groups = "drop"
  ) |>
  filter(has_0mo & has_6mo)

cat(sprintf("  Subjects with both time points: %d\n", nrow(paired_check)))

valid_subjects <- paired_check$subject_id
paired_mapping <- paired_mapping |>
  filter(subject_id %in% valid_subjects)

# Get column indices for paired samples
paired_0mo <- paired_mapping |> filter(group == "0-mo") |> arrange(subject_id)
paired_6mo <- paired_mapping |> filter(group == "6-mo") |> arrange(subject_id)

cat(sprintf("  Paired 0-mo samples: %d\n", nrow(paired_0mo)))
cat(sprintf("  Paired 6-mo samples: %d\n", nrow(paired_6mo)))

# Extract paired data from normalized matrices
paired_cols <- c(paired_0mo$col_name, paired_6mo$col_name)
paired_groups <- factor(c(rep("0-mo", nrow(paired_0mo)),
                          rep("6-mo", nrow(paired_6mo))),
                        levels = c("0-mo", "6-mo"))
paired_subject_ids <- c(paired_0mo$subject_id, paired_6mo$subject_id)

# Pareto-scaled data for PLS-DA
paired_scaled <- norm_data$scaled[, paired_cols, drop = FALSE]
# Log2 data for univariate
paired_log2 <- norm_data$log2_data[, paired_cols, drop = FALSE]

# ============================================================
# Step 11: Paired PLS-DA
# ============================================================
cat("\nStep 11: Paired PLS-DA...\n")

plsda_paired <- run_plsda(paired_scaled, paired_groups, ncomp = 2)

cat(sprintf("  R2 (2 comp): %.4f\n", plsda_paired$R2[3]))
cat(sprintf("  Q2 (2 comp): %.4f\n", plsda_paired$Q2[3]))
cat(sprintf("  Features with VIP > 1.0: %d\n", sum(plsda_paired$vip > 1.0)))

# Score plot with paired lines
plot_plsda_scores(plsda_paired, "Paired Samples",
                  file.path(output_dir, "10_plsda_scoreplot_paired.pdf"),
                  paired_lines = TRUE,
                  subject_ids = paired_subject_ids)

# VIP table
vip_paired_df <- data.frame(
  Alignment_ID    = metadata$`Alignment ID`,
  Average_Rt      = metadata$`Average Rt(min)`,
  Average_Mz      = metadata$`Average Mz`,
  Metabolite_name = metadata$`Metabolite name`,
  Ontology        = metadata$Ontology,
  VIP             = plsda_paired$vip
) |>
  arrange(desc(VIP))

write.csv(vip_paired_df, file.path(output_dir, "12_vip_scores_paired.csv"),
          row.names = FALSE, fileEncoding = "UTF-8")
cat(sprintf("  Saved: %s\n", file.path(output_dir, "12_vip_scores_paired.csv")))

# Permutation test
cat("  Running permutation test (n=1000)...\n")
perm_paired <- run_permutation_test(paired_scaled, paired_groups,
                                    ncomp = 2, n_perm = 1000)
cat(sprintf("  Permutation p-value: %.4f\n", perm_paired$p_value))

plot_permutation(perm_paired,
                 file.path(output_dir, "11_plsda_permutation_paired.pdf"))

# ============================================================
# Step 12: Paired Univariate Analysis
# ============================================================
cat("\nStep 12: Paired univariate analysis...\n")

# Paired t-test: for each feature, compare within-subject differences
n_features <- nrow(paired_log2)
n_pairs <- nrow(paired_0mo)

pvals_paired  <- numeric(n_features)
log2fc_paired <- numeric(n_features)

for (i in seq_len(n_features)) {
  vals_0 <- as.numeric(paired_log2[i, paired_0mo$col_name])
  vals_6 <- as.numeric(paired_log2[i, paired_6mo$col_name])

  # Paired t-test (two-sided)
  tt <- tryCatch(
    t.test(vals_6, vals_0, paired = TRUE),
    error = function(e) list(p.value = NA)
  )
  pvals_paired[i] <- tt$p.value

  # Fold change for paired subset
  log2fc_paired[i] <- mean(vals_6) - mean(vals_0)
}

# FDR correction
fdr_paired <- p.adjust(pvals_paired, method = "BH")

# Combine results
stats_paired <- data.frame(
  Alignment_ID    = metadata$`Alignment ID`,
  Average_Rt      = metadata$`Average Rt(min)`,
  Average_Mz      = metadata$`Average Mz`,
  Metabolite_name = metadata$`Metabolite name`,
  Ontology        = metadata$Ontology,
  FC              = 2^log2fc_paired,
  log2FC          = log2fc_paired,
  p_value         = pvals_paired,
  FDR             = fdr_paired,
  VIP             = plsda_paired$vip
)

# Summary
n_sig <- sum(fdr_paired < 0.05, na.rm = TRUE)
n_sig_fc <- sum(fdr_paired < 0.05 & abs(log2fc_paired) > 1, na.rm = TRUE)
n_up <- sum(fdr_paired < 0.05 & log2fc_paired > 1, na.rm = TRUE)
n_down <- sum(fdr_paired < 0.05 & log2fc_paired < -1, na.rm = TRUE)

cat(sprintf("  Significant (FDR<0.05): %d\n", n_sig))
cat(sprintf("  Significant + |log2FC|>1: %d (Up=%d, Down=%d)\n",
            n_sig_fc, n_up, n_down))

# Save stats table
write.csv(stats_paired, file.path(output_dir, "14_stats_results_paired.csv"),
          row.names = FALSE, fileEncoding = "UTF-8")
cat(sprintf("  Saved: %s\n", file.path(output_dir, "14_stats_results_paired.csv")))

# Volcano plot
plot_volcano(stats_paired, file.path(output_dir, "13_volcano_paired.pdf"),
             title_suffix = "Paired Samples")

cat("\n============================================================\n")
cat(sprintf("Phase 3 complete for %s!\n", dataset_name))
cat("============================================================\n")
