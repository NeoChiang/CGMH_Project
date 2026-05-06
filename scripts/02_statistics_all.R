# ============================================================
# 02_statistics_all.R — Phase 2: All-Sample Statistical Analysis
# PST Untargeted Metabolomics & Lipidomics Analysis
# ============================================================

# ============================================================
# CONFIGURATION — Change this block to switch datasets
# ============================================================
dataset_name   <- "Meta_Pos"   # Options: "Meta_Pos", "Meta_Neg", "Lip_Pos", "Lip_Neg"

output_dir     <- file.path("output", dataset_name)

set.seed(42)

# ============================================================
# Load shared functions and Phase 1 data
# ============================================================
source("scripts/00_functions.R")

cat("============================================================\n")
cat(sprintf("Phase 2: All-Sample Statistics — %s\n", dataset_name))
cat("============================================================\n\n")

phase1 <- readRDS(file.path(output_dir, "phase1_data.rds"))
metadata    <- phase1$metadata
intensity   <- phase1$intensity
sample_info <- phase1$sample_info

cat(sprintf("  Loaded %d features x %d samples\n", nrow(intensity), ncol(intensity)))

# ============================================================
# Step 7: Data Normalization
# ============================================================
cat("\nStep 7: Data normalization (sum-norm -> log2 -> Pareto)...\n")

norm_data <- normalize_data(intensity, sample_info)

cat(sprintf("  Sample groups: 0-mo=%d, 6-mo=%d\n",
            sum(norm_data$groups == "0-mo"),
            sum(norm_data$groups == "6-mo")))

# ============================================================
# Step 8: PLS-DA (All Samples)
# ============================================================
cat("\nStep 8: PLS-DA (all samples)...\n")

# Run PLS-DA on Pareto-scaled data
plsda_all <- run_plsda(norm_data$scaled, norm_data$groups, ncomp = 2)

cat(sprintf("  R2 (2 comp): %.4f\n", plsda_all$R2[3]))
cat(sprintf("  Q2 (2 comp): %.4f\n", plsda_all$Q2[3]))
cat(sprintf("  Features with VIP > 1.0: %d\n", sum(plsda_all$vip > 1.0)))

# Score plot
plot_plsda_scores(plsda_all, "All Samples",
                  file.path(output_dir, "05_plsda_scoreplot_all.pdf"))

# VIP table
vip_df <- data.frame(
  Alignment_ID    = metadata$`Alignment ID`,
  Average_Rt      = metadata$`Average Rt(min)`,
  Average_Mz      = metadata$`Average Mz`,
  Metabolite_name = metadata$`Metabolite name`,
  Ontology        = metadata$Ontology,
  VIP             = plsda_all$vip
) |>
  arrange(desc(VIP))

write.csv(vip_df, file.path(output_dir, "07_vip_scores_all.csv"),
          row.names = FALSE, fileEncoding = "UTF-8")
cat(sprintf("  Saved: %s\n", file.path(output_dir, "07_vip_scores_all.csv")))

# Permutation test
cat("  Running permutation test (n=1000)...\n")
perm_all <- run_permutation_test(norm_data$scaled, norm_data$groups,
                                 ncomp = 2, n_perm = 1000)
cat(sprintf("  Permutation p-value: %.4f\n", perm_all$p_value))

plot_permutation(perm_all,
                 file.path(output_dir, "06_plsda_permutation_all.pdf"))

# ============================================================
# Step 9: Univariate Analysis (All Samples)
# ============================================================
cat("\nStep 9: Univariate analysis (all samples)...\n")

# Use sum-normalized + log2 data (before Pareto) for t-test and FC
log2_mat <- norm_data$log2_data
groups   <- norm_data$groups

idx_0mo <- which(groups == "0-mo")
idx_6mo <- which(groups == "6-mo")

n_features <- nrow(log2_mat)
pvals  <- numeric(n_features)
log2fc <- numeric(n_features)

for (i in seq_len(n_features)) {
  vals_0 <- log2_mat[i, idx_0mo]
  vals_6 <- log2_mat[i, idx_6mo]

  # Welch's t-test (unpaired, two-sided)
  tt <- tryCatch(
    t.test(vals_6, vals_0, var.equal = FALSE),
    error = function(e) list(p.value = NA)
  )
  pvals[i] <- tt$p.value

  # Fold change: mean(6-mo) / mean(0-mo) on log2 scale = difference
  log2fc[i] <- mean(vals_6) - mean(vals_0)
}

# FDR correction
fdr <- p.adjust(pvals, method = "BH")

# Combine results
stats_all <- data.frame(
  Alignment_ID    = metadata$`Alignment ID`,
  Average_Rt      = metadata$`Average Rt(min)`,
  Average_Mz      = metadata$`Average Mz`,
  Metabolite_name = metadata$`Metabolite name`,
  Ontology        = metadata$Ontology,
  FC              = 2^log2fc,
  log2FC          = log2fc,
  p_value         = pvals,
  FDR             = fdr,
  VIP             = plsda_all$vip
)

# Summary
n_sig <- sum(fdr < 0.05, na.rm = TRUE)
n_sig_fc <- sum(fdr < 0.05 & abs(log2fc) > 1, na.rm = TRUE)
n_up <- sum(fdr < 0.05 & log2fc > 1, na.rm = TRUE)
n_down <- sum(fdr < 0.05 & log2fc < -1, na.rm = TRUE)

cat(sprintf("  Significant (FDR<0.05): %d\n", n_sig))
cat(sprintf("  Significant + |log2FC|>1: %d (Up=%d, Down=%d)\n",
            n_sig_fc, n_up, n_down))

# Save stats table
write.csv(stats_all, file.path(output_dir, "09_stats_results_all.csv"),
          row.names = FALSE, fileEncoding = "UTF-8")
cat(sprintf("  Saved: %s\n", file.path(output_dir, "09_stats_results_all.csv")))

# Volcano plot
plot_volcano(stats_all, file.path(output_dir, "08_volcano_all.pdf"),
             title_suffix = "All Samples")

# Save normalized data for Phase 3
saveRDS(norm_data, file.path(output_dir, "phase2_norm_data.rds"))

cat("\n============================================================\n")
cat(sprintf("Phase 2 complete for %s!\n", dataset_name))
cat("============================================================\n")
