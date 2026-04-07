#!/usr/bin/env Rscript
# ============================================================
# run_meta_noqc.R — Meta_Pos & Meta_Neg: skip QC filtering, run stats
# ============================================================

source("scripts/00_functions.R")
set.seed(42)

datasets <- c("Meta_Pos", "Meta_Neg")

input_files <- list(
  Meta_Pos = "Meta_Pos_Area_2_2026_03_17_15_52_18.xlsx",
  Meta_Neg = "Meta_Neg_Area_1_2026_03_17_14_48_29.xlsx"
)

grouping_file <- "Sample grouping.xlsx"

for (dataset_name in datasets) {
  cat("\n\n########################################################\n")
  cat(sprintf("### %s — No QC Filtering\n", dataset_name))
  cat("########################################################\n\n")

  output_dir <- file.path("output", paste0(dataset_name, "_noQC"))
  dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)

  # ============================================================
  # Step 1: Parse
  # ============================================================
  cat("Step 1: Parsing MS-DIAL export...\n")
  parsed <- parse_msdial_export(input_files[[dataset_name]])
  metadata  <- parsed$metadata
  intensity <- parsed$intensity
  class_info <- parsed$class_info

  grouping <- load_grouping(grouping_file)
  sample_info <- match_samples(colnames(intensity), class_info, grouping)

  valid_groups <- c("0-mo", "6-mo", "QC", "Blank")
  keep_cols <- sample_info$group %in% valid_groups
  intensity <- intensity[, keep_cols, drop = FALSE]
  sample_info <- sample_info[keep_cols, , drop = FALSE]

  cat(sprintf("  Features: %d, Samples: %d\n", nrow(intensity), ncol(intensity)))
  cat("  Group summary:\n")
  print(table(sample_info$group))

  # ============================================================
  # Skip QC filtering — only do imputation
  # ============================================================
  cat("\nSkipping blank/MV/RSD filtering. Imputing missing values...\n")
  intensity <- impute_half_min(intensity)

  filtering_log <- data.frame(
    Step = c("Raw (after parsing)", "After imputation (no QC filter)"),
    Features = c(nrow(intensity), nrow(intensity))
  )
  write.csv(filtering_log, file.path(output_dir, "01_filtering_summary.csv"),
            row.names = FALSE, fileEncoding = "UTF-8")

  # Save Phase 1 data
  saveRDS(list(metadata = metadata, intensity = intensity, sample_info = sample_info),
          file.path(output_dir, "phase1_data.rds"))

  # ============================================================
  # Phase 2: All-sample statistics
  # ============================================================
  cat("\n--- Phase 2: All-Sample Statistics ---\n")

  norm_data <- normalize_data(intensity, sample_info)
  cat(sprintf("  0-mo=%d, 6-mo=%d\n",
              sum(norm_data$groups == "0-mo"), sum(norm_data$groups == "6-mo")))

  # PLS-DA
  plsda_all <- run_plsda(norm_data$scaled, norm_data$groups, ncomp = 2)
  cat(sprintf("  R2=%.4f, Q2=%.4f, VIP>1: %d\n",
              plsda_all$R2[3], plsda_all$Q2[3], sum(plsda_all$vip > 1.0)))

  plot_plsda_scores(plsda_all, paste(dataset_name, "All (no QC filter)"),
                    file.path(output_dir, "05_plsda_scoreplot_all.pdf"))

  vip_df <- data.frame(
    Alignment_ID = metadata$`Alignment ID`,
    Average_Rt = metadata$`Average Rt(min)`,
    Average_Mz = metadata$`Average Mz`,
    Metabolite_name = metadata$`Metabolite name`,
    Ontology = metadata$Ontology,
    VIP = plsda_all$vip
  ) |> arrange(desc(VIP))
  write.csv(vip_df, file.path(output_dir, "07_vip_scores_all.csv"),
            row.names = FALSE, fileEncoding = "UTF-8")

  # Univariate
  log2_mat <- norm_data$log2_data
  groups <- norm_data$groups
  idx_0 <- which(groups == "0-mo")
  idx_6 <- which(groups == "6-mo")
  n_feat <- nrow(log2_mat)
  pvals <- log2fc <- numeric(n_feat)
  for (i in seq_len(n_feat)) {
    tt <- tryCatch(t.test(log2_mat[i, idx_6], log2_mat[i, idx_0], var.equal = FALSE),
                   error = function(e) list(p.value = NA))
    pvals[i] <- tt$p.value
    log2fc[i] <- mean(log2_mat[i, idx_6]) - mean(log2_mat[i, idx_0])
  }
  fdr <- p.adjust(pvals, method = "BH")

  stats_all <- data.frame(
    Alignment_ID = metadata$`Alignment ID`,
    Average_Rt = metadata$`Average Rt(min)`,
    Average_Mz = metadata$`Average Mz`,
    Metabolite_name = metadata$`Metabolite name`,
    Ontology = metadata$Ontology,
    FC = 2^log2fc, log2FC = log2fc,
    p_value = pvals, FDR = fdr, VIP = plsda_all$vip
  )
  n_sig <- sum(fdr < 0.05, na.rm = TRUE)
  n_up <- sum(fdr < 0.05 & log2fc > 1, na.rm = TRUE)
  n_down <- sum(fdr < 0.05 & log2fc < -1, na.rm = TRUE)
  cat(sprintf("  FDR<0.05: %d | |log2FC|>1: Up=%d Down=%d\n", n_sig, n_up, n_down))

  write.csv(stats_all, file.path(output_dir, "09_stats_results_all.csv"),
            row.names = FALSE, fileEncoding = "UTF-8")
  plot_volcano(stats_all, file.path(output_dir, "08_volcano_all.pdf"),
               title_suffix = paste(dataset_name, "All (no QC filter)"))

  saveRDS(norm_data, file.path(output_dir, "phase2_norm_data.rds"))

  # ============================================================
  # Phase 3: Paired analysis
  # ============================================================
  cat("\n--- Phase 3: Paired Analysis ---\n")

  paired_subjects <- c(63,64,65,66,67,71,76,77,81,88,91,92,94,95,96,97,99,
                       101,103,105,107,108,111,113,114,116,119,120,121,122,
                       123,124,125,126,127,128,131,133,135,138)

  sample_mapping <- sample_info |>
    filter(col_name %in% norm_data$col_names) |>
    mutate(subject_id = extract_subject_id(sample_name))

  paired_mapping <- sample_mapping |> filter(subject_id %in% paired_subjects)
  paired_check <- paired_mapping |>
    group_by(subject_id) |>
    summarize(has_0 = any(group == "0-mo"), has_6 = any(group == "6-mo"), .groups = "drop") |>
    filter(has_0 & has_6)
  paired_mapping <- paired_mapping |> filter(subject_id %in% paired_check$subject_id)

  paired_0mo <- paired_mapping |> filter(group == "0-mo") |> arrange(subject_id)
  paired_6mo <- paired_mapping |> filter(group == "6-mo") |> arrange(subject_id)
  cat(sprintf("  Paired subjects: %d\n", nrow(paired_0mo)))

  paired_cols <- c(paired_0mo$col_name, paired_6mo$col_name)
  paired_groups <- factor(c(rep("0-mo", nrow(paired_0mo)), rep("6-mo", nrow(paired_6mo))),
                          levels = c("0-mo", "6-mo"))
  paired_subject_ids <- c(paired_0mo$subject_id, paired_6mo$subject_id)
  paired_scaled <- norm_data$scaled[, paired_cols, drop = FALSE]
  paired_log2 <- norm_data$log2_data[, paired_cols, drop = FALSE]

  plsda_p <- run_plsda(paired_scaled, paired_groups, ncomp = 2)
  cat(sprintf("  R2=%.4f, Q2=%.4f, VIP>1: %d\n",
              plsda_p$R2[3], plsda_p$Q2[3], sum(plsda_p$vip > 1.0)))

  plot_plsda_scores(plsda_p, paste(dataset_name, "Paired (no QC filter)"),
                    file.path(output_dir, "10_plsda_scoreplot_paired.pdf"),
                    paired_lines = TRUE, subject_ids = paired_subject_ids)

  vip_p_df <- data.frame(
    Alignment_ID = metadata$`Alignment ID`,
    Average_Rt = metadata$`Average Rt(min)`,
    Average_Mz = metadata$`Average Mz`,
    Metabolite_name = metadata$`Metabolite name`,
    Ontology = metadata$Ontology,
    VIP = plsda_p$vip
  ) |> arrange(desc(VIP))
  write.csv(vip_p_df, file.path(output_dir, "12_vip_scores_paired.csv"),
            row.names = FALSE, fileEncoding = "UTF-8")

  # Paired t-test
  pvals_p <- log2fc_p <- numeric(n_feat)
  for (i in seq_len(n_feat)) {
    v0 <- as.numeric(paired_log2[i, paired_0mo$col_name])
    v6 <- as.numeric(paired_log2[i, paired_6mo$col_name])
    tt <- tryCatch(t.test(v6, v0, paired = TRUE), error = function(e) list(p.value = NA))
    pvals_p[i] <- tt$p.value
    log2fc_p[i] <- mean(v6) - mean(v0)
  }
  fdr_p <- p.adjust(pvals_p, method = "BH")

  stats_paired <- data.frame(
    Alignment_ID = metadata$`Alignment ID`,
    Average_Rt = metadata$`Average Rt(min)`,
    Average_Mz = metadata$`Average Mz`,
    Metabolite_name = metadata$`Metabolite name`,
    Ontology = metadata$Ontology,
    FC = 2^log2fc_p, log2FC = log2fc_p,
    p_value = pvals_p, FDR = fdr_p, VIP = plsda_p$vip
  )
  n_sig <- sum(fdr_p < 0.05, na.rm = TRUE)
  n_up <- sum(fdr_p < 0.05 & log2fc_p > 1, na.rm = TRUE)
  n_down <- sum(fdr_p < 0.05 & log2fc_p < -1, na.rm = TRUE)
  cat(sprintf("  FDR<0.05: %d | |log2FC|>1: Up=%d Down=%d\n", n_sig, n_up, n_down))

  write.csv(stats_paired, file.path(output_dir, "14_stats_results_paired.csv"),
            row.names = FALSE, fileEncoding = "UTF-8")
  plot_volcano(stats_paired, file.path(output_dir, "13_volcano_paired.pdf"),
               title_suffix = paste(dataset_name, "Paired (no QC filter)"))

  cat(sprintf("\n### %s (no QC) — ALL PHASES COMPLETE ###\n", dataset_name))
}

cat("\n========================================\n")
cat("Meta_Pos & Meta_Neg (no QC) DONE!\n")
cat("========================================\n")
