#!/usr/bin/env Rscript
# ============================================================
# run_all_remaining.R — Run Phase 1/2/3 for Meta_Neg, Lip_Pos, Lip_Neg
# ============================================================

datasets <- c("Meta_Neg", "Lip_Pos", "Lip_Neg")

for (ds in datasets) {
  cat("\n\n########################################################\n")
  cat(sprintf("### DATASET: %s\n", ds))
  cat("########################################################\n\n")

  # --- Phase 1 ---
  env1 <- new.env(parent = globalenv())
  env1$dataset_name <- ds
  tryCatch({
    # Override dataset_name before sourcing
    tmp1 <- tempfile(fileext = ".R")
    code1 <- readLines("scripts/01_data_parsing_qc.R")
    code1 <- sub('^dataset_name\\s*<-.*',
                 sprintf('dataset_name <- "%s"', ds), code1)
    writeLines(code1, tmp1)
    source(tmp1, local = env1)
    unlink(tmp1)
  }, error = function(e) {
    cat(sprintf("ERROR in Phase 1 for %s: %s\n", ds, e$message))
  })

  # --- Phase 2 ---
  env2 <- new.env(parent = globalenv())
  env2$dataset_name <- ds
  tryCatch({
    tmp2 <- tempfile(fileext = ".R")
    code2 <- readLines("scripts/02_statistics_all.R")
    code2 <- sub('^dataset_name\\s*<-.*',
                 sprintf('dataset_name <- "%s"', ds), code2)
    writeLines(code2, tmp2)
    source(tmp2, local = env2)
    unlink(tmp2)
  }, error = function(e) {
    cat(sprintf("ERROR in Phase 2 for %s: %s\n", ds, e$message))
  })

  # --- Phase 3 ---
  env3 <- new.env(parent = globalenv())
  env3$dataset_name <- ds
  tryCatch({
    tmp3 <- tempfile(fileext = ".R")
    code3 <- readLines("scripts/03_statistics_paired.R")
    code3 <- sub('^dataset_name\\s*<-.*',
                 sprintf('dataset_name <- "%s"', ds), code3)
    writeLines(code3, tmp3)
    source(tmp3, local = env3)
    unlink(tmp3)
  }, error = function(e) {
    cat(sprintf("ERROR in Phase 3 for %s: %s\n", ds, e$message))
  })

  cat(sprintf("\n### %s — ALL PHASES COMPLETE ###\n", ds))
}

cat("\n\n========================================\n")
cat("ALL DATASETS PROCESSED!\n")
cat("========================================\n")
