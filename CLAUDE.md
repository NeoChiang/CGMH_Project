# PST Untargeted Metabolomics & Lipidomics Analysis

## Project Overview

This project processes **MS-DIAL** exported Area data from an untargeted metabolomics and lipidomics experiment (LC-MS, Waters Xevo G2-XS QToF). The study is a longitudinal cohort comparing baseline (0-mo) vs follow-up (6-mo) time points.

There are **4 datasets** to analyze independently with the same workflow:

| Dataset | Mode | File |
|---------|------|------|
| Metabolomics Positive | ESI+ | `Meta_Pos_Area_2_2026_03_17_15_52_18.xlsx` |
| Metabolomics Negative | ESI- | `Meta_Neg_Area_1_2026_03_17_14_48_29.xlsx` |
| Lipidomics Positive | ESI+ | `Lip_Pos_Area_1_2026_03_17_14_42_20.xlsx` |
| Lipidomics Negative | ESI- | `Lip_Neg_Area_1_2025_10_08_21_17_38_raw.xlsx` |

All analysis scripts should be written in **R**. Statistical analysis and visualization follow standard metabolomics workflows compatible with **MetaboAnalyst**.

---

## Repository Structure

```
NeoChiang/PST/
├── CLAUDE.md                                              ← This file
├── Meta_Pos_Area_2_2026_03_17_15_52_18.xlsx               ← MS-DIAL export (Meta Pos)
├── Meta_Neg_Area_1_2026_03_17_14_48_29.xlsx               ← MS-DIAL export (Meta Neg)
├── Lip_Pos_Area_1_2026_03_17_14_42_20.xlsx                ← MS-DIAL export (Lip Pos)
├── Lip_Neg_Area_1_2025_10_08_21_17_38_raw.xlsx            ← MS-DIAL export (Lip Neg)
├── Sample grouping.xlsx                                   ← Sample grouping table
├── scripts/                                               ← R scripts (to be created)
│   ├── 00_functions.R                                     ← Shared utility functions
│   ├── 01_data_parsing_qc.R                               ← Phase 1: Parsing & QC
│   ├── 02_statistics_all.R                                 ← Phase 2: All-sample stats
│   └── 03_statistics_paired.R                              ← Phase 3: Paired analysis
└── output/                                                ← Results organized by dataset
    ├── Meta_Pos/
    ├── Meta_Neg/
    ├── Lip_Pos/
    └── Lip_Neg/
```

**GitHub repo:** `https://github.com/NeoChiang/PST`

**Claude Code usage:**
```bash
git clone https://github.com/NeoChiang/PST.git
cd PST
claude
# Then: "按照 CLAUDE.md 的 workflow，先處理 Meta_Pos 的 Phase 1"
```

---

## Input File Format (MS-DIAL Area Export)

All 4 xlsx files share the same MS-DIAL alignment result format:

- **Rows 1–4:** Header rows (Class, File type, Injection order, Batch ID)
- **Row 5:** Column names
- **Columns 1–34 (A–AH):** Feature metadata (Alignment ID, Average Rt, Average Mz, Metabolite name, Adduct type, Post curation result, Fill %, MS/MS assigned, Reference RT, Reference m/z, Formula, Ontology, INCHIKEY, SMILES, Annotation tag, RT matched, m/z matched, MS/MS matched, Comment, Manually modified for quantification, Manually modified for annotation, Isotope tracking parent ID, Isotope tracking weight number, RT similarity, m/z similarity, Simple dot product, Weighted dot product, Reverse dot product, Matched peaks count, Matched peaks percentage, Total score, S/N average, Spectrum reference file name, MS1 isotopic spectrum, MS/MS spectrum)
- **Columns 35 onward:** Sample intensity data
- **Last 12 columns:** Average and Stdev summary columns → **exclude from analysis**
- **Data type:** Normalized intensity (Area)

**Important parsing note:** The number of metadata columns (34) is consistent across all MS-DIAL exports, but the number of sample columns may vary between datasets. Always detect the boundary dynamically by reading the header row (row 5) and identifying where sample names begin.

---

## Sample Grouping

- **File:** `Sample grouping.xlsx` (note: filename has a space)
- **Columns:** Prefix, Sample number (B001–B168), Sample name, Group

### Matching MS-DIAL columns to grouping:
- MS-DIAL column names contain a prefix + Sample number, e.g. `20251203_PST_HS_Meta_Pos_B001`
- The prefix differs between datasets (Meta_Pos, Meta_Neg, Lip_Pos, Lip_Neg) and may also differ from the grouping file prefix
- **Match by Sample number only** (extract `B001`, `B002`, etc. from both sources)

### Groups:

| Group | N | Description |
|-------|---|-------------|
| 0-mo | 114 | Baseline samples |
| 6-mo | 39 | 6-month follow-up samples |
| 12-mo | 1 | Subject 121 — **reassign to 6-mo group** |
| QC | 10 | Pooled QC samples |
| Blank | 16 | Blank injections |
| Test | 2 | Test injections → **exclude** |

**Critical:** Sample `121_12-mo` (B168) must be reassigned to the **6-mo** group in all analyses.

---

## Paired Samples

40 subjects have both 0-mo and 6-mo/12-mo samples. Subject ID = number before underscore in sample name.

**Paired subjects (n=40):** 63, 64, 65, 66, 67, 71, 76, 77, 81, 88, 91, 92, 94, 95, 96, 97, 99, 101, 103, 105, 107, 108, 111, 113, 114, 116, 119, 120, 121, 122, 123, 124, 125, 126, 127, 128, 131, 133, 135, 138

Each paired subject has one 0-mo sample and one 6-mo sample (Subject 121: 0-mo ↔ 12-mo treated as 6-mo).

---

## Analysis Workflow

**Run the identical workflow for each of the 4 datasets independently.** Output to separate subdirectories under `output/`.

### Phase 1: Data Parsing & QC Assessment (on RAW intensity)

#### Step 1: Parse MS-DIAL Export
- Read xlsx, separate metadata (first 34 cols) from intensity matrix
- Remove last 12 columns (Average/Stdev summaries)
- Map column names to Sample number → Sample name → Group via grouping file
- Reassign 121_12-mo to 6-mo group
- Remove Test and NA columns; keep only: 0-mo, 6-mo, QC, Blank

#### Step 2: Blank Filtering
- For each feature, compare mean intensity in Blank vs mean in QC
- Remove features where `mean(Blank) / mean(QC) > 0.3`
- Also remove features where `mean(Blank) > mean(all Samples)`

#### Step 3: Missing Value Assessment
- Calculate detection rate (non-zero proportion) per feature in QC
- Remove features with QC detection rate < 80%
- Calculate detection rate per feature in each sample group (0-mo, 6-mo)
- Remove features with detection rate < 50% in **ALL** sample groups

#### Step 4: QC RSD Filtering (RAW intensity — NO normalization/scaling)
- Calculate RSD (CV%) = SD/mean × 100 for each feature across QC samples
- **Do NOT apply sum normalization, log transformation, or Pareto scaling**
- Remove features with QC RSD > 30%
- Report: total features before/after each filtering step

#### Step 5: Missing Value Imputation
- Replace remaining zeros/NAs with half-minimum value per feature

#### Step 6: QC PCA (raw intensity after filtering)
- PCA score plot: all samples (0-mo, 6-mo) + QC
- Color by group, different shapes for QC vs samples
- QC should cluster tightly near center
- Export as PDF

---

### Phase 2: Statistical Analysis — All Samples (114 × 0-mo vs 40 × 6-mo)

#### Step 7: Data Normalization
Apply to **sample data only** (exclude QC and Blank), in this order:
1. Sum normalization (normalize each sample to its total intensity)
2. Log2 transformation
3. Pareto scaling (mean-center, divide by √SD)

#### Step 8: PLS-DA
- PLS-DA: 0-mo vs 6-mo (2 groups)
- Score plot (PC1 vs PC2)
- VIP scores: extract features with VIP > 1.0
- Permutation test (n = 1000) for model validation
- Cross-validation (report Q2, R2)
- Export VIP table as CSV

#### Step 9: Univariate Analysis
- **Welch's t-test** (unpaired, two-sided): 0-mo vs 6-mo
- FDR correction (Benjamini-Hochberg)
- **Fold change:** FC = mean(6-mo) / mean(0-mo), calculated on **sum-normalized + log2-transformed** data (before Pareto scaling)
- Log2(FC)
- **Volcano plot:** x = log2(FC), y = −log10(FDR-adjusted p-value)
  - Thresholds: FDR < 0.05, |log2(FC)| > 1
  - Color: red = up in 6-mo, blue = down in 6-mo, grey = not significant
- Export full results table as CSV (Alignment ID, RT, m/z, Metabolite name, Ontology, FC, log2FC, p-value, FDR, VIP)

---

### Phase 3: Paired Analysis (40 paired subjects only)

#### Step 10: Extract Paired Subset
- Select only 40 paired subjects (80 samples)
- Use the same filtered & normalized data from Step 7

#### Step 11: Paired PLS-DA
- PLS-DA on paired subset: 0-mo vs 6-mo
- Score plot, VIP scores (VIP > 1.0)
- Permutation test (n = 1000), cross-validation
- Optional: connect paired samples with lines on score plot

#### Step 12: Paired Univariate Analysis
- **Paired t-test** (two-sided): within-subject 0-mo vs 6-mo
- FDR correction (Benjamini-Hochberg)
- Fold change: FC = mean(6-mo) / mean(0-mo) for paired subset
- Volcano plot with same thresholds
- Export paired results table as CSV

---

## Output Files (per dataset)

All outputs go to `./output/{Dataset}/` (e.g., `./output/Meta_Pos/`):

| File | Description |
|------|-------------|
| `01_filtering_summary.csv` | Feature count after each QC filtering step |
| `02_qc_rsd_distribution.pdf` | Histogram of QC RSD values |
| `03_pca_qc_scoreplot.pdf` | PCA score plot (QC assessment) |
| `04_filtered_data_matrix.csv` | Cleaned intensity matrix |
| `05_plsda_scoreplot_all.pdf` | PLS-DA score plot (all samples) |
| `06_plsda_permutation_all.pdf` | Permutation test result |
| `07_vip_scores_all.csv` | VIP scores for all features |
| `08_volcano_all.pdf` | Volcano plot (all samples) |
| `09_stats_results_all.csv` | Full statistical results table |
| `10_plsda_scoreplot_paired.pdf` | PLS-DA score plot (paired) |
| `11_plsda_permutation_paired.pdf` | Permutation test (paired) |
| `12_vip_scores_paired.csv` | VIP scores (paired) |
| `13_volcano_paired.pdf` | Volcano plot (paired) |
| `14_stats_results_paired.csv` | Full statistical results (paired) |

---

## Dataset Configuration

The R scripts should use a configuration block at the top to switch between datasets:

```r
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
```

---

## R Packages Required

```r
# Data handling
library(readxl)
library(dplyr)
library(tidyr)
library(tibble)
library(stringr)

# Statistics
library(ropls)          # PLS-DA (Bioconductor: BiocManager::install("ropls"))
# Alternative: library(mixOmics)

# Visualization
library(ggplot2)
library(ggrepel)
library(patchwork)
library(RColorBrewer)
```

---

## Coding Conventions

- Use tidyverse style (pipe `|>` or `%>%`)
- All plots with `ggplot2`, use `theme_bw()` or `theme_classic()`
- Color palette (colorblind-friendly):
  - 0-mo: `"#E64B35"` (red)
  - 6-mo: `"#4DBBD5"` (blue)
  - QC: `"#00A087"` (green)
  - Blank: `"#B09C85"` (grey-brown)
- PDF output: `width = 8, height = 6`
- CSV output: UTF-8 encoding
- `set.seed(42)` for all randomized procedures
- Comments in English
- Script should be modular: shared functions in `00_functions.R`, sourced by other scripts

---

## Important Notes

1. **QC RSD is calculated on raw intensity** — never after normalization/scaling
2. **Fold change is calculated on sum-normalized + log2 data** — before Pareto scaling
3. **Pareto scaling is only for multivariate analysis** (PCA, PLS-DA)
4. **The 4 datasets are analyzed independently** — do not merge them
5. **Subject 121 (12-mo) is always treated as 6-mo**
6. **Lip_Neg file has `_raw` suffix** — same MS-DIAL format, just different naming
