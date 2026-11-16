# Dataset C audit and exploratory analysis

## Scope and method

This report analyzes only `archive (8).zip` as Dataset C. It was extracted into `dataset_c_extracted`; no other breast cancer project or dataset was included. No model was trained, no source data were modified, and no feature/sample was removed. The archive contains two tabular datasets, a README, a notebook, and two existing QLattice model files. The audit calculations are reproducible using the included Python script and Python standard library. This environment did not have pandas, NumPy, SciPy, scikit-learn, matplotlib, seaborn, or nbformat installed; plots were generated as SVG directly.

## Dimensions and files

| File | Rows | Columns | Numeric | Categorical | Missing cells | Duplicate rows | Identical-value column groups |
|---|---:|---:|---:|---:|---:|---:|---:|
| `data.csv` | 705 | 1937 | 1937 | 0 | 0 (0.000000%) | 0 | 62 groups, 99 extra columns |
| `brca_data_w_subtypes.csv` | 705 | 1941 | 1937 | 4 | 389 (0.028427%) | 0 | 62 groups, 99 extra columns |

All columns and descriptive statistics are in `dataset_structure.csv`. First/last five rows were inspected computationally; to avoid exposing sample-level records in the report, the structure CSV contains per-column SHA-256 fingerprints for each five-row edge instead of raw values.

## Target and task

The only target explicitly used by notebook code and existing QLattice artifacts is `vital.status`. It is binary-coded: 0 = 611/705 (86.67%); 1 = 94/705 (13.33%); majority/minority ratio = 6.50:1. **The actual clinical meaning of 0 and 1 is not verified.** The README lists subtype or survival prediction as general possibilities, but does not define these labels. There is no survival duration, censoring, or named molecular subtype column. A provisional binary classification task is the only task evidenced by the implementation; it is not ready until the code mapping and prediction time are documented. Full counts for `PR.Status`, `ER.Status`, `HER2.Final.Status`, and `histological.type` are in `target_analysis.txt` and `CLASS` summaries in the figures.

## Modalities and integration

The table has 604 `rs_`, 249 `mu_`, 860 `cn_`, and 223 `pp_` prefixed numeric columns. These are observed feature families, not confirmed modalities. Names and README context are consistent with some biological measurements, but no data dictionary establishes their assays or values. In particular, `rs_` must not be assumed to mean SNP or expression. Four categorical clinical/phenotype-like columns are available only in `brca_data_w_subtypes.csv`: PR, ER, HER2 status, and histological type. No demographic features or independent clinical metadata file are present.

Both tables have 705 rows, share 1,937 columns exactly in the same order, and match exactly row by row over those columns. The enriched file adds the four categorical fields. No explicit patient/sample identifiers exist, so patient-level matching across omics or with external cohorts cannot be measured. Feature families are technically co-stored in one numeric matrix, but shared-patient provenance is unverified.

## Missingness and data quality

`data.csv` has 0 missing cells. `brca_data_w_subtypes.csv` has 389 blanks out of 1,368,405 cells (0.0284%), all in `PR.Status` (122), `ER.Status` (122), and `HER2.Final.Status` (145). No exact duplicate rows were found. Each file has 62 groups of columns with identical value vectors, totaling 99 redundant extra columns. Constant and near-constant columns, all descriptive statistics, and IQR outlier screening counts are in `dataset_structure.csv` and `data_quality_report.txt`. The 40-column reproducible correlation screen is explicitly labeled non-exhaustive; target correlations cover every numeric candidate column and are descriptive only.

## Provenance, target leakage, limitations

The README links the archive to a Kaggle dataset titled “Multi-Omics Breast Cancer Prediction QLattice.” It says the source is Kaggle/Computational Oncology and says the data are “likely derived from TCGA or similar.” That wording is uncertain, and no accession, cohort citation, codebook, or IDs are included. TCGA provenance is therefore not established.

`vital.status` is included as a column in the candidate predictor table and must be excluded from predictors if it is the target. If receptor or histology labels become targets, related columns may leak the outcome; without temporal metadata, leakage risk cannot be resolved. Other risks include p > n, only 94 examples in the minority-coded class, 99 redundant columns, unknown sample unit, uncertain modality meanings, no external validation, and possible optimistic selection in the existing notebook. Its narrative metrics are not reproduced or independently validated here.

## Readiness and next experiments

**Not ready for model development.** First obtain authoritative provenance, a codebook defining `vital.status`, sample unit and prediction horizon, feature/assay dictionary, timing, and patient/sample IDs. Preserve both source tables; investigate duplicate-valued columns with the publisher. Once clarified: (1) define one clinical estimand; (2) specify patient-level nested/repeated stratified CV and an external holdout; (3) keep all imputation/scaling/selection inside training folds; (4) compare prevalence/regularized logistic and linear SVM baselines; (5) evaluate modality groups separately only after labels are verified; (6) test integration only with documented patient alignment; (7) consider tree boosting and QLattice against baselines with nested tuning; (8) report PR-AUC, ROC-AUC, sensitivity, specificity, calibration, confidence intervals, and permutation/statistical checks; (9) assess explanation stability; (10) validate externally. Deep learning is not justified by the current 705-row tabular sample. No experiments above were run.

## Figures

See `figures/` for target distribution, missingness, feature family counts, example numerical distributions, a feature correlation heatmap, and row-alignment diagram. The alignment visualization depicts positional equality, not patient-ID matching.

## Reproducibility

Run from the workspace root: `python dataset_c_audit/run_audit.py`. It reads only `dataset_c_extracted/data.csv` and `dataset_c_extracted/brca_data_w_subtypes.csv` and rewrites audit outputs in `dataset_c_audit/`. No row-level values or identifiers are written to reports.
