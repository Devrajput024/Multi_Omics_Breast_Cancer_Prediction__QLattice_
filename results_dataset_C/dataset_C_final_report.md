# Dataset C Phase 2 — exploratory modelling report

## Scope and scientific limits

These experiments classify the dataset’s recorded `vital.status` code (0 alive, 1 deceased). They are **not survival-time prediction**: no survival duration, follow-up, event time, or censoring variable is available. No patient IDs/barcodes are present; the split is row-based, and patient-level independence cannot be independently verified. Results are exploratory and do not establish clinical usefulness or generalization to new patients.

## Data and split

- Input: `C:\Users\asus\OneDrive\Desktop\Breast part 2\dataset_c_extracted\data.csv`; 705 rows, 1,936 predictors and target; target counts 0=611, 1=94.
- Stratified row split: train=564 (0=489, 1=75), test=141 (0=122, 1=19); seed=20261110; test fraction=0.2.
- Stratified 5-fold CV on training split; folds shuffled with seed 20261110.
- Prefix counts: cn_=860, rs_=604, mu_=249, pp_=223; no unexpected predictors.
- CNV values treated as categorical and one-hot encoded; RNA/protein median-imputed and standardized; mutation features median-imputed and training-fold constant columns removed. All transforms and selectors live inside scikit-learn pipelines and are fit separately within CV folds.

## Baseline models and results

Compared class-weighted Logistic Regression, Random Forest (100 trees), and Extra Trees (100 trees) for CNV, RNA, mutation, protein, and integrated data. Same held-out rows used throughout. Metrics include accuracy, deceased precision/recall, alive specificity, F1, macro-F1, ROC-AUC, average precision (PR-AUC), class-specific precision/recall, and confusion counts; full values are in `metrics_summary.csv` and fold statistics in `cross_validation_results.csv`.

| Scope | Model | CV ROC-AUC mean ± SD | CV PR-AUC mean ± SD | Test ROC-AUC | Test PR-AUC | Test deceased recall | Test specificity | Test macro-F1 | Confusion TN/FP/FN/TP |
|---|---|---:|---:|---:|---:|---:|---:|---:|---|
| cnv | logistic_regression | 0.510 ± 0.071 | 0.147 ± 0.027 | 0.535 | 0.150 | 0.158 | 0.811 | 0.484 | 99/23/16/3 |
| cnv | random_forest | 0.531 ± 0.034 | 0.158 ± 0.016 | 0.448 | 0.151 | 0.000 | 0.984 | 0.460 | 120/2/19/0 |
| cnv | extra_trees | 0.523 ± 0.061 | 0.164 ± 0.034 | 0.403 | 0.111 | 0.000 | 0.992 | 0.462 | 121/1/19/0 |
| rna | logistic_regression | 0.648 ± 0.055 | 0.343 ± 0.049 | 0.626 | 0.310 | 0.474 | 0.754 | 0.566 | 92/30/10/9 |
| rna | random_forest | 0.690 ± 0.035 | 0.446 ± 0.048 | 0.579 | 0.348 | 0.316 | 0.967 | 0.673 | 118/4/13/6 |
| rna | extra_trees | 0.683 ± 0.066 | 0.453 ± 0.083 | 0.546 | 0.358 | 0.263 | 0.984 | 0.661 | 120/2/14/5 |
| mutation | logistic_regression | 0.596 ± 0.095 | 0.237 ± 0.081 | 0.390 | 0.122 | 0.158 | 0.770 | 0.465 | 94/28/16/3 |
| mutation | random_forest | 0.591 ± 0.065 | 0.221 ± 0.055 | 0.466 | 0.150 | 0.211 | 0.795 | 0.498 | 97/25/15/4 |
| mutation | extra_trees | 0.588 ± 0.066 | 0.213 ± 0.068 | 0.459 | 0.125 | 0.105 | 0.836 | 0.472 | 102/20/17/2 |
| protein | logistic_regression | 0.566 ± 0.049 | 0.241 ± 0.050 | 0.494 | 0.196 | 0.263 | 0.713 | 0.475 | 87/35/14/5 |
| protein | random_forest | 0.605 ± 0.070 | 0.296 ± 0.100 | 0.471 | 0.168 | 0.053 | 0.984 | 0.507 | 120/2/18/1 |
| protein | extra_trees | 0.564 ± 0.068 | 0.281 ± 0.039 | 0.438 | 0.174 | 0.053 | 0.992 | 0.511 | 121/1/18/1 |
| integrated | logistic_regression | 0.651 ± 0.060 | 0.342 ± 0.069 | 0.604 | 0.338 | 0.263 | 0.902 | 0.586 | 110/12/14/5 |
| integrated | random_forest | 0.665 ± 0.082 | 0.447 ± 0.082 | 0.605 | 0.356 | 0.316 | 0.975 | 0.683 | 119/3/13/6 |
| integrated | extra_trees | 0.684 ± 0.088 | 0.378 ± 0.121 | 0.538 | 0.324 | 0.158 | 0.984 | 0.590 | 120/2/16/3 |

Best individual baseline selected by training CV ROC-AUC: **rna / random_forest**, CV ROC-AUC 0.690 ± 0.035.
Best integrated baseline selected by training CV ROC-AUC: **extra_trees**, CV ROC-AUC 0.684 ± 0.088.
Highest integrated held-out ROC-AUC among the three prespecified integrated baselines was random_forest (0.605); highest individual held-out ROC-AUC was rna/logistic_regression (0.626). These test comparisons are descriptive and were not used to tune features. Using the models selected by training CV, integration did not improve the point estimates: best individual RNA/Random Forest had CV ROC-AUC 0.690, test ROC-AUC 0.579, and test PR-AUC 0.348; best integrated Extra Trees had CV ROC-AUC 0.684, test ROC-AUC 0.538, and test PR-AUC 0.324.

## CNV duplicate-vector comparison

The existing alias map’s 99 aliases were verified equal to representatives on all 705 observed rows. Keeping the first representative removed 99 columns: 860 → 761.

- cnv_dedup_original: CV ROC-AUC 0.510 ± 0.071; held-out ROC-AUC 0.535; held-out PR-AUC 0.150.
- cnv_dedup_deduplicated: CV ROC-AUC 0.513 ± 0.068; held-out ROC-AUC 0.524; held-out PR-AUC 0.149.
The equality is only of observed vectors in this dataset; it does not mean the gene aliases are biologically identical. See `duplicate_alias_map.csv`.

## Feature selection

Compared no selection with fold-contained mutual information, L1 logistic selection, and Extra-Trees model-based selection on CV-selected best individual scope (rna) and integrated data. Selector comparison is in `feature_selection/feature_selection_results.csv`; selected names/modalities are in `feature_selection/selected_features.csv`. No test labels were used for feature choice.

- rna / mutual_information: CV ROC-AUC 0.612 ± 0.038; held-out ROC-AUC 0.661; PR-AUC 0.348.
- rna / l1_logistic: CV ROC-AUC 0.652 ± 0.060; held-out ROC-AUC 0.630; PR-AUC 0.242.
- rna / tree_importance: CV ROC-AUC 0.662 ± 0.068; held-out ROC-AUC 0.551; PR-AUC 0.230.
- integrated / mutual_information: CV ROC-AUC 0.594 ± 0.032; held-out ROC-AUC 0.621; PR-AUC 0.309.
- integrated / l1_logistic: CV ROC-AUC 0.658 ± 0.041; held-out ROC-AUC 0.560; PR-AUC 0.261.
- integrated / tree_importance: CV ROC-AUC 0.630 ± 0.072; held-out ROC-AUC 0.578; PR-AUC 0.249.

Selected features from the fitted training-only selector pipelines (each row includes modality; these are model-dependent exploratory selections):
- rs_ADIPOQ (rna), rs_ADH1B (rna), rs_ANKRD30A (rna), rs_TAT (rna), rs_SLC30A8 (rna), rs_TUSC5 (rna), rs_LEP (rna), rs_FABP4 (rna), rs_CIDEC (rna), rs_CYP4F8 (rna), rs_PGR (rna), rs_PLIN1 (rna), rs_PI16 (rna), rs_NLRP2 (rna), rs_COL17A1 (rna), rs_PROL1 (rna), rs_SCGB2A1 (rna), rs_HBA1 (rna), rs_FOLR1 (rna), rs_MIA (rna)
- rs_CLEC3A (rna), rs_CPB1 (rna), rs_SCGB2A2 (rna), rs_SCGB1D2 (rna), rs_TFF1 (rna), rs_MUCL1 (rna), rs_GSTM1 (rna), rs_PIP (rna), rs_ADIPOQ (rna), rs_ADH1B (rna), rs_S100A7 (rna), rs_HMGCS2 (rna), rs_CYP2B7P1 (rna), rs_ANKRD30A (rna), rs_PRAME (rna), rs_TAT (rna), rs_SERPINA6 (rna), rs_AGR3 (rna), rs_TFAP2B (rna), rs_CYP4Z1 (rna)
- rs_CLEC3A (rna), rs_TAT (rna), rs_AGR3 (rna), rs_KCNJ3 (rna), rs_TUSC5 (rna), rs_CST1 (rna), rs_CIDEC (rna), rs_BMPR1B (rna), rs_TNNT1 (rna), rs_GPD1 (rna), rs_PLIN1 (rna), rs_CIDEA (rna), rs_MSLN (rna), rs_MUC15 (rna), rs_CALML3 (rna), rs_PLIN4 (rna), rs_ORM1 (rna), rs_HBA1 (rna), rs_PCK1 (rna), rs_COL11A1 (rna)
- cn_ARTN_0 (cnv), cn_DNASE2B_-1 (cnv), cn_NBPF4_0 (cnv), cn_AREG_2 (cnv), cn_ART3_0 (cnv), cn_ESR1_-1 (cnv), cn_PIP_-2 (cnv), cn_C8orf4_2 (cnv), cn_PREX2_0 (cnv), cn_MATN2_-1 (cnv), cn_IL33_2 (cnv), cn_SYNPO2L_2 (cnv), cn_KNDC1_1 (cnv), cn_EHF_1 (cnv), cn_DRD2_1 (cnv), cn_ADAMTS8_-1 (cnv), cn_POU4F1_2 (cnv), cn_SLC15A1_0 (cnv), cn_KIAA0125_2 (cnv), cn_CRYM_2 (cnv)
- cn_ISG15_-2 (cnv), cn_ISG15_-1 (cnv), cn_ISG15_0 (cnv), cn_ISG15_1 (cnv), cn_ISG15_2 (cnv), cn_PLCH2_-2 (cnv), cn_PLCH2_-1 (cnv), cn_PLCH2_0 (cnv), cn_PLCH2_1 (cnv), cn_PLCH2_2 (cnv), cn_SAMD11_-2 (cnv), cn_SAMD11_-1 (cnv), cn_SAMD11_0 (cnv), cn_SAMD11_1 (cnv), cn_SAMD11_2 (cnv), cn_TNFRSF18_-2 (cnv), cn_TNFRSF18_0 (cnv), cn_BAI2_0 (cnv), cn_CYP4Z2P_0 (cnv), cn_IL12RB2_0 (cnv)
- cn_BAI2_0 (cnv), cn_CYP4Z2P_0 (cnv), cn_NBPF4_0 (cnv), cn_NBPF6_0 (cnv), cn_FAM19A3_0 (cnv), cn_ZNF695_1 (cnv), cn_ROBO2_0 (cnv), cn_AQPEP_1 (cnv), cn_PCDHA10_0 (cnv), cn_PCDHA5_0 (cnv), cn_WISP3_2 (cnv), cn_GATA4_1 (cnv), cn_PEBP4_0 (cnv), cn_CTSW_1 (cnv), cn_ALDH3A1_0 (cnv), cn_GREB1L_0 (cnv), cn_TDRD12_1 (cnv), cn_MAG_1 (cnv), cn_FFAR2_1 (cnv), cn_SBSN_1 (cnv)

## Explainability

**Best Individual top features** (model-based descriptive importance; not causal or proof of biomarkers): rs_SCUBE2 [rna]; rs_CHRDL1 [rna]; rs_ALOX15 [rna]; rs_PLIN1 [rna]; rs_RBP4 [rna]; rs_ACVR1C [rna]; rs_LOC283392 [rna]; rs_ADH1B [rna]; rs_APOB [rna]; rs_ANO3 [rna]; rs_ABCA8 [rna]; rs_COL11A1 [rna]; rs_FOSB [rna]; rs_PPAPDC1A [rna]; rs_GLYAT [rna]; rs_MMP10 [rna]; rs_KRT17 [rna]; rs_HEPN1 [rna]; rs_SLC7A2 [rna]; rs_PAX7 [rna].
**Best Integrated top features** (model-based descriptive importance; not causal or proof of biomarkers): rs_ACVR1C [rna]; rs_ADH1C [rna]; rs_AQP7P1 [rna]; rs_GPD1 [rna]; rs_MRAP [rna]; rs_APOB [rna]; rs_C6 [rna]; rs_FIGF [rna]; rs_SMYD1 [rna]; rs_GLYAT [rna]; rs_ALOX15 [rna]; rs_HEPACAM [rna]; rs_FXYD1 [rna]; rs_IL33 [rna]; rs_CIDEC [rna]; rs_C14orf180 [rna]; rs_LOC283392 [rna]; cn_MAG_1 [cnv]; pp_PI3K.p85 [protein]; rs_SLC19A3 [rna].
Permutation importance for a limited set of top raw features is in `explainability/top_feature_importance.csv`; plots are in `figures/feature_importance/`.

SHAP status: unavailable (ModuleNotFoundError: No module named 'shap').

## QLattice and subtype analysis

QLattice/Feyn: unavailable; exact error: `ModuleNotFoundError: No module named 'feyn'`. Historical serialized models were not treated as validated.
Subtype columns were inspected separately and excluded from the main model. Their prediction-time availability is undocumented, so leakage/confounding cannot be assessed; no subtype model was run.

## Execution issues and recovery

An initial parallel tree fit failed with `PermissionError: [WinError 5] Access is denied` when joblib attempted worker creation; the run was repeated with one in-process worker and all baselines completed. An initial CNV dedup indexing attempt failed with `KeyError: "None of [Index([...], dtype=int64, length=564)] are in the [columns]"` because row indices were applied as DataFrame column labels; this was corrected to `.iloc`, and the CNV comparison completed. Full current-run traceback is in `experiment_log.txt`; the first error is recorded here.

## Saved outputs

- `metrics_summary.csv`, `cross_validation_results.csv`, `experiment_log.txt`.
- Final training-split pipelines and JSON configurations under `models/`; test predictions under `predictions/`.
- Confusion matrices, ROC and precision-recall curves, comparison plot, feature importance plots under `figures/`.
- Feature-selection results/selected features under `feature_selection/`; explainability under `explainability/`; optional-tool and subtype audits under `reports/`.

## Limitations

Small cohort (705 rows), 6.5:1 imbalance, p≫n, unknown sample independence, no identifiers, no follow-up/time-to-event/censoring, undocumented measurement/prediction timing, and source-derived feature semantics constrain interpretation. A single random row split is not external validation. All reported estimates are exploratory recorded-vital-status classification only; no clinical survival prediction, clinical utility, or new-patient generalization is established.
