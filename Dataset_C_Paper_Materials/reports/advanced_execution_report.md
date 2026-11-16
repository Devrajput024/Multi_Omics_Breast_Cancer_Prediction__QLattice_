# Final optimization addendum


## Final optimization and methodology verification (2026-11-10)


This addendum reports the later optimization run; earlier-phase results and tables were preserved. The new artifacts are in `../final_optimization/`. `model_comparison_with_dummy.csv` is the versioned comparison including the majority baseline.


### Methodology audit and new evaluation


The original modelling CSV has 705 rows and 1,936 predictors. The target is absent from all predictor matrices. All features belong to the four expected prefixes (CNV 860, RNA 604, mutation 249, protein 223). The previous 23 pipelines and 23 prediction files were audited: saved preprocessing and feature selection are pipeline-contained; all predictions use the same 141-row holdout; labels and confusion counts reconcile in class order `[0,1]` (TN, FP, FN, TP). The earlier report displayed top test scores descriptively; this can create optimistic post-hoc rankings, so the new model choice uses training-only nested CV.


The new evaluation used a stratified row split of 564 training rows and 141 holdout rows, seed 20261110. Nested stratified CV used 3 outer × 3 inner folds within training only. Imputation, one-hot encoding, scaling, mutation variance filtering and selectors were fitted within their respective training folds. The final holdout was evaluated after selection; it did not choose features, algorithms, or thresholds. Because no patient IDs/barcodes are supplied, independence between rows cannot be verified.


### Quality findings


The modelling matrix had 0 missing cells, 0 duplicate rows, 0 constant features, and 230 features with a single value in at least 95% of rows. The screen found 1,116 pairs with absolute Pearson correlation at least 0.95. It retained 62 groups of identical observed CNV feature vectors with 99 additional alias columns; this is not evidence of biological identity. Mutations were zero in 96.9% of cells; all observed mutation values were binary. IQR flags for RNA/protein are descriptive screens, not reasons to remove observations.


### Nested CV comparison and model choice


The CV-selected global candidate was **RNA Random Forest** by the prespecified average rank over ROC-AUC, PR-AUC, balanced accuracy, macro-F1, and deceased recall. The scope-level nested outer-fold estimates below evaluate the full inner-selection procedure; folds may select different algorithms. The CV table is therefore the safer estimate of the selection procedure, while candidate means in `model_selection_summary.csv` are inner-CV selection values and are not independent validation.


| Scope | Outer CV ROC-AUC | PR-AUC | Balanced accuracy | Macro-F1 | Deceased recall |
|---|---:|---:|---:|---:|---:|
| cnv | 0.543 ± 0.019 | 0.150 ± 0.006 | 0.480 ± 0.012 | 0.472 ± 0.018 | 0.093 ± 0.083 |
| rna | 0.632 ± 0.073 | 0.405 ± 0.053 | 0.640 ± 0.008 | 0.667 ± 0.017 | 0.320 ± 0.000 |
| mutation | 0.578 ± 0.075 | 0.200 ± 0.015 | 0.556 ± 0.047 | 0.535 ± 0.033 | 0.320 ± 0.120 |
| protein | 0.505 ± 0.138 | 0.208 ± 0.112 | 0.565 ± 0.023 | 0.576 ± 0.026 | 0.173 ± 0.061 |
| integrated | 0.671 ± 0.068 | 0.431 ± 0.076 | 0.630 ± 0.016 | 0.660 ± 0.015 | 0.293 ± 0.046 |


The SD values are across three outer folds; they are descriptive fold variability, not confidence intervals (outer-fold training sets overlap). Candidate means below average each algorithm’s inner CV across the three outer training partitions and are used for selection, so they are not unbiased generalization estimates.


| Scope | Candidate | Inner ROC-AUC | PR-AUC | Balanced accuracy | Macro-F1 | Deceased recall | Mean rank | Selected in scope |
|---|---|---:|---:|---:|---:|---:|---:|:---:|
| cnv | logistic_regression | 0.474 | 0.168 | 0.516 | 0.514 | 0.188 | 1.80 | yes |
| cnv | rbf_svm | 0.511 | 0.180 | 0.507 | 0.495 | 0.054 | 2.00 |  |
| cnv | hist_gradient_boosting | 0.470 | 0.151 | 0.504 | 0.503 | 0.106 | 3.00 |  |
| cnv | extra_trees | 0.511 | 0.151 | 0.495 | 0.462 | 0.000 | 4.00 |  |
| cnv | random_forest | 0.500 | 0.149 | 0.486 | 0.469 | 0.027 | 4.20 |  |
| integrated | random_forest | 0.659 | 0.405 | 0.627 | 0.659 | 0.280 | 1.40 | yes |
| integrated | hist_gradient_boosting | 0.668 | 0.377 | 0.601 | 0.625 | 0.234 | 2.40 |  |
| integrated | rbf_svm | 0.661 | 0.380 | 0.601 | 0.625 | 0.228 | 2.60 |  |
| integrated | logistic_regression | 0.641 | 0.320 | 0.596 | 0.607 | 0.261 | 3.80 |  |
| integrated | extra_trees | 0.636 | 0.343 | 0.542 | 0.539 | 0.094 | 4.80 |  |
| mutation | logistic_regression | 0.577 | 0.208 | 0.537 | 0.527 | 0.252 | 1.40 | yes |
| mutation | rbf_svm | 0.566 | 0.206 | 0.534 | 0.533 | 0.107 | 2.60 |  |
| mutation | random_forest | 0.560 | 0.187 | 0.537 | 0.526 | 0.252 | 3.00 |  |
| mutation | extra_trees | 0.560 | 0.191 | 0.514 | 0.513 | 0.166 | 3.80 |  |
| mutation | hist_gradient_boosting | 0.513 | 0.158 | 0.512 | 0.479 | 0.328 | 4.20 |  |
| protein | random_forest | 0.534 | 0.215 | 0.532 | 0.526 | 0.088 | 1.80 | yes |
| protein | rbf_svm | 0.532 | 0.197 | 0.510 | 0.504 | 0.095 | 2.60 |  |
| protein | logistic_regression | 0.509 | 0.191 | 0.513 | 0.496 | 0.281 | 3.00 |  |
| protein | extra_trees | 0.549 | 0.218 | 0.510 | 0.487 | 0.027 | 3.20 |  |
| protein | hist_gradient_boosting | 0.514 | 0.177 | 0.507 | 0.495 | 0.054 | 4.40 |  |
| rna | random_forest | 0.688 | 0.430 | 0.646 | 0.681 | 0.319 | 2.00 | yes |
| rna | rbf_svm | 0.671 | 0.424 | 0.646 | 0.669 | 0.339 | 2.40 |  |
| rna | extra_trees | 0.681 | 0.436 | 0.619 | 0.655 | 0.254 | 3.40 |  |
| rna | hist_gradient_boosting | 0.713 | 0.412 | 0.623 | 0.653 | 0.273 | 3.40 |  |
| rna | logistic_regression | 0.648 | 0.357 | 0.626 | 0.595 | 0.427 | 3.80 |  |


DummyClassifier predicted every sample alive. Its holdout accuracy was 0.865, balanced accuracy 0.500, ROC-AUC 0.500, PR-AUC 0.135, and deceased recall 0. This exposes why accuracy alone is misleading for the 611/94 imbalance. Exact baseline metrics and confusion counts are saved in `../final_optimization/tables/dummy_baseline_comparison.csv`.


### Descriptive holdout performance


The training-CV-selected model for each scope was fitted on all 564 training rows and evaluated on the fixed 141-row holdout (122 alive, 19 deceased). These test scores are descriptive, not model-selection criteria.


| Scope and selected model | ROC-AUC | PR-AUC | Accuracy | Balanced accuracy | Deceased precision | Deceased recall | Specificity | Macro-F1 | TN/FP/FN/TP |
|---|---:|---:|---:|---:|---:|---:|---:|---:|---:|
| cnv / logistic_regression | 0.535 | 0.150 | 0.723 | 0.485 | 0.115 | 0.158 | 0.811 | 0.484 | 99/23/16/3 |
| rna / random_forest | 0.579 | 0.348 | 0.879 | 0.642 | 0.600 | 0.316 | 0.967 | 0.673 | 118/4/13/6 |
| mutation / logistic_regression | 0.390 | 0.122 | 0.688 | 0.464 | 0.097 | 0.158 | 0.770 | 0.465 | 94/28/16/3 |
| protein / random_forest | 0.471 | 0.168 | 0.858 | 0.518 | 0.333 | 0.053 | 0.984 | 0.507 | 120/2/18/1 |
| integrated / random_forest | 0.605 | 0.356 | 0.887 | 0.646 | 0.667 | 0.316 | 0.975 | 0.683 | 119/3/13/6 |


Integrated Random Forest had higher descriptive holdout ROC-AUC and PR-AUC than RNA Random Forest, while deceased recall was equal (6/19). The nested CV integrated procedure also had higher ROC-AUC and PR-AUC point estimates than RNA, but slightly lower balanced accuracy, macro-F1 and deceased recall. The original phase comparison favored RNA over integrated Extra Trees, and prior holdout permutation tests were nonsignificant. Across these small samples and varying model selection, the evidence does not establish a reliable or meaningful integration gain.


### Feature selection and explainability


Mutual information, ANOVA SelectKBest, L1, and tree-importance selection were each fitted in five stratified training folds for RNA and integrated data. Fold-level selections and frequencies are saved under `../final_optimization/feature_selection/`. A feature selected often is only stable under that selector and this cohort; it is not a validated biomarker. Training-fold permutation importance was computed only for each fold’s native-importance top 30 variables on that fold’s validation data (3 repeats per variable). These are predictive associations, not causal evidence. SHAP was unavailable.


Top fold-validation permutation features (mean ROC-AUC decrease; fold support):


**Rna**

| Feature | Modality | Mean ROC-AUC decrease | Folds represented |
|---|---|---:|---:|

| rs_CYP4B1 | rna | 0.0103 | 1/3 |
| rs_TNMD | rna | 0.0072 | 1/3 |
| rs_ZNF385B | rna | 0.0072 | 1/3 |
| rs_GPM6A | rna | 0.0071 | 1/3 |
| rs_ALOX15 | rna | 0.0067 | 2/3 |
| rs_PCK1 | rna | 0.0058 | 1/3 |
| rs_GPR98 | rna | 0.0057 | 1/3 |
| rs_WFDC2 | rna | 0.0048 | 1/3 |
| rs_PSCA | rna | 0.0045 | 1/3 |
| rs_PI3 | rna | 0.0041 | 1/3 |

**Integrated**

| Feature | Modality | Mean ROC-AUC decrease | Folds represented |
|---|---|---:|---:|

| rs_SCUBE2 | rna | 0.0122 | 1/3 |
| rs_HBB | rna | 0.0111 | 1/3 |
| rs_MMP10 | rna | 0.0077 | 1/3 |
| rs_CIDEC | rna | 0.0070 | 1/3 |
| rs_LGALS12 | rna | 0.0066 | 1/3 |
| rs_ADH1C | rna | 0.0052 | 1/3 |
| rs_LOC283392 | rna | 0.0050 | 2/3 |
| rs_GPM6A | rna | 0.0047 | 1/3 |
| rs_TMEM132C | rna | 0.0042 | 2/3 |
| rs_GPD1 | rna | 0.0039 | 1/3 |


Examples selected in at least four of five folds by non-ANOVA selectors (feature stability is method-specific):


| Scope | Selector | Feature | Modality | Selected folds |
|---|---|---|---|---:|
| integrated | l1 | cn_ROBO2 | cnv | 5/5 |
| integrated | l1 | rs_ACTL8 | rna | 5/5 |
| integrated | l1 | rs_ANO3 | rna | 5/5 |
| rna | l1 | rs_ACTL8 | rna | 5/5 |
| rna | l1 | rs_ANO3 | rna | 5/5 |
| rna | l1 | rs_C6 | rna | 5/5 |
| integrated | mutual_information | rs_CHRDL1 | rna | 5/5 |
| integrated | mutual_information | rs_LEP | rna | 5/5 |
| integrated | mutual_information | rs_MMRN1 | rna | 5/5 |
| rna | mutual_information | rs_C14orf180 | rna | 5/5 |
| rna | mutual_information | rs_CHRDL1 | rna | 5/5 |
| rna | mutual_information | rs_FIGF | rna | 5/5 |
| integrated | tree_importance | rs_HEPN1 | rna | 5/5 |
| integrated | tree_importance | rs_ANGPTL7 | rna | 4/5 |
| integrated | tree_importance | rs_APOB | rna | 4/5 |
| rna | tree_importance | rs_ALOX15 | rna | 5/5 |
| rna | tree_importance | rs_ANO3 | rna | 5/5 |
| rna | tree_importance | rs_AQP7 | rna | 5/5 |


Selector comparison, five-fold training-only summary:


| Scope | Selector | ROC-AUC mean ± SD | PR-AUC | Balanced accuracy | Deceased recall |
|---|---|---:|---:|---:|---:|
| integrated | anova | 0.669 ± 0.063 | 0.343 | 0.618 | 0.493 |
| integrated | l1 | 0.663 ± 0.039 | 0.293 | 0.619 | 0.427 |
| integrated | mutual_information | 0.618 ± 0.043 | 0.313 | 0.597 | 0.467 |
| integrated | tree_importance | 0.640 ± 0.065 | 0.318 | 0.611 | 0.453 |
| rna | anova | 0.669 ± 0.064 | 0.343 | 0.615 | 0.493 |
| rna | l1 | 0.653 ± 0.056 | 0.313 | 0.635 | 0.493 |
| rna | mutual_information | 0.587 ± 0.022 | 0.290 | 0.582 | 0.480 |
| rna | tree_importance | 0.652 ± 0.047 | 0.292 | 0.636 | 0.547 |


Holdout Brier scores (no post-hoc calibrator was fitted):


| Model | Brier score |
|---|---:|
| cnv__logistic_regression | 0.2019 |
| integrated__random_forest | 0.1222 |
| mutation__logistic_regression | 0.2203 |
| protein__random_forest | 0.1428 |
| rna__random_forest | 0.1203 |
| dummy_prior | 0.1166 |

The Dummy prior Brier score is lower than the RNA and integrated model scores in this holdout; the model probability estimates therefore do not show better squared-probability accuracy here. Calibration curves are noisy with only 19 deceased cases.



### Calibration, uncertainty, and algorithms not run


No probability calibrator was fitted. Holdout calibration curves and Brier scores are descriptive only; the holdout has 19 deceased cases. Earlier stratified bootstrap intervals for RNA Random Forest were broad: test ROC-AUC 0.579, 95% percentile CI [0.413, 0.744], and deceased sensitivity 0.316, CI [0.105, 0.526]. Earlier one-sided permutation tests were not statistically significant for RNA (p=0.133) or integrated (p=0.305). These results do not support a stable performance claim.


Random Forest, Extra Trees, class-weighted Logistic Regression, RBF SVM, and HistGradientBoosting were compared. HistGradientBoosting required running outside the restricted sandbox because Windows worker-pipe creation was denied within it; all requested nested fits completed in the approved execution and the final run has 0 candidate errors. XGBoost, LightGBM, imbalanced-learn (Balanced Random Forest/SMOTE), and SHAP were unavailable and were not installed. Candidate-level inner-CV results for all estimators are in `../final_optimization/tables/model_selection_summary.csv`; they are model-selection estimates, not independent performance estimates.


### Model artifacts and limitations


CV-selected final scope pipelines, configs, predictions, fold indices, package versions, and preprocessing settings are under `../final_optimization/models/`, `predictions/`, and `reproducibility/`. Figures (PNG and PDF) include ROC, precision-recall, confusion matrices, calibration, nested-CV model comparison, feature-quality distributions/correlations, fold stability, and permutation importance.


The dataset contains only 94 deceased cases and the holdout only 19. There are no patient identifiers or follow-up durations, so row-level independence and survival endpoints cannot be assessed. There is no external validation. Recorded vital status is not survival-time prediction. Models are exploratory and not clinically validated; they should not be used to diagnose breast cancer or predict individual outcomes.
