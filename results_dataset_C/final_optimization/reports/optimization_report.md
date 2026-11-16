# Final Dataset C optimization report

## Execution and scientific scope

Work was limited to Dataset C. Original CSVs and previous result files were preserved. New models were trained only on the original 564-row training partition. The original 141-row holdout was evaluated after model selection and was not used to tune models, features, or thresholds. This is exploratory classification of recorded vital status (0 alive, 1 deceased), not survival-time prediction. There are no patient IDs or follow-up durations, so patient-level independence cannot be verified.

## Methodology audit

`methodology_audit.md` records review of 23 saved pipelines and prediction files. All 23 prediction files use identical test rows; recomputed confusion matrices match saved metrics in class order [0,1] (TN,FP,FN,TP). Target leakage was not found. Preprocessing and previous feature selectors are pipeline-contained. Earlier post-hoc display of highest test scores is acknowledged as optimistic ranking; this new selection uses nested training-only CV.

## Feature quality

There are 705 samples, 1936 predictors, 611 alive and 94 deceased. Quality summaries and exhaustive |Pearson r|≥0.95 pairs are in `feature_quality_summary.csv` and `correlation_analysis.csv`; see `feature_quality_report.md`. Mutation sparsity is in `reports/mutation_sparsity.json`; CNV value categories and RNA/protein distributions/outlier screens are summarized separately. A total of 62 exact duplicate-vector groups with 99 redundant CNV aliases remain documented; alias equality does not mean biological identity.

## Model optimization and validation

Compared class-weighted Logistic Regression, Random Forest, Extra Trees, RBF SVM, and HistGradientBoosting across CNV, RNA, mutation, protein, and integrated inputs using nested stratified 3-fold outer / 3-fold inner CV on training rows. Within inner CV, candidates were selected by average rank across ROC-AUC, PR-AUC, balanced accuracy, macro-F1, and deceased recall. Each outer fold evaluated its inner-selected model. `nested_cv_summary.csv`, `model_selection_summary.csv`, `nested_outer_results.csv`, and `model_fit_errors.csv` contain exact results.

The best overall training-CV scope/algorithm selected by the composite priority ranking was **rna / random_forest**. Its nested outer metrics were ROC-AUC 0.632 ± 0.073, PR-AUC 0.405 ± 0.053, balanced accuracy 0.640 ± 0.008, macro-F1 0.667 ± 0.017, deceased recall 0.320 ± 0.000. These are nested CV estimates on this small training cohort, not independent validation.

Per-scope final models were then refit on all 564 training rows and evaluated once on the original holdout. Exact held-out values are in `selected_scope_models_test_metrics.csv`; these descriptive values were not used to select the global winner. The holdout has only 19 deceased observations. The phase-3 majority baseline remains the comparator.

## Class imbalance and feature selection

Class weighting was used for supported estimators. Balanced Random Forest and SMOTE were not run because imbalanced-learn is not installed; no dependencies were installed. No validation/test samples were resampled. Mutual information, ANOVA SelectKBest, L1, and Extra-Trees importance selectors were compared within five training folds on RNA and integrated matrices; results and fold support are in `feature_selection/`. Deceased recall, PR-AUC, balanced accuracy, macro-F1, and specificity accompany ROC-AUC. Accuracy was not the selector criterion.

## Explainability and modality contribution

Outer-fold training importances and selector frequencies are in `feature_selection/`. These are model-dependent associations, not causal findings or validated biomarkers. SHAP/XGBoost/LightGBM were unavailable and skipped without installation. Initial Dataset C repeated-CV comparisons found RNA RF ROC-AUC 0.684 ± 0.008 versus integrated Extra Trees 0.671 ± 0.016. Prior held-out permutation tests were not significant (RNA p=0.133; integrated p=0.305). No integration benefit is presumed; updated nested results are reported directly.

## Limitations and conclusion

Only 94 of 705 labels are deceased, including just 19 in the original test set. Patient IDs/barcodes and follow-up duration/censoring data are unavailable; patient-level independence and survival outcomes cannot be verified. There is no independent external validation. These experiments cannot diagnose breast cancer or predict an individual's survival. Results are exploratory recorded-vital-status classification only and are not clinically validated.

## Reproducibility

See `reproducibility/` for Python/package versions, seeds, split indices, fold configuration, exact model and preprocessing settings, and fitted model paths. New final models are under `models/`; per-scope test predictions are under `predictions/`. The full run log is `logs/optimization.log`.


## Explicit baseline and final execution status

The DummyClassifier prior baseline predicted all cases alive: test accuracy 0.865, balanced accuracy 0.500, ROC-AUC 0.500, PR-AUC 0.135, and deceased recall 0.000. Its confusion counts are TN=122, FP=0, FN=19, TP=0. This is the accuracy-only trap under class imbalance.

HistGradientBoosting completed nested CV in the approved execution after the sandbox prevented worker-pipe creation. All 25 scope/model candidates completed with 0 fitting errors. The candidate metrics are in `../tables/model_selection_summary.csv`; original holdout scores remain descriptive only.
