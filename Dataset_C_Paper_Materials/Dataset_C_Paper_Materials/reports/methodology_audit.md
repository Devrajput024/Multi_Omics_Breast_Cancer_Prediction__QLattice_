# Methodology audit

Audited the previous Phase 2/3 scripts, all 23 saved central pipelines, 23 prediction CSVs, metric tables, logs, and saved split metadata. No source CSV or previous output was changed.

- `vital.status` is binary 0 alive / 1 deceased and is dropped from X before fitting. None of the 23 saved preprocessors includes the target.
- All 23 prediction files use the same exact 141 test row indices and true labels (122 alive, 19 deceased). Recomputed confusion matrices use class order `[0,1]`, so flattened cells are TN, FP, FN, TP; all counts match saved metrics.
- Imputation, scaling, one-hot encoding, mutation variance filtering, and selectors are pipeline steps. Phase 2 CV ran on the 564-row training partition; selectors are inside CV pipelines. Phase 3 audit independently agrees.
- Test scores were used descriptively to list highest held-out results in the previous report, creating post-hoc winner ranking. No selector, model hyperparameter, or threshold was chosen from test results. Treat all maximum test comparisons as optimistic/descriptive.
- Original CNV dedup indexing error and initial worker error were fixed and documented in earlier logs. Saved final results reconcile.
- No IDs or follow-up time exist; row-based splits cannot verify patient-level independence.

Audit summary: `{"saved_prediction_files": 23, "saved_model_files": 23, "reconciled_predictions": 23, "target_in_saved_pipeline_preprocessor": false, "phase2_target_dropped_in_source": true, "transforms_and_selection_inside_pipeline": true, "phase3_audit_confirmed": true, "phase2_cv_used_training_only": true, "all_saved_models_use_same_holdout": true, "all_saved_confusion_counts_match": true, "holdout_winner_selection_caveat": "Previous report ranked some models by held-out ROC-AUC descriptively. This is post-hoc and optimistic if interpreted as model selection. No Phase 2 model parameters, selectors, or thresholds were chosen from test scores; this optimization uses training-only inner CV and outer CV, with holdout metrics reported only after selection."}`

## Final optimization run checks

The 705-row source matrix has 1,936 predictors with zero columns outside the four declared modality prefixes; the only non-prefixed field is the target. No patient/sample/barcode or follow-up/survival-duration column is present. Target counts remain 611 alive and 94 deceased. The new run reconstructs and verifies the exact prior split, then fits all transforms and feature selectors only within nested training folds. All 25 scope/algorithm candidates (five model families × five scopes) completed with zero fit errors after HistGradientBoosting was run through the approved non-sandbox execution path. The final report uses nested outer results and the original fixed holdout; test metrics do not participate in candidate selection.

The new feature-selection comparisons use five stratified folds within the 564 training rows. Fold-validation permutation importance for RNA/integrated models uses only each outer-training model's native-importance top 30 features; no holdout labels were used for this analysis. SHAP and resampling packages were unavailable. The complete new CV, test, and fit-error tables are under `tables/`.
