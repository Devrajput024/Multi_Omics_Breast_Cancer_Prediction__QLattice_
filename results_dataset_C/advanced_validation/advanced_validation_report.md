# Dataset C — Phase 3 advanced validation report

## Purpose and limits

This phase audits and extends the existing exploratory classifiers for recorded `vital.status` (0 alive, 1 deceased). It is not survival-time prediction: there are no durations, follow-up, event times, censoring fields, or patient IDs/TCGA barcodes. Splits remain row-based; patient-level independence cannot be independently verified. Nothing here validates clinical utility or generalization to new patients. Original phase-2 results remain unchanged.

## Audit of Phase 2 implementation

- Found 23 saved pipelines, 23 saved prediction files, and 23 predictions reconciled against `metrics_summary.csv`. Every prediction file has 141 rows and exactly the same saved test row indices and true labels. Recomputed metrics match the summary within maximum absolute error 0.
- Target excluded: all saved preprocessor input columns exclude `vital.status`; the phase-2 runner built X by dropping the target before fitting.
- Leakage control: imputers, scaling, CNV one-hot encoding, mutation variance filtering, and selectors are Pipeline steps. `cross_validate` received training-partition rows only. Feature-selection selectors are nested inside the pipeline and therefore refit per training fold.
- The same stratified 564/141 row split (seed 20261110) is used by every saved prediction; class order is [0,1], so confusion counts are TN, FP, FN, TP.
- Test-ranking caveat: the previous report also listed the highest held-out ROC-AUC models after comparing many candidates. This is posthoc winner ranking and can be optimistic. No test-derived threshold or feature selection was found; Phase 3 freezes model families by Phase-2 training-CV results and treats all test metrics as descriptive.
- Phase-2 source has no patient identifiers and a single random split; no grouped or external validation is possible. See `reports/implementation_audit.json` and `reports/prediction_reconciliation.json`.

## Majority-class baseline and model comparison

The majority DummyClassifier predicts alive for all held-out rows (141 rows; 122 alive, 19 deceased). Full comparison against all 23 existing prediction files is in `baseline_comparison.csv`; it includes accuracy, balanced accuracy, macro-F1, ROC-AUC, PR-AUC, class-specific metrics, and confusion counts. The dummy holdout baseline is accuracy 0.865, balanced accuracy 0.500, macro-F1 0.464, ROC-AUC 0.500, PR-AUC 0.135, with zero deceased recall (TN/FP/FN/TP 122/0/19/0).

## Repeated stratified CV

Used 5 repeats of stratified 5-fold CV on the 564-row training partition. For modality comparisons, architectures were fixed in advance from the previous training-CV ranking: CNV/RF, RNA/RF, mutation/Logistic Regression, protein/RF, integrated/Extra Trees. Dummy majority was scored on the exact same folds. `repeated_cv_results.csv` contains all 150 fold-level rows (including dummy) and repeat-summary rows. Mean/SD and empirical 2.5–97.5% quantiles are computed over five repeat means; these are descriptive resampling ranges, not formal independent-sample confidence intervals because folds/repeats overlap.

| Scope | Fixed model | ROC-AUC mean ± SD of repeat means | PR-AUC mean ± SD of repeat means | empirical repeat-mean 95% range |
|---|---|---:|---:|---:|
| cnv | random_forest | 0.515 ± 0.016 | 0.156 ± 0.011 | 0.494–0.531 |
| dummy_majority | dummy | 0.500 ± 0.000 | 0.133 ± 0.000 | 0.500–0.500 |
| integrated | extra_trees | 0.671 ± 0.016 | 0.388 ± 0.014 | 0.647–0.684 |
| mutation | logistic_regression | 0.572 ± 0.024 | 0.228 ± 0.009 | 0.539–0.596 |
| protein | random_forest | 0.573 ± 0.018 | 0.270 ± 0.018 | 0.561–0.601 |
| rna | random_forest | 0.684 ± 0.008 | 0.436 ± 0.008 | 0.677–0.693 |

### Same held-out split results

| Scope | Fixed model | Test accuracy | Balanced accuracy | Macro-F1 | ROC-AUC | PR-AUC | Deceased sensitivity | Specificity | TN/FP/FN/TP |
|---|---|---:|---:|---:|---:|---:|---:|---:|---|
| cnv | random_forest | 0.851 | 0.492 | 0.460 | 0.448 | 0.151 | 0.000 | 0.984 | 120/2/19/0 |
| rna | random_forest | 0.879 | 0.642 | 0.673 | 0.579 | 0.348 | 0.316 | 0.967 | 118/4/13/6 |
| mutation | logistic_regression | 0.688 | 0.464 | 0.465 | 0.390 | 0.122 | 0.158 | 0.770 | 94/28/16/3 |
| protein | random_forest | 0.858 | 0.518 | 0.507 | 0.471 | 0.168 | 0.053 | 0.984 | 120/2/18/1 |
| integrated | extra_trees | 0.872 | 0.571 | 0.590 | 0.538 | 0.324 | 0.158 | 0.984 | 120/2/16/3 |

## Threshold analysis

Thresholds were selected from the mean of repeated out-of-fold training probabilities for RNA Random Forest and integrated Extra Trees by maximizing Youden J (sensitivity + specificity − 1). The test set did not choose thresholds. Both 0.5 and the selected threshold are reported in `threshold_analysis.csv`, with OOF and held-out sensitivity, specificity, deceased precision/F1, and confusion counts. The threshold chosen on OOF data can still be optimistic; test rows are a final descriptive check.

| Scope | OOF-selected threshold | Test sensitivity | Test specificity | Test deceased precision | Test F1 | Test TN/FP/FN/TP |
|---|---:|---:|---:|---:|---:|---|
| rna | 0.312 | 0.316 | 0.877 | 0.286 | 0.300 | 107/15/13/6 |
| integrated | 0.159 | 0.421 | 0.639 | 0.154 | 0.225 | 78/44/11/8 |

## Calibration

Sigmoid calibration used 5-fold CV within the training partition only. Test Brier score/log loss and calibration curves for uncalibrated and calibrated probabilities are in `calibration_results.csv` and `figures/`. Test calibration estimates are noisy with only 19 deceased cases. Calibration was not selected using test performance.

## Held-out uncertainty and permutation tests

Stratified bootstrap percentile intervals use 3000 resamples of frozen test predictions, resampling within alive/deceased strata; this quantifies uncertainty conditional on this particular holdout and fitted models, not training/model-selection variability. See `bootstrap_confidence_intervals.csv`. A fixed-model, held-out label permutation test (5000 permutations) is reported for the prespecified RNA and integrated models; test labels are used for inference only, not model/threshold selection. Results are in `reports/heldout_permutation_tests.csv`.


Test Brier scores: rna uncalibrated=0.120; rna sigmoid_cv_calibrated=0.101; integrated uncalibrated=0.113; integrated sigmoid_cv_calibrated=0.109. Sigmoid CV calibration reduced Brier score on this holdout for both candidates; this is not proof of stable calibration given the small test sample.


Fixed-prediction permutation p-values: rna ROC-AUC=0.579, one-sided p=0.133; integrated ROC-AUC=0.538, one-sided p=0.305. These are conditional tests on the frozen holdout predictions, not external validation.

## Modality contribution

All five modalities use the same training rows, repeated CV folds, and held-out rows. Integrated minus RNA paired-bootstrap differences for ROC-AUC, PR-AUC, balanced accuracy, and macro-F1 are in `reports/integration_minus_rna_paired_bootstrap.json`. Point estimates do not show an integration gain over RNA. Paired held-out bootstrap deltas (integrated minus RNA) were roc_auc -0.042 (95% percentile interval -0.161 to 0.073); pr_auc_average_precision -0.028 (95% percentile interval -0.136 to 0.077); balanced_accuracy -0.070 (95% percentile interval -0.158 to 0.004); macro_f1 -0.082 (95% percentile interval -0.196 to 0.009). Intervals include zero, so this small held-out comparison is inconclusive.

## Feature stability and interpretation

Repeated-fold tree importance was aggregated across 25 training folds for RNA Random Forest and integrated Extra Trees.
Top stable INTEGRATED variables (frequency in top-30 across 25 folds): rs_APOB [rna, 56%], rs_HEPACAM [rna, 48%], rs_ANO3 [rna, 48%], rs_FIGF [rna, 44%], cn_ABCA4 [cnv, 40%], rs_KIAA0408 [rna, 40%], rs_PLIN1 [rna, 40%], rs_SLC7A10 [rna, 40%], rs_ADH1C [rna, 36%], rs_RBP4 [rna, 32%].
Top stable RNA variables (frequency in top-30 across 25 folds): rs_CHRDL1 [rna, 96%], rs_PLIN1 [rna, 96%], rs_SCUBE2 [rna, 88%], rs_ALOX15 [rna, 88%], rs_APOB [rna, 88%], rs_ACVR1C [rna, 84%], rs_ANO3 [rna, 84%], rs_HEPN1 [rna, 76%], rs_LOC283392 [rna, 72%], rs_ADH1B [rna, 72%]. `feature_stability.csv` reports modality, top-30 frequency, mean rank, and overlap with previously selected features. It is a stability screen, not causal evidence or a biomarker claim. Figures are in `figures/`; raw fold ranks are under `feature_importance/`.

## Optional methods

SHAP is not installed in the existing environment; it was not installed because permutation/model-based explanations are feasible and changing the project environment is unnecessary for this phase. QLattice/Feyn was unavailable in Phase 2 (`ModuleNotFoundError: No module named feyn`) and was not pursued further.

## Key limitations

There are 705 observations, only 94 deceased labels, no patient IDs/barcodes, no survival times or follow-up, no censoring information, and no external cohort. Repeated CV and bootstrap quantify only parts of resampling uncertainty and cannot repair unknown patient dependence or selection bias. All outputs concern exploratory classification of recorded vital status only; they are not clinical survival predictions.
