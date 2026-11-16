# Dataset C Phase 2 implementation audit

## Scope

Audited the Phase 2 runner, all saved pipelines, all test prediction files, and metrics. Phase 2 outputs were not overwritten. Advanced methods are in this directory.

## Verified checks

- Target `vital.status` was dropped before predictor construction. No saved fitted preprocessor contains the target: **True**.
- Preprocessing and feature selectors are scikit-learn Pipeline steps; fold transforms/selection are fit by each training fold: **True / True**.
- CV was applied only to the 564-row training partition: **True**. Holdout: 141 rows (122 alive, 19 deceased).
- Every prediction file used identical saved holdout indices and labels: **True**.
- Reconciled predictions and summary metrics: **23/23 pass**, maximum absolute metric difference **0.0**.
- Confusion matrices use labels `[0,1]`; the flattened matrix order is TN, FP, FN, TP. Recomputed confusion counts agree exactly with the saved prediction labels.

## Methodological caveats and corrections

- The prior final report included ?highest held-out test ROC-AUC? model comparisons. This is posthoc winner ranking over multiple candidates and is optimistic if interpreted as unbiased model selection. The saved final pipelines were selected by training CV; advanced comparisons keep those architectures frozen and do not use test metrics to select thresholds/features.
- Phase 2?s CNV dedup code initially used DataFrame bracket indexing with row indices, causing a KeyError. It was corrected to `.iloc` and the original final results were regenerated. The full correction history remains in the Phase 2 experiment log.
- An initial parallel fitting attempt hit a Windows joblib worker PermissionError; it was rerun single-process. No result artifacts from the failed attempts are used as final results.
- Phase 2 test confusion order is correct, class weighting does not alter the target classes, and the saved probabilities reproduce the reported ROC-AUC/PR-AUC.

## Scope limits

There are no patient IDs/barcodes, no follow-up duration, and no censoring information. Splits are row-based; patient independence cannot be verified. These models classify recorded vital status only.
