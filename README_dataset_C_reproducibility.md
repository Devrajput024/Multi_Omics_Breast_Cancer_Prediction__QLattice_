# Dataset C reproducibility handover

This package reruns the final verified Dataset C experiments from the project CSV. It does not use saved metrics to create results and does not modify the source CSV.

## Input and scope

The default input is `dataset_c_extracted/data.csv`, resolved relative to this script. The script stops with an explicit error if that exact file is missing or if its target, modality prefixes, feature counts, or validated shape differ. It does not search for or substitute another dataset. An explicit input may be supplied with `--data`, but it must match the verified Dataset C schema.

The target is `vital.status`: `0 = Alive`, `1 = Deceased`. The script reports exploratory classification of recorded vital status. It does not predict survival time, has no follow-up/censoring data, and cannot establish patient-level independence because patient IDs/barcodes are absent. Results are not clinically validated.

## Environment and installation

The verified execution environment used Python 3.11.9 and the versions pinned in `requirements.txt` (NumPy, pandas, scikit-learn, Matplotlib, joblib and threadpoolctl). Use Python 3.11 for closest reproduction.

From PowerShell in the project root:

```powershell
py -3.11 -m venv .venv-dataset-c
.\.venv-dataset-c\Scripts\Activate.ps1
python -m pip install --upgrade pip
python -m pip install -r requirements.txt
python dataset_C_reproducibility.py
```

Or call the installed Python directly:

```powershell
py -3.11 dataset_C_reproducibility.py
```

The optional arguments are `--data <exact CSV path>` and `--output <results folder>`. By default outputs go to `dataset_C_reproducibility_results/` beside the script.

## Reproduced methodology

- Stratified 80/20 row split with random seed `20261110`; the split has 564 training and 141 holdout rows for the verified source.
- Nested stratified cross-validation: 3 outer folds × 3 inner folds. Preprocessing is fitted inside each training fold. The inner folds choose among class-weighted Logistic Regression, Random Forest, Extra Trees, RBF SVM, and HistGradientBoosting. These five model families were used in the verified final optimization.
- Selection uses average rank across ROC-AUC, PR-AUC (average precision), balanced accuracy, macro-F1, and recall for deceased cases. Equal scores use stable model-order tie-breaking, matching the verified project. Holdout scores are descriptive and never select a model.
- Class-specific predictions treat deceased (`1`) as positive. Confusion matrices use label order `[0, 1]`, corresponding to TN, FP, FN, TP.
- CNV values use most-frequent imputation and one-hot encoding with unseen-category handling. RNA/protein use median imputation and standard scaling. Mutation uses median imputation and a training-fitted constant-feature filter. Class weights address imbalance. No synthetic resampling is used.
- Original CNV columns with identical observed vectors are retained, matching the verified primary modelling procedure; duplicate observed vectors do not establish biological identity.
- Feature selection is separately evaluated in five stratified training-only folds using mutual information, ANOVA SelectKBest, L1 regularization, and tree importance. The first three requested methods are included, and ANOVA is retained because it was part of the verified project experiment.
- Native numeric-library threads are limited to one during model fits for stable execution on Windows. On scikit-learn 1.9, RBF SVM may print a deprecation warning for `SVC(probability=True)`; the archived final experiment used this estimator setting, which is preserved for methodological reproduction.

If a candidate fails on a different machine, the script prints the exception, writes it to `model_fit_errors.csv` and `execution_log.txt`, and excludes incomplete folds from model ranking. It does not silently substitute an estimator. Missing required packages cause a setup message and exit.

## Outputs

`dataset_C_reproducibility_results/` contains:

- `dataset_summary.csv`, `feature_names_by_modality.csv`, `highly_correlated_pairs.csv`, and `cnv_duplicate_feature_vectors.csv`.
- `candidate_cv_scores.csv`, `candidate_inner_fold_scores.csv`, `nested_fold_results.csv`, `nested_cv_scope_summary.csv`, `model_comparison.csv`, and `dummy_baseline.csv`.
- `selected_features.csv`, `selected_features_by_fold.csv`, `feature_selection_summary.csv`, `feature_selection_fold_metrics.csv`, and `classification_reports.csv`.
- `confusion_matrices/`, `evaluation_plots/`, `trained_models/`, and saved holdout/OOF predictions.
- `execution_log.txt`, `run_metadata.json`, and `final_summary.txt`.

The complete feature-name printout, CV candidate scores, classification reports, final comparison and run details are mirrored in `execution_log.txt`.

## Interpreting differences from earlier results

The script computes new metrics directly from the specified CSV, split, folds and estimators. It does not copy values from prior CSVs. Minor numerical differences can occur across Python/library versions or hardware math libraries. A different source hash, schema, package versions, missing candidates or exceptions will be recorded in the log and metadata. The `model_comparison.csv` separates nested outer-CV estimates from held-out metrics.
