# Dataset C final research package

This directory consolidates existing Dataset C results. The assembler does not train or refit any model.

## Main files

- `dataset_C_research_report.md`: final scientific report.
- `master_metrics_comparison.csv`: all 23 Phase 2 evaluation rows plus DummyClassifier, with Phase 2 five-fold CV and Phase 3 repeated-CV values in separate fields.
- `final_selected_features.csv`: Phase 2 selected features and Phase 3 repeated-fold top-30 importance/stability features. Phase 2 per-fold selector support counts were not retained and are blank; do not interpret top-30 importance fold counts as selector support.
- `figures/`: ten publication-style figure families, each exported as 300-dpi PNG and vector PDF; captions are in igures/figure_captions.md.
- `reproducibility/`: package versions, seeds/split indices, feature manifest, model configs and model inventory.

## Reproduction

Use the original inputs (read only):

- `C:\Users\asus\OneDrive\Desktop\Breast part 2\dataset_c_extracted\data.csv`
- `C:\Users\asus\OneDrive\Desktop\Breast part 2\dataset_c_extracted\brca_data_w_subtypes.csv` (used only for separate subtype-label audit; not the primary classifier)

The original Phase 2 training/evaluation script is `C:\Users\asus\OneDrive\Desktop\Breast part 2\results_dataset_C\run_experiments.py`. The Phase 3 validation script is `C:\Users\asus\OneDrive\Desktop\Breast part 2\results_dataset_C\advanced_validation\run_advanced_validation.py`. Exact copies are included in `reproducibility/scripts/` for reference. Existing package versions are recorded in `reproducibility/package_versions.json`; split and CV configuration is in `reproducibility/validation_config.json`. Phase 2 seed: 20261110; stratified 80/20 row split; 5-fold stratified training CV. Phase 3: 5 repeats × 5 folds on training data.

Run the package assembler from the workspace root with `python results_dataset_C/final_research/build_final_package.py`. This rebuilds tables, report, metadata, and figures from the existing results; it does not refit models. To rerun model experiments, use the canonical scripts above and the recorded environment. Do not use a different split if comparing directly with these holdout metrics.

## Trained model locations

No estimators are copied into this package. Phase 2 model artifacts remain under `C:\Users\asus\OneDrive\Desktop\Breast part 2\results_dataset_C\models`; modality-scoped copies are in `C:\Users\asus\OneDrive\Desktop\Breast part 2\results_dataset_C\individual_omics` and `C:\Users\asus\OneDrive\Desktop\Breast part 2\results_dataset_C\integrated`. Exact absolute paths are in `reproducibility/model_inventory.csv` and `reproducibility/central_model_paths.csv`. Training-partition calibrated variants from Phase 3 remain in `C:\Users\asus\OneDrive\Desktop\Breast part 2\results_dataset_C\advanced_validation\calibration_models`.

## Interpretation

The target is recorded vital status, not survival time. Patient IDs and follow-up durations are absent, so patient independence and time-to-event outcomes cannot be checked. Results are exploratory, non-clinical, and not validated for diagnosis or predicting an individual's survival.
