# Dataset C complete experiment audit

**Scope:** Dataset C breast-cancer multi-omics project files present in this workspace. This is a read-only audit of saved code, notebooks, reports, models, predictions, and result tables. No experiment was run, retrained, or changed for this audit. Dataset A and Dataset B are excluded.

**Audit date:** 2026-11-11 (workspace date). Paths below are relative to the project root unless stated otherwise.

## Executive summary

The project contains three distinguishable generations of conventional machine-learning experiments: (1) an initial fixed holdout / training-only 5-fold CV suite with 23 saved pipelines, (2) repeated cross-validation and secondary validation analyses, and (3) a later five-scope nested-CV selection study with 25 candidates per scope. The later final-optimization reports select RNA Random Forest as the overall model by a training-only composite rank; its holdout ROC-AUC is 0.5794 and average precision is 0.3484. The integrated Random Forest has the highest selected-scope holdout ROC-AUC (0.6053), but the study's stated selection rule is CV composite rank, not heldout ROC-AUC. Neither number is evidence of clinical utility.

The ordinary scikit-learn experiments have substantial saved artifacts and prediction-level checks. The final-optimization output verification records 25 candidate models, zero fit errors, five saved final pipelines, 141 predictions per selected scope, and confusion matrices consistent with the saved metrics. The earlier 23-model suite has an explicit reconciliation report showing all 23 prediction sets use the same holdout and match summary metrics exactly.

QLattice/Feyn is a separate exploratory line. The notebooks contain executed cells and narrative metric claims, and two serialized `.model` graphs are present, but their training code depends on a hosted QLattice/Feyn service or an unavailable historical runtime, has no prediction/metric table for independent checking, and does not preserve the split indices and full execution environment needed to reproduce the stated scores. Treat the Feyn claims as historical, unverified notebook results, not as a verified final experiment. No SHAP results were found. XGBoost and LightGBM are documented as skipped because they were unavailable; there are no completed runs or metrics for them.

## Dataset identity and limits

The primary model input is `dataset_c_extracted/data.csv`: **705 rows Ã— 1,937 columns**, comprising 1,936 molecular predictors and `vital.status`. The target mapping in the scripts and reports is **0 = Alive (611; 86.67%)** and **1 = Deceased (94; 13.33%)**. Predictor prefixes/configuration are CNV `cn_` (860), RNA `rs_` (604), mutation `mu_` (249), and protein `pp_` (223). These sum to 1,936. The primary matrix is reported as complete (no missing cells) and without duplicate rows. The final quality analysis reports 230 near-constant predictors at its chosen 95% dominant-value criterion, 1,116 feature pairs with absolute Pearson correlation at least 0.95, and 62 exact duplicate-vector CNV groups (99 redundant columns beyond one representative per group). The initial experiment also tested original versus deduplicated CNV predictors; the later primary final-optimization methodology retained aliases.

`dataset_c_extracted/brca_data_w_subtypes.csv` and recovered `brca_data_w_meta.csv` are different augmented tables, not the primary 705 Ã— 1,937 model matrix. The audit/recovery documents say they contain additional clinical/subtype fields and missingness and did not substitute them into the principal predictive experiments. Do not confuse them with `data.csv`.

The available records do not include usable patient identifiers, survival duration, or censoring information. Accordingly, the saved target supports only the recorded vital-status classification task; it does not establish survival prediction, breast-cancer-specific mortality, independent patient-level validation, or clinical validation. Provenance investigations characterize the source as a likely TCGA-BRCA-derived subset, but exact sample linkage/release provenance is not established from the files available here.

## Canonical experiment generations and methodology

### A. Initial fixed-holdout suite (23 models)

**Code:** `results_dataset_C/run_experiments.py`. **Primary outputs:** `results_dataset_C/metrics_summary.csv`, `cross_validation_results.csv`, `predictions/`, `models/`, and modality subfolders. The script's documented split is stratified 80/20, seed 20261110: 564 train rows (489 Alive / 75 Deceased) and 141 test rows (122 / 19). Training-only validation uses stratified 5-fold CV. Per-modality and integrated pipelines fit preprocessing inside the model/CV pipeline; class weights are balanced. Baselines are Logistic Regression (`C=0.1`, liblinear in the saved configurations), Random Forest (100 trees), and Extra Trees (100 trees). The suite contains 15 modality/model combinations, two CNV duplicate-handling variants, and six Logistic Regression feature-selection variants (MI, L1, tree importance on RNA and integrated), for 23 saved prediction/model outputs. RNA and integrated feature selection is pipeline-contained; it is not fitted on the test set.

The experiment log records an initial failed CNV deduplication attempt (`KeyError` from incorrect column indexing) and a Windows worker permission failure. The log then documents corrected execution. The final saved artifacts are complete; the initial failures are not the status of the final outputs.

Selected result summary below reports training CV ROC-AUC mean Â± SD, then heldout ROC-AUC / average precision, and heldout confusion tuple TN/FP/FN/TP. Values are read from the saved result summaries and reconciled predictions; these are not fresh results.

| Configuration | CV ROC-AUC mean Â± SD | Test ROC-AUC / AP | TN/FP/FN/TP |
|---|---:|---:|---:|
| CNV Logistic Regression | .510 Â± .071 | .535 / .150 | 99/23/16/3 |
| CNV Random Forest | .531 Â± .034 | .448 / .151 | 120/2/19/0 |
| CNV Extra Trees | .523 Â± .061 | .403 / .111 | 121/1/19/0 |
| RNA Logistic Regression | .648 Â± .055 | .626 / .310 | 92/30/10/9 |
| RNA Random Forest | .690 Â± .035 | .579 / .348 | 118/4/13/6 |
| RNA Extra Trees | .683 Â± .066 | .546 / .358 | 120/2/14/5 |
| Mutation Logistic Regression | .596 Â± .095 | .390 / .122 | 94/28/16/3 |
| Mutation Random Forest | .591 Â± .065 | .466 / .150 | 97/25/15/4 |
| Mutation Extra Trees | .588 Â± .066 | .459 / .125 | 102/20/17/2 |
| Protein Logistic Regression | .566 Â± .049 | .494 / .196 | 87/35/14/5 |
| Protein Random Forest | .605 Â± .070 | .471 / .168 | 120/2/18/1 |
| Protein Extra Trees | .564 Â± .068 | .438 / .174 | 121/1/18/1 |
| Integrated Logistic Regression | .651 Â± .060 | .604 / .338 | 110/12/14/5 |
| Integrated Random Forest | .665 Â± .082 | .605 / .356 | 119/3/13/6 |
| Integrated Extra Trees | .684 Â± .088 | .538 / .324 | 120/2/16/3 |
| CNV dedup: original-feature LR | .510 Â± .071 | .535 / .150 | 99/23/16/3 |
| CNV dedup: deduplicated LR | .513 Â± .068 | .524 / .149 | 100/22/16/3 |
| RNA LR + mutual information | .612 Â± .038 | .661 / .348 | 75/47/9/10 |
| RNA LR + L1 | .652 Â± .060 | .630 / .242 | 90/32/10/9 |
| RNA LR + tree importance | .662 Â± .068 | .551 / .230 | 74/48/10/9 |
| Integrated LR + mutual information | .594 Â± .032 | .621 / .309 | 86/36/8/11 |
| Integrated LR + L1 | .658 Â± .041 | .560 / .261 | 92/30/13/6 |
| Integrated LR + tree importance | .630 Â± .072 | .578 / .249 | 90/32/12/7 |

The detailed `metrics_summary.csv` also contains accuracy, deceased-class precision/recall/F1, specificity, macro-F1, balanced accuracy, ROC-AUC and PR-AUC. The prediction reconciliation (`advanced_validation/reports/prediction_reconciliation.json`) verifies all 23 prediction CSVs have 141 rows, the same test indices and labels, `[0,1]` confusion order, and exact metric-summary matches (maximum absolute metric difference 0.0). `final_optimization/reports/previous_methodology_audit.json` similarly records all 23 saved models/predictions, matching confusion counts, preprocessing and selectors inside pipelines, and CV restricted to the training split. One caveat is explicitly recorded there: a previous report described some heldout-ROC rankings after looking at test results; such rankings are post-hoc and optimistic if treated as model selection. The later optimization uses training-only CV selection.

### B. Repeated CV and secondary validation (Phase 3)

**Code:** `results_dataset_C/advanced_validation/run_advanced_validation.py`. It performs repeated stratified 5-fold CV (5 repeats) on the 564-row training split. Candidate architectures were frozen from the preceding scope selections: CNV Random Forest, RNA Random Forest, mutation Logistic Regression, protein Random Forest, and integrated Extra Trees. It saves repeated CV summaries and out-of-fold predictions, then separate holdout threshold, calibration, bootstrap and permutation analyses. It is not the same estimate as the later 3Ã—3 nested CV below.

| Scope/model | Repeated-CV ROC-AUC, mean Â± SD across repeat means | Holdout ROC-AUC / AP | Holdout TN/FP/FN/TP |
|---|---:|---:|---:|
| CNV / Random Forest | .5145 Â± .0160 | .4482 / .1508 | 120/2/19/0 |
| RNA / Random Forest | .6837 Â± .0078 | .5794 / .3484 | 118/4/13/6 |
| Mutation / Logistic Regression | .5721 Â± .0239 | .3904 / .1222 | 94/28/16/3 |
| Protein / Random Forest | .5732 Â± .0180 | .4707 / .1684 | 120/2/18/1 |
| Integrated / Extra Trees | .6712 Â± .0161 | .5382 / .3243 | 120/2/16/3 |

The same outputs report a threshold analysis selected on training OOF data with Youden's J: RNA threshold 0.312397 (OOF sensitivity .4133, specificity .9018, F1 .4026); integrated Extra Trees threshold 0.158532 (OOF sensitivity .6267, specificity .6912, F1 .3443). At those thresholds, holdout RNA has sensitivity .3158, specificity .8770, F1 .3000, TN/FP/FN/TP 107/15/13/6; integrated has sensitivity .4211, specificity .6393, F1 .2254, 78/44/11/8. The ordinary 0.5 holdout scores are shown in the table above. Thresholds are not selected on test labels.

Sigmoid calibration is reported as training-split 5-fold CV fitted. On holdout RNA, Brier score / log loss change from .120344 / .409673 uncalibrated to .100893 / .356597 calibrated. Integrated Extra Trees changes from .113345 / .400543 to .108571 / .379231. These are internal holdout calibration results, not clinical calibration.

The saved stratified bootstrap (3,000 resamples) gives RNA RF holdout ROC-AUC .5794, 95% interval [.4129, .7442], and deceased sensitivity .3158, interval [.1053, .5263]. Integrated Extra Trees ROC-AUC is .5382, interval [.3723, .6980]. Paired bootstrap mean differences (integrated minus RNA) are ROC-AUC âˆ’.04205 [âˆ’.16135, .07293], AP âˆ’.02804 [âˆ’.13595, .07671], balanced accuracy âˆ’.06956 [âˆ’.15789, .00410], and macro-F1 âˆ’.08184 [âˆ’.19617, .00930]. Intervals condition on the fixed heldout predictions and do not represent new independent cohorts.

The 5,000-permutation holdout-label inference reports RNA RF p=.132773 (null 95th percentile .615207) and integrated Extra Trees p=.305139 (null 95th percentile .616695); the reports mark test labels as inference-only and say the model/threshold were not selected using test. These procedures and their outputs are present, but are secondary analyses on the same small fixed holdout, not external validation.

### C. Final nested-CV model-selection experiment

**Code:** `results_dataset_C/final_optimization/run_optimization.py`; detailed report `final_optimization/reports/final_optimization_report.md`; main tables under `final_optimization/tables/`. It uses a stratified 80/20 outer holdout with the same seed/data split and 3 outer Ã— 3 inner stratified CV on the 564-row training portion. For each of five scopes, it compares five algorithm families: balanced Logistic Regression, balanced Random Forest, balanced Extra Trees, RBF SVC, and HistGradientBoosting (25 scope/model candidates). Selection uses a mean composite rank across ROC-AUC, average precision, balanced accuracy, macro-F1, and deceased recall from inner CV. The holdout is evaluated after selection and is not the selection metric.

Parameters documented in the script/config: LR `C=0.1`, liblinear, balanced, `max_iter=2500`; RF and Extra Trees `n_estimators=100`, `class_weight=balanced`, `max_features=sqrt`, `min_samples_leaf=2`, `n_jobs=1`; SVC RBF `C=1`, `gamma=scale`, balanced, probability enabled; HistGradientBoosting `max_iter=100`, learning rate `.08`, `max_leaf_nodes=15`, `l2_regularization=1`, balanced, early stopping disabled. CNV uses most-frequent imputation and one-hot encoding with unknown categories ignored; mutation uses median imputation and variance filtering; RNA/protein use median imputation and standardization; integrated processing is modality-aware. Transformations and feature selectors are fitted within folds. The main final-optimization pipeline retains exact duplicate CNV aliases, consistent with its declared configuration.

The per-scope CV-selected models and final test scores are:

| Scope / selected algorithm | Nested outer ROC-AUC mean Â± SD | Outer AP mean Â± SD | Test accuracy | Test deceased precision / recall | Specificity | F1 / macro-F1 | Test ROC-AUC / AP | TN/FP/FN/TP |
|---|---:|---:|---:|---:|---:|---:|---:|---:|
| CNV / LR | .5428 Â± .0190 | .1500 Â± .0055 | .7234 | .1154 / .1579 | .8115 | .1333 / .4844 | .5345 / .1503 | 99/23/16/3 |
| RNA / RF | .6318 Â± .0734 | .4050 Â± .0528 | .8794 | .6000 / .3158 | .9672 | .4138 / .6733 | .5794 / .3484 | 118/4/13/6 |
| Mutation / LR | .5782 Â± .0752 | .1996 Â± .0147 | .6879 | .0968 / .1579 | .7705 | .1200 / .4652 | .3904 / .1222 | 94/28/16/3 |
| Protein / RF | .5050 Â± .1382 | .2082 Â± .1123 | .8582 | .3333 / .0526 | .9836 | .0909 / .5070 | .4707 / .1684 | 120/2/18/1 |
| Integrated / RF | .6708 Â± .0677 | .4310 Â± .0756 | .8865 | .6667 / .3158 | .9754 | .4286 / .6828 | .6053 / .3561 | 119/3/13/6 |

The overall best by the predefined training-CV composite rank is **RNA Random Forest** (rank 1.0), not the model with best heldout ROC-AUC. Nested CV estimates scope/model selection stability; test estimates are one fixed internal holdout. CV mean deceased recall/specificity/macro-F1 are in `nested_cv_summary.csv`; the table above includes ROC/AP and other fold means. The final selected-scope comparison table includes all listed metrics. As demonstrated by the low deceased recall, high overall accuracy here is partly driven by class imbalance.

The nested candidate model summary contains all 25 candidates and their CV metrics/ranks. Winners by scope: CNV LR (mean rank 1.8), RNA RF (2.0), mutation LR (1.4), protein RF (1.8), integrated RF (1.4). HistGradientBoosting was attempted as a candidate family and completed; verification records five such candidates and zero fit errors. This does not mean XGBoost or LightGBM ran: they are different libraries and were skipped.

Feature selection was evaluated in RNA and integrated scopes using mutual information, ANOVA, L1, and Extra Trees feature importance, nested within training folds. Per-fold selectors targeted 50 encoded features; the saved unique source-feature counts are RNA ANOVA 70, L1 126, MI 124, tree importance 138; integrated ANOVA 71, L1 133, MI 187, tree importance 179. The saved CSVs include feature names, selection frequencies/scores and modality labels where available. Fold-native tree importance and validation-fold permutation importance outputs are available. They are predictive feature rankings, not validated biomarkers.

`final_optimization/reports/output_verification.json` reports PASS, source shape 705Ã—1,937, class counts 611/94, scope feature counts, 564/141 split, 25 candidates, zero candidate fit errors, five saved selected pipelines, 12 PNG and 12 PDF figures, and test confusion matrices that agree with the table above. These are saved verification results, not a rerun for this audit.

### D. Separate `dataset_C_reproducibility.py` execution outputs

The root script and `dataset_C_reproducibility_results/` are a later portability/reproduction deliverable. Its saved metadata declares the same dataset/split, all five selected scope models, candidate and nested-fold outputs, selected features, reports, figures and joblib pipelines. Its final selected metrics reproduce the final-optimization holdout ROC-AUC/AP exactly (0.0 absolute differences on all scopes), and the five scope winners reproduce. The script's `reproduction_verification.json`, however, also records a non-match on protein nested outer folds: one selected-algorithm mismatch and maximum ROC-AUC difference 0.2638 (all outer folds do not match), while the other four scopes' nested outer summary differences are zero. Thus test scores/winners match the saved canonical table, but the stored verification does **not** support a claim that every protein outer-fold result was independently reproduced. Pipeline prediction checks in that verification are positive for all five scope models. Report this discrepancy rather than hiding it.

## Feature importance and explainability

Saved conventional explainability outputs include `results_dataset_C/explainability/top_feature_importance.csv`, phase-2 selected-feature tables, advanced-validation feature stability and repeated-fold top-30 importance, and final-optimization selected-by-fold, selector frequency, outer-fold importance/stability and permutation-importance CSVs. They are based on MI/L1/tree importance and permutation methods in the saved scikit-learn workflows. These CSVs carry the actual selected names/scores and should be read directly for feature-level detail; this audit does not reproduce hundreds of individual rows. No SHAP values, SHAP tables, SHAP summary/dependence plots, or SHAP-generating completed experiment were found. `optional_tools.json` records SHAP unavailable/not installed. Therefore SHAP analysis is **not completed** in the saved project.

## QLattice / Feyn review (high scrutiny)

Relevant files: `dataset_c_extracted/Interpretable_Multi_Omics_Integration_for_Breast_Cancer_Severity_Prediction_using_QLattice_Symbolic_Regression.ipynb`, its exact inspection copy under `.multiomics_audit_inspect/`, `dataset_c_extracted/model1_200its_mc7.model`, `dataset_c_extracted/model_200its_md2_bic_0.model`, and `dataset_c_recovery/notebooks__bc_case_multi.ipynb`.

The first notebook has executed cells and uses `data.csv` (705 Ã— 1,937), stratifies target `vital.status` with `test_size=.34`, `random_state=42`, then calls Feyn/QLattice classification `auto_run` with 25 epochs and `max_complexity=3`. The notebook narrative states a two-input expression model (`rs_CA3`, `rs_CALB2`) with training accuracy .886, ROC-AUC .629, precision .737, recall .226; test accuracy .871, ROC-AUC .656, precision .556, recall .156. It later describes a second four-input graph (`mu_TNXB`, `rs_KRT23`, `rs_APOB`, `cn_ANKRD30B`) with addition, Gaussian, and multiplication operators and reports training/test accuracy .895/.875, AUC .648/.660, precision .810/.667, recall .274/.125. The notebook loads/saves `/content/model1_200its_mc7.model`; the second local serialized model is not clearly tied to a reproducible training cell. These claims are narrative notebook output, not saved score tables or prediction files; no independent metric recomputation is possible from stored predictions.

The notebook records a failed attempt to install Feyn 2.1.4 (version unavailable) and a successful notebook-runtime install of Feyn 3.4.1. The current project environment audit says Feyn is not installed. It relies on Colab `/content` paths and a remote QLattice/community service; the saved artifacts do not include portable training data references, fold/split indices, prediction files, complete environment lock, or a full local replay workflow. The `.model` artifacts were inspected as serialized graphs, not executed. Consequently, the cells did execute historically, but the metrics/models are **not independently verified or reproducible from this workspace as supplied**. Do not describe these scores as a validated or canonical final result.

The companion notebook `dataset_c_recovery/notebooks__bc_case_multi.ipynb` is a separate recovered experiment. It imports helper functions from an unavailable `functions` module and contains historical notebook outputs for 5-fold validation. Those outputs report QLattice ROC-AUC validation mean .63823 Â± .07407 (train .67124 Â± .03333), and a second configuration validation .633219 Â± .041365 (train .765106 Â± .033119). A notebook benchmark table shows QLattice .633, LASSO .568, ElasticNet .611, Random Forest all-features .604 (MI top-10 .628, F-test top-10 .558, MI top-5 .594, F-test top-5 .541, LASSO-selected .604), and Gradient Boosting .570 (MI top-10 .621, F-test top-10 .541, MI top-5 .591, F-test top-5 .590, LASSO-selected .594). These are captured notebook outputs, not regenerated scores. Missing imported helpers, service dependency, and exact runtime prevent reproduction; treat as **historical/unverified** and separate from the 705-row primary final pipeline. No patient-level holdout linkage or prediction artifacts were found for this notebook's claimed estimates.

## XGBoost and LightGBM

The saved optional-tool/environment status and final-optimization documentation say XGBoost (`xgboost`) and LightGBM (`lightgbm`) were unavailable and skipped. No XGBoost/LightGBM model files, predictions, parameter/result tables, or plots were found. Status for both is **not run / skipped**, not failed performance and not completed. There are no numerical results to report. HistGradientBoosting in the final nested study did run, but is not XGBoost or LightGBM.

## Calibration, permutation, bootstrap, and CV completeness

| Analysis | Status and saved evidence | Reproducibility assessment |
|---|---|---|
| Stratified 5-fold train CV (initial suite) | Completed; `cross_validation_results.csv` and per-model metrics | Reproducible from `run_experiments.py`, source matrix, fixed seed and dependencies; saved pipelines/preds cross-check metrics. |
| Repeated 5Ã—5-fold train CV | Completed; `advanced_validation/repeated_cv_results.csv` and five OOF files | Code and outputs present; same train split and seeds are recorded. |
| Nested 3Ã—3 CV and candidate selection | Completed; `final_optimization/tables/nested_outer_results.csv`, `nested_cv_summary.csv`, `model_selection_summary.csv` | Code and outputs present; 25 candidates, no recorded fit errors. Separate reproduction output has the protein outer-fold mismatch described above. |
| MI/L1/tree/ANOVA feature selection | Completed within training folds for RNA/integrated final optimization; selected-by-fold and frequency tables | Scores/names are saved; folds and code available. Initial phase-2 had MI/L1/tree LR variants. |
| Threshold selection | Completed on training OOF, RNA and integrated | Thresholds and OOF/test summaries saved; no test-set threshold selection recorded. |
| Sigmoid calibration | Completed for RNA RF/integrated ET using training-only CV calibration | Metrics and curves saved; external calibration not available. |
| Bootstrap uncertainty | Completed, 3,000 resamples of fixed holdout predictions | Intervals saved; conditional on one holdout. |
| Permutation inference | Completed, 5,000 permutations, two selected models | Results saved; test labels used for inference only per report. |
| SHAP | Not found / unavailable | No reproducible SHAP execution or outputs. |
| XGBoost, LightGBM | Skipped/unavailable | No code outputs or metrics. |
| Feyn/QLattice | Historical notebook run and serialized graphs found; independent verification incomplete | Not reproducible from included artifacts/runtime as supplied. |

## Discrepancies and caveats

1. The fixed-holdout Phase-2 winners and final nested-CV winners differ by design: Phase 2 compared preselected LR/RF/ET configurations and included feature-selection variants; final optimization added SVC and HistGradientBoosting and used a composite nested-CV rank. Their CV estimates are not interchangeable.
2. Repeated-CV Phase 3 compares frozen architectures and reports RNA RF repeated mean ROC-AUC .6837; final nested CV reports RNA RF scope nested outer mean .6318. These use different resampling/selection protocols, so differences do not alone indicate corruption.
3. The reproducibility-run check matches test metrics and winners but reports protein nested-fold mismatch (max ROC-AUC absolute difference .2638; one algorithm mismatch). This is a genuine stored inconsistency and should be resolved before claiming exact nested-fold reproduction.
4. Existing reports sometimes foreground accuracy despite strong class imbalance. For the selected RNA RF, 87.9% accuracy coexists with deceased recall 31.6%; protein recall is only 5.3%. Report minority-class metrics alongside accuracy.
5. A 141-row holdout has only 19 deceased cases. Bootstrap/permutation intervals and observed scores are correspondingly uncertain. There is no independent cohort/patient-level validation.
6. The main model study has no evidence of fitting preprocessing or selecting features on the complete dataset; available audits state pipeline-contained transforms and training-only fitting. The earlier report's post-hoc heldout ranking caveat remains.
7. Feyn narrative scores have no saved prediction probabilities/labels to check and rely on artifacts/services unavailable locally. They must remain explicitly unverified.

## Complete relevant file inventory

Inventory includes Dataset C source data, runnable code/notebooks, saved experiment outputs, reports, handover/reference packages, and Dataset C data-quality/provenance investigations. Directories named `__pycache__` and transient cache files are not research outputs and are excluded. `Dataset_C_Paper_Materials` contains a nested copied package and an extracted ZIP; these are derivative copies, not separate experiments. `Dataset_C_Supervisor_Package` is a portable handover copy, not a new experiment.

### Root and source materials

- `dataset_C_reproducibility.py`; `README_dataset_C_reproducibility.md`; `requirements.txt`.
- `dataset_c_extracted/data.csv`; `dataset_c_extracted/brca_data_w_subtypes.csv`; `dataset_c_extracted/Interpretable_Multi_Omics_Integration_for_Breast_Cancer_Severity_Prediction_using_QLattice_Symbolic_Regression.ipynb`; `dataset_c_extracted/model1_200its_mc7.model`; `dataset_c_extracted/model_200its_md2_bic_0.model`; `dataset_c_extracted/Multi_Omics_Breast_Cancer_Prediction__QLattice readme.md`.
- Exact inspection-copy set under `.multiomics_audit_inspect/`: `data.csv`, `brca_data_w_subtypes.csv`, the QLattice notebook, both `.model` files, and the QLattice readme.
- `dataset_c_recovery/data__brca_data_w_meta.csv`; `data__brca_data_w_subtypes.csv`; `notebooks__bc_case_multi.ipynb`; `README.md`; `assay_preprocessing_provenance.md`; `clinical_followup_investigation.md`; `histology_discrepancy.md`; `identifier_recovery.md`; `recovery_feasibility_conclusion.md`; `repository_inspection.md`.
- `results_dataset_C/run_experiments.py`; `results_dataset_C/experiment_log.txt`; `dataset_C_final_report.md`; `dataset_validation.md`; `metrics_summary.csv`; `cross_validation_results.csv`; `feature_summary.csv`; `duplicate_alias_map.csv`.

### `results_dataset_C/` output files (full path inventory by subfolder)

- `advanced_validation/`: `run_advanced_validation.py`; `baseline_comparison.csv`; `bootstrap_confidence_intervals.csv`; `calibration_results.csv`; `feature_stability.csv`; `modality_comparison.csv`; `repeated_cv_results.csv`; `threshold_analysis.csv`; `README.md`; OOF predictions for `cnv`, `integrated`, `mutation`, `protein`, `rna`; `predictions/dummy_majority_test_predictions.csv`; `reports/heldout_permutation_tests.csv`; `reports/prediction_reconciliation.json`; `reports/validation_audit.md`; `feature_importance/repeated_fold_top30_importance.csv`; plots under `figures/`.
- `explainability/top_feature_importance.csv`.
- `feature_selection/feature_selection_results.csv`; `selected_features.csv`.
- `figures/`: 72 saved figures from the initial experiment and validation suite.
- `final_optimization/`: `run_optimization.py`; `reports/final_optimization_report.md`; `reports/methodology_audit.md`; `reports/output_verification.json`; `reports/previous_methodology_audit.json`; `reports/optional_tools.json`; `reports/optimization_execution.log`; tables `cnv_category_distribution.csv`, `dummy_baseline_comparison.csv`, `modality_quality_summary.csv`, `model_comparison_master.csv`, `model_comparison_with_dummy.csv`, `model_fit_errors.csv`, `model_selection_summary.csv`, `nested_cv_summary.csv`, `nested_inner_and_importance.csv`, `nested_outer_results.csv`, `optimized_calibration_metrics.csv`, `selected_scope_models_test_metrics.csv`; predictions `dummy_baseline_test_predictions.csv`, `nested_oof_predictions.csv`, `selected_scope_test_predictions.csv`; reproducibility `train_test_indices.csv`; feature-selection files `best_cv_scope_feature_frequency.csv`, `feature_selection_frequency.csv`, `outer_fold_feature_importance.csv`, `outer_fold_feature_stability.csv`, `outer_fold_permutation_importance.csv`, `permutation_importance_summary.csv`, `selected_features_by_fold.csv`, `selector_fold_metrics.csv`, `selector_summary.csv`; `correlation_analysis.csv`, `feature_quality_summary.csv`; 12 PNG + 12 PDF figures under `figures/`.
- `final_research/`: `dataset_C_research_report.md`; `advanced_execution_report.md`; `methodology_audit.md`; `master_metrics_comparison.csv`; `final_selected_features.csv`; reproducibility `central_model_paths.csv`, `feature_manifest.csv`, `model_inventory.csv`, `train_test_indices.csv`; tables `bootstrap_confidence_intervals.csv`, `calibration_results.csv`, `heldout_permutation_tests.csv`, `modality_winner_comparison.csv`, `repeated_cv_results.csv`, `threshold_analysis.csv`; figures under `figures/`.
- `individual_omics/cnv/`, `rna/`, `mutation/`, `protein/`: model prediction CSVs for the algorithms listed in the Phase-2 table, with confusion-matrix tables and associated plots in the scope subfolders. RNA includes the MI/L1/tree Logistic Regression variants.
- `integrated/`: integrated LR/RF/Extra Trees predictions, the MI/L1/tree Logistic Regression variants, confusion matrices and plots.
- `models/`: 23 Phase-2 `.joblib` pipelines and their associated model manifests/metadata.
- `predictions/`: 23 Phase-2 heldout prediction CSVs, corresponding to all 23 configurations listed in the Phase-2 table.
- `reports/cnv_dedup_comparison.csv`; other report/summary files are listed above.

### Reproduction-run files

`dataset_C_reproducibility_results/` contains `candidate_cv_scores.csv`, `candidate_inner_fold_scores.csv`, `classification_reports.csv`, `cnv_duplicate_feature_vectors.csv`, `dataset_summary.csv`, `dummy_baseline.csv`, `dummy_baseline_predictions.csv`, `execution_log.txt`, `feature_names_by_modality.csv`, `feature_selection_fold_metrics.csv`, `feature_selection_summary.csv`, `final_summary.txt`, `heldout_predictions.csv`, `highly_correlated_pairs.csv`, `model_comparison.csv`, `model_fit_errors.csv`, `nested_cv_scope_summary.csv`, `nested_fold_results.csv`, `nested_oof_predictions.csv`, `reproduction_verification.json`, `run_metadata.json`, `selected_features.csv`, `selected_features_by_fold.csv`, `selected_scope_models.csv`, `train_test_indices.csv`, five confusion-matrix CSV/PNG pairs, three evaluation plot PDF/PNG pairs, and five selected-scope trained `.joblib` pipelines with metadata JSON.

### Dataset C audits and handover/report copies

- `dataset_c_audit/`: `dataset_c_audit_complete.md`, `dataset_inventory.txt`, `dataset_structure.csv`, `data_quality_report.txt`, `duplicate_feature_analysis.md`, `feature_category_analysis.txt`, `feature_dictionary.md`, `missing_values.csv`, `provenance_investigation.md`, `qlattice_model_inspection.md`, `research_readiness_report.txt`, `run_audit.py`, `sample_alignment_investigation.md`, `sample_overlap_report.txt`, `second_stage_conclusions.md`, `target_analysis.txt`, `target_definition.md`, and six SVG figures under `figures/`.
- `multiomics_dataset_audit/`: `CLASS_DISTRIBUTION.csv`, `DATASET_ISSUES_AND_RECOMMENDATIONS.md`, `DATASET_OVERVIEW.md`, `FEATURE_SUMMARY.csv`, `MISSING_VALUES.csv`, `PATIENT_ID_MATCHING.md`.
- `new_dataset_audit/`: `clinical_audit.md`, `final_conclusion.md`, `image_audit.md`, `inventory.md`, `molecular_audit.md`, `patient_id_matching.md`, `survival_variable_audit.md`.
- `Dataset_C_Paper_Materials/`: `CONTENTS.md`, root confusion-matrix CSV/PNG copies, and a nested copy `Dataset_C_Paper_Materials/` of final reports, comparison/feature tables, figures (PNG/PDF), and confusion matrices; it also contains `Dataset_C_Paper_Materials_Extracted.zip.zip`. These are packaging copies of existing verified outputs, not new experiment runs.
- `Dataset_C_Supervisor_Package/`: `dataset_C_reproducibility.py`, `requirements.txt`, `README.md`, `dataset_c_extracted/data.csv`, `reference_outputs/REFERENCE_ONLY_README.md`, `dataset_summary.csv`, `model_comparison.csv`, and `final_summary.txt`.
- Root archives: `Dataset_C_Paper_Materials.zip`; `Dataset_C_Supervisor_Package.zip`.

## Final completed / attempted status

**Completed with saved code and auditable outputs:** initial 23-model fixed-holdout suite; duplicate-CNV sensitivity comparison; MI/L1/tree selection variants; repeated stratified CV with OOF predictions; later nested 3Ã—3 model selection over five scopes and five algorithms; fold-contained feature selection/importance; training-OOF threshold analysis; sigmoid calibration; fixed-holdout bootstrap and permutation analyses; dataset QC and final figures/reports.

**Attempted historically but not currently reproducible/independently verified:** Feyn/QLattice notebook runs and serialized symbolic models; the recovered companion QLattice notebook depends on missing code and service/runtime. Their narrative metrics are recorded above with their provenance caveat.

**Skipped / no results:** SHAP, XGBoost and LightGBM. SHAP unavailable/no artifacts; XGBoost and LightGBM unavailable and explicitly skipped. No numeric scores should be inferred for any of these.

**Reproduction discrepancy:** the saved portability verification confirms the five holdout winners/metrics but flags one protein outer-fold algorithm mismatch and max ROC-AUC difference 0.2638; nested protein exact reproduction therefore remains unresolved in the artifacts.

### Full Phase-2 per-model heldout metrics

The preceding phase table gives CV ROC-AUC and key discrimination summaries. For completeness, the saved test metrics for all 23 models are reproduced here (positive class = Deceased; confusion values are TN/FP/FN/TP). Values are copied from metrics_summary.csv, not recalculated.

| Configuration | Accuracy | Precision | Recall | Specificity | F1 | Macro-F1 | ROC-AUC | PR-AUC | TN/FP/FN/TP |
|---|---:|---:|---:|---:|---:|---:|---:|---:|---:|
| cnv / logistic_regression | 0.723 | 0.115 | 0.158 | 0.811 | 0.133 | 0.484 | 0.535 | 0.150 | 99/23/16/3 |
| cnv / random_forest | 0.851 | 0.000 | 0.000 | 0.984 | 0.000 | 0.460 | 0.448 | 0.151 | 120/2/19/0 |
| cnv / extra_trees | 0.858 | 0.000 | 0.000 | 0.992 | 0.000 | 0.462 | 0.403 | 0.111 | 121/1/19/0 |
| rna / logistic_regression | 0.716 | 0.231 | 0.474 | 0.754 | 0.310 | 0.566 | 0.626 | 0.310 | 92/30/10/9 |
| rna / random_forest | 0.879 | 0.600 | 0.316 | 0.967 | 0.414 | 0.673 | 0.579 | 0.348 | 118/4/13/6 |
| rna / extra_trees | 0.887 | 0.714 | 0.263 | 0.984 | 0.385 | 0.661 | 0.546 | 0.358 | 120/2/14/5 |
| mutation / logistic_regression | 0.688 | 0.097 | 0.158 | 0.770 | 0.120 | 0.465 | 0.390 | 0.122 | 94/28/16/3 |
| mutation / random_forest | 0.716 | 0.138 | 0.211 | 0.795 | 0.167 | 0.498 | 0.466 | 0.150 | 97/25/15/4 |
| mutation / extra_trees | 0.738 | 0.091 | 0.105 | 0.836 | 0.098 | 0.472 | 0.459 | 0.125 | 102/20/17/2 |
| protein / logistic_regression | 0.652 | 0.125 | 0.263 | 0.713 | 0.169 | 0.475 | 0.494 | 0.196 | 87/35/14/5 |
| protein / random_forest | 0.858 | 0.333 | 0.053 | 0.984 | 0.091 | 0.507 | 0.471 | 0.168 | 120/2/18/1 |
| protein / extra_trees | 0.865 | 0.500 | 0.053 | 0.992 | 0.095 | 0.511 | 0.438 | 0.174 | 121/1/18/1 |
| integrated / logistic_regression | 0.816 | 0.294 | 0.263 | 0.902 | 0.278 | 0.586 | 0.604 | 0.338 | 110/12/14/5 |
| integrated / random_forest | 0.887 | 0.667 | 0.316 | 0.975 | 0.429 | 0.683 | 0.605 | 0.356 | 119/3/13/6 |
| integrated / extra_trees | 0.872 | 0.600 | 0.158 | 0.984 | 0.250 | 0.590 | 0.538 | 0.324 | 120/2/16/3 |
| cnv_dedup_original / logistic_regression | 0.723 | 0.115 | 0.158 | 0.811 | 0.133 | 0.484 | 0.535 | 0.150 | 99/23/16/3 |
| cnv_dedup_deduplicated / logistic_regression | 0.730 | 0.120 | 0.158 | 0.820 | 0.136 | 0.488 | 0.524 | 0.149 | 100/22/16/3 |
| rna / logistic_regression + mutual_information | 0.603 | 0.175 | 0.526 | 0.615 | 0.263 | 0.496 | 0.661 | 0.348 | 75/47/9/10 |
| rna / logistic_regression + l1_logistic | 0.702 | 0.220 | 0.474 | 0.738 | 0.300 | 0.555 | 0.630 | 0.242 | 90/32/10/9 |
| rna / logistic_regression + tree_importance | 0.589 | 0.158 | 0.474 | 0.607 | 0.237 | 0.478 | 0.551 | 0.230 | 74/48/10/9 |
| integrated / logistic_regression + mutual_information | 0.688 | 0.234 | 0.579 | 0.705 | 0.333 | 0.565 | 0.621 | 0.309 | 86/36/8/11 |
| integrated / logistic_regression + l1_logistic | 0.695 | 0.167 | 0.316 | 0.754 | 0.218 | 0.514 | 0.560 | 0.261 | 92/30/13/6 |
| integrated / logistic_regression + tree_importance | 0.688 | 0.179 | 0.368 | 0.738 | 0.241 | 0.522 | 0.578 | 0.249 | 90/32/12/7 |
