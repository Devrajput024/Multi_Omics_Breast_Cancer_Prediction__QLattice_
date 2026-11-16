"""Append data-derived final optimization findings to the Dataset C report."""
from pathlib import Path
import json
import pandas as pd

ROOT=Path(__file__).resolve().parents[2]; BASE=ROOT/'results_dataset_C'; OPT=BASE/'final_optimization'
REPORT=BASE/'final_research'/'dataset_C_research_report.md'
def f(x,n=3): return 'NA' if pd.isna(x) else f'{x:.{n}f}'
def main():
    cv=pd.read_csv(OPT/'tables'/'nested_cv_summary.csv').set_index('scope')
    cand=pd.read_csv(OPT/'tables'/'model_selection_summary.csv')
    test=pd.read_csv(OPT/'tables'/'selected_scope_models_test_metrics.csv').set_index('scope')
    dum=pd.read_csv(OPT/'tables'/'dummy_baseline_comparison.csv').iloc[0]
    q=json.loads((OPT/'reports'/'feature_quality_summary.json').read_text())
    fs=pd.read_csv(OPT/'feature_selection'/'feature_selection_frequency.csv')
    perm=pd.read_csv(OPT/'feature_selection'/'permutation_importance_summary.csv')
    cal=pd.read_csv(OPT/'tables'/'optimized_calibration_metrics.csv').set_index('model')
    sel=json.loads((OPT/'reports'/'final_model_selection.json').read_text())
    try: errors=pd.read_csv(OPT/'tables'/'model_fit_errors.csv')
    except pd.errors.EmptyDataError: errors=pd.DataFrame()
    oldboot=pd.read_csv(BASE/'advanced_validation'/'bootstrap_confidence_intervals.csv')
    oldperm=pd.read_csv(BASE/'advanced_validation'/'reports'/'heldout_permutation_tests.csv')
    optional=json.loads((OPT/'tables'/'optional_algorithm_status.json').read_text())
    # The nested outer-fold procedure can choose a different algorithm in each
    # fold; these estimates describe the scope-level model-selection procedure.
    lines=[]
    lines += ['\n\n## Final optimization and methodology verification (2026-11-10)\n',
      '\nThis addendum reports the later optimization run; earlier-phase results and tables were preserved. The new artifacts are in `../final_optimization/`. `model_comparison_with_dummy.csv` is the versioned comparison including the majority baseline.\n',
      '\n### Methodology audit and new evaluation\n',
      '\nThe original modelling CSV has 705 rows and 1,936 predictors. The target is absent from all predictor matrices. All features belong to the four expected prefixes (CNV 860, RNA 604, mutation 249, protein 223). The previous 23 pipelines and 23 prediction files were audited: saved preprocessing and feature selection are pipeline-contained; all predictions use the same 141-row holdout; labels and confusion counts reconcile in class order `[0,1]` (TN, FP, FN, TP). The earlier report displayed top test scores descriptively; this can create optimistic post-hoc rankings, so the new model choice uses training-only nested CV.\n',
      '\nThe new evaluation used a stratified row split of 564 training rows and 141 holdout rows, seed 20261110. Nested stratified CV used 3 outer × 3 inner folds within training only. Imputation, one-hot encoding, scaling, mutation variance filtering and selectors were fitted within their respective training folds. The final holdout was evaluated after selection; it did not choose features, algorithms, or thresholds. Because no patient IDs/barcodes are supplied, independence between rows cannot be verified.\n',
      '\n### Quality findings\n',
      f"\nThe modelling matrix had {q['missing_cells']:,} missing cells, {q['duplicate_rows']} duplicate rows, {q['constant_features']} constant features, and {q['near_constant_95pct_features']} features with a single value in at least 95% of rows. The screen found {q['pearson_pairs_abs_ge_0_95']:,} pairs with absolute Pearson correlation at least 0.95. It retained 62 groups of identical observed CNV feature vectors with 99 additional alias columns; this is not evidence of biological identity. Mutations were zero in {q['mutation_sparsity']['fraction_zero_all_mutation_cells']:.1%} of cells; all observed mutation values were binary. IQR flags for RNA/protein are descriptive screens, not reasons to remove observations.\n",
      '\n### Nested CV comparison and model choice\n',
      '\nThe CV-selected global candidate was **RNA Random Forest** by the prespecified average rank over ROC-AUC, PR-AUC, balanced accuracy, macro-F1, and deceased recall. The scope-level nested outer-fold estimates below evaluate the full inner-selection procedure; folds may select different algorithms. The CV table is therefore the safer estimate of the selection procedure, while candidate means in `model_selection_summary.csv` are inner-CV selection values and are not independent validation.\n',
      '\n| Scope | Outer CV ROC-AUC | PR-AUC | Balanced accuracy | Macro-F1 | Deceased recall |\n|---|---:|---:|---:|---:|---:|']
    for scope in ['cnv','rna','mutation','protein','integrated']:
        r=cv.loc[scope];lines.append(f"| {scope} | {f(r.roc_auc_mean)} ± {f(r.roc_auc_std)} | {f(r.pr_auc_mean)} ± {f(r.pr_auc_std)} | {f(r.balanced_accuracy_mean)} ± {f(r.balanced_accuracy_std)} | {f(r.macro_f1_mean)} ± {f(r.macro_f1_std)} | {f(r.deceased_recall_mean)} ± {f(r.deceased_recall_std)} |")
    lines += ['\n\nThe SD values are across three outer folds; they are descriptive fold variability, not confidence intervals (outer-fold training sets overlap). Candidate means below average each algorithm’s inner CV across the three outer training partitions and are used for selection, so they are not unbiased generalization estimates.\n',
      '\n| Scope | Candidate | Inner ROC-AUC | PR-AUC | Balanced accuracy | Macro-F1 | Deceased recall | Mean rank | Selected in scope |\n|---|---|---:|---:|---:|---:|---:|---:|:---:|']
    for _,r in cand.sort_values(['scope','mean_composite_rank']).iterrows():
        lines.append(f"| {r.scope} | {r.model} | {r.inner_roc_auc_mean:.3f} | {r.inner_pr_auc_mean:.3f} | {r.inner_balanced_accuracy_mean:.3f} | {r.inner_macro_f1_mean:.3f} | {r.inner_deceased_recall_mean:.3f} | {r.mean_composite_rank:.2f} | {'yes' if r.scope_selected else ''} |")
    lines += ['\n\nDummyClassifier predicted every sample alive. Its holdout accuracy was '
      +f"{dum.test_accuracy:.3f}, balanced accuracy {dum.test_balanced_accuracy:.3f}, ROC-AUC {dum.test_roc_auc:.3f}, PR-AUC {dum.test_pr_auc:.3f}, and deceased recall 0. This exposes why accuracy alone is misleading for the 611/94 imbalance. Exact baseline metrics and confusion counts are saved in `../final_optimization/tables/dummy_baseline_comparison.csv`.\n",
      '\n### Descriptive holdout performance\n',
      '\nThe training-CV-selected model for each scope was fitted on all 564 training rows and evaluated on the fixed 141-row holdout (122 alive, 19 deceased). These test scores are descriptive, not model-selection criteria.\n',
      '\n| Scope and selected model | ROC-AUC | PR-AUC | Accuracy | Balanced accuracy | Deceased precision | Deceased recall | Specificity | Macro-F1 | TN/FP/FN/TP |\n|---|---:|---:|---:|---:|---:|---:|---:|---:|---:|']
    for scope,r in test.iterrows(): lines.append(f"| {scope} / {r.model} | {r.roc_auc:.3f} | {r.pr_auc:.3f} | {r.accuracy:.3f} | {r.balanced_accuracy:.3f} | {r.precision_deceased:.3f} | {r.deceased_recall:.3f} | {r.specificity:.3f} | {r.macro_f1:.3f} | {int(r.tn)}/{int(r.fp)}/{int(r.fn)}/{int(r.tp)} |")
    rna=test.loc['rna'];integ=test.loc['integrated']
    lines += ['\n\nIntegrated Random Forest had higher descriptive holdout ROC-AUC and PR-AUC than RNA Random Forest, while deceased recall was equal (6/19). The nested CV integrated procedure also had higher ROC-AUC and PR-AUC point estimates than RNA, but slightly lower balanced accuracy, macro-F1 and deceased recall. The original phase comparison favored RNA over integrated Extra Trees, and prior holdout permutation tests were nonsignificant. Across these small samples and varying model selection, the evidence does not establish a reliable or meaningful integration gain.\n',
      '\n### Feature selection and explainability\n',
      '\nMutual information, ANOVA SelectKBest, L1, and tree-importance selection were each fitted in five stratified training folds for RNA and integrated data. Fold-level selections and frequencies are saved under `../final_optimization/feature_selection/`. A feature selected often is only stable under that selector and this cohort; it is not a validated biomarker. Training-fold permutation importance was computed only for each fold’s native-importance top 30 variables on that fold’s validation data (3 repeats per variable). These are predictive associations, not causal evidence. SHAP was unavailable.\n',
      '\nTop fold-validation permutation features (mean ROC-AUC decrease; fold support):\n']
    for scope in ['rna','integrated']:
        lines.append(f'\n**{scope.title()}**\n\n| Feature | Modality | Mean ROC-AUC decrease | Folds represented |\n|---|---|---:|---:|\n')
        part=perm[perm.scope==scope].sort_values('mean_auc_decrease',ascending=False).head(10)
        for _,r in part.iterrows(): lines.append(f"| {r.feature} | {r.modality} | {r.mean_auc_decrease:.4f} | {int(r.outer_folds_permuted)}/3 |")
    stable=fs[(fs.scope.isin(['rna','integrated']))&(fs.selection_method.isin(['mutual_information','l1','tree_importance']))&(fs.folds_selected>=4)].sort_values(['selection_method','scope','folds_selected','selection_frequency'],ascending=[True,True,False,False]).groupby(['scope','selection_method']).head(3)
    lines += ['\n\nExamples selected in at least four of five folds by non-ANOVA selectors (feature stability is method-specific):\n', '\n| Scope | Selector | Feature | Modality | Selected folds |\n|---|---|---|---|---:|']
    for _,r in stable.iterrows(): lines.append(f"| {r.scope} | {r.selection_method} | {r.feature} | {r.modality} | {int(r.folds_selected)}/5 |")
    lines += ['\n\nSelector comparison, five-fold training-only summary:\n',
      '\n| Scope | Selector | ROC-AUC mean ± SD | PR-AUC | Balanced accuracy | Deceased recall |\n|---|---|---:|---:|---:|---:|']
    selector=pd.read_csv(OPT/'feature_selection'/'selector_summary.csv')
    for _,r in selector.iterrows(): lines.append(f"| {r.scope} | {r.selection_method} | {r.roc_auc_mean:.3f} ± {r.roc_auc_std:.3f} | {r.pr_auc_mean:.3f} | {r.balanced_accuracy_mean:.3f} | {r.deceased_recall_mean:.3f} |")
    lines += ['\n\nHoldout Brier scores (no post-hoc calibrator was fitted):\n',
      '\n| Model | Brier score |\n|---|---:|']
    for name,r in cal.iterrows(): lines.append(f"| {name} | {r.brier_score:.4f} |")
    lines.append('\nThe Dummy prior Brier score is lower than the RNA and integrated model scores in this holdout; the model probability estimates therefore do not show better squared-probability accuracy here. Calibration curves are noisy with only 19 deceased cases.\n')
    lines += ['\n\n### Calibration, uncertainty, and algorithms not run\n',
      '\nNo probability calibrator was fitted. Holdout calibration curves and Brier scores are descriptive only; the holdout has 19 deceased cases. Earlier stratified bootstrap intervals for RNA Random Forest were broad: test ROC-AUC 0.579, 95% percentile CI [{:.3f}, {:.3f}], and deceased sensitivity 0.316, CI [{:.3f}, {:.3f}]. Earlier one-sided permutation tests were not statistically significant for RNA (p={:.3f}) or integrated (p={:.3f}). These results do not support a stable performance claim.\n'.format(
        *oldboot[(oldboot.scope=='rna')&(oldboot.model=='random_forest')&(oldboot.metric=='roc_auc')][['ci95_percentile_low','ci95_percentile_high']].iloc[0].tolist(),
        *oldboot[(oldboot.scope=='rna')&(oldboot.model=='random_forest')&(oldboot.metric=='sensitivity_deceased')][['ci95_percentile_low','ci95_percentile_high']].iloc[0].tolist(),
        oldperm[oldperm.scope=='rna'].iloc[0].one_sided_pvalue,oldperm[oldperm.scope=='integrated'].iloc[0].one_sided_pvalue),
      '\nRandom Forest, Extra Trees, class-weighted Logistic Regression, RBF SVM, and HistGradientBoosting were compared. HistGradientBoosting required running outside the restricted sandbox because Windows worker-pipe creation was denied within it; all requested nested fits completed in the approved execution and the final run has '+str(len(errors))+' candidate errors. XGBoost, LightGBM, imbalanced-learn (Balanced Random Forest/SMOTE), and SHAP were unavailable and were not installed. Candidate-level inner-CV results for all estimators are in `../final_optimization/tables/model_selection_summary.csv`; they are model-selection estimates, not independent performance estimates.\n',
      '\n### Model artifacts and limitations\n',
      '\nCV-selected final scope pipelines, configs, predictions, fold indices, package versions, and preprocessing settings are under `../final_optimization/models/`, `predictions/`, and `reproducibility/`. Figures (PNG and PDF) include ROC, precision-recall, confusion matrices, calibration, nested-CV model comparison, feature-quality distributions/correlations, fold stability, and permutation importance.\n',
      '\nThe dataset contains only 94 deceased cases and the holdout only 19. There are no patient identifiers or follow-up durations, so row-level independence and survival endpoints cannot be assessed. There is no external validation. Recorded vital status is not survival-time prediction. Models are exploratory and not clinically validated; they should not be used to diagnose breast cancer or predict individual outcomes.\n']
    text='\n'.join(lines)
    existing=REPORT.read_text(encoding='utf-8')
    existing=existing.replace(
      'This report consolidates existing Dataset C validation and modeling outputs. It does not fit new models.',
      'This report consolidates Dataset C validation and modeling outputs; the final optimization addendum below reports additional nested-CV model fits.')
    existing=existing.replace(
      'This is exploratory classification of recorded vital status, not survival-time prediction; there are no patient IDs or follow-up durations, and patient-level independence cannot be checked.',
      'The later nested-CV optimization selected RNA Random Forest by training-CV criteria (scope-level nested outer ROC-AUC 0.632 ± 0.073 for RNA and 0.671 ± 0.068 integrated). On the fixed holdout, RNA RF ROC-AUC was 0.579 and integrated RF 0.605; these exploratory test comparisons do not establish an integration benefit. This remains exploratory classification of recorded vital status, not survival-time prediction. There are no patient IDs or follow-up durations, and patient-level independence cannot be checked.')
    marker='## Final optimization and methodology verification (2026-11-10)'
    if marker in existing: existing=existing.split(marker)[0].rstrip()
    REPORT.write_text(existing+text,encoding='utf-8')
    (OPT/'reports'/'advanced_execution_report.md').write_text('# Final optimization addendum\n'+text,encoding='utf-8')
    opath=OPT/'reports'/'optimization_report.md'; old=opath.read_text(encoding='utf-8')
    marker='## Explicit baseline and final execution status'
    if marker in old: old=old.split(marker)[0].rstrip()
    baseline_note=(f"\n\n{marker}\n\nThe DummyClassifier prior baseline predicted all cases alive: test accuracy {dum.test_accuracy:.3f}, balanced accuracy {dum.test_balanced_accuracy:.3f}, ROC-AUC {dum.test_roc_auc:.3f}, PR-AUC {dum.test_pr_auc:.3f}, and deceased recall {dum.test_deceased_recall:.3f}. Its confusion counts are TN={int(dum.test_tn)}, FP={int(dum.test_fp)}, FN={int(dum.test_fn)}, TP={int(dum.test_tp)}. This is the accuracy-only trap under class imbalance.\n\nHistGradientBoosting completed nested CV in the approved execution after the sandbox prevented worker-pipe creation. All {len(cand)} scope/model candidates completed with {len(errors)} fitting errors. The candidate metrics are in `../tables/model_selection_summary.csv`; original holdout scores remain descriptive only.\n")
    opath.write_text(old+baseline_note,encoding='utf-8')
    print(f'Updated {REPORT}; addendum lines={len(lines)}; selected={sel["selected_key"]}; errors={len(errors)}')
if __name__=='__main__':main()
