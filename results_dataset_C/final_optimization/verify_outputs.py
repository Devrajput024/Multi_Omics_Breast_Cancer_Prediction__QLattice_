"""Validate saved Dataset C optimization artifacts against the source CSV."""
from pathlib import Path
import json
import numpy as np
import pandas as pd
import joblib
from sklearn.metrics import (accuracy_score, balanced_accuracy_score, precision_score,
                             recall_score, f1_score, roc_auc_score,
                             average_precision_score, confusion_matrix)

ROOT=Path(__file__).resolve().parents[2]; OPT=ROOT/'results_dataset_C'/'final_optimization'
def main():
    errors=[]; data=pd.read_csv(ROOT/'dataset_c_extracted'/'data.csv');
    y=data['vital.status'].astype(int);x=data.drop(columns='vital.status')
    prefixes={'cnv':'cn_','rna':'rs_','mutation':'mu_','protein':'pp_'}
    pcounts={m:sum(c.startswith(p) for c in x.columns) for m,p in prefixes.items()}
    unexpected=[c for c in x.columns if not c.startswith(tuple(prefixes.values()))]
    split=pd.read_csv(OPT/'reproducibility'/'train_test_indices.csv'); testidx=split.loc[split.split=='test','original_row_index_zero_based'].to_numpy(int)
    pred=pd.read_csv(OPT/'predictions'/'selected_scope_test_predictions.csv')
    testmetrics=pd.read_csv(OPT/'tables'/'selected_scope_models_test_metrics.csv')
    predchecks={}
    for name,g in pred.groupby('model'):
        if len(g)!=len(testidx) or set(g.original_row_index_zero_based)!=set(testidx):errors.append(f'{name}: holdout indices mismatch')
        tru=y.iloc[g.original_row_index_zero_based.to_numpy(int)].to_numpy()
        if not np.array_equal(tru,g.true_vital_status.to_numpy(int)):errors.append(f'{name}: saved truth mismatch')
        cm=confusion_matrix(tru,g.predicted_vital_status,labels=[0,1]);tn,fp,fn,tp=map(int,cm.ravel())
        modelpath=OPT/'models'/f'{name}.joblib'
        try:
            fitted=joblib.load(modelpath)
            ordered_rows=g.original_row_index_zero_based.to_numpy(int)
            rebuilt_prob=fitted.predict_proba(x.iloc[ordered_rows])[:,1]
            rebuilt_pred=fitted.predict(x.iloc[ordered_rows])
            if not np.allclose(rebuilt_prob,g.probability_deceased.to_numpy(float),rtol=0,atol=1e-12):errors.append(f'{name}: saved pipeline probabilities mismatch predictions CSV')
            if not np.array_equal(rebuilt_pred,g.predicted_vital_status.to_numpy(int)):errors.append(f'{name}: saved pipeline labels mismatch predictions CSV')
        except Exception as exc:
            errors.append(f'{name}: saved pipeline load/predict failed: {type(exc).__name__}: {exc}')
        row=testmetrics[testmetrics.scope==name.split('__')[0]]
        if len(row)!=1: errors.append(f'{name}: missing/duplicate metric row')
        else:
            r=row.iloc[0]
            checks={'accuracy':accuracy_score(tru,g.predicted_vital_status),'balanced_accuracy':balanced_accuracy_score(tru,g.predicted_vital_status),
                    'precision_deceased':precision_score(tru,g.predicted_vital_status,pos_label=1,zero_division=0),
                    'deceased_recall':recall_score(tru,g.predicted_vital_status,pos_label=1,zero_division=0),
                    'specificity':tn/(tn+fp),'f1_deceased':f1_score(tru,g.predicted_vital_status,pos_label=1,zero_division=0),
                    'macro_f1':f1_score(tru,g.predicted_vital_status,average='macro',zero_division=0),
                    'roc_auc':roc_auc_score(tru,g.probability_deceased),'pr_auc':average_precision_score(tru,g.probability_deceased)}
            for k,v in checks.items():
                if not np.isclose(v,float(r[k]),rtol=0,atol=1e-12):errors.append(f'{name}: {k} does not match predictions')
            if (tn,fp,fn,tp)!=(int(r.tn),int(r.fp),int(r.fn),int(r.tp)):errors.append(f'{name}: confusion counts mismatch')
        predchecks[name]={'rows':len(g),'unique_test_indices':int(g.original_row_index_zero_based.nunique()),'tn_fp_fn_tp':[tn,fp,fn,tp]}
    cand=pd.read_csv(OPT/'tables'/'model_selection_summary.csv');fiterr=pd.read_csv(OPT/'tables'/'model_fit_errors.csv')
    if cand.shape[0]!=25:errors.append(f'expected 25 model candidates; found {cand.shape[0]}')
    hgb=int((cand.model=='hist_gradient_boosting').sum())
    if hgb!=5:errors.append(f'expected 5 HGB scope candidates; found {hgb}')
    if len(fiterr):errors.append(f'{len(fiterr)} model fit errors remain')
    required=['methodology_audit.md','feature_quality_report.md','feature_quality_summary.csv','correlation_analysis.csv',
      'tables/nested_cv_summary.csv','tables/model_selection_summary.csv','tables/dummy_baseline_comparison.csv',
      'tables/model_comparison_with_dummy.csv','tables/selected_scope_models_test_metrics.csv',
      'feature_selection/feature_selection_frequency.csv','feature_selection/outer_fold_permutation_importance.csv',
      'reports/optimization_report.md','reports/advanced_execution_report.md',
      'reproducibility/package_versions.json','reproducibility/train_test_indices.csv']
    missing=[f for f in required if not (OPT/f).is_file()]
    if missing:errors.append(f'missing required artifacts: {missing}')
    model_files=sorted(p.name for p in (OPT/'models').glob('*.joblib'))
    if len(model_files)!=5:errors.append(f'expected five final selected pipelines, found {len(model_files)}')
    figures=list((OPT/'figures').glob('*.png'));pdfs=list((OPT/'figures').glob('*.pdf'))
    report=ROOT/'results_dataset_C'/'final_research'/'dataset_C_research_report.md'
    if not report.is_file():errors.append('updated final research report missing')
    result={'status':'PASS' if not errors else 'FAIL','errors':errors,'source_shape':list(data.shape),
      'target_counts':{str(k):int(v) for k,v in y.value_counts().sort_index().items()},'modality_feature_counts':pcounts,
      'unexpected_predictors':unexpected,'train_rows':int((split.split=='train').sum()),'test_rows':len(testidx),
      'saved_test_models':predchecks,'candidate_models':int(len(cand)),'hist_gradient_boosting_candidates':hgb,
      'model_fit_errors':len(fiterr),'final_pipeline_paths':[f'models/{p}' for p in model_files],
      'figure_png_count':len(figures),'figure_pdf_count':len(pdfs),'final_research_report_exists':report.is_file()}
    (OPT/'reports'/'output_verification.json').write_text(json.dumps(result,indent=2),encoding='utf-8')
    print(json.dumps(result,indent=2))
    if errors:raise SystemExit(1)
if __name__=='__main__':main()
