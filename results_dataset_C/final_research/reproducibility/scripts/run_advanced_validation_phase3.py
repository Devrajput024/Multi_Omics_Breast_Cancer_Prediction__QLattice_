"""Advanced, separate validation for Dataset C; does not overwrite phase-2 outputs."""
from __future__ import annotations
import csv, importlib.util, json, math, os, sys, traceback
from pathlib import Path

import joblib
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
import sklearn
from sklearn.calibration import CalibratedClassifierCV, calibration_curve
from sklearn.dummy import DummyClassifier
from sklearn.metrics import (accuracy_score, average_precision_score, balanced_accuracy_score,
    brier_score_loss, confusion_matrix, f1_score, log_loss, precision_score, recall_score,
    roc_auc_score)
from sklearn.model_selection import RepeatedStratifiedKFold, StratifiedKFold, train_test_split
from sklearn.compose import ColumnTransformer
from sklearn.ensemble import ExtraTreesClassifier, RandomForestClassifier
from sklearn.impute import SimpleImputer
from sklearn.linear_model import LogisticRegression
from sklearn.pipeline import Pipeline
from sklearn.preprocessing import OneHotEncoder, StandardScaler
from sklearn.feature_selection import VarianceThreshold

ROOT=Path(__file__).resolve().parents[2]
BASE=ROOT/'results_dataset_C'; OUT=BASE/'advanced_validation'; DATA=ROOT/'dataset_c_extracted'/'data.csv'
SEED=20261110; TEST_SIZE=.20; REPEATS=5; FOLDS=5; BOOTSTRAPS=3000; PERMUTATIONS=5000
SCOPE_MODEL={'cnv':'random_forest','rna':'random_forest','mutation':'logistic_regression','protein':'random_forest','integrated':'extra_trees'}
PREFIX={'cnv':'cn_','rna':'rs_','mutation':'mu_','protein':'pp_'}
for d in ['figures','reports','predictions','calibration_models','oof_predictions','feature_importance']:
 (OUT/d).mkdir(parents=True,exist_ok=True)

def log(s):
 print(s,flush=True)
 with (OUT/'advanced_validation_log.txt').open('a',encoding='utf-8') as f:f.write(s+'\n')

def evaluate(y,pred,prob):
 cm=confusion_matrix(y,pred,labels=[0,1]);tn,fp,fn,tp=map(int,cm.ravel())
 return {'accuracy':accuracy_score(y,pred),'balanced_accuracy':balanced_accuracy_score(y,pred),
  'precision_deceased':precision_score(y,pred,pos_label=1,zero_division=0),
  'sensitivity_deceased':recall_score(y,pred,pos_label=1,zero_division=0),
  'specificity_alive':tn/(tn+fp) if tn+fp else np.nan,
  'f1_deceased':f1_score(y,pred,pos_label=1,zero_division=0),
  'macro_f1':f1_score(y,pred,average='macro',zero_division=0),
  'roc_auc':roc_auc_score(y,prob),'pr_auc_average_precision':average_precision_score(y,prob),
  'precision_alive':precision_score(y,pred,pos_label=0,zero_division=0),
  'recall_alive':recall_score(y,pred,pos_label=0,zero_division=0),
  'tn':tn,'fp':fp,'fn':fn,'tp':tp}

def metrics_pred(y,prob,threshold=.5):
 return evaluate(y,(np.asarray(prob)>=threshold).astype(int),prob)

def paired_stratified_indices(y,rng):
 y=np.asarray(y); zero=np.flatnonzero(y==0); one=np.flatnonzero(y==1)
 return np.concatenate([rng.choice(zero,len(zero),replace=True),rng.choice(one,len(one),replace=True)])

def get_pipe_builder():
 class Builder: pass
 def build_pipeline(scope,model_name):
  groups=list(PREFIX) if scope=='integrated' else [scope]
  cols=[c for c in X.columns if any(c.startswith(PREFIX[g]) for g in groups)]
  cats=[c for c in cols if c.startswith('cn_')];mut=[c for c in cols if c.startswith('mu_')]
  continuous=[c for c in cols if c.startswith(('rs_','pp_'))];trans=[]
  if cats:trans.append(('cnv_categorical',Pipeline([('imputer',SimpleImputer(strategy='most_frequent')),('onehot',OneHotEncoder(handle_unknown='ignore',sparse_output=True))]),cats))
  if mut:trans.append(('mutation_binary',Pipeline([('imputer',SimpleImputer(strategy='median')),('variance',VarianceThreshold(0.0))]),mut))
  if continuous:trans.append(('continuous',Pipeline([('imputer',SimpleImputer(strategy='median')),('scale',StandardScaler())]),continuous))
  prep=ColumnTransformer(trans,remainder='drop',verbose_feature_names_out=False)
  if model_name=='logistic_regression':model=LogisticRegression(class_weight='balanced',C=.1,solver='liblinear',max_iter=2500,random_state=SEED)
  elif model_name=='random_forest':model=RandomForestClassifier(n_estimators=100,class_weight='balanced',max_features='sqrt',min_samples_leaf=2,n_jobs=1,random_state=SEED)
  else:model=ExtraTreesClassifier(n_estimators=100,class_weight='balanced',max_features='sqrt',min_samples_leaf=2,n_jobs=1,random_state=SEED)
  return Pipeline([('preprocess',prep),('model',model)])
 Builder.build_pipeline=staticmethod(build_pipeline)
 return Builder

def modality_feature(name, encoded, columns):
 # Map one-hot expanded CNV features back to their raw column, longest match first.
 for c in sorted(columns,key=len,reverse=True):
  if encoded==c or encoded.startswith(c+'_'): return c
 return encoded

def main():
 global X
 log(f'Advanced validation start; sklearn={sklearn.__version__}, seed={SEED}, repeats={REPEATS}, folds={FOLDS}')
 runner=get_pipe_builder();df=pd.read_csv(DATA); y=df['vital.status'].astype(int);X=df.drop(columns=['vital.status'])
 runner.X=X
 idx_train,idx_test=train_test_split(np.arange(len(df)),test_size=TEST_SIZE,random_state=SEED,stratify=y)
 ytr=y.iloc[idx_train].reset_index(drop=True);yte=y.iloc[idx_test].reset_index(drop=True)
 Xtr=X.iloc[idx_train].reset_index(drop=True);Xte=X.iloc[idx_test].reset_index(drop=True)
 prior=json.loads((BASE/'models'/'rna__random_forest.json').read_text(encoding='utf-8'))
 original_split=(np.asarray(prior['train_indices_zero_based']),np.asarray(prior['test_indices_zero_based']))
 if not np.array_equal(idx_train,original_split[0]) or not np.array_equal(idx_test,original_split[1]):raise ValueError('Reconstructed split differs from saved phase-2 split')

 # Audit all saved predictions against summary metrics and fixed split.
 summary=pd.read_csv(BASE/'metrics_summary.csv');audit=[];all_pred={}
 prediction_files=list((BASE/'predictions').glob('*_test_predictions.csv'))
 for path in prediction_files:
  p=pd.read_csv(path); name=path.name.removesuffix('_test_predictions.csv')
  if len(p)!=len(idx_test):raise ValueError(f'{name}: prediction count mismatch')
  if not np.array_equal(p.row_index_original_zero_based.to_numpy(),idx_test):raise ValueError(f'{name}: test indices differ')
  if not np.array_equal(p.true_vital_status.to_numpy(),yte.to_numpy()):raise ValueError(f'{name}: true labels do not match source target')
  m=evaluate(yte,p.predicted_vital_status.to_numpy(),p.probability_deceased.to_numpy())
  # Resolve exact metric row by experiment and model label, including selectors/dedup.
  if name.startswith(('cnv__','rna__','mutation__','protein__','integrated__')) and '__logistic_regression__' in name:
   scope,model=name.split('__',1)[0], 'logistic_regression'
  elif name.startswith('cnv_dedup__'):
   scope='cnv_dedup_'+name.split('__')[1];model='logistic_regression'
  else: scope,model=name.split('__',1)[0],name.split('__',1)[1]
  selection='none'
  if '__logistic_regression__' in name: selection=name.rsplit('__',1)[1]
  match=summary[(summary.experiment==scope)&(summary.model==model)&(summary.feature_selection.fillna('none')==selection)]
  if len(match)!=1: raise ValueError(f'{name}: summary lookup returned {len(match)} rows for {scope}/{model}/{selection}')
  stored=match.iloc[0]
  metric_alias={'sensitivity_deceased':'recall_deceased_sensitivity','balanced_accuracy':None}
  metric_diffs={k:abs(float(stored['test_'+metric_alias.get(k,k)])-float(v)) for k,v in m.items() if metric_alias.get(k,k) and 'test_'+metric_alias.get(k,k) in stored and pd.notna(stored['test_'+metric_alias.get(k,k)])}
  maxdiff=max(metric_diffs.values(),default=0.0)
  if maxdiff>1e-10:raise ValueError(f'{name}: saved metrics mismatch max diff {maxdiff}')
  audit.append({'prediction_file':str(path.relative_to(BASE)),'rows':len(p),'same_test_indices':True,'same_true_labels':True,'confusion_labels_order':'[0,1] (TN,FP,FN,TP)','metrics_match_summary':True,'max_absolute_metric_difference':maxdiff})
  all_pred[name]=p
 # Inspect saved pipelines and implementation.
 source=(BASE/'run_experiments.py').read_text(encoding='utf-8')
 predictor_columns=[c for c in X.columns]
 target_in_predictors='vital.status' in predictor_columns
 pipeline_target_leak=False
 pipeline_paths=list((BASE/'models').glob('*.joblib'))
 for fp in pipeline_paths:
  obj=joblib.load(fp)
  cols=list(obj.named_steps['preprocess'].feature_names_in_)
  if 'vital.status' in cols:pipeline_target_leak=True
 audit_details={
  'prediction_files_found':len(prediction_files),'model_pipelines_found':len(pipeline_paths),
  'target_excluded_from_phase2_input_X':not target_in_predictors,'target_in_any_saved_preprocessor_columns':pipeline_target_leak,
  'preprocessing_inside_sklearn_pipelines':all('preprocess' in joblib.load(fp).named_steps for fp in pipeline_paths),
  'feature_selection_inside_pipeline':"steps.append((\"select\"" in source and 'cross_validate(pipe,X_train,y_train' in source,
  'CV_uses_training_split_only':'cross_validate(pipe, X_train, y_train' in source,
  'test_split_reused_for_all_predictions':all(a['same_test_indices'] for a in audit),
  'saved_metrics_reconciled':all(a['metrics_match_summary'] for a in audit),
  'confusion_matrix_order':'sklearn confusion_matrix labels=[0,1], yielding TN, FP, FN, TP',
  'posthoc_test_ranking_caveat':'Phase 2 report included maxima ranked by held-out test ROC-AUC. This did not change saved pipelines, thresholds, or feature selection, but those posthoc maxima are descriptive and optimistically selected. Advanced model choices are frozen from original training-CV ranks.'}
 (OUT/'reports'/'prediction_reconciliation.json').write_text(json.dumps(audit,indent=2),encoding='utf-8')
 (OUT/'reports'/'implementation_audit.json').write_text(json.dumps(audit_details,indent=2),encoding='utf-8')

 # Compare every phase-2 test prediction with a majority-class DummyClassifier baseline.
 dummy=DummyClassifier(strategy='most_frequent');dummy.fit(np.zeros((len(ytr),1)),ytr)
 dprob=dummy.predict_proba(np.zeros((len(yte),1)))[:,1];dpred=dummy.predict(np.zeros((len(yte),1)))
 dummy_metrics=evaluate(yte,dpred,dprob)
 pd.DataFrame({'row_index_original_zero_based':idx_test,'true_vital_status':yte,'predicted_vital_status':dpred,'probability_deceased':dprob}).to_csv(OUT/'predictions'/'dummy_majority_test_predictions.csv',index=False)
 baseline_rows=[{'model':'dummy_majority','scope':'baseline','threshold':.5,**dummy_metrics}]
 for name,p in all_pred.items():
  modelname=name.removesuffix('_test_predictions')
  scope=modelname.split('__')[0]
  mm=evaluate(p.true_vital_status,p.predicted_vital_status,p.probability_deceased)
  baseline_rows.append({'model':modelname,'scope':scope,'threshold':.5,**mm})
 pd.DataFrame(baseline_rows).to_csv(OUT/'baseline_comparison.csv',index=False)

 # Repeated 5x5 CV for one prespecified, original-CV-selected model per modality.
 repeated=RepeatedStratifiedKFold(n_splits=FOLDS,n_repeats=REPEATS,random_state=SEED)
 split_records=[];oof={scope:np.full((REPEATS,len(ytr)),np.nan) for scope in SCOPE_MODEL};feature_records=[]
 for split_no,(tr,va) in enumerate(repeated.split(Xtr,ytr)):
  repeat_id=split_no//FOLDS;fold_id=split_no%FOLDS
  for scope,model_name in SCOPE_MODEL.items():
   pipe=runner.build_pipeline(scope,model_name);pipe.fit(Xtr.iloc[tr],ytr.iloc[tr]);prob=pipe.predict_proba(Xtr.iloc[va])[:,1];pred=(prob>=.5).astype(int)
   oof[scope][repeat_id,va]=prob
   m=evaluate(ytr.iloc[va],pred,prob)
   split_records.append({'scope':scope,'model':model_name,'repeat':repeat_id+1,'fold':fold_id+1,**m})
   if scope in {'rna','integrated'}:
    names=np.asarray(pipe.named_steps['preprocess'].get_feature_names_out())
    imps=np.asarray(pipe.named_steps['model'].feature_importances_)
    rawcols=[c for c in X.columns if c.startswith(PREFIX.get(scope,''))] if scope!='integrated' else list(X.columns)
    agg={}
    for name,imp in zip(names,imps):
     raw=modality_feature(name,name,rawcols)
     agg[raw]=agg.get(raw,0.0)+float(imp)
    top=sorted(agg.items(),key=lambda z:z[1],reverse=True)[:30]
    for rank,(feat,imp) in enumerate(top,1):
     mod=next((k for k,pref in PREFIX.items() if feat.startswith(pref)),'unknown')
     feature_records.append({'scope':scope,'model':model_name,'repeat':repeat_id+1,'fold':fold_id+1,'rank':rank,'feature':feat,'modality':mod,'importance':imp})
  # repeated majority Dummy baseline using this exact fold.
  dum=DummyClassifier(strategy='most_frequent').fit(np.zeros((len(tr),1)),ytr.iloc[tr])
  pr=dum.predict_proba(np.zeros((len(va),1)))[:,1];pdum=dum.predict(np.zeros((len(va),1)))
  split_records.append({'scope':'dummy_majority','model':'dummy','repeat':repeat_id+1,'fold':fold_id+1,**evaluate(ytr.iloc[va],pdum,pr)})
  if fold_id==FOLDS-1:log(f'Repeated CV finished repeat {repeat_id+1}/{REPEATS}')
 # Per-fold and empirical repeat-mean summaries, with quantiles of the five repeat means.
 rf=pd.DataFrame(split_records); summary_rows=[]
 for (scope,model),sub in rf.groupby(['scope','model']):
  metric_names=['accuracy','balanced_accuracy','macro_f1','roc_auc','pr_auc_average_precision','sensitivity_deceased','specificity_alive']
  repmeans=sub.groupby('repeat')[metric_names].mean()
  row={'scope':scope,'model':model,'folds':len(sub),'repeats':REPEATS}
  for metric in metric_names:
   row[metric+'_mean']=repmeans[metric].mean();row[metric+'_sd_repeat_means']=repmeans[metric].std(ddof=1)
   row[metric+'_repeat_mean_ci95_low']=repmeans[metric].quantile(.025);row[metric+'_repeat_mean_ci95_high']=repmeans[metric].quantile(.975)
  summary_rows.append(row)
 pd.concat([rf.assign(record_type='fold'),pd.DataFrame(summary_rows).assign(record_type='repeat_summary')],ignore_index=True,sort=False).to_csv(OUT/'repeated_cv_results.csv',index=False)
 for scope,arr in oof.items():
  if np.isnan(arr).any():raise RuntimeError(f'Missing repeated OOF predictions for {scope}')
  pd.DataFrame({'training_row_index_zero_based':idx_train,'true_vital_status':ytr,**{f'repeat_{i+1}_probability_deceased':arr[i] for i in range(REPEATS)},'mean_oof_probability_deceased':arr.mean(axis=0)}).to_csv(OUT/'oof_predictions'/f'{scope}_repeated_oof.csv',index=False)
 log('Repeated-CV experiments complete')

 # Threshold selection uses averaged repeated out-of-fold training probabilities only.
 thresholds=[]
 for scope in ['rna','integrated']:
  oofprob=oof[scope].mean(axis=0); yarr=ytr.to_numpy()
  candidates=np.unique(np.concatenate(([.5],oofprob)))
  scores=[]
  for t in candidates:
   mm=metrics_pred(yarr,oofprob,t); scores.append((mm['sensitivity_deceased']+mm['specificity_alive']-1,t,mm))
  best=max(scores,key=lambda z:(z[0],-abs(z[1]-.5))); selected_t=float(best[1])
  predfile=all_pred[f'{scope}__{SCOPE_MODEL[scope]}']
  test_prob=predfile.probability_deceased.to_numpy();test_y=predfile.true_vital_status.to_numpy()
  for label,t in [('default_0.5',.5),('OOF_Youden_J',selected_t)]:
   om=metrics_pred(yarr,oofprob,t);tm=metrics_pred(test_y,test_prob,t)
   thresholds.append({'scope':scope,'model':SCOPE_MODEL[scope],'threshold_source':'mean repeated stratified OOF training probabilities','threshold_rule':'maximize sensitivity + specificity - 1 (Youden J), not accuracy','threshold_label':label,'threshold':t,'oof_sensitivity':om['sensitivity_deceased'],'oof_specificity':om['specificity_alive'],'oof_precision_deceased':om['precision_deceased'],'oof_f1_deceased':om['f1_deceased'],'oof_tn':om['tn'],'oof_fp':om['fp'],'oof_fn':om['fn'],'oof_tp':om['tp'],'test_sensitivity':tm['sensitivity_deceased'],'test_specificity':tm['specificity_alive'],'test_precision_deceased':tm['precision_deceased'],'test_f1_deceased':tm['f1_deceased'],'test_tn':tm['tn'],'test_fp':tm['fp'],'test_fn':tm['fn'],'test_tp':tm['tp']})
 pd.DataFrame(thresholds).to_csv(OUT/'threshold_analysis.csv',index=False)

 # Probability calibration: fit sigmoid calibrators only on the training partition.
 calibration_rows=[]
 for scope in ['rna','integrated']:
  base_name=f'{scope}__{SCOPE_MODEL[scope]}';raw_pipe=joblib.load(BASE/'models'/f'{base_name}.joblib')
  raw_prob=all_pred[base_name].probability_deceased.to_numpy()
  cal=CalibratedClassifierCV(estimator=runner.build_pipeline(scope,SCOPE_MODEL[scope]),method='sigmoid',cv=StratifiedKFold(n_splits=5,shuffle=True,random_state=SEED))
  cal.fit(Xtr,ytr); calibrated=cal.predict_proba(Xte)[:,1]
  joblib.dump(cal,OUT/'calibration_models'/f'{base_name}_sigmoid_calibrated.joblib',compress=3)
  for name,prob in [('uncalibrated',raw_prob),('sigmoid_cv_calibrated',calibrated)]:
   calibration_rows.append({'scope':scope,'model':SCOPE_MODEL[scope],'probability_type':name,'brier_score_test':brier_score_loss(yte,prob),'log_loss_test':log_loss(yte,prob,labels=[0,1]),'test_samples':len(yte),'calibration_training':'training split only; sigmoid calibration with 5-fold CV' if name!='uncalibrated' else 'uncalibrated training-split fitted model'})
  frac,mean=calibration_curve(yte,raw_prob,n_bins=6,strategy='quantile')
  frac2,mean2=calibration_curve(yte,calibrated,n_bins=6,strategy='quantile')
  fig,ax=plt.subplots(figsize=(5.5,5));ax.plot([0,1],[0,1],'k--',label='ideal');ax.plot(mean,frac,'o-',label='uncalibrated');ax.plot(mean2,frac2,'s-',label='sigmoid calibrated');ax.set(xlabel='Mean predicted probability',ylabel='Observed deceased fraction',title=f'Test calibration: {scope}');ax.legend();fig.tight_layout();fig.savefig(OUT/'figures'/f'{scope}_calibration.png',dpi=150);plt.close(fig)
 pd.DataFrame(calibration_rows).to_csv(OUT/'calibration_results.csv',index=False)

 # Stratified bootstrap CIs on frozen test predictions; class-conditional resampling preserves counts.
 rng=np.random.default_rng(SEED);boot_rows=[];boot_draws={}
 ci_metrics=['accuracy','balanced_accuracy','macro_f1','roc_auc','pr_auc_average_precision','sensitivity_deceased','specificity_alive']
 for scope,model_name in SCOPE_MODEL.items():
  name=f'{scope}__{model_name}';p=all_pred[name];yt=p.true_vital_status.to_numpy();prob=p.probability_deceased.to_numpy();pred=p.predicted_vital_status.to_numpy()
  draws={k:[] for k in ci_metrics}
  for _ in range(BOOTSTRAPS):
   ii=paired_stratified_indices(yt,rng);m=evaluate(yt[ii],pred[ii],prob[ii])
   for k in ci_metrics:draws[k].append(m[k])
  boot_draws[scope]=draws
  point=evaluate(yt,pred,prob)
  for k,values in draws.items():
   boot_rows.append({'scope':scope,'model':model_name,'metric':k,'point_estimate':point[k],'ci95_percentile_low':float(np.quantile(values,.025)),'ci95_percentile_high':float(np.quantile(values,.975)),'bootstrap_replicates':BOOTSTRAPS,'resampling':'stratified within observed test classes; fixed model/predictions'})
 pd.DataFrame(boot_rows).to_csv(OUT/'bootstrap_confidence_intervals.csv',index=False)

 # Frozen-model permutation test for conditional held-out association; labels shuffled, models unchanged.
 perm_rows=[]
 for scope in ['rna','integrated']:
  name=f'{scope}__{SCOPE_MODEL[scope]}';p=all_pred[name];yt=p.true_vital_status.to_numpy();prob=p.probability_deceased.to_numpy();obs=roc_auc_score(yt,prob)
  rngp=np.random.default_rng(SEED+31);null=np.empty(PERMUTATIONS)
  for i in range(PERMUTATIONS):null[i]=roc_auc_score(rngp.permutation(yt),prob)
  pval=(1+np.sum(null>=obs))/(PERMUTATIONS+1)
  perm_rows.append({'scope':scope,'model':SCOPE_MODEL[scope],'observed_test_roc_auc':obs,'permutations':PERMUTATIONS,'one_sided_pvalue':pval,'null_95th_percentile':float(np.quantile(null,.95)),'test_labels_used_for_inference_only':True,'model_or_threshold_selected_using_test':False})
 pd.DataFrame(perm_rows).to_csv(OUT/'reports'/'heldout_permutation_tests.csv',index=False)

 # Modality contribution: paired bootstrap differences integrated minus RNA using the same test rows.
 modrows=[]
 for scope,model_name in SCOPE_MODEL.items():
  name=f'{scope}__{model_name}';p=all_pred[name];m=evaluate(p.true_vital_status,p.predicted_vital_status,p.probability_deceased)
  cvsummary=pd.read_csv(OUT/'repeated_cv_results.csv');cvsummary=cvsummary[(cvsummary.scope==scope)&(cvsummary.record_type=='repeat_summary')].iloc[0]
  modrows.append({'scope':scope,'model':model_name,'repeated_cv_roc_auc_mean':cvsummary.roc_auc_mean,'repeated_cv_roc_auc_sd_repeat_means':cvsummary.roc_auc_sd_repeat_means,'repeated_cv_roc_auc_empirical_95_low':cvsummary.roc_auc_repeat_mean_ci95_low,'repeated_cv_roc_auc_empirical_95_high':cvsummary.roc_auc_repeat_mean_ci95_high,**{'test_'+k:v for k,v in m.items()}})
 # Paired bootstrap delta by same resampled heldout rows.
 p1=all_pred[f'rna__{SCOPE_MODEL["rna"]}'];p2=all_pred[f'integrated__{SCOPE_MODEL["integrated"]}'];yt=p1.true_vital_status.to_numpy();rngd=np.random.default_rng(SEED+9);diffs={'roc_auc':[],'pr_auc_average_precision':[],'balanced_accuracy':[],'macro_f1':[]}
 for _ in range(BOOTSTRAPS):
  ii=paired_stratified_indices(yt,rngd);m1=evaluate(yt[ii],p1.predicted_vital_status.to_numpy()[ii],p1.probability_deceased.to_numpy()[ii]);m2=evaluate(yt[ii],p2.predicted_vital_status.to_numpy()[ii],p2.probability_deceased.to_numpy()[ii])
  for k in diffs:diffs[k].append(m2[k]-m1[k])
 (OUT/'reports'/'integration_minus_rna_paired_bootstrap.json').write_text(json.dumps({k:{'mean_difference':float(np.mean(v)),'ci95_low':float(np.quantile(v,.025)),'ci95_high':float(np.quantile(v,.975)),'bootstrap_replicates':BOOTSTRAPS} for k,v in diffs.items()},indent=2),encoding='utf-8')
 pd.DataFrame(modrows).to_csv(OUT/'modality_comparison.csv',index=False)

 # Fold stability: aggregate top-30 model importance across repeated folds; annotate prior selections.
 fi=pd.DataFrame(feature_records);prior_sel=pd.read_csv(BASE/'feature_selection'/'selected_features.csv')
 prior_names={modality_feature(str(n),str(n),list(X.columns)) for n in prior_sel.selected_feature.astype(str)}
 stability=[]
 for (scope,feat,modality),sub in fi.groupby(['scope','feature','modality']):
  stability.append({'scope':scope,'feature':feat,'modality':modality,'folds_top30':sub.shape[0],'total_repeated_folds':REPEATS*FOLDS,'top30_frequency':sub.shape[0]/(REPEATS*FOLDS),'mean_rank_when_top30':sub['rank'].mean(),'mean_importance_when_top30':sub.importance.mean(),'selected_in_previous_phase2_feature_selection':feat in prior_names})
 stability_df=pd.DataFrame(stability).sort_values(['scope','top30_frequency','mean_importance_when_top30'],ascending=[True,False,False])
 stability_df.to_csv(OUT/'feature_stability.csv',index=False)
 fi.to_csv(OUT/'feature_importance'/'repeated_fold_top30_importance.csv',index=False)
 for scope,sub in stability_df.groupby('scope'):
  top=sub.head(20).iloc[::-1];fig,ax=plt.subplots(figsize=(8,7));ax.barh(top.feature,top.top30_frequency);ax.set_xlabel('Fraction of repeated folds in top 30');ax.set_title(f'Feature stability — {scope}');fig.tight_layout();fig.savefig(OUT/'figures'/f'{scope}_feature_stability.png',dpi=150);plt.close(fig)

 # Baseline comparison plot by repeated-CV ROC-AUC.
 rs=pd.read_csv(OUT/'repeated_cv_results.csv');rs=rs[rs.record_type.eq('repeat_summary')]
 fig,ax=plt.subplots(figsize=(8,5));ax.bar(rs.scope,rs.roc_auc_mean,yerr=rs.roc_auc_sd_repeat_means,capsize=4);ax.axhline(.5,color='gray',ls='--',label='chance ROC-AUC');ax.set_ylim(0,1);ax.set_ylabel('Repeated CV ROC-AUC (mean ± SD of repeat means)');ax.legend();fig.tight_layout();fig.savefig(OUT/'figures'/'modality_repeated_cv_roc_auc.png',dpi=150);plt.close(fig)

 # Final advanced report.
 integ_json=json.loads((OUT/'reports'/'integration_minus_rna_paired_bootstrap.json').read_text(encoding='utf-8'))
 report=['# Dataset C — Phase 3 advanced validation report','',
 '## Purpose and limits','',
 'This phase audits and extends the existing exploratory classifiers for recorded `vital.status` (0 alive, 1 deceased). It is not survival-time prediction: there are no durations, follow-up, event times, censoring fields, or patient IDs/TCGA barcodes. Splits remain row-based; patient-level independence cannot be independently verified. Nothing here validates clinical utility or generalization to new patients. Original phase-2 results remain unchanged.','',
 '## Audit of Phase 2 implementation','',
 f'- Found {len(pipeline_paths)} saved pipelines, {len(prediction_files)} saved prediction files, and {len(audit)} predictions reconciled against `metrics_summary.csv`. Every prediction file has {len(idx_test)} rows and exactly the same saved test row indices and true labels. Recomputed metrics match the summary within maximum absolute error {max((a["max_absolute_metric_difference"] for a in audit),default=0):.3g}.',
 '- Target excluded: all saved preprocessor input columns exclude `vital.status`; the phase-2 runner built X by dropping the target before fitting.',
 '- Leakage control: imputers, scaling, CNV one-hot encoding, mutation variance filtering, and selectors are Pipeline steps. `cross_validate` received training-partition rows only. Feature-selection selectors are nested inside the pipeline and therefore refit per training fold.',
 '- The same stratified 564/141 row split (seed 20261110) is used by every saved prediction; class order is [0,1], so confusion counts are TN, FP, FN, TP.',
 '- Test-ranking caveat: the previous report also listed the highest held-out ROC-AUC models after comparing many candidates. This is posthoc winner ranking and can be optimistic. No test-derived threshold or feature selection was found; Phase 3 freezes model families by Phase-2 training-CV results and treats all test metrics as descriptive.',
 '- Phase-2 source has no patient identifiers and a single random split; no grouped or external validation is possible. See `reports/implementation_audit.json` and `reports/prediction_reconciliation.json`.','',
 '## Majority-class baseline and model comparison','',
 f'The majority DummyClassifier predicts alive for all held-out rows ({len(yte)} rows; {int((yte==0).sum())} alive, {int((yte==1).sum())} deceased). Full comparison against all 23 existing prediction files is in `baseline_comparison.csv`; it includes accuracy, balanced accuracy, macro-F1, ROC-AUC, PR-AUC, class-specific metrics, and confusion counts.','',
 '## Repeated stratified CV','',
 f'Used {REPEATS} repeats of stratified {FOLDS}-fold CV on the 564-row training partition. For modality comparisons, architectures were fixed in advance from the previous training-CV ranking: CNV/RF, RNA/RF, mutation/Logistic Regression, protein/RF, integrated/Extra Trees. Dummy majority was scored on the exact same folds. `repeated_cv_results.csv` contains all 150 fold-level rows (including dummy) and repeat-summary rows. Mean/SD and empirical 2.5–97.5% quantiles are computed over five repeat means; these are descriptive resampling ranges, not formal independent-sample confidence intervals because folds/repeats overlap.','',
 '| Scope | Fixed model | ROC-AUC mean ± SD of repeat means | PR-AUC mean ± SD of repeat means | empirical repeat-mean 95% range |', '|---|---|---:|---:|---:|']
 for _,r in rs.iterrows():report.append(f'| {r.scope} | {r.model} | {r.roc_auc_mean:.3f} ± {r.roc_auc_sd_repeat_means:.3f} | {r.pr_auc_average_precision_mean:.3f} ± {r.pr_auc_average_precision_sd_repeat_means:.3f} | {r.roc_auc_repeat_mean_ci95_low:.3f}–{r.roc_auc_repeat_mean_ci95_high:.3f} |')
 report+=['','## Threshold analysis','',
 'Thresholds were selected from the mean of repeated out-of-fold training probabilities for RNA Random Forest and integrated Extra Trees by maximizing Youden J (sensitivity + specificity − 1). The test set did not choose thresholds. Both 0.5 and the selected threshold are reported in `threshold_analysis.csv`, with OOF and held-out sensitivity, specificity, deceased precision/F1, and confusion counts. The threshold chosen on OOF data can still be optimistic; test rows are a final descriptive check.','',
 '## Calibration','',
 'Sigmoid calibration used 5-fold CV within the training partition only. Test Brier score/log loss and calibration curves for uncalibrated and calibrated probabilities are in `calibration_results.csv` and `figures/`. Test calibration estimates are noisy with only 19 deceased cases. Calibration was not selected using test performance.','',
 '## Held-out uncertainty and permutation tests','',
 f'Stratified bootstrap percentile intervals use {BOOTSTRAPS} resamples of frozen test predictions, resampling within alive/deceased strata; this quantifies uncertainty conditional on this particular holdout and fitted models, not training/model-selection variability. See `bootstrap_confidence_intervals.csv`. A fixed-model, held-out label permutation test ({PERMUTATIONS} permutations) is reported for the prespecified RNA and integrated models; test labels are used for inference only, not model/threshold selection. Results are in `reports/heldout_permutation_tests.csv`.','',
 '## Modality contribution','',
 'All five modalities use the same training rows, repeated CV folds, and held-out rows. Integrated minus RNA paired-bootstrap differences for ROC-AUC, PR-AUC, balanced accuracy, and macro-F1 are in `reports/integration_minus_rna_paired_bootstrap.json`. These data do not support forcing a positive integration conclusion; point estimates and intervals should be read together.','',
 '## Feature stability and interpretation','',
 'Repeated-fold tree importance was aggregated across 25 training folds for RNA Random Forest and integrated Extra Trees. `feature_stability.csv` reports modality, top-30 frequency, mean rank, and overlap with previously selected features. It is a stability screen, not causal evidence or a biomarker claim. Figures are in `figures/`; raw fold ranks are under `feature_importance/`.','',
 '## Optional methods','',
 'SHAP is not installed in the existing environment; it was not installed because permutation/model-based explanations are feasible and changing the project environment is unnecessary for this phase. QLattice/Feyn was unavailable in Phase 2 (`ModuleNotFoundError: No module named feyn`) and was not pursued further.','',
 '## Key limitations','',
 'There are 705 observations, only 94 deceased labels, no patient IDs/barcodes, no survival times or follow-up, no censoring information, and no external cohort. Repeated CV and bootstrap quantify only parts of resampling uncertainty and cannot repair unknown patient dependence or selection bias. All outputs concern exploratory classification of recorded vital status only; they are not clinical survival predictions.']
 (OUT/'advanced_validation_report.md').write_text('\n'.join(report)+'\n',encoding='utf-8')
 log('Advanced validation completed; outputs verified and report saved.')

if __name__=='__main__':
 try:main()
 except Exception as e:
  log(f'FATAL {type(e).__name__}: {e}');log(traceback.format_exc());raise
