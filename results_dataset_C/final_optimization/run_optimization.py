"""Dataset C final optimization; all selection/CV uses training rows only."""
from __future__ import annotations
from pathlib import Path
import json,sys,platform,traceback,importlib.util,importlib.metadata,hashlib
import numpy as np
import pandas as pd
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
import sklearn
from sklearn.compose import ColumnTransformer
from sklearn.dummy import DummyClassifier
from sklearn.ensemble import ExtraTreesClassifier,RandomForestClassifier,HistGradientBoostingClassifier
from sklearn.feature_selection import SelectFromModel,SelectKBest,mutual_info_classif,f_classif,VarianceThreshold
from sklearn.impute import SimpleImputer
from sklearn.linear_model import LogisticRegression
from sklearn.metrics import (accuracy_score,average_precision_score,balanced_accuracy_score,confusion_matrix,
 f1_score,make_scorer,precision_score,recall_score,roc_auc_score,roc_curve,precision_recall_curve,brier_score_loss)
from sklearn.model_selection import StratifiedKFold,cross_validate,train_test_split
from sklearn.pipeline import Pipeline
from sklearn.preprocessing import OneHotEncoder,StandardScaler
from sklearn.svm import SVC
import joblib

ROOT=Path(__file__).resolve().parents[2];BASE=ROOT/'results_dataset_C';ADV=BASE/'advanced_validation';OUT=BASE/'final_optimization'
DATA=ROOT/'dataset_c_extracted'/'data.csv';ENR=ROOT/'dataset_c_extracted'/'brca_data_w_subtypes.csv'
SEED=20261110;TEST_SIZE=.20;OUTER_FOLDS=3;INNER_FOLDS=3;SELECTOR_FOLDS=5;TREES=100;N_JOBS=1
PREFIX={'cnv':'cn_','rna':'rs_','mutation':'mu_','protein':'pp_'}
for d in ['reports','tables','models','predictions','figures','feature_selection','logs','reproducibility']:(OUT/d).mkdir(parents=True,exist_ok=True)
LOG=OUT/'logs'/'optimization.log'

def log(s):
 print(s,flush=True)
 with LOG.open('a',encoding='utf-8') as f:f.write(s+'\n')

def metric_dict(y,pred,prob):
 cm=confusion_matrix(y,pred,labels=[0,1]);tn,fp,fn,tp=map(int,cm.ravel())
 return {'accuracy':accuracy_score(y,pred),'balanced_accuracy':balanced_accuracy_score(y,pred),
  'precision_deceased':precision_score(y,pred,pos_label=1,zero_division=0),
  'deceased_recall':recall_score(y,pred,pos_label=1,zero_division=0),
  'specificity':tn/(tn+fp) if tn+fp else np.nan,'f1_deceased':f1_score(y,pred,pos_label=1,zero_division=0),
  'macro_f1':f1_score(y,pred,average='macro',zero_division=0),'roc_auc':roc_auc_score(y,prob),
  'pr_auc':average_precision_score(y,prob),'tn':tn,'fp':fp,'fn':fn,'tp':tp}

def specificity(y,pred):
 c=confusion_matrix(y,pred,labels=[0,1]);return c[0,0]/max(1,c[0].sum())
SCORING={'accuracy':'accuracy','balanced_accuracy':'balanced_accuracy','macro_f1':'f1_macro',
 'roc_auc':'roc_auc','pr_auc':'average_precision','deceased_recall':make_scorer(recall_score,pos_label=1,zero_division=0),
 'specificity':make_scorer(specificity)}
PRIORITY=['roc_auc','pr_auc','balanced_accuracy','macro_f1','deceased_recall']

def make_pipe(scope,model,columns,seed=SEED,selector=None):
 groups=list(PREFIX) if scope=='integrated' else [scope]
 cols=[c for c in columns if any(c.startswith(PREFIX[g]) for g in groups)]
 cn=[c for c in cols if c.startswith('cn_')];mu=[c for c in cols if c.startswith('mu_')]
 cont=[c for c in cols if c.startswith(('rs_','pp_'))];trs=[]
 if cn:trs.append(('cnv',Pipeline([('impute',SimpleImputer(strategy='most_frequent')),('onehot',OneHotEncoder(handle_unknown='ignore',sparse_output=False))]),cn))
 if mu:trs.append(('mutation',Pipeline([('impute',SimpleImputer(strategy='median')),('variance',VarianceThreshold(0.0))]),mu))
 if cont:trs.append(('continuous',Pipeline([('impute',SimpleImputer(strategy='median')),('scale',StandardScaler())]),cont))
 prep=ColumnTransformer(trs,remainder='drop',sparse_threshold=0,verbose_feature_names_out=False)
 steps=[('preprocess',prep)]
 if selector:
  kind,k=selector
  if kind=='mutual_information':steps.append(('select',SelectKBest(lambda a,b:mutual_info_classif(a,b,random_state=seed),k=k)))
  elif kind=='anova':steps.append(('select',SelectKBest(f_classif,k=k)))
  elif kind=='l1':steps.append(('select',SelectFromModel(LogisticRegression(l1_ratio=1.0,solver='liblinear',C=.08,class_weight='balanced',max_iter=2000,random_state=seed),threshold=-np.inf,max_features=k)))
  elif kind=='tree_importance':steps.append(('select',SelectFromModel(ExtraTreesClassifier(n_estimators=60,class_weight='balanced',n_jobs=N_JOBS,random_state=seed,min_samples_leaf=2),threshold=-np.inf,max_features=k)))
 if model=='logistic_regression':clf=LogisticRegression(class_weight='balanced',C=.1,solver='liblinear',max_iter=2500,random_state=seed)
 elif model=='random_forest':clf=RandomForestClassifier(n_estimators=TREES,class_weight='balanced',max_features='sqrt',min_samples_leaf=2,n_jobs=N_JOBS,random_state=seed)
 elif model=='extra_trees':clf=ExtraTreesClassifier(n_estimators=TREES,class_weight='balanced',max_features='sqrt',min_samples_leaf=2,n_jobs=N_JOBS,random_state=seed)
 elif model=='rbf_svm':clf=SVC(kernel='rbf',C=1.0,gamma='scale',class_weight='balanced',probability=True,random_state=seed)
 elif model=='hist_gradient_boosting':clf=HistGradientBoostingClassifier(max_iter=100,learning_rate=.08,max_leaf_nodes=15,l2_regularization=1.0,class_weight='balanced',early_stopping=False,random_state=seed)
 else:raise ValueError(model)
 steps.append(('model',clf));return Pipeline(steps)

def composite_rank(means):
 # rank by each prespecified priority metric, then average ranks (lower is better).
 ranks={}
 for m in PRIORITY:
  vals={k:v[m] for k,v in means.items()};ordered=sorted(vals,key=lambda k:vals[k],reverse=True)
  for rank,k in enumerate(ordered,1):ranks[k]=ranks.get(k,0)+rank/len(PRIORITY)
 return ranks

def feature_modality(c):
 return next((m for m,p in PREFIX.items() if c.startswith(p)),'unknown')

def audit_previous(df,y,train_idx,test_idx):
 source=(BASE/'run_experiments.py').read_text(encoding='utf-8')
 phase3=json.loads((ADV/'reports'/'implementation_audit.json').read_text(encoding='utf-8'))
 metrics=pd.read_csv(BASE/'metrics_summary.csv');pred_dir=BASE/'predictions';pfiles=list(pred_dir.glob('*_test_predictions.csv'))
 audit=[];target_found=False
 expected=set(test_idx.tolist())
 for modelp in list((BASE/'models').glob('*.joblib')):
  pipe=joblib.load(modelp)
  if 'preprocess' not in pipe.named_steps:raise ValueError(f'No preprocess step in {modelp}')
  if 'vital.status' in list(pipe.named_steps['preprocess'].feature_names_in_):target_found=True
 all_pred={}
 for p in pfiles:
  pred=pd.read_csv(p);name=p.stem.removesuffix('_test_predictions')
  if len(pred)!=len(test_idx) or set(pred.row_index_original_zero_based)!=expected:raise ValueError(f'Inconsistent test indices in {p}')
  if not np.array_equal(pred.true_vital_status.to_numpy(),y.iloc[test_idx].to_numpy()):raise ValueError(f'True labels mismatch: {p}')
  all_pred[name]=pred
  if name.startswith('cnv_dedup__'):exp='cnv_dedup_'+name.split('__')[1];mod='logistic_regression';fs='none'
  elif '__logistic_regression__' in name:exp=name.split('__')[0];mod='logistic_regression';fs=name.rsplit('__',1)[1]
  else:exp,mod=name.split('__',1);fs='none'
  row=metrics[(metrics.experiment==exp)&(metrics.model==mod)&(metrics.feature_selection.fillna('none')==fs)]
  if len(row)!=1:raise ValueError(f'Metric row lookup failure {name}')
  cm=confusion_matrix(pred.true_vital_status,pred.predicted_vital_status,labels=[0,1]);tn,fp,fn,tp=cm.ravel()
  if (tn,fp,fn,tp)!=(int(row.iloc[0].test_tn),int(row.iloc[0].test_fp),int(row.iloc[0].test_fn),int(row.iloc[0].test_tp)):raise ValueError(f'Confusion matrix mismatch {name}')
  audit.append({'model':name,'test_rows':len(pred),'same_test_index_set':True,'true_labels_match':True,'confusion_order':'[0,1], TN/FP/FN/TP','confusion_matches_summary':True})
 findings={'saved_prediction_files':len(pfiles),'saved_model_files':len(list((BASE/'models').glob('*.joblib'))),
  'reconciled_predictions':len(audit),'target_in_saved_pipeline_preprocessor':target_found,
  'phase2_target_dropped_in_source':"drop(columns=[target_col])" in source,
  'transforms_and_selection_inside_pipeline':'steps = [("preprocess", prep)]' in source and 'cross_validate(pipe,X_train,y_train' in source,
  'phase3_audit_confirmed':phase3['preprocessing_inside_sklearn_pipelines'] and phase3['feature_selection_inside_pipeline'] and phase3['CV_uses_training_split_only'],
  'phase2_cv_used_training_only':'cross_validate(pipe, X_train, y_train' in source,
  'all_saved_models_use_same_holdout':len(audit)==len(pfiles) and all(x['same_test_index_set'] for x in audit),
  'all_saved_confusion_counts_match':all(x['confusion_matches_summary'] for x in audit),
  'holdout_winner_selection_caveat':'Previous report ranked some models by held-out ROC-AUC descriptively. This is post-hoc and optimistic if interpreted as model selection. No Phase 2 model parameters, selectors, or thresholds were chosen from test scores; this optimization uses training-only inner CV and outer CV, with holdout metrics reported only after selection.'}
 (OUT/'reports'/'previous_methodology_audit.json').write_text(json.dumps({'findings':findings,'per_model':audit},indent=2),encoding='utf-8')
 return findings,all_pred

def feature_qc(df):
 cols=[c for c in df.columns if c!='vital.status'];arr=df[cols].to_numpy(dtype=float);n,p=arr.shape
 summary=[];cnv_dist=[];global_cat={}
 alias=pd.read_csv(BASE/'duplicate_alias_map.csv');alias_count=len(alias)
 for j,c in enumerate(cols):
  x=arr[:,j];counts=pd.Series(x).value_counts(dropna=False);uniq=len(counts);miss=int(np.isnan(x).sum());is_cnv=c.startswith('cn_');is_mu=c.startswith('mu_');
  q1,q3=np.nanpercentile(x,[25,75]);iqr=q3-q1;out=int(np.sum((x<q1-1.5*iqr)|(x>q3+1.5*iqr))) if not is_cnv and not is_mu else 0
  top_share=float(counts.iloc[0]/n) if uniq else np.nan
  summary.append({'feature':c,'modality':feature_modality(c),'n':n,'missing':miss,'missing_pct':miss/n,
   'n_unique':uniq,'constant':uniq<=1,'near_constant_95pct':top_share>=.95,'dominant_value_fraction':top_share,
   'zero_fraction':float(np.mean(x==0)),'mean':float(np.nanmean(x)),'sd':float(np.nanstd(x,ddof=1)),
   'min':float(np.nanmin(x)),'q01':float(np.nanpercentile(x,1)),'q25':float(q1),'median':float(np.nanmedian(x)),
   'q75':float(q3),'q99':float(np.nanpercentile(x,99)),'max':float(np.nanmax(x)),
   'iqr_outlier_count':out,'iqr_outlier_pct':out/n if not is_cnv and not is_mu else np.nan})
  if is_cnv:
   for val,cnt in counts.items():cnv_dist.append({'feature':c,'category':val,'count':int(cnt),'fraction':cnt/n})
   global_cat.update({})
 qc=pd.DataFrame(summary);qc.to_csv(OUT/'feature_quality_summary.csv',index=False);pd.DataFrame(cnv_dist).to_csv(OUT/'tables'/'cnv_category_distribution.csv',index=False)
 # Exhaustive Pearson correlation screen over all numeric feature columns.
 mu=np.nanmean(arr,axis=0);sd=np.nanstd(arr,axis=0,ddof=1);z=(arr-mu)/sd
 z=np.nan_to_num(z,nan=0,posinf=0,neginf=0)
 corr=(z.T@z)/(n-1);np.fill_diagonal(corr,0)
 ii,jj=np.where(np.triu(np.abs(corr)>=.95,k=1));pairs=pd.DataFrame({'feature_1':[cols[i] for i in ii],'feature_2':[cols[j] for j in jj],'pearson_r':corr[ii,jj]})
 pairs['modality_1']=pairs.feature_1.map(feature_modality);pairs['modality_2']=pairs.feature_2.map(feature_modality)
 pairs['exact_duplicate_vector']=pairs.apply(lambda r:r.feature_1 in set(alias.representative_feature)|set(alias.alias_feature) and r.feature_2 in set(alias.representative_feature)|set(alias.alias_feature) and abs(r.pearson_r-1)<1e-12,axis=1)
 pairs.to_csv(OUT/'correlation_analysis.csv',index=False)
 # Mutation sparsity and category summaries.
 muts=[c for c in cols if c.startswith('mu_')];mus=arr[:,[cols.index(c) for c in muts]]
 sparsity={'mutation_feature_count':len(muts),'fraction_zero_all_mutation_cells':float(np.mean(mus==0)),
  'fraction_nonzero_all_mutation_cells':float(np.mean(mus!=0)),'mean_nonzero_features_per_sample':float(np.mean(np.sum(mus!=0,axis=1))),
  'median_nonzero_features_per_sample':float(np.median(np.sum(mus!=0,axis=1))),
  'nonbinary_mutation_features':{c:sorted(set(df[c].dropna().unique().tolist())) for c in muts if not set(df[c].dropna().unique()).issubset({0,1})}}
 (OUT/'reports'/'mutation_sparsity.json').write_text(json.dumps(sparsity,indent=2),encoding='utf-8')
 # Modality-level distribution summary.
 modrows=[]
 for m,pref in PREFIX.items():
  cc=[c for c in cols if c.startswith(pref)];aa=df[cc].to_numpy(float)
  modrows.append({'modality':m,'features':len(cc),'missing_cells':int(np.isnan(aa).sum()),'constant_features':int(qc[qc.modality==m].constant.sum()),
   'near_constant_95pct_features':int(qc[qc.modality==m].near_constant_95pct.sum()),
   'median_feature_zero_fraction':float(np.median(np.mean(aa==0,axis=0))),
   'median_feature_unique_values':float(np.median([len(np.unique(aa[:,i])) for i in range(aa.shape[1])])),
   'iqr_outlier_cells_continuous':int(qc[qc.modality==m].iqr_outlier_count.sum())})
 pd.DataFrame(modrows).to_csv(OUT/'tables'/'modality_quality_summary.csv',index=False)
 # Figures for QC.
 fig,ax=plt.subplots(figsize=(7,4));ax.bar([x['modality'] for x in modrows],[x['near_constant_95pct_features'] for x in modrows],color='#577590');ax.set_ylabel('Features with dominant value ≥95%');ax.set_title('Near-constant feature counts by modality');fig.tight_layout();fig.savefig(OUT/'figures'/'qc_near_constant.png',dpi=300);fig.savefig(OUT/'figures'/'qc_near_constant.pdf');plt.close(fig)
 fig,axs=plt.subplots(2,2,figsize=(10,7))
 for ax,m in zip(axs.flat,['cnv','rna','mutation','protein']):
  c=[c for c in cols if c.startswith(PREFIX[m])];vals=arr[:,[cols.index(x) for x in c]].ravel()
  if m=='cnv':
   uu,nn=np.unique(vals,return_counts=True);ax.bar([str(x) for x in uu],nn/len(vals),color='#577590');ax.set_ylabel('Cell fraction')
  else:
   vals=vals[np.isfinite(vals)];ax.hist(vals,bins=35,color='#43aa8b',alpha=.85,density=True);ax.set_ylabel('Density')
  ax.set_title(m.upper());ax.set_xlabel('Observed value')
 fig.suptitle('Observed feature-value distributions by modality');fig.tight_layout();fig.savefig(OUT/'figures'/'qc_modality_distributions.png',dpi=300);fig.savefig(OUT/'figures'/'qc_modality_distributions.pdf');plt.close(fig)
 # Correlation heatmap from 40 highly variable features (plot subset only; CSV screen was exhaustive).
 vari=np.nanvar(arr,axis=0);sel=np.argsort(vari)[-40:];sub=corr[np.ix_(sel,sel)];fig,ax=plt.subplots(figsize=(10,9));im=ax.imshow(sub,vmin=-1,vmax=1,cmap='coolwarm');labs=[cols[i] for i in sel];ax.set_xticks(range(40),labs,rotation=90,fontsize=6);ax.set_yticks(range(40),labs,fontsize=6);fig.colorbar(im,ax=ax,label='Pearson correlation');ax.set_title('Correlation heatmap: 40 highest-variance features');fig.tight_layout();fig.savefig(OUT/'figures'/'qc_correlation_heatmap.png',dpi=300);fig.savefig(OUT/'figures'/'qc_correlation_heatmap.pdf');plt.close(fig)
 dups=[];g=alias.groupby('representative_feature').alias_feature.apply(list)
 dup_summary={'features':p,'missing_cells':int(np.isnan(arr).sum()),'duplicate_rows':int(len(df)-len(df.drop_duplicates())),
  'exact_duplicate_vector_groups':int(g.shape[0]),'redundant_extra_columns':alias_count,'constant_features':int(qc.constant.sum()),
  'near_constant_95pct_features':int(qc.near_constant_95pct.sum()),'pearson_pairs_abs_ge_0_95':int(len(pairs)),
  'correlation_cutoff':.95,'iqr_outlier_counts_continuous_features_only':int(qc.iqr_outlier_count.sum()),'mutation_sparsity':sparsity}
 (OUT/'reports'/'feature_quality_summary.json').write_text(json.dumps(dup_summary,indent=2),encoding='utf-8')
 report=f'''# Dataset C feature-quality report\n\nThe primary matrix has {n} rows and {p} predictors. No source values were modified. The complete per-feature table is `feature_quality_summary.csv`; high absolute correlations (|Pearson r| ≥ 0.95) are in `correlation_analysis.csv`.\n\n## Missingness and constants\n\nMissing cells in the primary modelling matrix: {int(np.isnan(arr).sum())}. Constant features: {int(qc.constant.sum())}; features with one value in at least 95% of rows: {int(qc.near_constant_95pct.sum())}. Feature-level dominance, zero fractions, unique counts and distribution quantiles are in the summary CSV.\n\n## Duplicate and redundant features\n\nThe previous audit found {int(g.shape[0])} groups of exactly identical vectors and {alias_count} extra CNV aliases. These are identical observed data vectors, not biologically identical genes. The exhaustive numeric correlation screen reports {len(pairs)} pairs at |r|≥0.95, including exact CNV aliases; see CSV for identifiers and correlations.\n\n## Modality distributions\n\n- CNV observed categories and per-feature counts are in `tables/cnv_category_distribution.csv`; the modelling pipeline treats them categorically.\n- Mutation matrix zero fraction is {sparsity['fraction_zero_all_mutation_cells']:.4f}; the mean nonzero mutation-feature count per sample is {sparsity['mean_nonzero_features_per_sample']:.2f}. Observed values outside binary 0/1 are listed in `reports/mutation_sparsity.json` if any.\n- RNA and protein feature distributions, quantiles and IQR-screened outlier counts are in `feature_quality_summary.csv`. IQR flags are descriptive, not automatic errors; no rows or values were removed.\n\n## Correlation and outlier caveats\n\nPearson correlations were computed across all numeric columns, including categorical GISTIC-style CNV codes treated numerically for this screening only. This highlights redundancy but is not a biological association analysis. IQR outlier counts were computed only for RNA/protein continuous features; CNV states and sparse mutation indicators were not treated as continuous outlier candidates.\n'''
 (OUT/'feature_quality_report.md').write_text(report,encoding='utf-8')
 return qc,pairs,modrows,sparsity

def run_nested(Xtr,ytr,columns):
 models=['logistic_regression','random_forest','extra_trees','rbf_svm','hist_gradient_boosting']
 scopes=['cnv','rna','mutation','protein','integrated'];outer=StratifiedKFold(n_splits=OUTER_FOLDS,shuffle=True,random_state=SEED)
 outerrows=[];innerr=[];oof=[];errorrows=[];scope_inner={}
 for scope in scopes:
  log('NESTED start '+scope)
  for ofold,(oi,ov) in enumerate(outer.split(Xtr,ytr),1):
   Xi,yi=Xtr.iloc[oi],ytr.iloc[oi];Xv,yv=Xtr.iloc[ov],ytr.iloc[ov]
   inner=StratifiedKFold(n_splits=INNER_FOLDS,shuffle=True,random_state=SEED+ofold)
   summary={}
   for model in models:
    try:
     pipe=make_pipe(scope,model,columns)
     cv=cross_validate(pipe,Xi,yi,cv=inner,scoring=SCORING,n_jobs=1,error_score='raise',return_train_score=False)
     means={m:float(np.mean(cv['test_'+m])) for m in SCORING};summary[model]=means
     for fold in range(INNER_FOLDS):
      innerr.append({'scope':scope,'outer_fold':ofold,'model':model,'inner_fold':fold+1,**{m:float(cv['test_'+m][fold]) for m in SCORING},'status':'ok'})
    except Exception as e:
     errorrows.append({'scope':scope,'outer_fold':ofold,'model':model,'error':f'{type(e).__name__}: {e}'})
     log(f'ERROR nested {scope}/{model}/outer{ofold}: {type(e).__name__}: {e}')
   if not summary:continue
   ranks=composite_rank(summary);best=sorted(summary,key=lambda m:(ranks[m],-summary[m]['roc_auc'],-summary[m]['pr_auc']))[0]
   pipe=make_pipe(scope,best,columns);pipe.fit(Xi,yi);prob=pipe.predict_proba(Xv)[:,1];pred=pipe.predict(Xv);mm=metric_dict(yv,pred,prob)
   outerrows.append({'scope':scope,'outer_fold':ofold,'selected_algorithm_by_inner_cv':best,'inner_composite_mean_rank':ranks[best],**mm})
   for rowix,yt,predicted,pr in zip(ov,yv,pred,prob):oof.append({'scope':scope,'training_row_index':int(rowix),'true_vital_status':int(yt),'predicted_vital_status':int(predicted),'probability_deceased':float(pr),'selected_algorithm_by_inner_cv':best,'outer_fold':ofold})
   # Store outer-training-fold model importances for stable feature reporting where available.
   estimator=pipe.named_steps['model'];names=np.asarray(pipe.named_steps['preprocess'].get_feature_names_out())
   importance=getattr(estimator,'feature_importances_',None)
   if importance is None and hasattr(estimator,'coef_'):importance=np.mean(np.abs(estimator.coef_),axis=0)
   if importance is not None and len(importance)==len(names):
    rawcols=[c for c in columns if c.startswith(PREFIX['cnv'])] if scope=='cnv' else [c for c in columns if c.startswith(PREFIX[scope])] if scope!='integrated' else list(columns)
    agg={}
    for nm,im in zip(names,importance):
     raw=nm if nm in rawcols else max([c for c in rawcols if nm.startswith(c+'_')],key=len,default=nm)
     agg[raw]=agg.get(raw,0)+float(im)
    for rank,(feat,imp) in enumerate(sorted(agg.items(),key=lambda z:z[1],reverse=True)[:30],1):
     innerr.append({'scope':scope,'outer_fold':ofold,'model':best,'inner_fold':0,'feature':feat,'modality':feature_modality(feat),'importance':imp,'importance_rank':rank,'status':'outer_train_importance'})
   scope_inner.setdefault(scope,[]).append(summary)
  log('NESTED done '+scope)
 # Aggregate inner candidate performance across the 3 outer-train partitions; this uses training data only.
 cand=[]
 for scope,fold_summaries in scope_inner.items():
  for model in models:
   vals=[f[model] for f in fold_summaries if model in f]
   if not vals:continue
   mean={m:float(np.mean([v[m] for v in vals])) for m in PRIORITY};cand.append({'scope':scope,'model':model,'outer_training_inner_cv_replicates':len(vals),**{f'inner_{m}_mean':v for m,v in mean.items()}})
 byscope={s:[r for r in cand if r['scope']==s] for s in scopes};scope_winners={}
 for scope,rr in byscope.items():
  if rr:
   # Rank all algorithms by the five prespecified average-CV metrics; ROC/AP/balanced/F1/recall.
   means={r['model']:{m:r[f'inner_{m}_mean'] for m in PRIORITY} for r in rr};ranks=composite_rank(means)
   for r in rr:r['mean_composite_rank']=ranks[r['model']];r['scope_selected']=False
   selected=sorted(rr,key=lambda r:(r['mean_composite_rank'],-r['inner_roc_auc_mean'],-r['inner_pr_auc_mean']))[0]
   selected['scope_selected']=True;scope_winners[scope]=selected['model']
 outerdf=pd.DataFrame(outerrows);innerdf=pd.DataFrame(innerr);canddf=pd.DataFrame(cand)
 outerdf.to_csv(OUT/'tables'/'nested_outer_results.csv',index=False);innerdf.to_csv(OUT/'tables'/'nested_inner_and_importance.csv',index=False);canddf.to_csv(OUT/'tables'/'model_selection_summary.csv',index=False)
 pd.DataFrame(oof).to_csv(OUT/'predictions'/'nested_oof_predictions.csv',index=False)
 pd.DataFrame(errorrows,columns=['scope','outer_fold','model','error']).to_csv(OUT/'tables'/'model_fit_errors.csv',index=False)
 # Outer-fold and candidate summary.
 agg=outerdf.groupby('scope').agg({m:['mean','std'] for m in ['roc_auc','pr_auc','balanced_accuracy','macro_f1','deceased_recall','specificity']})
 agg.columns=['_'.join(c) for c in agg.columns];agg.reset_index().to_csv(OUT/'tables'/'nested_cv_summary.csv',index=False)
 return models,scopes,canddf,outerdf,scope_winners,errorrows

def run_selectors(Xtr,ytr,columns):
 records=[];chosen=[];fold=StratifiedKFold(n_splits=SELECTOR_FOLDS,shuffle=True,random_state=SEED)
 for scope in ['rna','integrated']:
  for method in ['mutual_information','anova','l1','tree_importance']:
   pipe_sel=(method,50)
   for fi,(itr,iva) in enumerate(fold.split(Xtr,ytr),1):
    pipe=make_pipe(scope,'logistic_regression',columns,selector=pipe_sel);pipe.fit(Xtr.iloc[itr],ytr.iloc[itr])
    prob=pipe.predict_proba(Xtr.iloc[iva])[:,1];pred=pipe.predict(Xtr.iloc[iva]);metrics=metric_dict(ytr.iloc[iva],pred,prob)
    records.append({'scope':scope,'selection_method':method,'fold':fi,**metrics})
    names=np.asarray(pipe.named_steps['preprocess'].get_feature_names_out());support=pipe.named_steps['select'].get_support();selnames=names[support]
    rawcols=[c for c in columns if c.startswith(PREFIX[scope])] if scope!='integrated' else list(columns)
    for en in selnames:
     raw=en if en in rawcols else max([c for c in rawcols if en.startswith(c+'_')],key=len,default=en)
     chosen.append({'scope':scope,'selection_method':method,'fold':fi,'selected_encoded_feature':en,'feature':raw,'modality':feature_modality(raw)})
   log(f'SELECTOR {scope}/{method} complete')
 res=pd.DataFrame(records);res.to_csv(OUT/'feature_selection'/'selector_fold_metrics.csv',index=False)
 sel=pd.DataFrame(chosen);sel.to_csv(OUT/'feature_selection'/'selected_features_by_fold.csv',index=False)
 freq=sel.groupby(['scope','selection_method','feature','modality']).fold.nunique().reset_index(name='folds_selected')
 freq['selection_frequency']=freq.folds_selected/SELECTOR_FOLDS
 freq['stable_4_of_5_or_more']=freq.folds_selected>=4
 freq.to_csv(OUT/'feature_selection'/'feature_selection_frequency.csv',index=False)
 summary=res.groupby(['scope','selection_method']).agg({m:['mean','std'] for m in ['roc_auc','pr_auc','balanced_accuracy','macro_f1','deceased_recall','specificity']})
 summary.columns=['_'.join(c) for c in summary.columns];summary.reset_index().to_csv(OUT/'feature_selection'/'selector_summary.csv',index=False)
 return res,sel,freq

def run_dummy_baselines(Xtr,ytr,Xte,yte,test_indices,columns):
 rows=[];predrows=[];scopes=['cnv','rna','mutation','protein','integrated']
 fold=StratifiedKFold(n_splits=OUTER_FOLDS,shuffle=True,random_state=SEED)
 for scope in scopes:
  cols=[c for c in columns if any(c.startswith(PREFIX[s]) for s in (list(PREFIX) if scope=='integrated' else [scope]))]
  Xscope=Xtr[cols];oofprob=np.zeros(len(Xtr));oofpred=np.zeros(len(Xtr),dtype=int)
  for itr,iva in fold.split(Xscope,ytr):
   clf=DummyClassifier(strategy='prior').fit(Xscope.iloc[itr],ytr.iloc[itr]);oofprob[iva]=clf.predict_proba(Xscope.iloc[iva])[:,1];oofpred[iva]=clf.predict(Xscope.iloc[iva])
  cvm=metric_dict(ytr,oofpred,oofprob)
  clf=DummyClassifier(strategy='prior').fit(Xscope,ytr);prob=clf.predict_proba(Xte[cols])[:,1];pred=clf.predict(Xte[cols]);tm=metric_dict(yte,pred,prob)
  rows.append({'scope':scope,'model':'dummy_prior',**{f'cv_{k}':v for k,v in cvm.items()},**{f'test_{k}':v for k,v in tm.items()}})
  for i,(yt,pdv,pr) in enumerate(zip(yte,pred,prob)):
   predrows.append({'scope':scope,'model':'dummy_prior','original_row_index_zero_based':int(test_indices[i]),'true_vital_status':int(yt),'predicted_vital_status':int(pdv),'probability_deceased':float(pr)})
 pd.DataFrame(rows).to_csv(OUT/'tables'/'dummy_baseline_comparison.csv',index=False)
 pd.DataFrame(predrows).to_csv(OUT/'predictions'/'dummy_baseline_test_predictions.csv',index=False)
 return pd.DataFrame(rows)

def main():
 if LOG.exists():LOG.unlink()
 log(f'Optimization run; Python={sys.version.split()[0]}, sklearn={sklearn.__version__}, seed={SEED}; outer={OUTER_FOLDS}, inner={INNER_FOLDS}')
 df=pd.read_csv(DATA);enr=pd.read_csv(ENR);y=df['vital.status'].astype(int);X=df.drop(columns=['vital.status']);columns=list(X.columns)
 if set(y.unique())!={0,1} or X.shape!=(705,1936):raise ValueError(f'Unexpected dataset shape/target: {X.shape}, {y.value_counts().to_dict()}')
 unexpected=[c for c in columns if not c.startswith(tuple(PREFIX.values()))]
 if unexpected:raise ValueError(f'Unexpected non-modality predictors: {unexpected}')
 train_idx,test_idx=train_test_split(np.arange(len(df)),test_size=TEST_SIZE,random_state=SEED,stratify=y)
 old=json.loads((BASE/'models'/'rna__random_forest.json').read_text(encoding='utf-8'))
 if train_idx.tolist()!=old['train_indices_zero_based'] or test_idx.tolist()!=old['test_indices_zero_based']:raise ValueError('Reconstructed phase 2 split mismatch')
 findings,oldpred=audit_previous(df,y,train_idx,test_idx)
 auditmd=f'''# Methodology audit\n\nAudited the previous Phase 2/3 scripts, all 23 saved central pipelines, 23 prediction CSVs, metric tables, logs, and saved split metadata. No source CSV or previous output was changed.\n\n- `vital.status` is binary 0 alive / 1 deceased and is dropped from X before fitting. None of the 23 saved preprocessors includes the target.\n- All 23 prediction files use the same exact 141 test row indices and true labels (122 alive, 19 deceased). Recomputed confusion matrices use class order `[0,1]`, so flattened cells are TN, FP, FN, TP; all counts match saved metrics.\n- Imputation, scaling, one-hot encoding, mutation variance filtering, and selectors are pipeline steps. Phase 2 CV ran on the 564-row training partition; selectors are inside CV pipelines. Phase 3 audit independently agrees.\n- Test scores were used descriptively to list highest held-out results in the previous report, creating post-hoc winner ranking. No selector, model hyperparameter, or threshold was chosen from test results. Treat all maximum test comparisons as optimistic/descriptive.\n- Original CNV dedup indexing error and initial worker error were fixed and documented in earlier logs. Saved final results reconcile.\n- No IDs or follow-up time exist; row-based splits cannot verify patient-level independence.\n\nAudit summary: `{json.dumps(findings,ensure_ascii=False)}`\n'''
 (OUT/'methodology_audit.md').write_text(auditmd,encoding='utf-8')
 # Data quality, correlations, category and outlier analysis.
 qc,pairs,modrows,sparsity=feature_qc(df)
 # Dataset and dependency manifest, plus optional availability.
 optional={p:bool(importlib.util.find_spec(p)) for p in ['xgboost','lightgbm','imblearn','shap']}
 (OUT/'reproducibility'/'optional_package_availability.json').write_text(json.dumps(optional,indent=2),encoding='utf-8')
 with (OUT/'reproducibility'/'package_versions.json').open('w',encoding='utf-8') as f:json.dump({'python':sys.version,'platform':platform.platform(),'scikit_learn':sklearn.__version__,'numpy':np.__version__,'pandas':pd.__version__,'seed':SEED,'outer_folds':OUTER_FOLDS,'inner_folds':INNER_FOLDS,'selector_folds':SELECTOR_FOLDS,'optional_packages':optional},f,indent=2)
 # Save split and preprocessing manifest.
 pd.DataFrame({'original_row_index_zero_based':np.arange(len(df)),'split':['train' if i in set(train_idx) else 'test' for i in range(len(df))],'vital.status':y}).to_csv(OUT/'reproducibility'/'train_test_indices.csv',index=False)
 (OUT/'reproducibility'/'preprocessing_config.json').write_text(json.dumps({'CNV':'most-frequent imputation, OneHotEncoder(handle_unknown=ignore)','RNA':'median imputation + StandardScaler','mutation':'median imputation + VarianceThreshold fitted per fold; no resampling','protein':'median imputation + StandardScaler','all':'ColumnTransformer inside Pipeline; fit only on corresponding training fold','target_excluded':True},indent=2),encoding='utf-8')
 (OUT/'reproducibility'/'cv_config.json').write_text(json.dumps({'selection':'nested stratified CV; 3 outer folds × 3 inner folds on original training partition','outer_shuffle':True,'inner_shuffle':True,'seed':SEED,'priority_metrics':PRIORITY,'selection_method':'mean rank over ROC-AUC, average precision, balanced accuracy, macro-F1, deceased recall','final_test':'original 141-row holdout evaluated only after model selection; not used for selection'},indent=2),encoding='utf-8')
 # Nested comparison; handles expected failures per model internally.
 Xtr=X.iloc[train_idx].reset_index(drop=True);ytr=y.iloc[train_idx].reset_index(drop=True);Xte=X.iloc[test_idx].reset_index(drop=True);yte=y.iloc[test_idx].reset_index(drop=True)
 models,scopes,cand,outer,scope_winners,errors=run_nested(Xtr,ytr,columns)
 baseline=run_dummy_baselines(Xtr,ytr,Xte,yte,test_idx,columns)
 # Selector comparison on RNA/integrated with all selector fitting inside each of five training folds.
 sr,selected,frequency=run_selectors(Xtr,ytr,columns)
 # Select final scope/model solely from training CV: best algorithm per scope by five-metric rank,
 # then best scope/model by the same training-CV composite (not final test score).
 valid=cand[cand.scope_selected==True].copy()
 # Use dict rows because model identifiers and metric names are dynamic.
 global_means={f"{r['scope']}__{r['model']}":{m:r[f'inner_{m}_mean'] for m in PRIORITY} for r in valid.to_dict('records')}
 global_ranks=composite_rank(global_means);bestkey=sorted(global_means,key=lambda k:(global_ranks[k],-global_means[k]['roc_auc'],-global_means[k]['pr_auc']))[0]
 bestscope,bestmodel=bestkey.split('__',1)
 # Fit one final winner per scope on all training rows; evaluate each once on untouched test for descriptive reporting.
 finalrows=[];predrows=[];finalmodels=[]
 for scope,model in scope_winners.items():
  pipe=make_pipe(scope,model,columns);pipe.fit(Xtr,ytr);prob=pipe.predict_proba(Xte)[:,1];pred=pipe.predict(Xte);m=metric_dict(yte,pred,prob)
  name=f'{scope}__{model}';joblib.dump(pipe,OUT/'models'/f'{name}.joblib',compress=3)
  conf={'scope':scope,'model':model,'features': [c for c in columns if any(c.startswith(PREFIX[s]) for s in (list(PREFIX) if scope=='integrated' else [scope]))],
   'seed':SEED,'train_rows':len(train_idx),'test_rows':len(test_idx),'cv_selection':'nested 3x3 CV, training only','selection_metric_priority':PRIORITY,
   'preprocessing':json.loads((OUT/'reproducibility'/'preprocessing_config.json').read_text(encoding='utf-8')),
   'model_params':{k:str(v) for k,v in pipe.named_steps['model'].get_params().items()},
   'train_indices':list(map(int,train_idx)),'test_indices':list(map(int,test_idx))}
  (OUT/'models'/f'{name}.json').write_text(json.dumps(conf,indent=2),encoding='utf-8')
  finalrows.append({'scope':scope,'model':model,'selected_by_training_CV':True,'global_best_by_CV':name==bestkey,'selection_composite_rank':global_ranks[name],**m})
  for ix,yt,predicted,pr in zip(test_idx,yte,pred,prob):predrows.append({'model':name,'original_row_index_zero_based':int(ix),'true_vital_status':int(yt),'predicted_vital_status':int(predicted),'probability_deceased':float(pr)})
  finalmodels.append((name,pipe))
  log(f'FINAL train/evaluate {name}: CV-inner-ranked; test ROC-AUC={m["roc_auc"]:.4f}, PR-AUC={m["pr_auc"]:.4f}')
 pd.DataFrame(finalrows).to_csv(OUT/'tables'/'selected_scope_models_test_metrics.csv',index=False)
 pd.DataFrame(predrows).to_csv(OUT/'predictions'/'selected_scope_test_predictions.csv',index=False)
 (OUT/'reports'/'final_model_selection.json').write_text(json.dumps({'selected_key':bestkey,'scope':bestscope,'model':bestmodel,'selection_source':'training-only nested CV composite rank; never test metrics','composite_rank':global_ranks[bestkey],'scope_winners':scope_winners},indent=2),encoding='utf-8')
 # CV-only selected feature list with support count from fold-contained selectors.
 bestfeat=frequency[(frequency.scope==bestscope)].sort_values(['selection_frequency','folds_selected'],ascending=False)
 bestfeat.to_csv(OUT/'feature_selection'/'best_cv_scope_feature_frequency.csv',index=False)
 # Fold importance stability plot based only on outer training fits in nested CV.
 imp= pd.read_csv(OUT/'tables'/'nested_inner_and_importance.csv'); imp=imp[imp.status=='outer_train_importance'] if 'status' in imp else pd.DataFrame()
 imp.to_csv(OUT/'feature_selection'/'outer_fold_feature_importance.csv',index=False)
 if len(imp):
  st=imp.groupby(['scope','feature','modality']).agg(folds_top30=('outer_fold','nunique'),mean_importance=('importance','mean'),mean_rank=('importance_rank','mean')).reset_index()
  st['fold_frequency']=st.folds_top30/OUTER_FOLDS;st.to_csv(OUT/'feature_selection'/'outer_fold_feature_stability.csv',index=False)
  show=st[st.scope==bestscope].sort_values('mean_importance').tail(20);fig,ax=plt.subplots(figsize=(8,7));ax.barh(show.feature,show.mean_importance,color='#577590');ax.set_title(f'{bestscope} outer-training fold importance');ax.set_xlabel('Mean model importance');fig.tight_layout();fig.savefig(OUT/'figures'/'best_model_feature_importance.png',dpi=300);fig.savefig(OUT/'figures'/'best_model_feature_importance.pdf');plt.close(fig)
 # Final model holdout plots for CV-selected per-scope configurations.
 finalpred=pd.read_csv(OUT/'predictions'/'selected_scope_test_predictions.csv');fig,ax=plt.subplots(figsize=(7,6))
 for scope,grp in finalpred.groupby('model'):
  fpr,tpr,_=roc_curve(grp.true_vital_status,grp.probability_deceased);ax.plot(fpr,tpr,label=f'{scope} AUC={roc_auc_score(grp.true_vital_status,grp.probability_deceased):.3f}')
 ax.plot([0,1],[0,1],'k--');ax.set(xlabel='False-positive rate',ylabel='True-positive rate',title='Final training-CV-selected models — exploratory holdout');ax.legend();fig.tight_layout();fig.savefig(OUT/'figures'/'optimized_roc_curves.png',dpi=300);fig.savefig(OUT/'figures'/'optimized_roc_curves.pdf');plt.close(fig)
 fig,axs=plt.subplots(2,3,figsize=(12,7))
 for ax,(model,grp) in zip(axs.flat,finalpred.groupby('model')):
  cm=confusion_matrix(grp.true_vital_status,grp.predicted_vital_status,labels=[0,1]);ax.imshow(cm,cmap='Blues');ax.set_title(model);ax.set_xticks([0,1],['Alive','Deceased']);ax.set_yticks([0,1],['Alive','Deceased']);
  for (i,j),v in np.ndenumerate(cm):ax.text(j,i,str(v),ha='center',va='center')
 fig.suptitle('Optimized final-scope confusion matrices; threshold 0.5');fig.tight_layout();fig.savefig(OUT/'figures'/'optimized_confusion_matrices.png',dpi=300);fig.savefig(OUT/'figures'/'optimized_confusion_matrices.pdf');plt.close(fig)
 # Candidate availability and failure log.
 (OUT/'tables'/'optional_algorithm_status.json').write_text(json.dumps({'hist_gradient_boosting':'available in scikit-learn; nested CV completed','xgboost':'not installed; skipped without installation','lightgbm':'not installed; skipped without installation','balanced_random_forest':'imbalanced-learn unavailable; skipped','SMOTE':'imbalanced-learn unavailable; no synthetic resampling performed','SHAP':'not installed; skipped'},indent=2),encoding='utf-8')
 # Statistical report.
 nested=pd.read_csv(OUT/'tables'/'nested_cv_summary.csv');cvrow=nested[nested.scope==bestscope].iloc[0]
 oldrep=pd.read_csv(ADV/'reports'/'heldout_permutation_tests.csv')
 report=f'''# Final Dataset C optimization report\n\n## Execution and scientific scope\n\nWork was limited to Dataset C. Original CSVs and previous result files were preserved. New models were trained only on the original 564-row training partition. The original 141-row holdout was evaluated after model selection and was not used to tune models, features, or thresholds. This is exploratory classification of recorded vital status (0 alive, 1 deceased), not survival-time prediction. There are no patient IDs or follow-up durations, so patient-level independence cannot be verified.\n\n## Methodology audit\n\n`methodology_audit.md` records review of 23 saved pipelines and prediction files. All 23 prediction files use identical test rows; recomputed confusion matrices match saved metrics in class order [0,1] (TN,FP,FN,TP). Target leakage was not found. Preprocessing and previous feature selectors are pipeline-contained. Earlier post-hoc display of highest test scores is acknowledged as optimistic ranking; this new selection uses nested training-only CV.\n\n## Feature quality\n\nThere are {len(y)} samples, {X.shape[1]} predictors, 611 alive and 94 deceased. Quality summaries and exhaustive |Pearson r|≥0.95 pairs are in `feature_quality_summary.csv` and `correlation_analysis.csv`; see `feature_quality_report.md`. Mutation sparsity is in `reports/mutation_sparsity.json`; CNV value categories and RNA/protein distributions/outlier screens are summarized separately. A total of 62 exact duplicate-vector groups with 99 redundant CNV aliases remain documented; alias equality does not mean biological identity.\n\n## Model optimization and validation\n\nCompared class-weighted Logistic Regression, Random Forest, Extra Trees, RBF SVM, and HistGradientBoosting across CNV, RNA, mutation, protein, and integrated inputs using nested stratified 3-fold outer / 3-fold inner CV on training rows. Within inner CV, candidates were selected by average rank across ROC-AUC, PR-AUC, balanced accuracy, macro-F1, and deceased recall. Each outer fold evaluated its inner-selected model. `nested_cv_summary.csv`, `model_selection_summary.csv`, `nested_outer_results.csv`, and `model_fit_errors.csv` contain exact results.\n\nThe best overall training-CV scope/algorithm selected by the composite priority ranking was **{bestscope} / {bestmodel}**. Its nested outer metrics were ROC-AUC {cvrow.roc_auc_mean:.3f} ± {cvrow.roc_auc_std:.3f}, PR-AUC {cvrow.pr_auc_mean:.3f} ± {cvrow.pr_auc_std:.3f}, balanced accuracy {cvrow.balanced_accuracy_mean:.3f} ± {cvrow.balanced_accuracy_std:.3f}, macro-F1 {cvrow.macro_f1_mean:.3f} ± {cvrow.macro_f1_std:.3f}, deceased recall {cvrow.deceased_recall_mean:.3f} ± {cvrow.deceased_recall_std:.3f}. These are nested CV estimates on this small training cohort, not independent validation.\n\nPer-scope final models were then refit on all 564 training rows and evaluated once on the original holdout. Exact held-out values are in `selected_scope_models_test_metrics.csv`; these descriptive values were not used to select the global winner. The holdout has only 19 deceased observations. The phase-3 majority baseline remains the comparator.\n\n## Class imbalance and feature selection\n\nClass weighting was used for supported estimators. Balanced Random Forest and SMOTE were not run because imbalanced-learn is not installed; no dependencies were installed. No validation/test samples were resampled. Mutual information, ANOVA SelectKBest, L1, and Extra-Trees importance selectors were compared within five training folds on RNA and integrated matrices; results and fold support are in `feature_selection/`. Deceased recall, PR-AUC, balanced accuracy, macro-F1, and specificity accompany ROC-AUC. Accuracy was not the selector criterion.\n\n## Explainability and modality contribution\n\nOuter-fold training importances and selector frequencies are in `feature_selection/`. These are model-dependent associations, not causal findings or validated biomarkers. SHAP/XGBoost/LightGBM were unavailable and skipped without installation. Initial Dataset C repeated-CV comparisons found RNA RF ROC-AUC {pd.read_csv(ADV/'modality_comparison.csv').query("scope=='rna'").iloc[0].repeated_cv_roc_auc_mean:.3f} ± {pd.read_csv(ADV/'modality_comparison.csv').query("scope=='rna'").iloc[0].repeated_cv_roc_auc_sd_repeat_means:.3f} versus integrated Extra Trees {pd.read_csv(ADV/'modality_comparison.csv').query("scope=='integrated'").iloc[0].repeated_cv_roc_auc_mean:.3f} ± {pd.read_csv(ADV/'modality_comparison.csv').query("scope=='integrated'").iloc[0].repeated_cv_roc_auc_sd_repeat_means:.3f}. Prior held-out permutation tests were not significant (RNA p={oldrep.query("scope=='rna'").iloc[0].one_sided_pvalue:.3f}; integrated p={oldrep.query("scope=='integrated'").iloc[0].one_sided_pvalue:.3f}). No integration benefit is presumed; updated nested results are reported directly.\n\n## Limitations and conclusion\n\nOnly 94 of 705 labels are deceased, including just 19 in the original test set. Patient IDs/barcodes and follow-up duration/censoring data are unavailable; patient-level independence and survival outcomes cannot be verified. There is no independent external validation. These experiments cannot diagnose breast cancer or predict an individual's survival. Results are exploratory recorded-vital-status classification only and are not clinically validated.\n\n## Reproducibility\n\nSee `reproducibility/` for Python/package versions, seeds, split indices, fold configuration, exact model and preprocessing settings, and fitted model paths. New final models are under `models/`; per-scope test predictions are under `predictions/`. The full run log is `logs/optimization.log`.\n'''
 (OUT/'reports'/'optimization_report.md').write_text(report,encoding='utf-8')
 log('Optimization complete; final models/report written.')
 log(f'Candidate errors recorded={len(errors)}; model scope winners={scope_winners}; selected global={bestkey}')

if __name__=='__main__':
 try:main()
 except Exception as e:log(f'FATAL {type(e).__name__}: {e}\n{traceback.format_exc()}');raise
