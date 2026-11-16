"""Assemble a final research/reproducibility package from existing Dataset C outputs.

This script does not fit models. It reads original CSVs and existing result artifacts.
"""
from pathlib import Path
import json, shutil, sys, platform, importlib.metadata
import numpy as np
import pandas as pd
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from sklearn.calibration import calibration_curve
from sklearn.metrics import roc_curve, precision_recall_curve, auc, confusion_matrix
import joblib, sklearn

ROOT=Path(__file__).resolve().parents[2]
BASE=ROOT/'results_dataset_C'; ADV=BASE/'advanced_validation'; OUT=BASE/'final_research'
DATA=ROOT/'dataset_c_extracted'/'data.csv'; ENRICHED=ROOT/'dataset_c_extracted'/'brca_data_w_subtypes.csv'
SEED=20261110
for d in ['figures','tables','reproducibility/model_configs','reproducibility/scripts','reproducibility/feature_lists']:
 (OUT/d).mkdir(parents=True,exist_ok=True)
plt.rcParams.update({'font.family':'DejaVu Sans','font.size':10,'axes.titlesize':13,'axes.labelsize':11,'legend.fontsize':9,'figure.dpi':150,'savefig.dpi':300})

def savefig(fig,name):
 fig.tight_layout();fig.savefig(OUT/'figures'/f'{name}.png',dpi=300,bbox_inches='tight');fig.savefig(OUT/'figures'/f'{name}.pdf',bbox_inches='tight');plt.close(fig)

def raw_feature(encoded, columns):
 if encoded in columns:return encoded
 candidates=[c for c in columns if encoded.startswith(c+'_')]
 return max(candidates,key=len) if candidates else encoded

def main():
 df=pd.read_csv(DATA); enriched=pd.read_csv(ENRICHED)
 phase2=pd.read_csv(BASE/'metrics_summary.csv')
 baseline=pd.read_csv(ADV/'baseline_comparison.csv')
 repeated=pd.read_csv(ADV/'repeated_cv_results.csv'); repeated=repeated[repeated.record_type=='repeat_summary']
 selected=pd.read_csv(BASE/'feature_selection'/'selected_features.csv')
 importance=pd.read_csv(BASE/'explainability'/'top_feature_importance.csv')
 stability=pd.read_csv(ADV/'feature_stability.csv')
 pred_dir=BASE/'predictions'; pred={p.stem.removesuffix('_test_predictions'):pd.read_csv(p) for p in pred_dir.glob('*_test_predictions.csv')}
 y=df['vital.status'].astype(int);features=[c for c in df.columns if c!='vital.status']
 test_idx=json.loads((BASE/'models'/'rna__random_forest.json').read_text(encoding='utf-8'))['test_indices_zero_based']
 test_idx=np.asarray(test_idx,dtype=int); ytest=y.iloc[test_idx].to_numpy(); Xtest=df.iloc[test_idx][features]
 if len(test_idx)!=141 or sorted(pd.Series(ytest).value_counts().to_dict().items())!=[(0,122),(1,19)]: raise ValueError('Expected original held-out split not found')

 # One master table: original 5-fold CV for every trained model, Phase-2 heldout metrics,
 # plus repeated 5x5 CV metrics for the prespecified model per modality.
 rows=[]
 for _,r in phase2.iterrows():
  exp=str(r.experiment);model=str(r.model);selmethod=str(r.feature_selection) if pd.notna(r.feature_selection) else 'none'
  if exp.startswith('cnv_dedup_'):key='cnv_dedup__'+exp.removeprefix('cnv_dedup_')
  elif selmethod!='none':key=f'{exp}__{model}__{selmethod}'
  else:key=f'{exp}__{model}'
  if key not in pred:raise ValueError(f'Missing saved predictions for {key}')
  b=baseline[baseline.model==key]
  if len(b)!=1:raise ValueError(f'Expected one baseline comparison row for {key}, got {len(b)}')
  br=b.iloc[0]
  rep=repeated[(repeated.scope==exp)&(repeated.model==model)]
  rr=rep.iloc[0] if len(rep)==1 else None
  rows.append({'scope':exp,'model':model,'feature_selection':selmethod,'baseline':False,
   'phase2_cv_roc_auc_mean':r.roc_auc_mean,'phase2_cv_roc_auc_sd':r.roc_auc_std,
   'repeated_cv_roc_auc_mean':rr.roc_auc_mean if rr is not None else np.nan,
   'repeated_cv_roc_auc_sd_repeat_means':rr.roc_auc_sd_repeat_means if rr is not None else np.nan,
   'test_accuracy':br.accuracy,'test_balanced_accuracy':br.balanced_accuracy,
   'test_precision_deceased':br.precision_deceased,'test_deceased_recall_sensitivity':br.sensitivity_deceased,
   'test_specificity_alive':br.specificity_alive,'test_f1_deceased':br.f1_deceased,'test_macro_f1':br.macro_f1,
   'test_roc_auc':br.roc_auc,'test_pr_auc_average_precision':br.pr_auc_average_precision,
   'tn':br.tn,'fp':br.fp,'fn':br.fn,'tp':br.tp,
   'cv_stat_note':'Phase 2 5-fold CV on training split; repeated CV fields populated only for prespecified modality winners',
   'test_stat_note':'Exploratory fixed holdout; not used for model/feature/threshold selection'})
 # Dummy comparison row
 dum=baseline[baseline.model=='dummy_majority'].iloc[0]
 rd=repeated[repeated.scope=='dummy_majority'].iloc[0]
 rows.append({'scope':'baseline','model':'DummyClassifier(most_frequent)','feature_selection':'none','baseline':True,
  'phase2_cv_roc_auc_mean':np.nan,'phase2_cv_roc_auc_sd':np.nan,'repeated_cv_roc_auc_mean':rd.roc_auc_mean,
  'repeated_cv_roc_auc_sd_repeat_means':rd.roc_auc_sd_repeat_means,'test_accuracy':dum.accuracy,
  'test_balanced_accuracy':dum.balanced_accuracy,'test_precision_deceased':dum.precision_deceased,
  'test_deceased_recall_sensitivity':dum.sensitivity_deceased,'test_specificity_alive':dum.specificity_alive,
  'test_f1_deceased':dum.f1_deceased,'test_macro_f1':dum.macro_f1,'test_roc_auc':dum.roc_auc,
  'test_pr_auc_average_precision':dum.pr_auc_average_precision,'tn':dum.tn,'fp':dum.fp,'fn':dum.fn,'tp':dum.tp,
  'cv_stat_note':'Repeated 5x5 CV on training split','test_stat_note':'Exploratory fixed holdout'})
 master=pd.DataFrame(rows)
 master['cv_selected_winner_within_scope']=False
 for scope,sub in master[(master.feature_selection=='none')&(~master.baseline)&(~master.scope.str.startswith('cnv_dedup_'))].groupby('scope'):
  ix=sub.phase2_cv_roc_auc_mean.idxmax();master.loc[ix,'cv_selected_winner_within_scope']=True
 master.to_csv(OUT/'master_metrics_comparison.csv',index=False)
 # Compact modality comparison with phase-3 repeated CV values and common test metrics.
 modality=master[(master.cv_selected_winner_within_scope)|master.baseline].copy()
 modality.to_csv(OUT/'tables'/'modality_winner_comparison.csv',index=False)
 for src,name in [(ADV/'repeated_cv_results.csv','repeated_cv_results.csv'),(ADV/'threshold_analysis.csv','threshold_analysis.csv'),(ADV/'calibration_results.csv','calibration_results.csv'),(ADV/'bootstrap_confidence_intervals.csv','bootstrap_confidence_intervals.csv'),(ADV/'reports'/'heldout_permutation_tests.csv','heldout_permutation_tests.csv'),(ADV/'reports'/'integration_minus_rna_paired_bootstrap.json','integration_minus_rna_paired_bootstrap.json')]:
  shutil.copy2(src,OUT/'tables'/name)

 # Final selected/stable feature table. Phase-2 saved selector artifact does not retain per-fold
 # support counts; report those as missing and provide separate repeated tree top-30 fold frequency.
 original_columns=features
 stabmap={(r.scope,r.feature):r for r in stability.itertuples()}
 impmap={}
 for r in importance.itertuples():
  if str(r.importance_method) not in {'test-set permutation ROC-AUC decrease; top-ranked raw variables only'}:
   impmap.setdefault(str(r.feature),float(r.importance))
 finalf=[]
 for r in selected.itertuples():
  scope=str(r.experiment).split('__')[0];encoded=str(r.selected_feature);raw=raw_feature(encoded,original_columns)
  st=stabmap.get((scope,raw))
  mod=str(r.modality)
  finalf.append({'feature_name':encoded,'raw_feature_name':raw,'modality':mod,
   'selection_method':str(r.experiment).rsplit('__',1)[-1],
   'importance_score':impmap.get(encoded,st.mean_importance_when_top30 if st else np.nan),
   'importance_score_source':'Phase 2 model importance if ranked; otherwise repeated-fold mean importance when in top 30' if encoded in impmap or st else 'not available in saved artifacts',
   'folds_selected_by_phase2_selector':np.nan,
   'folds_top30_importance_phase3':int(st.folds_top30) if st else 0,
   'repeated_fold_top30_frequency':float(st.top30_frequency) if st else 0.0,
   'stable_status':'stable (>=60% of 25 folds)' if st and st.top30_frequency>=.60 else ('intermittently top-30' if st and st.folds_top30>0 else 'not observed in repeated top-30'),
   'previous_selection_experiment':str(r.experiment)})
 # Include stable repeated-CV importance features not in phase-2 selected list.
 seen={(r['raw_feature_name'],r['selection_method']) for r in finalf}
 for r in stability.itertuples():
  key=(str(r.feature),'repeated_cv_tree_top30')
  if key in seen:continue
  finalf.append({'feature_name':str(r.feature),'raw_feature_name':str(r.feature),'modality':str(r.modality),
   'selection_method':'repeated_cv_tree_top30','importance_score':float(r.mean_importance_when_top30),
   'importance_score_source':'mean tree importance conditional on top-30 inclusion',
   'folds_selected_by_phase2_selector':np.nan,'folds_top30_importance_phase3':int(r.folds_top30),
   'repeated_fold_top30_frequency':float(r.top30_frequency),
   'stable_status':'stable (>=60% of 25 folds)' if r.top30_frequency>=.60 else 'intermittently top-30',
   'previous_selection_experiment':'not in Phase 2 selected-feature list'})
 final_features=pd.DataFrame(finalf)
 final_features.to_csv(OUT/'final_selected_features.csv',index=False)

 # Feature manifest and exact modality feature lists.
 manifest=pd.read_csv(BASE/'feature_summary.csv');manifest=manifest[manifest.dataset_file=='data.csv'].copy()
 manifest=manifest[manifest.column_name!='vital.status']
 manifest['modality']=manifest.column_name.map(lambda c: next((m for m,p in [('CNV','cn_'),('RNA','rs_'),('Mutation','mu_'),('Protein','pp_')] if c.startswith(p)),'unclassified'))
 manifest[['column_name','modality','inferred_type','missing_count','constant']].to_csv(OUT/'reproducibility'/'feature_manifest.csv',index=False)
 feature_lists={m:[c for c in features if c.startswith(p)] for m,p in [('CNV','cn_'),('RNA','rs_'),('Mutation','mu_'),('Protein','pp_')]}
 (OUT/'reproducibility'/'feature_lists'/'modality_feature_lists.json').write_text(json.dumps(feature_lists,indent=2),encoding='utf-8')

 # Reproducibility metadata and original artifacts inventory.
 (OUT/'reproducibility'/'package_versions.json').write_text(json.dumps({'python':sys.version,'platform':platform.platform(),'packages':{n:importlib.metadata.version(n) for n in ['numpy','pandas','scikit-learn','scipy','matplotlib','joblib']},'seed':SEED},indent=2),encoding='utf-8')
 cfg=next((BASE/'models').glob('rna__random_forest.json'))
 split=json.loads(cfg.read_text(encoding='utf-8'))
 pd.DataFrame({'original_row_index_zero_based':np.arange(len(df)),'split':['train' if i in set(split['train_indices_zero_based']) else 'test' for i in range(len(df))],'target_vital_status':y}).to_csv(OUT/'reproducibility'/'train_test_indices.csv',index=False)
 (OUT/'reproducibility'/'validation_config.json').write_text(json.dumps({'split_seed':SEED,'test_fraction':.2,'stratified':True,'row_based':True,'cv_phase2':{'type':'StratifiedKFold','folds':5,'shuffle':True,'seed':SEED},'cv_phase3':{'type':'RepeatedStratifiedKFold','folds':5,'repeats':5,'seed':SEED},'threshold_selection':'mean repeated OOF training probabilities; maximize Youden J','bootstrap_replicates':3000,'permutation_replicates':5000,'class_weight':'balanced for fitted classifiers','note':'No patient IDs: patient-level independence not verifiable.'},indent=2),encoding='utf-8')
 configs=list((BASE/'models').glob('*.json'))
 configs+=list((BASE/'individual_omics').rglob('*.json'))+list((BASE/'integrated').rglob('*.json'))
 copied=set()
 for c in configs:
  if c.name in copied:continue
  copied.add(c.name);shutil.copy2(c,OUT/'reproducibility'/'model_configs'/c.name)
 # Archive canonical training/validation scripts for review; actual paths remain at project root.
 shutil.copy2(BASE/'run_experiments.py',OUT/'reproducibility'/'scripts'/'run_experiments_phase2.py')
 shutil.copy2(ADV/'run_advanced_validation.py',OUT/'reproducibility'/'scripts'/'run_advanced_validation_phase3.py')
 inventory=[]
 for p in list((BASE/'models').glob('*.joblib'))+list((BASE/'individual_omics').rglob('*.joblib'))+list((BASE/'integrated').rglob('*.joblib')):
  inventory.append({'model_path':str(p.resolve()),'file_name':p.name,'size_bytes':p.stat().st_size,'location_type':'central model artifact' if p.parent==BASE/'models' else 'scoped copy'})
 pd.DataFrame(inventory).drop_duplicates('model_path').to_csv(OUT/'reproducibility'/'model_inventory.csv',index=False)

 # Figures: distributions.
 counts=[len(feature_lists[k]) for k in ['CNV','RNA','Mutation','Protein']]
 fig,ax=plt.subplots(figsize=(7,4));bars=ax.bar(['CNV','RNA','Mutation','Protein'],counts,color=['#577590','#43aa8b','#f9c74f','#f9844a']);ax.set_ylabel('Number of molecular features');ax.set_title('Dataset C modality feature counts')
 for b,v in zip(bars,counts):ax.text(b.get_x()+b.get_width()/2,v+12,str(v),ha='center')
 savefig(fig,'01_modality_distribution')
 vc=y.value_counts().sort_index();fig,ax=plt.subplots(figsize=(6,4));bars=ax.bar(['Alive (0)','Deceased (1)'],[vc[0],vc[1]],color=['#4c78a8','#e45756']);ax.set_ylabel('Samples');ax.set_title('Recorded vital-status classes (n=705)')
 for b,v in zip(bars,[vc[0],vc[1]]):ax.text(b.get_x()+b.get_width()/2,v+8,f'{v} ({v/len(y):.1%})',ha='center')
 savefig(fig,'02_target_distribution')

 # Repeated CV comparisons.
 winners=modality[~modality.baseline].copy();winners=winners.set_index('scope').loc[['cnv','rna','mutation','protein','integrated']].reset_index()
 fig,ax=plt.subplots(figsize=(8,4.8));xx=np.arange(len(winners));means=winners.repeated_cv_roc_auc_mean.astype(float);sd=winners.repeated_cv_roc_auc_sd_repeat_means.astype(float)
 ax.bar(xx,means,yerr=sd,capsize=4,color=['#577590','#43aa8b','#f9c74f','#f9844a','#9b5de5']);ax.axhline(.5,color='gray',ls='--',label='chance ROC-AUC');ax.set_xticks(xx,['CNV','RNA','Mutation','Protein','Integrated']);ax.set_ylim(0,1);ax.set_ylabel('Repeated-CV ROC-AUC (mean ± SD of repeat means)');ax.set_title('Individual and integrated modality comparison');ax.legend();savefig(fig,'03_modality_comparison')
 fig,ax=plt.subplots(figsize=(6,4.5));labels=['RNA only','Integrated'];vals=[float(winners[winners.scope=='rna'].repeated_cv_roc_auc_mean.iloc[0]),float(winners[winners.scope=='integrated'].repeated_cv_roc_auc_mean.iloc[0])];errs=[float(winners[winners.scope=='rna'].repeated_cv_roc_auc_sd_repeat_means.iloc[0]),float(winners[winners.scope=='integrated'].repeated_cv_roc_auc_sd_repeat_means.iloc[0])]
 ax.bar(labels,vals,yerr=errs,capsize=5,color=['#43aa8b','#9b5de5']);ax.set_ylim(0,1);ax.set_ylabel('Repeated-CV ROC-AUC');ax.set_title('Integrated model versus strongest individual modality');savefig(fig,'04_integrated_vs_individual')

 # ROC/PR curves from frozen, saved test probabilities for one CV-selected model per scope.
 choices={'Dummy baseline':'dummy_majority','CNV':'cnv__random_forest','RNA':'rna__random_forest','Mutation':'mutation__logistic_regression','Protein':'protein__random_forest','Integrated':'integrated__extra_trees'}
 colors=['#555555','#577590','#43aa8b','#f9c74f','#f9844a','#9b5de5']
 fig,ax=plt.subplots(figsize=(7,6))
 for (label,key),color in zip(choices.items(),colors):
  p=pd.read_csv(ADV/'predictions'/'dummy_majority_test_predictions.csv') if key=='dummy_majority' else pred[key]
  fpr,tpr,_=roc_curve(p.true_vital_status,p.probability_deceased);from sklearn.metrics import roc_auc_score
  ax.plot(fpr,tpr,label=f'{label} (AUC={roc_auc_score(p.true_vital_status,p.probability_deceased):.3f})',color=color)
 ax.plot([0,1],[0,1],'k--',lw=1,label='chance');ax.set(xlabel='False-positive rate',ylabel='True-positive rate',title='ROC curves on the common exploratory holdout');ax.legend(loc='lower right');savefig(fig,'05_roc_curves')
 fig,ax=plt.subplots(figsize=(7,6))
 for (label,key),color in zip(choices.items(),colors):
  p=pd.read_csv(ADV/'predictions'/'dummy_majority_test_predictions.csv') if key=='dummy_majority' else pred[key]
  precision,recall,_=precision_recall_curve(p.true_vital_status,p.probability_deceased);score=np.mean(p.true_vital_status)
  ax.plot(recall,precision,label=f'{label} (AP={score:.3f})' if key=='dummy_majority' else f'{label} (AP={__import__("sklearn.metrics",fromlist=["average_precision_score"]).average_precision_score(p.true_vital_status,p.probability_deceased):.3f})',color=color)
 ax.axhline(float(np.mean(ytest)),color='gray',ls='--',lw=1,label='test prevalence');ax.set(xlabel='Recall',ylabel='Precision',title='Precision–recall curves on the common exploratory holdout');ax.legend(loc='best');savefig(fig,'06_precision_recall_curves')

 # Confusion matrices in explicit [0,1] order.
 fig,axs=plt.subplots(2,3,figsize=(12,7));
 for ax,(label,key) in zip(axs.flat,choices.items()):
  p=pd.read_csv(ADV/'predictions'/'dummy_majority_test_predictions.csv') if key=='dummy_majority' else pred[key]
  cm=confusion_matrix(p.true_vital_status,p.predicted_vital_status,labels=[0,1]);ax.imshow(cm,cmap='Blues');ax.set_title(label);ax.set_xticks([0,1],['Alive','Deceased']);ax.set_yticks([0,1],['Alive','Deceased']);ax.set_xlabel('Predicted');ax.set_ylabel('True')
  for (i,j),v in np.ndenumerate(cm):ax.text(j,i,str(v),ha='center',va='center',color='black')
 fig.suptitle('Confusion matrices at threshold 0.5 (class order [alive, deceased])');savefig(fig,'07_confusion_matrices')

 # Feature importance figures by stable repeated folds.
 for scope in ['rna','integrated']:
  sub=stability[stability.scope==scope].sort_values(['mean_importance_when_top30','top30_frequency'],ascending=False).head(20).iloc[::-1]
  fig,ax=plt.subplots(figsize=(9,7));ax.barh(sub.feature,sub.mean_importance_when_top30,color='#577590');ax.set_xlabel('Mean importance when in top 30');ax.set_title(f'Top repeated-fold importance — {scope.upper()}')
  savefig(fig,f'08_feature_importance_{scope}')
 top=stability[stability.scope.isin(['rna','integrated'])].sort_values('top30_frequency',ascending=False).groupby('scope').head(15)
 fig,ax=plt.subplots(figsize=(9,8));top=top.sort_values(['scope','top30_frequency']);labels=[f'{r.feature} ({r.scope})' for r in top.itertuples()];ax.barh(labels,top.top30_frequency,color=['#43aa8b' if r=='rna' else '#9b5de5' for r in top.scope]);ax.axvline(.6,color='black',ls='--',lw=1,label='stable threshold: 60%');ax.set_xlim(0,1);ax.set_xlabel('Fraction of 25 repeated-CV folds in top 30');ax.set_title('Feature-selection / importance stability');ax.legend();savefig(fig,'09_feature_stability')

 # Calibration curves: saved uncalibrated probabilities plus saved training-fitted calibrators.
 fig,axs=plt.subplots(1,2,figsize=(11,5))
 for ax,scope in zip(axs,['rna','integrated']):
  model='random_forest' if scope=='rna' else 'extra_trees';name=f'{scope}__{model}'
  raw=pred[name].probability_deceased.to_numpy();cal=joblib.load(ADV/'calibration_models'/f'{name}_sigmoid_calibrated.joblib');calprob=cal.predict_proba(Xtest)[:,1]
  for probv,label,marker in [(raw,'Uncalibrated','o'),(calprob,'Sigmoid (training CV)','s')]:
   frac,mean=calibration_curve(ytest,probv,n_bins=6,strategy='quantile');ax.plot(mean,frac,marker+'-',label=label)
  ax.plot([0,1],[0,1],'k--',label='ideal');ax.set(xlim=(0,1),ylim=(0,1),xlabel='Mean predicted probability',ylabel='Observed deceased fraction',title=scope.upper());ax.legend()
 fig.suptitle('Calibration on held-out samples (n=141; 19 deceased)');savefig(fig,'10_calibration_curves')

 captions='''# Final figure captions\n\n1. **Modality feature distribution.** Counts of the 1,936 prefixed feature columns in the primary matrix: CNV 860, RNA 604, mutation 249, protein 223.\n2. **Recorded vital-status distribution.** Alive code 0: 611; deceased code 1: 94. The target is status only, with no follow-up time.\n3. **Individual and integrated modality comparison.** Repeated 5×5 training CV ROC-AUC means and standard deviations across five repeat means for the prespecified architecture per scope. Error bars describe resampling variability, not independent-sample CIs.\n4. **Integrated versus strongest individual modality.** Repeated-CV ROC-AUC comparison for RNA Random Forest and integrated Extra Trees, selected by prior training-CV results.\n5. **ROC curves.** Frozen predictions on the common 141-row exploratory holdout for Dummy and the CV-selected model architecture in each modality. Test curves are descriptive.\n6. **Precision–recall curves.** Frozen holdout precision/recall curves; the dashed horizontal line marks deceased prevalence in the test set.\n7. **Confusion matrices.** Held-out confusion counts at threshold 0.5; rows/columns use class order alive (0), deceased (1).\n8. **Feature importance.** Top repeated-fold tree importance summaries for RNA and integrated models; importance is associative and model-dependent.\n9. **Feature stability.** Frequency with which variables ranked in the top 30 across 25 training folds; the dashed threshold marks the report’s descriptive 60% stability cutoff.\n10. **Calibration curves.** Uncalibrated probabilities from the saved holdout predictions and probabilities from sigmoid calibrators fitted using training-only cross-validation; curves use the same held-out cases (n=141, 19 deceased).\n\nEach figure is exported as a 300-dpi PNG and a vector PDF.\n'''
 (OUT/'figures'/'figure_captions.md').write_text(captions,encoding='utf-8')

 # Export exact trained-model paths without copying/retraining estimators.
 pd.DataFrame([{'model_file':str(p.resolve()),'size_bytes':p.stat().st_size} for p in sorted((BASE/'models').glob('*.joblib'))]).to_csv(OUT/'reproducibility'/'central_model_paths.csv',index=False)
 # Copy phase-2 / phase-3 evaluation scripts exactly as provenance artifacts.
 shutil.copy2(BASE/'run_experiments.py',OUT/'reproducibility'/'scripts'/'run_experiments_phase2.py')
 shutil.copy2(ADV/'run_advanced_validation.py',OUT/'reproducibility'/'scripts'/'run_advanced_validation_phase3.py')

 # Final report values derived from files above.
 winner=master[(master.cv_selected_winner_within_scope)&(master.scope=='rna')].iloc[0]
 inte=master[(master.cv_selected_winner_within_scope)&(master.scope=='integrated')].iloc[0]
 dummyrow=master[master.baseline].iloc[0]
 rcv=repeated[repeated.scope=='rna'].iloc[0];icv=repeated[repeated.scope=='integrated'].iloc[0]
 caldf=pd.read_csv(ADV/'calibration_results.csv')
 permdf=pd.read_csv(ADV/'reports'/'heldout_permutation_tests.csv')
 threshdf=pd.read_csv(ADV/'threshold_analysis.csv')
 bootdf=pd.read_csv(ADV/'bootstrap_confidence_intervals.csv')
 diff=json.loads((ADV/'reports'/'integration_minus_rna_paired_bootstrap.json').read_text(encoding='utf-8'))
 thr_rna=threshdf.query("scope=='rna' and threshold_label=='OOF_Youden_J'").iloc[0]
 thr_int=threshdf.query("scope=='integrated' and threshold_label=='OOF_Youden_J'").iloc[0]
 def boot(scope,metric):return bootdf[(bootdf.scope==scope)&(bootdf.metric==metric)].iloc[0]
 ref_lines=[
 '[1] Christensen NJ, et al. Identifying interactions in omics data for clinical biomarker discovery using symbolic regression. *Bioinformatics*. 2022;38(15):3749–3758. https://doi.org/10.1093/bioinformatics/btac405',
 '[2] Ciriello G, et al.; TCGA Research Network. Comprehensive Molecular Portraits of Invasive Lobular Breast Cancer. *Cell*. 2015;163(2):506–519. https://doi.org/10.1016/j.cell.2015.09.033',
 '[3] QLattice clinical-omics companion repository. https://github.com/abzu-ai/QLattice-clinical-omics',
 '[4] GDC TCGA metadata documentation. https://docs.cancergenomicscloud.org/docs/tcga-metadata',
 '[5] Dataset source identified in prior audit: Kaggle, BRCA Multi-Omics (TCGA). https://www.kaggle.com/datasets/samdemharter/brca-multiomics-tcga']
 report=f'''# Dataset C multi-omics research report: exploratory recorded vital-status classification

## Abstract

We evaluated machine-learning models on a 705-row, 1,936-feature breast-cancer multi-omics matrix with a binary recorded vital-status target (611 alive, 94 deceased). CNV, RNA, mutation, protein and integrated classifiers were compared using training-only stratified cross-validation and a fixed 20% stratified row holdout. The strongest prespecified individual modality by Phase 2 training CV was RNA Random Forest (5-fold CV ROC-AUC {winner.phase2_cv_roc_auc_mean:.3f} ± {winner.phase2_cv_roc_auc_sd:.3f}; test ROC-AUC {winner.test_roc_auc:.3f}; PR-AUC {winner.test_pr_auc_average_precision:.3f}). The integrated Extra Trees model selected by training CV had 5-fold CV ROC-AUC {inte.phase2_cv_roc_auc_mean:.3f} ± {inte.phase2_cv_roc_auc_sd:.3f} and test ROC-AUC {inte.test_roc_auc:.3f}. Repeated 5×5 CV yielded RNA ROC-AUC {rcv.roc_auc_mean:.3f} ± {rcv.roc_auc_sd_repeat_means:.3f} and integrated ROC-AUC {icv.roc_auc_mean:.3f} ± {icv.roc_auc_sd_repeat_means:.3f}. Integration did not demonstrate better performance. This is exploratory classification of recorded vital status, not survival-time prediction; there are no patient IDs or follow-up durations, and patient-level independence cannot be checked.

## Introduction

Multi-omics classifiers can combine complementary measurement families, but apparent gains can reflect sample size, class imbalance, feature redundancy, or selection effects. This report consolidates existing Dataset C validation and modeling outputs. It does not fit new models. The underlying exported files are wide numerical matrices distributed with a Kaggle/QLattice analysis; prior provenance review links the source work to TCGA-BRCA and Ciriello et al. [1–3]. The precise TCGA barcodes and data-release mapping are absent from the supplied CSVs.

## Research objective

Compare individual CNV, RNA, mutation, and protein feature blocks with an integrated model for predicting the provided `vital.status` code. This target is recorded vital status (0=alive, 1=deceased), not survival time, cancer-specific mortality, or a defined prediction-horizon outcome. No clinical diagnostic or individual survival-prediction claim is made.

## Dataset description

The primary `data.csv` contains 705 rows, 1,936 prefixed candidate features, and `vital.status`. Counts are CNV (`cn_`) 860, RNA (`rs_`) 604, mutation (`mu_`) 249, and protein (`pp_`) 223. The class distribution is 611 alive and 94 deceased. The separate subtype-augmented CSV has four additional fields and 389 missing cells; these fields were excluded from the main models because prediction-time availability is undocumented. The CSVs contain no patient IDs, sample barcodes, or follow-up-duration/censoring fields. There are 62 exact duplicate-value column groups, with 99 extra CNV aliases; equality here means identical observed vectors only.

## Data preprocessing

The original pipelines excluded `vital.status` from predictors. CNV states were one-hot encoded as categorical with unknown categories ignored. RNA/protein columns used median imputation and standard scaling; mutation features used median imputation and training-fold variance filtering. Preprocessing and selectors were inside scikit-learn pipelines, so each training fold fitted its own transforms. Sources were not modified. Model training artifacts and predictions remain in the parent `results_dataset_C/` directory.

## Experimental design

The original split used stratified 80/20 row sampling, seed {SEED}: 564 training rows (489 alive, 75 deceased) and 141 held-out rows (122 alive, 19 deceased). The same holdout was used for all saved modality/model predictions. Five-fold stratified CV was applied only to training rows. Phase 3 added five repeats of stratified five-fold CV, with architecture fixed from the original training-CV winner for each modality (CNV/RF, RNA/RF, mutation/Logistic Regression, protein/RF, integrated/Extra Trees). Since sample IDs are unavailable, all splits are row-based; patient-level independence cannot be independently verified.

## Individual modality modelling

Class-weighted Logistic Regression, Random Forest, and Extra Trees were compared separately for each feature block. The complete table—including all feature-selection and CNV deduplication evaluations—is `master_metrics_comparison.csv`. By training CV, RNA Random Forest was the strongest individual unselected baseline. All held-out comparisons are exploratory and are not external validation.

## Integrated multi-omics modelling

The same three classifier families were trained on the combined four blocks. The training-CV winner was Extra Trees. In the Phase 2 5-fold CV, integrated Extra Trees ROC-AUC was {inte.phase2_cv_roc_auc_mean:.3f} ± {inte.phase2_cv_roc_auc_sd:.3f}, compared with {winner.phase2_cv_roc_auc_mean:.3f} ± {winner.phase2_cv_roc_auc_sd:.3f} for RNA Random Forest. In repeated CV, integrated ROC-AUC was {icv.roc_auc_mean:.3f} ± {icv.roc_auc_sd_repeat_means:.3f}, compared with RNA {rcv.roc_auc_mean:.3f} ± {rcv.roc_auc_sd_repeat_means:.3f}. Paired held-out bootstrap difference in ROC-AUC (integrated minus RNA) was {diff['roc_auc']['mean_difference']:.3f}, 95% percentile interval [{diff['roc_auc']['ci95_low']:.3f}, {diff['roc_auc']['ci95_high']:.3f}]. Integration did not demonstrate improved performance.

## Feature selection

Mutual information, L1 Logistic Regression, and Extra-Trees model-based selectors were fit inside CV pipelines for the best individual and integrated baselines. `final_selected_features.csv` combines the saved Phase 2 selector output with the repeated-fold top-30 importance screen. Phase 2 did not retain selector support counts for each fold; those counts are marked unavailable rather than inferred. The separate `folds_top30_importance_phase3` field records in how many of 25 tree-importance folds a raw feature ranked in the top 30. A descriptive stability threshold of 60% (15/25 folds) is used; none of these features is a validated biomarker.

## Model evaluation

Metrics include accuracy, balanced accuracy, class-specific precision/recall, deceased sensitivity, alive specificity, F1, macro-F1, ROC-AUC, average precision (PR-AUC), and confusion counts. The majority Dummy baseline on test had accuracy {dummyrow.test_accuracy:.3f}, balanced accuracy {dummyrow.test_balanced_accuracy:.3f}, macro-F1 {dummyrow.test_macro_f1:.3f}, ROC-AUC {dummyrow.test_roc_auc:.3f}, PR-AUC {dummyrow.test_pr_auc_average_precision:.3f}, and deceased recall {dummyrow.test_deceased_recall_sensitivity:.3f}. The RNA CV-selected model had test accuracy {winner.test_accuracy:.3f}, balanced accuracy {winner.test_balanced_accuracy:.3f}, deceased precision {winner.test_precision_deceased:.3f}, deceased recall {winner.test_deceased_recall_sensitivity:.3f}, specificity {winner.test_specificity_alive:.3f}, F1 {winner.test_f1_deceased:.3f}, macro-F1 {winner.test_macro_f1:.3f}, ROC-AUC {winner.test_roc_auc:.3f}, PR-AUC {winner.test_pr_auc_average_precision:.3f}, and confusion matrix TN/FP/FN/TP {int(winner.tn)}/{int(winner.fp)}/{int(winner.fn)}/{int(winner.tp)}. It identified 6/19 deceased test cases at the default threshold.

## Advanced validation

OOF threshold selection used training-only repeated OOF probabilities and Youden J, not test labels. At thresholds {thr_rna.threshold:.3f} (RNA) and {thr_int.threshold:.3f} (integrated), held-out deceased sensitivity/specificity/F1 were {thr_rna.test_sensitivity:.3f}/{thr_rna.test_specificity:.3f}/{thr_rna.test_f1_deceased:.3f} and {thr_int.test_sensitivity:.3f}/{thr_int.test_specificity:.3f}/{thr_int.test_f1_deceased:.3f}, respectively. These thresholds did not improve holdout F1 and are not proposed as validated operating points. Training-only sigmoid calibration gave holdout Brier scores RNA {caldf.query("scope=='rna' and probability_type=='uncalibrated'").iloc[0].brier_score_test:.3f} raw vs {caldf.query("scope=='rna' and probability_type=='sigmoid_cv_calibrated'").iloc[0].brier_score_test:.3f} calibrated; integrated {caldf.query("scope=='integrated' and probability_type=='uncalibrated'").iloc[0].brier_score_test:.3f} raw vs {caldf.query("scope=='integrated' and probability_type=='sigmoid_cv_calibrated'").iloc[0].brier_score_test:.3f} calibrated. Stratified bootstrap intervals use 3,000 resamples of frozen holdout predictions. RNA holdout ROC-AUC was {boot('rna','roc_auc').point_estimate:.3f} (95% percentile interval {boot('rna','roc_auc').ci95_percentile_low:.3f}?{boot('rna','roc_auc').ci95_percentile_high:.3f}); integrated was {boot('integrated','roc_auc').point_estimate:.3f} ({boot('integrated','roc_auc').ci95_percentile_low:.3f}?{boot('integrated','roc_auc').ci95_percentile_high:.3f}). Fixed-prediction label permutation tests used 5,000 shuffles: p={permdf.query("scope=='rna'").iloc[0].one_sided_pvalue:.3f} for RNA and p={permdf.query("scope=='integrated'").iloc[0].one_sided_pvalue:.3f} for integrated; neither was statistically significant at 0.05. These tests are conditional on the same small holdout and do not establish generalization.

## Explainability

Tree impurity importance, prior feature-selector outputs, and repeated-fold stability rankings are descriptive. The most recurrent RNA top-30 importance features included {', '.join(stability[(stability.scope=='rna')].sort_values('top30_frequency',ascending=False).head(5).feature)}. Integrated top-30 rankings were dominated by RNA features, with contributions from other modalities in some folds. Feature importance is not causal evidence and does not establish clinical biomarkers. SHAP was not installed; QLattice/Feyn was unavailable.

## Results

See `master_metrics_comparison.csv`, `final_selected_features.csv`, and the figure set under `figures/`. The master table preserves Phase 2 5-fold CV and Phase 3 repeated-CV statistics in separate columns. Detailed advanced metric files are copied under `tables/`. For repeated CV, the reported standard deviation and empirical range summarize five repeat means; they are not formal independent-sample confidence intervals because folds/repeats overlap. The holdout contains only 19 deceased cases, producing substantial uncertainty.

See `figures/figure_captions.md` for figure descriptions and interpretation notes. Figures are high-resolution PNG and vector PDF exports.

## Discussion

RNA alone had higher training CV and held-out point estimates than the CV-selected integrated model in these runs. The paired bootstrap interval includes zero, so the small holdout does not establish a reliable modality difference. The majority baseline illustrates why accuracy alone is inadequate under the 611:94 class imbalance. Threshold changes trade specificity for sensitivity and did not yield a clear F1 gain on holdout. Calibration scores improved on this holdout after sigmoid calibration, but this is not sufficient evidence of stable calibration. A single row-based holdout and repeated CV cannot resolve unknown patient overlap.

## Limitations

There are only 705 observations and 94 deceased samples total, including 19 in the test set. There are no patient IDs/barcodes, no survival follow-up durations, no censoring indicators, and no independently available sample manifest. Patient-level independence cannot be verified. Row-wise random splitting may therefore overstate effective sample independence. Target timing and prediction horizon are not defined. The model is not clinically validated, cannot diagnose breast cancer, and must not be presented as predicting an individual's survival. The findings are exploratory recorded-vital-status classification, not survival-time prediction.

## Conclusion

Within this dataset export, the training-CV-selected RNA model outperformed the integrated model on point estimates, and integration did not demonstrate added predictive value. Permutation tests were not statistically significant. These results are hypothesis-generating only and do not establish clinical usefulness, diagnostic performance, or generalization to new patients.

## Reproducibility details

Python {sys.version.split()[0]}; scikit-learn {sklearn.__version__}; pandas {pd.__version__}; NumPy {np.__version__}; matplotlib {matplotlib.__version__}; joblib {joblib.__version__}. Seed {SEED}; stratified 80/20 train/test; Phase 2 5-fold stratified CV; Phase 3 5 repeats × 5 folds; class weights balanced; ROC-AUC/AP computed from deceased-class probabilities. See `reproducibility/` for package versions, exact split indices, modality feature manifest/lists, model configs, and saved script copies. Canonical model files remain at `../models/` (and scoped copies under `../individual_omics/` and `../integrated/`); exact paths are listed in `reproducibility/model_inventory.csv` and `central_model_paths.csv`. Recreate this report/package with `python results_dataset_C/final_research/build_final_package.py` from the workspace root. This assembles outputs only and does not refit models.

## References

{chr(10).join(ref_lines)}
'''
 (OUT/'dataset_C_research_report.md').write_text(report,encoding='utf-8')
 # Final integrity manifest for expected package deliverables.
 expected=['dataset_C_research_report.md','master_metrics_comparison.csv','final_selected_features.csv','reproducibility/README.md']
 # README written below by this assembler.
 readme=f'''# Dataset C final research package\n\nThis directory consolidates existing Dataset C results. The assembler does not train or refit any model.\n\n## Main files\n\n- `dataset_C_research_report.md`: final scientific report.\n- `master_metrics_comparison.csv`: all 23 Phase 2 evaluation rows plus DummyClassifier, with Phase 2 five-fold CV and Phase 3 repeated-CV values in separate fields.\n- `final_selected_features.csv`: Phase 2 selected features and Phase 3 repeated-fold top-30 importance/stability features. Phase 2 per-fold selector support counts were not retained and are blank; do not interpret top-30 importance fold counts as selector support.\n- `figures/`: ten publication-style figure families, each exported as 300-dpi PNG and vector PDF; captions are in igures/figure_captions.md.\n- `reproducibility/`: package versions, seeds/split indices, feature manifest, model configs and model inventory.\n\n## Reproduction\n\nUse the original inputs (read only):\n\n- `{DATA}`\n- `{ENRICHED}` (used only for separate subtype-label audit; not the primary classifier)\n\nThe original Phase 2 training/evaluation script is `{BASE/'run_experiments.py'}`. The Phase 3 validation script is `{ADV/'run_advanced_validation.py'}`. Exact copies are included in `reproducibility/scripts/` for reference. Existing package versions are recorded in `reproducibility/package_versions.json`; split and CV configuration is in `reproducibility/validation_config.json`. Phase 2 seed: {SEED}; stratified 80/20 row split; 5-fold stratified training CV. Phase 3: 5 repeats × 5 folds on training data.\n\nRun the package assembler from the workspace root with `python results_dataset_C/final_research/build_final_package.py`. This rebuilds tables, report, metadata, and figures from the existing results; it does not refit models. To rerun model experiments, use the canonical scripts above and the recorded environment. Do not use a different split if comparing directly with these holdout metrics.\n\n## Trained model locations\n\nNo estimators are copied into this package. Phase 2 model artifacts remain under `{BASE/'models'}`; modality-scoped copies are in `{BASE/'individual_omics'}` and `{BASE/'integrated'}`. Exact absolute paths are in `reproducibility/model_inventory.csv` and `reproducibility/central_model_paths.csv`. Training-partition calibrated variants from Phase 3 remain in `{ADV/'calibration_models'}`.\n\n## Interpretation\n\nThe target is recorded vital status, not survival time. Patient IDs and follow-up durations are absent, so patient independence and time-to-event outcomes cannot be checked. Results are exploratory, non-clinical, and not validated for diagnosis or predicting an individual's survival.\n'''
 (OUT/'reproducibility'/'README.md').write_text(readme,encoding='utf-8')
 print('package created',OUT)
 print('master rows',len(master),'selected feature rows',len(final_features),'figures',len(list((OUT/'figures').glob('*.png'))),'pdfs',len(list((OUT/'figures').glob('*.pdf'))))

if __name__=='__main__':main()
