"""Create post-fit summaries/figures from saved Dataset C optimization outputs.

No model choices or thresholds are made from holdout predictions here.
"""
from pathlib import Path
import json
import numpy as np
import pandas as pd
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from sklearn.calibration import calibration_curve
from sklearn.metrics import (roc_curve, roc_auc_score, precision_recall_curve,
                             average_precision_score, confusion_matrix,
                             brier_score_loss)
from sklearn.model_selection import StratifiedKFold
import joblib

ROOT=Path(__file__).resolve().parents[2]
BASE=ROOT/'results_dataset_C'; OUT=BASE/'final_optimization'
DATA=ROOT/'dataset_c_extracted'/'data.csv'
SEED=20261110
PREF={'cnv':'cn_','rna':'rs_','mutation':'mu_','protein':'pp_'}

def main():
    figdir=OUT/'figures';tabs=OUT/'tables';pred_dir=OUT/'predictions'
    data=pd.read_csv(DATA); y=data['vital.status'].astype(int)
    tr=pd.read_csv(OUT/'reproducibility'/'train_test_indices.csv')
    train_idx=tr.loc[tr.split=='train','original_row_index_zero_based'].to_numpy(int)
    test_idx=tr.loc[tr.split=='test','original_row_index_zero_based'].to_numpy(int)
    test=pd.read_csv(pred_dir/'selected_scope_test_predictions.csv')
    dummy=pd.read_csv(pred_dir/'dummy_baseline_test_predictions.csv')
    # Probability curves and calibration are descriptive holdout evaluations only.
    fig,ax=plt.subplots(figsize=(7,6))
    for name,g in test.groupby('model'):
        ax.step(*precision_recall_curve(g.true_vital_status,g.probability_deceased)[1::-1],where='post',label=f"{name}: AP={average_precision_score(g.true_vital_status,g.probability_deceased):.3f}")
    dg=dummy[dummy.scope=='rna']
    pr,rc,_=precision_recall_curve(dg.true_vital_status,dg.probability_deceased)
    ax.step(rc,pr,where='post',linestyle='--',color='black',label=f"Dummy: AP={average_precision_score(dg.true_vital_status,dg.probability_deceased):.3f}")
    ax.set(xlabel='Recall (deceased)',ylabel='Precision (deceased)',title='Precision-recall curves on the original holdout; descriptive only',xlim=(0,1),ylim=(0,1));ax.legend(fontsize=7);fig.tight_layout();fig.savefig(figdir/'optimized_precision_recall_curves.png',dpi=300);fig.savefig(figdir/'optimized_precision_recall_curves.pdf');plt.close(fig)
    rows=[];fig,ax=plt.subplots(figsize=(7,6))
    for name,g in test.groupby('model'):
        frac,mean=calibration_curve(g.true_vital_status,g.probability_deceased,n_bins=6,strategy='quantile')
        ax.plot(mean,frac,marker='o',label=name)
        rows.append({'model':name,'brier_score':brier_score_loss(g.true_vital_status,g.probability_deceased),'test_rows':len(g),'calibration_bins':len(frac),'calibration_method':'none; evaluation only'})
    gd=dummy[dummy.scope=='rna'];frac,mean=calibration_curve(gd.true_vital_status,gd.probability_deceased,n_bins=6,strategy='quantile');ax.plot(mean,frac,'k--',marker='o',label='Dummy')
    rows.append({'model':'dummy_prior','brier_score':brier_score_loss(gd.true_vital_status,gd.probability_deceased),'test_rows':len(gd),'calibration_bins':len(frac),'calibration_method':'none; evaluation only'})
    ax.plot([0,1],[0,1],':',color='gray');ax.set(xlabel='Mean predicted probability',ylabel='Observed deceased fraction',title='Holdout calibration curves (19 deceased cases)');ax.legend(fontsize=8);fig.tight_layout();fig.savefig(figdir/'optimized_calibration_curves.png',dpi=300);fig.savefig(figdir/'optimized_calibration_curves.pdf');plt.close(fig)
    pd.DataFrame(rows).to_csv(tabs/'optimized_calibration_metrics.csv',index=False)
    # CV vs descriptive holdout model comparison.
    outer=pd.read_csv(tabs/'nested_cv_summary.csv'); cand=pd.read_csv(tabs/'model_selection_summary.csv'); testm=pd.read_csv(tabs/'selected_scope_models_test_metrics.csv')
    selected=cand[cand.scope_selected][['scope','model','inner_roc_auc_mean','inner_pr_auc_mean','inner_balanced_accuracy_mean','inner_macro_f1_mean','inner_deceased_recall_mean','mean_composite_rank']]
    comparison=outer.merge(selected,on='scope',how='left').merge(testm,on=['scope','model'],suffixes=('_nested','_test'),how='left')
    comparison.to_csv(tabs/'model_comparison_master.csv',index=False)
    # include dummy rows as an explicit baseline comparison
    baseline=pd.read_csv(tabs/'dummy_baseline_comparison.csv')
    b_rows=[]
    for _,r in baseline.iterrows():
        b_rows.append({'scope':r.scope,'model':'dummy_prior','roc_auc_mean':r.cv_roc_auc,'pr_auc_mean':r.cv_pr_auc,'balanced_accuracy_mean':r.cv_balanced_accuracy,'macro_f1_mean':r.cv_macro_f1,'deceased_recall_mean':r.cv_deceased_recall,'roc_auc':r.test_roc_auc,'pr_auc':r.test_pr_auc,'accuracy':r.test_accuracy,'balanced_accuracy':r.test_balanced_accuracy,'precision_deceased':r.test_precision_deceased,'deceased_recall':r.test_deceased_recall,'specificity':r.test_specificity,'f1_deceased':r.test_f1_deceased,'macro_f1':r.test_macro_f1,'tn':r.test_tn,'fp':r.test_fp,'fn':r.test_fn,'tp':r.test_tp})
    comparison.rename(columns={'roc_auc_mean':'roc_auc_mean','pr_auc_mean':'pr_auc_mean','balanced_accuracy_mean':'balanced_accuracy_mean','macro_f1_mean':'macro_f1_mean','deceased_recall_mean':'deceased_recall_mean','roc_auc':'roc_auc','pr_auc':'pr_auc'},inplace=True)
    pd.concat([comparison,pd.DataFrame(b_rows)],ignore_index=True,sort=False).to_csv(tabs/'model_comparison_with_dummy.csv',index=False)
    fig,ax=plt.subplots(figsize=(9,5));
    for _,r in comparison.iterrows(): ax.scatter(r.roc_auc_mean,r.pr_auc_mean,s=60,label=f"{r.scope} ({r.model})")
    ax.scatter(.5,611/705,marker='x',s=90,color='black',label='Dummy prior (training prevalence)')
    ax.set(xlabel='Nested outer CV ROC-AUC',ylabel='Nested outer CV PR-AUC',title='CV-selected model comparison; training rows only');ax.legend(fontsize=7);fig.tight_layout();fig.savefig(figdir/'model_comparison_nested_cv.png',dpi=300);fig.savefig(figdir/'model_comparison_nested_cv.pdf');plt.close(fig)
    # Training-fold permutation importance. Only the top 30 variables by the
    # corresponding outer-training model's native importance are permuted.
    imp_source=pd.read_csv(tabs/'nested_inner_and_importance.csv')
    folds=pd.read_csv(tabs/'nested_outer_results.csv')
    imp_rows=[]
    for scope in ['rna','integrated']:
        cols=[c for c in data.columns if c!='vital.status' and (scope=='integrated' or c.startswith(PREF[scope]))]
        fold=StratifiedKFold(n_splits=3,shuffle=True,random_state=SEED)
        for fi,(itr,iva) in enumerate(fold.split(data.iloc[train_idx],y.iloc[train_idx]),1):
            outer=folds[(folds.scope==scope)&(folds.outer_fold==fi)].iloc[0]
            model=outer.selected_algorithm_by_inner_cv
            if model not in ['random_forest','extra_trees']:
                continue
            base=imp_source[(imp_source.scope==scope)&(imp_source.outer_fold==fi)&(imp_source.status=='outer_train_importance')].sort_values('importance_rank').head(30)
            selected=base.feature.tolist()
            if not selected: continue
            # Reuse saved final pipeline configuration, rebuilding the estimator
            # identical to the selection script to preserve fold separation.
            from run_optimization import make_pipe
            pipe=make_pipe(scope,model,list(data.drop(columns='vital.status').columns));pipe.fit(data.iloc[train_idx[itr]],y.iloc[train_idx[itr]])
            Xv=data.iloc[train_idx[iva]][cols].copy();yv=y.iloc[train_idx[iva]].to_numpy();base=roc_auc_score(yv,pipe.predict_proba(Xv)[:,1]);rng=np.random.default_rng(SEED+fi)
            for feat in selected:
                drops=[]
                for _ in range(3):
                    xp=Xv.copy();xp[feat]=rng.permutation(xp[feat].to_numpy());drops.append(base-roc_auc_score(yv,pipe.predict_proba(xp)[:,1]))
                imp_rows.append({'scope':scope,'model':model,'outer_fold':fi,'feature':feat,'modality':next((m for m,p in PREF.items() if feat.startswith(p)),'unknown'),'permutation_roc_auc_decrease_mean':float(np.mean(drops)),'within_fold_permutation_sd':float(np.std(drops,ddof=1))})
    pi=pd.DataFrame(imp_rows);pi.to_csv(OUT/'feature_selection'/'outer_fold_permutation_importance.csv',index=False)
    if len(pi):
        stab=pi.groupby(['scope','feature','modality']).agg(outer_folds_permuted=('outer_fold','nunique'),mean_auc_decrease=('permutation_roc_auc_decrease_mean','mean'),sd_auc_decrease=('permutation_roc_auc_decrease_mean','std')).reset_index();stab['fold_frequency']=stab.outer_folds_permuted/3;stab.to_csv(OUT/'feature_selection'/'permutation_importance_summary.csv',index=False)
        bestscope=json.loads((OUT/'reports'/'final_model_selection.json').read_text(encoding='utf-8'))['scope']
        show=stab[stab.scope==bestscope].sort_values('mean_auc_decrease').tail(20);fig,ax=plt.subplots(figsize=(8,7));ax.barh(show.feature,show.mean_auc_decrease,color='#f8961e');ax.set(xlabel='Outer-validation ROC-AUC decrease after permutation',title=f'{bestscope} top-20 training-fold permutation importance');fig.tight_layout();fig.savefig(figdir/'permutation_importance.png',dpi=300);fig.savefig(figdir/'permutation_importance.pdf');plt.close(fig)
    # Feature selection stability plot from fold-contained selectors.
    freq=pd.read_csv(OUT/'feature_selection'/'feature_selection_frequency.csv')
    for sc in ['rna','integrated']:
        sub=freq[(freq.scope==sc)&(freq.selection_method.isin(['mutual_information','l1','tree_importance']))].sort_values('selection_frequency').tail(20)
        if len(sub):
            fig,ax=plt.subplots(figsize=(9,7));ax.barh(sub.feature,sub.selection_frequency,color='#277da1');ax.set(xlabel='Fraction of five training folds selected',xlim=(0,1),title=f'{sc} feature-selection stability (top frequencies)');fig.tight_layout();fig.savefig(figdir/f'{sc}_feature_selection_stability.png',dpi=300);fig.savefig(figdir/f'{sc}_feature_selection_stability.pdf');plt.close(fig)
    print('post-fit figures, combined comparison, calibration, and fold-validation permutation importance complete')

if __name__=='__main__':main()
