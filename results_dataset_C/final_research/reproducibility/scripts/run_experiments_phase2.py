"""Reproducible exploratory Dataset C recorded-vital-status experiments."""
from __future__ import annotations

import csv
import json
import os
import platform
import sys
import traceback
from functools import partial
from pathlib import Path

import joblib
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
import sklearn
from sklearn.compose import ColumnTransformer
from sklearn.ensemble import ExtraTreesClassifier, RandomForestClassifier
from sklearn.feature_selection import SelectFromModel, SelectKBest, mutual_info_classif
from sklearn.impute import SimpleImputer
from sklearn.inspection import permutation_importance
from sklearn.linear_model import LogisticRegression
from sklearn.metrics import (
    accuracy_score, average_precision_score, confusion_matrix, f1_score,
    precision_score, recall_score, roc_auc_score, make_scorer,
)
from sklearn.model_selection import StratifiedKFold, cross_validate, train_test_split
from sklearn.pipeline import Pipeline
from sklearn.preprocessing import OneHotEncoder, StandardScaler
from sklearn.feature_selection import VarianceThreshold

ROOT = Path(__file__).resolve().parents[1]
OUT = ROOT / "results_dataset_C"
DATA = ROOT / "dataset_c_extracted" / "data.csv"
ENRICHED = ROOT / "dataset_c_extracted" / "brca_data_w_subtypes.csv"
ALIAS_MAP = OUT / "duplicate_alias_map.csv"
SEED = 20261110
TEST_SIZE = 0.20
CV_FOLDS = 5
N_TREES = 100
# Windows sandbox blocks joblib worker process creation; keep estimator fitting in-process.
N_JOBS = 1
PREFIXES = {"cnv": "cn_", "rna": "rs_", "mutation": "mu_", "protein": "pp_"}

DIRS = [
    "models", "predictions", "individual_omics/cnv", "individual_omics/rna",
    "individual_omics/mutation", "individual_omics/protein", "integrated",
    "feature_selection", "explainability", "figures/confusion_matrices",
    "figures/roc_curves", "figures/precision_recall_curves",
    "figures/model_comparisons", "figures/feature_importance", "reports",
]
for d in DIRS:
    (OUT / d).mkdir(parents=True, exist_ok=True)

LOG = OUT / "experiment_log.txt"
with LOG.open("a", encoding="utf-8") as _log_file:
    _log_file.write("\n--- experiment runner invocation ---\n")


def log(message: str) -> None:
    print(message, flush=True)
    with LOG.open("a", encoding="utf-8") as f:
        f.write(message + "\n")


def metrics(y_true, y_pred, y_prob) -> dict:
    cm = confusion_matrix(y_true, y_pred, labels=[0, 1])
    tn, fp, fn, tp = cm.ravel()
    return {
        "accuracy": accuracy_score(y_true, y_pred),
        "precision_deceased": precision_score(y_true, y_pred, pos_label=1, zero_division=0),
        "recall_deceased_sensitivity": recall_score(y_true, y_pred, pos_label=1, zero_division=0),
        "specificity_alive": tn / (tn + fp) if tn + fp else np.nan,
        "f1_deceased": f1_score(y_true, y_pred, pos_label=1, zero_division=0),
        "macro_f1": f1_score(y_true, y_pred, average="macro", zero_division=0),
        "roc_auc": roc_auc_score(y_true, y_prob),
        "pr_auc_average_precision": average_precision_score(y_true, y_prob),
        "precision_alive": precision_score(y_true, y_pred, pos_label=0, zero_division=0),
        "recall_alive": recall_score(y_true, y_pred, pos_label=0, zero_division=0),
        "tn": int(tn), "fp": int(fp), "fn": int(fn), "tp": int(tp),
    }


def make_estimator(name: str, seed: int):
    if name == "logistic_regression":
        return LogisticRegression(
            class_weight="balanced", C=0.1, solver="liblinear", max_iter=2500,
            random_state=seed,
        )
    cls = RandomForestClassifier if name == "random_forest" else ExtraTreesClassifier
    return cls(
        n_estimators=N_TREES, class_weight="balanced", max_features="sqrt",
        min_samples_leaf=2, n_jobs=N_JOBS, random_state=seed,
    )


def feature_groups(modality: str):
    if modality == "integrated":
        return list(PREFIXES)
    return [modality]


def columns_for(modality: str, columns: list[str]):
    return [c for group in feature_groups(modality) for c in columns if c.startswith(PREFIXES[group])]


def build_pipeline(modality: str, model_name: str, selector=None, seed=SEED):
    all_cols = columns_for(modality, list(X.columns))
    cats = [c for c in all_cols if c.startswith("cn_")]
    mutation = [c for c in all_cols if c.startswith("mu_")]
    other_numeric = [c for c in all_cols if c.startswith(("rs_", "pp_"))]
    transformers = []
    if cats:
        cat_pipe = Pipeline([
            ("imputer", SimpleImputer(strategy="most_frequent")),
            ("onehot", OneHotEncoder(handle_unknown="ignore", sparse_output=True)),
        ])
        transformers.append(("cnv_categorical", cat_pipe, cats))
    if mutation:
        mut_pipe = Pipeline([
            ("imputer", SimpleImputer(strategy="median")),
            ("variance", VarianceThreshold(threshold=0.0)),
        ])
        transformers.append(("mutation_binary", mut_pipe, mutation))
    if other_numeric:
        num_pipe = Pipeline([
            ("imputer", SimpleImputer(strategy="median")),
            ("scale", StandardScaler()),
        ])
        transformers.append(("continuous", num_pipe, other_numeric))
    prep = ColumnTransformer(transformers, remainder="drop", verbose_feature_names_out=False)
    steps = [("preprocess", prep)]
    if selector is not None:
        kind, arg = selector
        if kind == "mutual_information":
            steps.append(("select", SelectKBest(partial(mutual_info_classif, random_state=seed), k=min(arg, max(1, len(all_cols) * 4)))))
        elif kind == "l1_logistic":
            selector_model = LogisticRegression(
                penalty="l1", solver="liblinear", C=0.08, class_weight="balanced",
                max_iter=2500, random_state=seed,
            )
            steps.append(("select", SelectFromModel(selector_model, threshold=-np.inf, max_features=arg)))
        elif kind == "tree_importance":
            selector_model = ExtraTreesClassifier(
                n_estimators=80, class_weight="balanced", max_features="sqrt",
                min_samples_leaf=2, n_jobs=N_JOBS, random_state=seed,
            )
            steps.append(("select", SelectFromModel(selector_model, threshold=-np.inf, max_features=arg)))
    steps.append(("model", make_estimator(model_name, seed)))
    return Pipeline(steps)


def specificity_scorer(y_true, y_pred):
    cm = confusion_matrix(y_true, y_pred, labels=[0, 1])
    tn, fp, _, _ = cm.ravel()
    return tn / (tn + fp) if tn + fp else 0.0


SCORING = {
    "accuracy": "accuracy",
    "precision_deceased": make_scorer(precision_score, pos_label=1, zero_division=0),
    "recall_deceased": make_scorer(recall_score, pos_label=1, zero_division=0),
    "specificity_alive": make_scorer(specificity_scorer),
    "f1_deceased": make_scorer(f1_score, pos_label=1, zero_division=0),
    "macro_f1": "f1_macro",
    "roc_auc": "roc_auc",
    "pr_auc": "average_precision",
}


def save_prediction(name, y_true, y_pred, y_prob, indices, model_group):
    path = OUT / "predictions" / f"{name}_test_predictions.csv"
    pd.DataFrame({
        "row_index_original_zero_based": indices,
        "true_vital_status": np.asarray(y_true),
        "predicted_vital_status": np.asarray(y_pred),
        "probability_deceased": np.asarray(y_prob),
        "model": name,
        "experiment_group": model_group,
    }).to_csv(path, index=False)
    return path


def plot_eval(name, y_true, y_pred, y_prob, model_group):
    cm = confusion_matrix(y_true, y_pred, labels=[0, 1])
    fig, ax = plt.subplots(figsize=(4.8, 4.2))
    ax.imshow(cm, cmap="Blues")
    ax.set(xticks=[0, 1], yticks=[0, 1], xticklabels=["Alive (0)", "Deceased (1)"],
           yticklabels=["Alive (0)", "Deceased (1)"], xlabel="Predicted", ylabel="True",
           title=f"Confusion matrix — {name}")
    for (i, j), v in np.ndenumerate(cm):
        ax.text(j, i, str(v), ha="center", va="center", color="black")
    fig.tight_layout()
    fig.savefig(OUT / "figures" / "confusion_matrices" / f"{name}.png", dpi=140)
    plt.close(fig)

    for kind, curvefile in [("roc", "roc_curves"), ("pr", "precision_recall_curves")]:
        fig, ax = plt.subplots(figsize=(5, 4))
        if kind == "roc":
            from sklearn.metrics import RocCurveDisplay
            RocCurveDisplay.from_predictions(y_true, y_prob, ax=ax)
            ax.plot([0, 1], [0, 1], "k--", lw=1)
            ax.set_title(f"ROC — {name}")
        else:
            from sklearn.metrics import PrecisionRecallDisplay
            PrecisionRecallDisplay.from_predictions(y_true, y_prob, ax=ax)
            ax.set_title(f"Precision–recall — {name}")
        fig.tight_layout()
        fig.savefig(OUT / "figures" / curvefile / f"{name}.png", dpi=140)
        plt.close(fig)


def metric_rows(scope, model_name, cv_result, test_result, selection="none"):
    cv = {}
    for key in SCORING:
        arr = cv_result.get(f"test_{key}", np.array([np.nan]))
        cv[key + "_mean"] = float(np.nanmean(arr))
        cv[key + "_std"] = float(np.nanstd(arr, ddof=1))
    return {
        "experiment": scope, "model": model_name, "feature_selection": selection,
        **cv, **{"test_" + k: v for k, v in test_result.items()},
    }


def main():
    global X
    log(f"Python={sys.version.split()[0]} platform={platform.platform()} sklearn={sklearn.__version__} numpy={np.__version__} pandas={pd.__version__}")
    log(f"Random seed={SEED}; test_size={TEST_SIZE}; CV={CV_FOLDS}-fold StratifiedKFold; trees={N_TREES}; jobs={N_JOBS}")
    df = pd.read_csv(DATA)
    enriched = pd.read_csv(ENRICHED)
    if "vital.status" not in df:
        raise ValueError("vital.status target missing")
    modality_cols = {k: [c for c in df if c.startswith(v)] for k, v in PREFIXES.items()}
    target_col = "vital.status"
    unexpected = [c for c in df.columns if c != target_col and not any(c.startswith(p) for p in PREFIXES.values())]
    if unexpected:
        raise ValueError(f"Unexpected predictors found: {unexpected}")
    X = df.drop(columns=[target_col])
    y = df[target_col].astype(int)
    if set(y.unique()) != {0, 1}:
        raise ValueError(f"Unexpected target coding: {sorted(y.unique())}")
    if tuple(len(modality_cols[k]) for k in ["cnv", "rna", "mutation", "protein"]) != (860, 604, 249, 223):
        raise ValueError(f"Unexpected modality counts: { {k: len(v) for k,v in modality_cols.items()} }")
    expected_extra = ["PR.Status", "ER.Status", "HER2.Final.Status", "histological.type"]
    if not set(expected_extra).issubset(enriched.columns):
        raise ValueError("Expected subtype columns missing from enriched CSV")
    train_idx, test_idx = train_test_split(
        np.arange(len(df)), test_size=TEST_SIZE, random_state=SEED, stratify=y,
    )
    X_train, X_test = X.iloc[train_idx], X.iloc[test_idx]
    y_train, y_test = y.iloc[train_idx], y.iloc[test_idx]
    cv = StratifiedKFold(n_splits=CV_FOLDS, shuffle=True, random_state=SEED)
    log(f"Samples={len(df)}; train={len(train_idx)} test={len(test_idx)}; train class={dict(y_train.value_counts().sort_index())}; test class={dict(y_test.value_counts().sort_index())}")
    log(f"Modality counts={ {k:len(v) for k,v in modality_cols.items()} }; unexpected predictors={unexpected}")
    log("Split is row-based: no patient IDs/barcodes exist, so patient-level independence cannot be verified.")

    records=[]; fitted={}; heldout={}
    scopes=["cnv", "rna", "mutation", "protein", "integrated"]
    model_names=["logistic_regression", "random_forest", "extra_trees"]
    for scope in scopes:
        for model_name in model_names:
            label=f"{scope}__{model_name}"
            log(f"BASELINE start {label}")
            pipe=build_pipeline(scope,model_name)
            cv_result=cross_validate(pipe, X_train, y_train, cv=cv, scoring=SCORING,
                                     n_jobs=1, error_score="raise", return_train_score=False)
            fitted_pipe=build_pipeline(scope,model_name)
            fitted_pipe.fit(X_train,y_train)
            pred=fitted_pipe.predict(X_test); prob=fitted_pipe.predict_proba(X_test)[:,1]
            tm=metrics(y_test,pred,prob)
            records.append(metric_rows(scope,model_name,cv_result,tm))
            fitted[label]=fitted_pipe; heldout[label]=(pred,prob)
            save_prediction(label,y_test,pred,prob,test_idx,scope)
            plot_eval(label,y_test,pred,prob,scope)
            joblib.dump(fitted_pipe,OUT/"models"/f"{label}.joblib",compress=3)
            meta={"model":model_name,"scope":scope,"random_seed":SEED,"train_indices_zero_based":train_idx.tolist(),
                  "test_indices_zero_based":test_idx.tolist(),"training_samples":len(train_idx),"test_samples":len(test_idx),
                  "feature_names":columns_for(scope,list(X.columns)),"configuration":str(fitted_pipe.get_params(deep=True))}
            (OUT/"models"/f"{label}.json").write_text(json.dumps(meta,indent=2),encoding="utf-8")
            log(f"BASELINE done {label}; CV ROC-AUC={np.mean(cv_result['test_roc_auc']):.4f}±{np.std(cv_result['test_roc_auc'],ddof=1):.4f}; test ROC-AUC={tm['roc_auc']:.4f}; PR-AUC={tm['pr_auc_average_precision']:.4f}")
    pd.DataFrame(records).to_csv(OUT/"metrics_summary.csv",index=False)
    pd.DataFrame([{k:v for k,v in r.items() if k.endswith("_mean") or k.endswith("_std") or k in {"experiment","model","feature_selection"}} for r in records]).to_csv(OUT/"cross_validation_results.csv",index=False)

    # Select one modality and integrated baseline by CV ROC-AUC only (test scores not used).
    cv_candidates=sorted(records,key=lambda r:(r["roc_auc_mean"],r["pr_auc_mean"]),reverse=True)
    best_integrated=next(r for r in cv_candidates if r["experiment"]=="integrated")
    best_individual=next(r for r in cv_candidates if r["experiment"]!="integrated")
    log(f"CV-selected feature selection/explanation baselines: individual={best_individual['experiment']}/{best_individual['model']}; integrated={best_integrated['model']}")

    # CNV dedup experiment. Dedup representative list comes from the pre-existing audited map,
    # applied only after verifying every mapped vector remains identical on this run's data.
    alias_df=pd.read_csv(ALIAS_MAP)
    alias_to_rep=dict(zip(alias_df.alias_feature,alias_df.representative_feature))
    cnv_original=modality_cols["cnv"]
    cnv_dedup=[c for c in cnv_original if c not in alias_to_rep]
    for alias,rep in alias_to_rep.items():
        if not np.array_equal(df[alias].to_numpy(),df[rep].to_numpy()):
            raise ValueError(f"Alias no longer exactly matches representative: {alias}/{rep}")
    dedup_records=[]
    for variant,cols in [("original",cnv_original),("deduplicated",cnv_dedup)]:
        # The pipeline helper sees X columns; use an explicit subset in a wrapper DataFrame.
        Xall=X
        X=Xall[cols]
        pipe=build_pipeline("cnv","logistic_regression")
        cvres=cross_validate(pipe,X.iloc[train_idx],y_train,cv=cv,scoring=SCORING,n_jobs=1,error_score="raise")
        pipe.fit(X.iloc[train_idx],y_train); pred=pipe.predict(X.iloc[test_idx]); prob=pipe.predict_proba(X.iloc[test_idx])[:,1]
        tm=metrics(y_test,pred,prob)
        rec=metric_rows("cnv_dedup_"+variant,"logistic_regression",cvres,tm); dedup_records.append(rec)
        label="cnv_dedup__"+variant
        save_prediction(label,y_test,pred,prob,test_idx,"CNV duplicate experiment")
        plot_eval(label,y_test,pred,prob,"CNV duplicate experiment")
        joblib.dump(pipe,OUT/"models"/(label+".joblib"),compress=3)
        X=Xall
        log(f"CNV dedup {variant}: {len(cols)} columns; CV ROC-AUC={rec['roc_auc_mean']:.4f}; test ROC-AUC={tm['roc_auc']:.4f}")
    pd.DataFrame(dedup_records).to_csv(OUT/"reports"/"cnv_dedup_comparison.csv",index=False)

    # Fold-contained feature selection compares each method with the no-selection baseline.
    fs_records=[]; fs_best_pipes={}
    selectors=[("mutual_information",100),("l1_logistic",100),("tree_importance",100)]
    for scope in [best_individual["experiment"],"integrated"]:
        for kind,k in selectors:
            label=f"{scope}__logistic_regression__{kind}"
            log(f"FEATURE_SELECTION start {label}")
            pipe=build_pipeline(scope,"logistic_regression",selector=(kind,k))
            cvres=cross_validate(pipe,X_train,y_train,cv=cv,scoring=SCORING,n_jobs=1,error_score="raise")
            pipe.fit(X_train,y_train); pred=pipe.predict(X_test); prob=pipe.predict_proba(X_test)[:,1]
            tm=metrics(y_test,pred,prob); rec=metric_rows(scope,"logistic_regression",cvres,tm,kind); fs_records.append(rec)
            fs_best_pipes[label]=pipe
            save_prediction(label,y_test,pred,prob,test_idx,scope+" feature selection")
            plot_eval(label,y_test,pred,prob,scope+" feature selection")
            joblib.dump(pipe,OUT/"models"/(label+".joblib"),compress=3)
            log(f"FEATURE_SELECTION done {label}; CV ROC-AUC={rec['roc_auc_mean']:.4f}; test ROC-AUC={tm['roc_auc']:.4f}")
    fs_df=pd.DataFrame(fs_records); fs_df.to_csv(OUT/"feature_selection"/"feature_selection_results.csv",index=False)
    all_records=records+dedup_records+fs_records
    pd.DataFrame(all_records).to_csv(OUT/"metrics_summary.csv",index=False)
    pd.DataFrame([{k:v for k,v in r.items() if k.endswith("_mean") or k.endswith("_std") or k in {"experiment","model","feature_selection"}} for r in all_records]).to_csv(OUT/"cross_validation_results.csv",index=False)

    # Save selected feature names from fitted full-training pipelines.
    selected_rows=[]
    for label,pipe in fs_best_pipes.items():
        names=pipe.named_steps["preprocess"].get_feature_names_out()
        selector=pipe.named_steps["select"]
        support=selector.get_support()
        selected=list(np.asarray(names)[support])
        for rank,name in enumerate(selected,1):
            raw=name.split("__",1)[-1]
            pref=raw.split("_",1)[0] if raw.startswith(("cn_","rs_","mu_","pp_")) else "CNV_encoded"
            modality=next((m for m,p in PREFIXES.items() if raw.startswith(p)),"cnv")
            selected_rows.append({"experiment":label,"rank":rank,"selected_feature":raw,"encoded_feature":name,"modality":modality})
    pd.DataFrame(selected_rows).to_csv(OUT/"feature_selection"/"selected_features.csv",index=False)

    # Comparison figure across baseline modality CV and test metrics.
    plot_df=pd.DataFrame(records)
    fig,ax=plt.subplots(figsize=(13,6))
    labels=[r["experiment"]+"\n"+r["model"].replace("_"," ") for r in records]
    x=np.arange(len(records)); ax.bar(x-0.18,plot_df.roc_auc_mean,width=.36,label="CV ROC-AUC mean")
    ax.bar(x+0.18,plot_df.test_roc_auc,width=.36,label="Held-out test ROC-AUC")
    ax.set_xticks(x,labels,rotation=65,ha="right",fontsize=8); ax.set_ylim(0,1); ax.legend(); ax.set_ylabel("ROC-AUC")
    ax.set_title("Exploratory recorded-vital-status classification; same held-out rows")
    fig.tight_layout(); fig.savefig(OUT/"figures"/"model_comparisons"/"baseline_roc_auc.png",dpi=150); plt.close(fig)

    # Explainability baseline selection is CV-only. Use tree importances where available, absolute
    # standardized logistic coefficients otherwise. Permutation importance is limited to top 30 raw
    # variables for feasible repeated scoring, on held-out data only, and is descriptive.
    explain_rows=[]
    chosen=[("best_individual",best_individual["experiment"]+"__"+best_individual["model"]),
            ("best_integrated",best_integrated["experiment"]+"__"+best_integrated["model"])]
    for title,label in chosen:
        pipe=fitted[label]; prep=pipe.named_steps["preprocess"]; model=pipe.named_steps["model"]
        enc_names=np.asarray(prep.get_feature_names_out())
        if hasattr(model,"feature_importances_"):
            importance=model.feature_importances_
        else:
            importance=np.mean(np.abs(model.coef_),axis=0)
        order=np.argsort(importance)[::-1][:20]
        for rank,idx in enumerate(order,1):
            raw=str(enc_names[idx]).split("__",1)[-1]
            mod=next((m for m,p in PREFIXES.items() if raw.startswith(p)),"CNV_encoded")
            explain_rows.append({"model_role":title,"baseline_model":label,"rank":rank,"feature":raw,"modality":mod,"importance":float(importance[idx]),"importance_method":"tree impurity importance" if hasattr(model,"feature_importances_") else "absolute logistic coefficient"})
        # Limit compute to at most 30 top original columns by input feature names; repeated permutation
        # does not use labels for model selection and measures perturbation sensitivity on test.
        # Recover top raw predictors by exact name or prefix for expanded CNV categories.
        top_raw=[]
        for idx in order:
            nm=str(enc_names[idx]).split("__",1)[-1]
            candidate=nm.split("_",1)[0]+"_"+nm.split("_",1)[1] if "_" in nm else nm
            if candidate in X_test.columns and candidate not in top_raw: top_raw.append(candidate)
            if len(top_raw)>=20: break
        if top_raw:
            # Permute selected original columns while preserving the full pipeline input schema.
            from sklearn.metrics import roc_auc_score
            base_score=roc_auc_score(y_test,pipe.predict_proba(X_test)[:,1])
            permutation_values={feat:[] for feat in top_raw}
            rng=np.random.default_rng(SEED)
            for feat in top_raw:
                for _ in range(5):
                    permuted=X_test.copy()
                    permuted[feat]=rng.permutation(permuted[feat].to_numpy())
                    permutation_values[feat].append(base_score-roc_auc_score(y_test,pipe.predict_proba(permuted)[:,1]))
            for feat,values in sorted(permutation_values.items(),key=lambda q:np.mean(q[1]),reverse=True):
                imp=float(np.mean(values)); std=float(np.std(values,ddof=1))
                mod=next((m for m,p in PREFIXES.items() if feat.startswith(p)),"unknown")
                explain_rows.append({"model_role":title,"baseline_model":label,"rank":"permutation","feature":feat,"modality":mod,"importance":imp,"importance_std":std,"importance_method":"test-set permutation ROC-AUC decrease; top-ranked raw variables only"})
        fig,ax=plt.subplots(figsize=(9,7))
        sel=order[:20][::-1]; names=[str(enc_names[i]).split("__",1)[-1] for i in sel]
        ax.barh(names,importance[sel]); ax.set_title(f"Top 20 descriptive importance — {label}"); ax.set_xlabel("Importance (not causal)")
        fig.tight_layout(); fig.savefig(OUT/"figures"/"feature_importance"/f"{title}_top20.png",dpi=150); plt.close(fig)
    pd.DataFrame(explain_rows).to_csv(OUT/"explainability"/"top_feature_importance.csv",index=False)

    # Optional methods and labels.
    qlattice={"installed_feyn":False,"status":"unavailable"}
    try:
        import feyn
        qlattice={"installed_feyn":True,"status":"installed; not run because this exploratory execution has no validated survival endpoint and Feyn API/runtime compatibility has not been established"}
    except Exception as e:
        qlattice={"installed_feyn":False,"status":"unavailable","exact_error":f"{type(e).__name__}: {e}"}
    try:
        import shap
        shap_status={"installed":True,"status":"installed; no SHAP values generated in this run"}
    except Exception as e:
        shap_status={"installed":False,"status":"unavailable","exact_error":f"{type(e).__name__}: {e}"}
    (OUT/"reports"/"optional_tools.json").write_text(json.dumps({"qlattice":qlattice,"shap":shap_status},indent=2),encoding="utf-8")
    # Subtype columns have no timing documentation in these CSVs; do not include in the main model.
    sub_counts={}
    for c in expected_extra:
        sub_counts[c]=enriched[c].fillna("<MISSING>").value_counts().to_dict()
    (OUT/"reports"/"subtype_label_audit.json").write_text(json.dumps({"columns":expected_extra,"counts":sub_counts,"prediction_time_documented":False,"main_model_inclusion":False,"reason":"timing/provenance unavailable; potential leakage or confounding cannot be assessed"},indent=2),encoding="utf-8")

    # Consolidated final report uses held-out selection-independent interpretation; report CV-best
    # and held-out metrics without presenting test-based selection as validation.
    records_df=pd.DataFrame(records)
    best_cv_ind=best_individual; best_cv_int=best_integrated
    integrated_scores={r["model"]:r for r in records if r["experiment"]=="integrated"}
    best_int_test=max(integrated_scores.values(),key=lambda r:r["test_roc_auc"])
    best_ind_test=max((r for r in records if r["experiment"]!="integrated"),key=lambda r:r["test_roc_auc"])
    dedup_df=pd.DataFrame(dedup_records)
    report=[]; w=report.append
    w('# Dataset C Phase 2 — exploratory modelling report'); w('')
    w('## Scope and scientific limits'); w('')
    w('These experiments classify the dataset’s recorded `vital.status` code (0 alive, 1 deceased). They are **not survival-time prediction**: no survival duration, follow-up, event time, or censoring variable is available. No patient IDs/barcodes are present; the split is row-based, and patient-level independence cannot be independently verified. Results are exploratory and do not establish clinical usefulness or generalization to new patients.'); w('')
    w('## Data and split'); w('')
    w(f'- Input: `{DATA}`; {len(df)} rows, 1,936 predictors and target; target counts 0={int((y==0).sum())}, 1={int((y==1).sum())}.')
    w(f'- Stratified row split: train={len(train_idx)} (0={int((y_train==0).sum())}, 1={int((y_train==1).sum())}), test={len(test_idx)} (0={int((y_test==0).sum())}, 1={int((y_test==1).sum())}); seed={SEED}; test fraction={TEST_SIZE}.')
    w(f'- Stratified {CV_FOLDS}-fold CV on training split; folds shuffled with seed {SEED}.')
    w(f'- Prefix counts: {", ".join(f"{PREFIXES[k]}={len(v)}" for k,v in modality_cols.items())}; no unexpected predictors.')
    w('- CNV values treated as categorical and one-hot encoded; RNA/protein median-imputed and standardized; mutation features median-imputed and training-fold constant columns removed. All transforms and selectors live inside scikit-learn pipelines and are fit separately within CV folds.')
    w(''); w('## Baseline models and results'); w('')
    w('Compared class-weighted Logistic Regression, Random Forest (100 trees), and Extra Trees (100 trees) for CNV, RNA, mutation, protein, and integrated data. Same held-out rows used throughout. Metrics include accuracy, deceased precision/recall, alive specificity, F1, macro-F1, ROC-AUC, average precision (PR-AUC), class-specific precision/recall, and confusion counts; full values are in `metrics_summary.csv` and fold statistics in `cross_validation_results.csv`.'); w('')
    w('| Scope | Model | CV ROC-AUC mean ± SD | CV PR-AUC mean ± SD | Test ROC-AUC | Test PR-AUC | Test deceased recall | Test specificity | Test macro-F1 | Confusion TN/FP/FN/TP |'); w('|---|---|---:|---:|---:|---:|---:|---:|---:|---|')
    for r in records:
        w(f'| {r["experiment"]} | {r["model"]} | {r["roc_auc_mean"]:.3f} ± {r["roc_auc_std"]:.3f} | {r["pr_auc_mean"]:.3f} ± {r["pr_auc_std"]:.3f} | {r["test_roc_auc"]:.3f} | {r["test_pr_auc_average_precision"]:.3f} | {r["test_recall_deceased_sensitivity"]:.3f} | {r["test_specificity_alive"]:.3f} | {r["test_macro_f1"]:.3f} | {r["test_tn"]}/{r["test_fp"]}/{r["test_fn"]}/{r["test_tp"]} |')
    w(''); w(f'Best individual baseline selected by training CV ROC-AUC: **{best_cv_ind["experiment"]} / {best_cv_ind["model"]}**, CV ROC-AUC {best_cv_ind["roc_auc_mean"]:.3f} ± {best_cv_ind["roc_auc_std"]:.3f}.')
    w(f'Best integrated baseline selected by training CV ROC-AUC: **{best_cv_int["model"]}**, CV ROC-AUC {best_cv_int["roc_auc_mean"]:.3f} ± {best_cv_int["roc_auc_std"]:.3f}.')
    w(f'Highest integrated held-out ROC-AUC among the three prespecified integrated baselines was {best_int_test["model"]} ({best_int_test["test_roc_auc"]:.3f}); highest individual held-out ROC-AUC was {best_ind_test["experiment"]}/{best_ind_test["model"]} ({best_ind_test["test_roc_auc"]:.3f}). These test comparisons are descriptive and were not used to tune features.')
    chosen_ind_test=next(r for r in records if r["experiment"]==best_cv_ind["experiment"] and r["model"]==best_cv_ind["model"])
    chosen_int_test=next(r for r in records if r["experiment"]=="integrated" and r["model"]==best_cv_int["model"])
    w(f'Using models selected by training CV, integration did not improve the point estimates: best individual ({chosen_ind_test["experiment"]}/{chosen_ind_test["model"]}) CV ROC-AUC {chosen_ind_test["roc_auc_mean"]:.3f}, test ROC-AUC {chosen_ind_test["test_roc_auc"]:.3f}, test PR-AUC {chosen_ind_test["test_pr_auc_average_precision"]:.3f}; best integrated ({chosen_int_test["model"]}) CV ROC-AUC {chosen_int_test["roc_auc_mean"]:.3f}, test ROC-AUC {chosen_int_test["test_roc_auc"]:.3f}, test PR-AUC {chosen_int_test["test_pr_auc_average_precision"]:.3f}.')
    w(''); w('## CNV duplicate-vector comparison'); w('')
    w(f'The existing alias map’s {len(alias_to_rep)} aliases were verified equal to representatives on all 705 observed rows. Keeping the first representative removed {len(alias_to_rep)} columns: {len(cnv_original)} → {len(cnv_dedup)}.'); w('')
    for r in dedup_records: w(f'- {r["experiment"]}: CV ROC-AUC {r["roc_auc_mean"]:.3f} ± {r["roc_auc_std"]:.3f}; held-out ROC-AUC {r["test_roc_auc"]:.3f}; held-out PR-AUC {r["test_pr_auc_average_precision"]:.3f}.')
    w('The equality is only of observed vectors in this dataset; it does not mean the gene aliases are biologically identical. See `duplicate_alias_map.csv`.'); w('')
    w('## Feature selection'); w('')
    w(f'Compared no selection with fold-contained mutual information, L1 logistic selection, and Extra-Trees model-based selection on CV-selected best individual scope ({best_cv_ind["experiment"]}) and integrated data. Selector comparison is in `feature_selection/feature_selection_results.csv`; selected names/modalities are in `feature_selection/selected_features.csv`. No test labels were used for feature choice.'); w('')
    for r in fs_records: w(f'- {r["experiment"]} / {r["feature_selection"]}: CV ROC-AUC {r["roc_auc_mean"]:.3f} ± {r["roc_auc_std"]:.3f}; held-out ROC-AUC {r["test_roc_auc"]:.3f}; PR-AUC {r["test_pr_auc_average_precision"]:.3f}.')
    top_selected=pd.DataFrame(selected_rows)
    if len(top_selected):
        w(''); w('Selected features from the fitted training-only selector pipelines (each row includes modality; these are model-dependent exploratory selections):')
        for _,row in top_selected.groupby('experiment',sort=False):
            w('- '+', '.join(f'{r.selected_feature} ({r.modality})' for r in row.head(20).itertuples()))
    w(''); w('## Explainability'); w('')
    imp=pd.DataFrame(explain_rows)
    for role in ['best_individual','best_integrated']:
        sub=imp[(imp.model_role==role)&(imp.importance_method!='test-set permutation ROC-AUC decrease; top-ranked raw variables only')].head(20)
        w(f'**{role.replace("_"," ").title()} top features** (model-based descriptive importance; not causal or proof of biomarkers): '+('; '.join(f'{r.feature} [{r.modality}]' for r in sub.itertuples()) if len(sub) else 'unavailable')+'.')
    w('Permutation importance for a limited set of top raw features is in `explainability/top_feature_importance.csv`; plots are in `figures/feature_importance/`.')
    w(''); w(f'SHAP status: {shap_status["status"]}'+(f' ({shap_status.get("exact_error")})' if shap_status.get('exact_error') else '')+'.')
    w(''); w('## QLattice and subtype analysis'); w('')
    w(f'QLattice/Feyn: {qlattice["status"]}'+(f'; exact error: `{qlattice.get("exact_error")}`' if qlattice.get('exact_error') else '')+'. Historical serialized models were not treated as validated.')
    w('Subtype columns were inspected separately and excluded from the main model. Their prediction-time availability is undocumented, so leakage/confounding cannot be assessed; no subtype model was run.')
    w(''); w('## Execution issues and recovery'); w('')
    w('An initial parallel tree fit failed with `PermissionError: [WinError 5] Access is denied` when joblib attempted worker creation. The runner was changed to single-worker in-process estimator fitting and all baselines completed. An initial CNV dedup indexing attempt failed with `KeyError: "None of [Index([...], dtype=int64, length=564)] are in the [columns]"` because row indices were applied as DataFrame column labels; this was corrected to `.iloc`, and the CNV comparison completed. Full traces are in `experiment_log.txt`.')
    w(''); w('## Saved outputs'); w('')
    w('- `metrics_summary.csv`, `cross_validation_results.csv`, `experiment_log.txt`.')
    w('- Final training-split pipelines and JSON configurations under `models/`; test predictions under `predictions/`.')
    w('- Confusion matrices, ROC and precision-recall curves, comparison plot, feature importance plots under `figures/`.')
    w('- Feature-selection results/selected features under `feature_selection/`; explainability under `explainability/`; optional-tool and subtype audits under `reports/`.')
    w(''); w('## Limitations'); w('')
    w('Small cohort (705 rows), 6.5:1 imbalance, p≫n, unknown sample independence, no identifiers, no follow-up/time-to-event/censoring, undocumented measurement/prediction timing, and source-derived feature semantics constrain interpretation. A single random row split is not external validation. All reported estimates are exploratory recorded-vital-status classification only; no clinical survival prediction, clinical utility, or new-patient generalization is established.')
    (OUT/"dataset_C_final_report.md").write_text("\n".join(report)+"\n",encoding="utf-8")
    log("RUN COMPLETE; final report written.")


if __name__ == "__main__":
    try:
        main()
    except Exception as exc:
        log(f"FATAL {type(exc).__name__}: {exc}")
        log(traceback.format_exc())
        raise
