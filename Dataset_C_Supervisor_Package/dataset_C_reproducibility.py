#!/usr/bin/env python
"""Reproduce the verified Dataset C multi-omics experiments.

This script uses only ``dataset_c_extracted/data.csv`` (unless an explicit
``--data`` path is supplied). It never reads saved metrics to construct
results. All preprocessing and model/feature selection is fitted within the
appropriate training fold. The held-out set is used only after training-CV
selection, for a final descriptive evaluation.

Scientific scope: binary classification of the supplied recorded
``vital.status`` (0=Alive, 1=Deceased). This is not survival-time prediction,
does not have a prediction horizon, and is not clinically validated. The
source file has no patient identifiers, so patient-level independence cannot
be verified.
"""
from __future__ import annotations

import argparse
import hashlib
import json
import os
import platform
import sys
import textwrap
import traceback
from datetime import datetime
from pathlib import Path
from typing import Any

try:
    import joblib
    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt
    import numpy as np
    import pandas as pd
    import sklearn
    from sklearn.compose import ColumnTransformer
    from sklearn.dummy import DummyClassifier
    from sklearn.ensemble import ExtraTreesClassifier, HistGradientBoostingClassifier, RandomForestClassifier
    from sklearn.feature_selection import SelectFromModel, SelectKBest, VarianceThreshold, f_classif, mutual_info_classif
    from sklearn.impute import SimpleImputer
    from sklearn.linear_model import LogisticRegression
    from sklearn.metrics import (
        accuracy_score, average_precision_score, balanced_accuracy_score,
        classification_report, confusion_matrix, f1_score,
        get_scorer, make_scorer, precision_score, precision_recall_curve, recall_score,
        roc_auc_score, roc_curve,
    )
    from sklearn.model_selection import StratifiedKFold, train_test_split
    from sklearn.pipeline import Pipeline
    from sklearn.preprocessing import OneHotEncoder, StandardScaler
    from threadpoolctl import threadpool_limits
    from sklearn.svm import SVC
except ModuleNotFoundError as exc:  # clear setup help for a supervisor's machine
    print(f"Missing required Python package: {exc.name}", file=sys.stderr)
    print("Install the pinned dependencies with: python -m pip install -r requirements.txt", file=sys.stderr)
    raise SystemExit(2) from exc


ROOT = Path(__file__).resolve().parent
DEFAULT_DATA = ROOT / "dataset_c_extracted" / "data.csv"
OUTPUT = ROOT / "dataset_C_reproducibility_results"
SEED = 20261110
TEST_SIZE = 0.20
OUTER_FOLDS = 3
INNER_FOLDS = 3
SELECTOR_FOLDS = 5
SELECTOR_K = 50
TREES = 100
N_JOBS = 1
TARGET = "vital.status"
TARGET_MAPPING = {0: "Alive", 1: "Deceased"}
MODALITIES = {"CNV": "cn_", "RNA": "rs_", "Mutation": "mu_", "Protein": "pp_"}
SCOPES = {"CNV": ["CNV"], "RNA": ["RNA"], "Mutation": ["Mutation"],
          "Protein": ["Protein"], "Integrated four-omics": list(MODALITIES)}
MODELS = ["logistic_regression", "random_forest", "extra_trees", "rbf_svm", "hist_gradient_boosting"]
PRIORITY_METRICS = ["roc_auc", "pr_auc", "balanced_accuracy", "macro_f1", "deceased_recall"]
MODEL_DISPLAY = {
    "logistic_regression": "Logistic Regression",
    "random_forest": "Random Forest",
    "extra_trees": "Extra Trees",
    "rbf_svm": "RBF SVM",
    "hist_gradient_boosting": "HistGradientBoosting",
}


class RunLog:
    """Mirror all progress, tables and errors to terminal and output log."""
    def __init__(self, path: Path):
        self.path = path
        self.path.parent.mkdir(parents=True, exist_ok=True)
        self.handle = self.path.open("w", encoding="utf-8", newline="\n")

    def write(self, message: str = "") -> None:
        print(message, flush=True)
        self.handle.write(str(message) + "\n")
        self.handle.flush()

    def close(self) -> None:
        self.handle.close()


def file_sha256(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as stream:
        for block in iter(lambda: stream.read(1024 * 1024), b""):
            digest.update(block)
    return digest.hexdigest()


def modality_of(feature: str) -> str:
    for modality, prefix in MODALITIES.items():
        if feature.startswith(prefix):
            return modality
    return "Unknown"


def selected_columns(scope: str, columns: list[str]) -> list[str]:
    wanted = SCOPES[scope]
    return [c for c in columns if any(c.startswith(MODALITIES[m]) for m in wanted)]


def build_pipeline(scope: str, model_name: str, columns: list[str],
                   seed: int = SEED, selector: tuple[str, int] | None = None) -> Pipeline:
    """Rebuild the verified preprocessing and estimator pipeline."""
    cols = selected_columns(scope, columns)
    cnv = [c for c in cols if c.startswith(MODALITIES["CNV"])]
    mutation = [c for c in cols if c.startswith(MODALITIES["Mutation"])]
    continuous = [c for c in cols if c.startswith((MODALITIES["RNA"], MODALITIES["Protein"]))]
    transformers: list[tuple[str, Any, list[str]]] = []
    if cnv:
        transformers.append(("cnv", Pipeline([
            ("impute", SimpleImputer(strategy="most_frequent")),
            ("onehot", OneHotEncoder(handle_unknown="ignore", sparse_output=False)),
        ]), cnv))
    if mutation:
        transformers.append(("mutation", Pipeline([
            ("impute", SimpleImputer(strategy="median")),
            ("variance", VarianceThreshold(threshold=0.0)),
        ]), mutation))
    if continuous:
        transformers.append(("continuous", Pipeline([
            ("impute", SimpleImputer(strategy="median")),
            ("scale", StandardScaler()),
        ]), continuous))
    preprocessor = ColumnTransformer(
        transformers, remainder="drop", sparse_threshold=0,
        verbose_feature_names_out=False,
    )
    steps: list[tuple[str, Any]] = [("preprocess", preprocessor)]
    if selector:
        kind, k = selector
        if kind == "mutual_information":
            selection = SelectKBest(lambda a, b: mutual_info_classif(a, b, random_state=seed), k=k)
        elif kind == "anova":
            selection = SelectKBest(f_classif, k=k)
        elif kind == "l1":
            selection = SelectFromModel(
                LogisticRegression(l1_ratio=1.0, solver="liblinear", C=0.08,
                                   class_weight="balanced", max_iter=2000, random_state=seed),
                threshold=-np.inf, max_features=k,
            )
        elif kind == "tree_importance":
            selection = SelectFromModel(
                ExtraTreesClassifier(n_estimators=60, class_weight="balanced", n_jobs=N_JOBS,
                                     random_state=seed, min_samples_leaf=2),
                threshold=-np.inf, max_features=k,
            )
        else:
            raise ValueError(f"Unknown selector: {kind}")
        steps.append(("select", selection))

    if model_name == "logistic_regression":
        estimator = LogisticRegression(class_weight="balanced", C=0.1, solver="liblinear",
                                       max_iter=2500, random_state=seed)
    elif model_name == "random_forest":
        estimator = RandomForestClassifier(n_estimators=TREES, class_weight="balanced",
                                           max_features="sqrt", min_samples_leaf=2,
                                           n_jobs=N_JOBS, random_state=seed)
    elif model_name == "extra_trees":
        estimator = ExtraTreesClassifier(n_estimators=TREES, class_weight="balanced",
                                         max_features="sqrt", min_samples_leaf=2,
                                         n_jobs=N_JOBS, random_state=seed)
    elif model_name == "rbf_svm":
        estimator = SVC(kernel="rbf", C=1.0, gamma="scale", class_weight="balanced",
                        probability=True, random_state=seed)
    elif model_name == "hist_gradient_boosting":
        estimator = HistGradientBoostingClassifier(max_iter=100, learning_rate=0.08,
                                                   max_leaf_nodes=15, l2_regularization=1.0,
                                                   class_weight="balanced", early_stopping=False,
                                                   random_state=seed)
    else:
        raise ValueError(f"Unknown model: {model_name}")
    steps.append(("model", estimator))
    return Pipeline(steps)


def metric_values(y_true: Any, y_pred: Any, y_probability: Any) -> dict[str, float | int]:
    """Metrics with deceased (1) as positive; confusion order is [0, 1]."""
    cm = confusion_matrix(y_true, y_pred, labels=[0, 1])
    tn, fp, fn, tp = (int(v) for v in cm.ravel())
    return {
        "accuracy": float(accuracy_score(y_true, y_pred)),
        "precision_deceased": float(precision_score(y_true, y_pred, pos_label=1, zero_division=0)),
        "deceased_recall": float(recall_score(y_true, y_pred, pos_label=1, zero_division=0)),
        "specificity_alive": float(tn / (tn + fp)) if tn + fp else float("nan"),
        "f1_deceased": float(f1_score(y_true, y_pred, pos_label=1, zero_division=0)),
        "macro_f1": float(f1_score(y_true, y_pred, average="macro", zero_division=0)),
        "balanced_accuracy": float(balanced_accuracy_score(y_true, y_pred)),
        "roc_auc": float(roc_auc_score(y_true, y_probability)),
        "pr_auc_average_precision": float(average_precision_score(y_true, y_probability)),
        "tn": tn, "fp": fp, "fn": fn, "tp": tp,
    }


def metric_scorers() -> dict[str, Any]:
    def specificity_score(y_true: Any, y_pred: Any) -> float:
        matrix = confusion_matrix(y_true, y_pred, labels=[0, 1])
        return float(matrix[0, 0] / matrix[0].sum()) if matrix[0].sum() else 0.0
    return {
        "accuracy": "accuracy",
        "balanced_accuracy": "balanced_accuracy",
        "macro_f1": "f1_macro",
        "roc_auc": "roc_auc",
        "pr_auc": "average_precision",
        "deceased_recall": make_scorer(recall_score, pos_label=1, zero_division=0),
        "specificity": make_scorer(specificity_score),
    }


def fit_pipeline(pipeline: Pipeline, X: pd.DataFrame, y: pd.Series) -> Pipeline:
    """Fit with a single native thread for deterministic Windows execution."""
    with threadpool_limits(limits=1):
        pipeline.fit(X, y)
    return pipeline


def average_rank(metric_means: dict[str, dict[str, float]]) -> dict[str, float]:
    """Average stable rank, matching the verified project's tie-breaking."""
    ranks = {name: 0.0 for name in metric_means}
    for metric in PRIORITY_METRICS:
        # Python's stable sort preserves model insertion order for equal scores,
        # matching the original final_optimization/run_optimization.py.
        ordered = sorted(metric_means, key=lambda name: metric_means[name][metric], reverse=True)
        for position, name in enumerate(ordered, start=1):
            ranks[name] += position / len(PRIORITY_METRICS)
    return ranks


def qc_dataset(df: pd.DataFrame, x: pd.DataFrame, target: pd.Series,
               source: Path, logger: RunLog) -> tuple[pd.DataFrame, pd.DataFrame, pd.DataFrame]:
    logger.write("\n=== Dataset and feature audit ===")
    counts = target.value_counts().reindex([0, 1], fill_value=0)
    predictors = list(x.columns)
    prefix_counts = {m: sum(c.startswith(p) for c in predictors) for m, p in MODALITIES.items()}
    missing_predictor_cells = int(x.isna().sum().sum())
    duplicate_rows = int(df.duplicated().sum())
    feature_unique = x.nunique(dropna=False)
    constant_features = feature_unique[feature_unique <= 1].index.tolist()
    near_constant = []
    for col in predictors:
        frequencies = x[col].value_counts(dropna=False, normalize=True)
        if len(frequencies) and float(frequencies.iloc[0]) >= 0.95:
            near_constant.append(col)

    # This is a descriptive QC screen, not model preprocessing. The current
    # verified source contains no missing predictors. If a changed input has
    # missing cells, median imputation is used only to calculate this screen;
    # modelling imputers remain fold-contained below.
    numeric = x.apply(pd.to_numeric, errors="coerce")
    qc_imputed = numeric.fillna(numeric.median())
    values = qc_imputed.to_numpy(dtype=np.float64)
    means = values.mean(axis=0)
    std = values.std(axis=0, ddof=1)
    std[std == 0] = 1.0
    standardized = (values - means) / std
    corr = (standardized.T @ standardized) / (len(values) - 1)
    np.fill_diagonal(corr, 0.0)
    upper_i, upper_j = np.triu_indices(len(predictors), k=1)
    keep = np.abs(corr[upper_i, upper_j]) >= 0.95
    high_corr = pd.DataFrame({
        "feature_1": np.asarray(predictors, dtype=object)[upper_i[keep]],
        "feature_2": np.asarray(predictors, dtype=object)[upper_j[keep]],
        "pearson_r": corr[upper_i[keep], upper_j[keep]],
    })
    high_corr["modality_1"] = high_corr.feature_1.map(modality_of)
    high_corr["modality_2"] = high_corr.feature_2.map(modality_of)
    high_corr.to_csv(OUTPUT / "highly_correlated_pairs.csv", index=False)

    # Hash each observed CNV vector, then confirm equality within hash groups.
    # The primary models retain these columns, matching the verified final run.
    cnv_columns = [c for c in predictors if c.startswith(MODALITIES["CNV"])]
    vector_groups: dict[str, list[str]] = {}
    for col in cnv_columns:
        vector = x[col].to_numpy(dtype=np.float64)
        canonical = np.nan_to_num(vector, nan=np.finfo(np.float64).max)
        signature = hashlib.sha256(np.ascontiguousarray(canonical).tobytes()).hexdigest()
        vector_groups.setdefault(signature, []).append(col)
    duplicate_group_rows = []
    identical_groups = 0
    redundant_aliases = 0
    for names in vector_groups.values():
        if len(names) < 2:
            continue
        reference = x[names[0]].to_numpy(dtype=np.float64)
        identical = all(np.array_equal(reference, x[name].to_numpy(dtype=np.float64), equal_nan=True)
                        for name in names[1:])
        if identical:
            identical_groups += 1
            redundant_aliases += len(names) - 1
            for alias in names[1:]:
                duplicate_group_rows.append({"representative_feature": names[0],
                                             "alias_feature": alias,
                                             "observed_vectors_identical": True,
                                             "biological_identity_inferred": False})
    duplicate_alias_table = pd.DataFrame(duplicate_group_rows, columns=[
        "representative_feature", "alias_feature", "observed_vectors_identical", "biological_identity_inferred"
    ])
    duplicate_alias_table.to_csv(OUTPUT / "cnv_duplicate_feature_vectors.csv", index=False)

    feature_rows = [{"feature_name": c, "modality": modality_of(c), "position": i + 1}
                    for i, c in enumerate(predictors)]
    feature_names = pd.DataFrame(feature_rows)
    feature_names.to_csv(OUTPUT / "feature_names_by_modality.csv", index=False)

    summary = {
        "source_file": str(source), "source_sha256": file_sha256(source),
        "rows": int(df.shape[0]), "columns_including_target": int(df.shape[1]),
        "predictor_count": int(x.shape[1]), "target_column": TARGET,
        "target_mapping": json.dumps({str(k): v for k, v in TARGET_MAPPING.items()}),
        "alive_count_0": int(counts.loc[0]), "alive_percent_0": float(counts.loc[0] / len(target) * 100),
        "deceased_count_1": int(counts.loc[1]), "deceased_percent_1": float(counts.loc[1] / len(target) * 100),
        "missing_values_target": int(target.isna().sum()),
        "missing_values_predictors": missing_predictor_cells,
        "duplicate_rows_extra_copies": duplicate_rows,
        "constant_predictor_count": len(constant_features),
        "constant_predictors": ";".join(constant_features),
        "near_constant_95pct_count": len(near_constant),
        "near_constant_95pct_features": ";".join(near_constant),
        "high_correlation_threshold": 0.95,
        "highly_correlated_pair_count_abs_pearson_ge_095": int(len(high_corr)),
        "cnv_identical_feature_vector_group_count": identical_groups,
        "cnv_redundant_extra_columns": redundant_aliases,
        "correlation_qc_missing_strategy": "median-filled for descriptive QC only" if missing_predictor_cells else "no missing values",
        "duplicate_CNV_aliases": "retained; verified primary final models did not deduplicate CNV aliases",
        **{f"{m.lower()}_feature_count": int(n) for m, n in prefix_counts.items()},
    }
    pd.DataFrame([summary]).to_csv(OUTPUT / "dataset_summary.csv", index=False)

    logger.write(f"Source CSV: {source}")
    logger.write(f"Source SHA-256: {summary['source_sha256']}")
    logger.write(f"Dataset shape: {df.shape}; samples={df.shape[0]}; predictors={x.shape[1]}")
    logger.write(f"Target: {TARGET}; mapping: 0=Alive, 1=Deceased")
    logger.write(f"Target counts: Alive={counts.loc[0]} ({summary['alive_percent_0']:.2f}%), Deceased={counts.loc[1]} ({summary['deceased_percent_1']:.2f}%)")
    logger.write(f"Missing: target={summary['missing_values_target']}, predictors={missing_predictor_cells}")
    logger.write(f"Duplicate rows (extra copies): {duplicate_rows}")
    logger.write(f"Constant features ({len(constant_features)}): {constant_features}")
    logger.write(f"Near-constant features, dominant value >=95% ({len(near_constant)}): {near_constant}")
    logger.write(f"Highly correlated feature pairs, |Pearson r| >= 0.95: {len(high_corr)}")
    logger.write(f"CNV identical observed-vector groups={identical_groups}; redundant extra columns={redundant_aliases} (retained in models)")
    for modality, n in prefix_counts.items():
        names = [c for c in predictors if c.startswith(MODALITIES[modality])]
        logger.write(f"\n{modality} features ({n}):")
        for chunk_start in range(0, len(names), 8):
            logger.write("  " + ", ".join(names[chunk_start:chunk_start + 8]))
    logger.write("\nCNV duplicate-alias policy: retain all verified columns; duplicate observed vectors do not imply biological identity.")
    return pd.DataFrame([summary]), feature_names, high_corr


def audit_schema(df: pd.DataFrame) -> tuple[pd.DataFrame, pd.Series]:
    if TARGET not in df.columns:
        raise ValueError(f"Required target {TARGET!r} is absent. Refusing to select a substitute dataset/target.")
    x = df.drop(columns=[TARGET]).copy()
    y_raw = df[TARGET]
    if y_raw.isna().any():
        raise ValueError(f"Target {TARGET} contains {int(y_raw.isna().sum())} missing values; inspect the input rather than silently dropping rows.")
    try:
        y = y_raw.astype(int)
    except (TypeError, ValueError) as exc:
        raise ValueError(f"Target {TARGET} must be encoded as integers 0/1; found non-integer values.") from exc
    if set(y.unique()) != {0, 1}:
        raise ValueError(f"Target must contain exactly codes 0 and 1; found {sorted(y.unique().tolist())}.")
    unexpected = [c for c in x.columns if not c.startswith(tuple(MODALITIES.values()))]
    expected_counts = {"CNV": 860, "RNA": 604, "Mutation": 249, "Protein": 223}
    actual_counts = {m: sum(c.startswith(p) for c in x.columns) for m, p in MODALITIES.items()}
    if unexpected or actual_counts != expected_counts:
        raise ValueError(f"Unexpected predictor schema. Non-modality columns={unexpected}; modality counts={actual_counts}; expected={expected_counts}.")
    non_numeric = [c for c in x if not pd.api.types.is_numeric_dtype(x[c])]
    if non_numeric:
        raise ValueError(f"Expected numeric molecular predictors, found non-numeric columns: {non_numeric[:20]}")
    if x.shape != (705, 1936) or y.value_counts().to_dict() != {0: 611, 1: 94}:
        raise ValueError(f"Dataset differs from verified Dataset C schema/labels. X shape={x.shape}; target counts={y.value_counts().to_dict()}.")
    return x, y


def nested_evaluation(X_train: pd.DataFrame, y_train: pd.Series,
                     columns: list[str], logger: RunLog) -> tuple[pd.DataFrame, pd.DataFrame, pd.DataFrame, pd.DataFrame, dict[str, str]]:
    # sklearn's cross_validate accepts names as strings; direct calls below
    # need callable scorer objects. Resolve the named scorers explicitly.
    scorers = {name: get_scorer(value) if isinstance(value, str) else value
               for name, value in metric_scorers().items()}
    outer = StratifiedKFold(n_splits=OUTER_FOLDS, shuffle=True, random_state=SEED)
    outer_rows: list[dict[str, Any]] = []
    inner_rows: list[dict[str, Any]] = []
    outer_predictions: list[dict[str, Any]] = []
    errors: list[dict[str, Any]] = []
    inner_by_scope: dict[str, list[dict[str, dict[str, float]]]] = {s: [] for s in SCOPES}

    for scope in SCOPES:
        logger.write(f"\n[NESTED CV] Scope: {scope}")
        for outer_fold, (outer_train_pos, outer_valid_pos) in enumerate(outer.split(X_train, y_train), start=1):
            X_inner_train = X_train.iloc[outer_train_pos]
            y_inner_train = y_train.iloc[outer_train_pos]
            X_outer_valid = X_train.iloc[outer_valid_pos]
            y_outer_valid = y_train.iloc[outer_valid_pos]
            inner = StratifiedKFold(n_splits=INNER_FOLDS, shuffle=True, random_state=SEED + outer_fold)
            inner_splits = list(inner.split(X_inner_train, y_inner_train))
            candidate_metrics: dict[str, dict[str, float]] = {}
            for model_name in MODELS:
                fold_metrics: list[dict[str, float]] = []
                fit_failed = False
                for inner_fold, (fit_pos, valid_pos) in enumerate(inner_splits, start=1):
                    try:
                        pipeline = build_pipeline(scope, model_name, columns, seed=SEED)
                        fit_pipeline(pipeline, X_inner_train.iloc[fit_pos], y_inner_train.iloc[fit_pos])
                        predicted = pipeline.predict(X_inner_train.iloc[valid_pos])
                        probabilities = pipeline.predict_proba(X_inner_train.iloc[valid_pos])[:, 1]
                        values = metric_values(y_inner_train.iloc[valid_pos], predicted, probabilities)
                        # Use the same sklearn scorer objects as the verified
                        # cross_validate implementation. For SVC, ROC-AUC and
                        # average precision prefer decision_function over
                        # calibrated predict_proba scores.
                        values_for_cv = {
                            name: float(scorer(pipeline, X_inner_train.iloc[valid_pos],
                                               y_inner_train.iloc[valid_pos]))
                            for name, scorer in scorers.items()
                        }
                        fold_metrics.append({k: float(values_for_cv[k]) for k in scorers})
                        inner_rows.append({"scope": scope, "outer_fold": outer_fold,
                                           "model": model_name, "inner_fold": inner_fold,
                                           **{k: float(values_for_cv[k]) for k in scorers}, "status": "ok"})
                    except Exception as exc:
                        fit_failed = True
                        message = f"{type(exc).__name__}: {exc}"
                        errors.append({"stage": "inner_cv", "scope": scope, "outer_fold": outer_fold,
                                       "model": model_name, "inner_fold": inner_fold, "error": message})
                        logger.write(f"ERROR inner CV {scope}/{model_name}/outer{outer_fold}/inner{inner_fold}: {message}")
                        logger.write(traceback.format_exc().rstrip())
                # Incomplete candidate folds are not used in ranking.
                if not fit_failed and len(fold_metrics) == INNER_FOLDS:
                    candidate_metrics[model_name] = {
                        metric: float(np.mean([entry[metric] for entry in fold_metrics])) for metric in scorers
                    }
            inner_by_scope[scope].append(candidate_metrics)
            if not candidate_metrics:
                errors.append({"stage": "outer_selection", "scope": scope, "outer_fold": outer_fold,
                               "model": "all", "inner_fold": "", "error": "No complete inner-CV candidate models."})
                continue
            ranks = average_rank(candidate_metrics)
            selected_model = sorted(candidate_metrics,
                                    key=lambda m: (ranks[m], -candidate_metrics[m]["roc_auc"],
                                                   -candidate_metrics[m]["pr_auc"]))[0]
            try:
                outer_pipeline = build_pipeline(scope, selected_model, columns, seed=SEED)
                fit_pipeline(outer_pipeline, X_inner_train, y_inner_train)
                outer_pred = outer_pipeline.predict(X_outer_valid)
                outer_prob = outer_pipeline.predict_proba(X_outer_valid)[:, 1]
                outer_values = metric_values(y_outer_valid, outer_pred, outer_prob)
                outer_rows.append({"scope": scope, "outer_fold": outer_fold,
                                   "selected_algorithm_by_inner_cv": selected_model,
                                   "inner_average_priority_rank": float(ranks[selected_model]),
                                   **outer_values})
                for local_pos, true_value, pred_value, probability in zip(
                    outer_valid_pos, y_outer_valid, outer_pred, outer_prob
                ):
                    outer_predictions.append({
                        "scope": scope, "training_row_position": int(local_pos),
                        "true_vital_status": int(true_value),
                        "predicted_vital_status": int(pred_value),
                        "probability_deceased": float(probability),
                        "outer_fold": outer_fold, "selected_algorithm": selected_model,
                    })
            except Exception as exc:
                message = f"{type(exc).__name__}: {exc}"
                errors.append({"stage": "outer_fit", "scope": scope, "outer_fold": outer_fold,
                               "model": selected_model, "inner_fold": "", "error": message})
                logger.write(f"ERROR outer fit {scope}/{selected_model}/outer{outer_fold}: {message}")
                logger.write(traceback.format_exc().rstrip())
        logger.write(f"Completed nested CV scope: {scope}")

    # Rank each algorithm per scope from inner CV means across outer training
    # partitions, then rank the scope winners across the same validation metrics.
    candidate_records: list[dict[str, Any]] = []
    scope_winners: dict[str, str] = {}
    for scope in SCOPES:
        outer_summaries = inner_by_scope[scope]
        candidates: dict[str, dict[str, float]] = {}
        for model_name in MODELS:
            available = [folds[model_name] for folds in outer_summaries if model_name in folds]
            if available:
                candidates[model_name] = {m: float(np.mean([a[m] for a in available])) for m in PRIORITY_METRICS}
        if not candidates:
            continue
        ranks = average_rank(candidates)
        for model_name, means in candidates.items():
            fold_values = [row for row in inner_rows
                           if row["scope"] == scope and row["model"] == model_name
                           and row["status"] == "ok"]
            candidate_records.append({"scope": scope, "model": model_name,
                                      "outer_training_inner_cv_replicates": len([f for f in outer_summaries if model_name in f]),
                                      **{f"inner_{m}_mean": v for m, v in means.items()},
                                      **{f"inner_{m}_sd_fold_scores": float(np.std([r[m] for r in fold_values], ddof=1))
                                         for m in PRIORITY_METRICS},
                                      "mean_composite_rank": float(ranks[model_name]),
                                      "scope_selected": False})
        winner = sorted(candidates, key=lambda m: (ranks[m], -candidates[m]["roc_auc"],
                                                    -candidates[m]["pr_auc"]))[0]
        scope_winners[scope] = winner
        for row in candidate_records:
            if row["scope"] == scope and row["model"] == winner:
                row["scope_selected"] = True

    outer_df = pd.DataFrame(outer_rows)
    inner_df = pd.DataFrame(inner_rows)
    candidate_df = pd.DataFrame(candidate_records)
    oof_df = pd.DataFrame(outer_predictions)
    outer_df.to_csv(OUTPUT / "nested_fold_results.csv", index=False)
    inner_df.to_csv(OUTPUT / "candidate_inner_fold_scores.csv", index=False)
    candidate_df.to_csv(OUTPUT / "candidate_cv_scores.csv", index=False)
    oof_df.to_csv(OUTPUT / "nested_oof_predictions.csv", index=False)
    pd.DataFrame(errors, columns=["stage", "scope", "outer_fold", "model", "inner_fold", "error"]).to_csv(
        OUTPUT / "model_fit_errors.csv", index=False)
    return outer_df, inner_df, candidate_df, oof_df, scope_winners


def run_dummy_baseline(X_train: pd.DataFrame, y_train: pd.Series,
                       X_test: pd.DataFrame, y_test: pd.Series,
                       test_positions: np.ndarray, logger: RunLog) -> dict[str, Any]:
    fold = StratifiedKFold(n_splits=OUTER_FOLDS, shuffle=True, random_state=SEED)
    oof_prob = np.zeros(len(y_train), dtype=float)
    oof_pred = np.zeros(len(y_train), dtype=int)
    for fit_pos, valid_pos in fold.split(X_train, y_train):
        dummy = DummyClassifier(strategy="prior").fit(X_train.iloc[fit_pos], y_train.iloc[fit_pos])
        oof_prob[valid_pos] = dummy.predict_proba(X_train.iloc[valid_pos])[:, 1]
        oof_pred[valid_pos] = dummy.predict(X_train.iloc[valid_pos])
    cv_metrics = metric_values(y_train, oof_pred, oof_prob)
    final_dummy = DummyClassifier(strategy="prior").fit(X_train, y_train)
    probability = final_dummy.predict_proba(X_test)[:, 1]
    prediction = final_dummy.predict(X_test)
    test_metrics = metric_values(y_test, prediction, probability)
    pd.DataFrame([{"model": "DummyClassifier prior", **{f"cv_{k}": v for k, v in cv_metrics.items()},
                   **{f"test_{k}": v for k, v in test_metrics.items()}}]).to_csv(OUTPUT / "dummy_baseline.csv", index=False)
    pd.DataFrame({"original_row_position": test_positions, "true_vital_status": y_test.to_numpy(),
                  "predicted_vital_status": prediction, "probability_deceased": probability}).to_csv(
        OUTPUT / "dummy_baseline_predictions.csv", index=False)
    logger.write("DummyClassifier baseline complete (prior strategy; threshold predictions from training prevalence).")
    logger.write("Dummy baseline CV metrics: " + json.dumps(cv_metrics, sort_keys=True))
    logger.write("Dummy baseline holdout metrics: " + json.dumps(test_metrics, sort_keys=True))
    return {"cv": cv_metrics, "test": test_metrics}


def selector_score_vector(selector: Any, feature_count: int) -> np.ndarray:
    if hasattr(selector, "scores_") and selector.scores_ is not None:
        return np.asarray(selector.scores_, dtype=float)
    estimator = getattr(selector, "estimator_", None)
    if estimator is not None and hasattr(estimator, "coef_"):
        return np.asarray(np.abs(estimator.coef_).mean(axis=0), dtype=float)
    if estimator is not None and hasattr(estimator, "feature_importances_"):
        return np.asarray(estimator.feature_importances_, dtype=float)
    return np.full(feature_count, np.nan)


def map_encoded_name(encoded: str, raw_columns: list[str]) -> str:
    if encoded in raw_columns:
        return encoded
    matches = [c for c in raw_columns if encoded.startswith(c + "_")]
    return max(matches, key=len) if matches else encoded


def run_feature_selection(X_train: pd.DataFrame, y_train: pd.Series,
                          columns: list[str], logger: RunLog) -> tuple[pd.DataFrame, pd.DataFrame]:
    records: list[dict[str, Any]] = []
    fold_metrics: list[dict[str, Any]] = []
    splitter = StratifiedKFold(n_splits=SELECTOR_FOLDS, shuffle=True, random_state=SEED)
    methods = ["mutual_information", "l1", "tree_importance", "anova"]
    for scope in ["RNA", "Integrated four-omics"]:
        raw_scope_columns = selected_columns(scope, columns)
        for method in methods:
            for fold_number, (fit_pos, valid_pos) in enumerate(splitter.split(X_train, y_train), start=1):
                pipeline = build_pipeline(scope, "logistic_regression", columns, seed=SEED,
                                          selector=(method, SELECTOR_K))
                pipeline.fit(X_train.iloc[fit_pos], y_train.iloc[fit_pos])
                predicted = pipeline.predict(X_train.iloc[valid_pos])
                probability = pipeline.predict_proba(X_train.iloc[valid_pos])[:, 1]
                fold_metrics.append({"scope": scope, "selection_method": method, "fold": fold_number,
                                     **metric_values(y_train.iloc[valid_pos], predicted, probability)})
                encoded_names = np.asarray(pipeline.named_steps["preprocess"].get_feature_names_out(), dtype=object)
                selector = pipeline.named_steps["select"]
                support = selector.get_support()
                scores = selector_score_vector(selector, len(encoded_names))
                for encoded_index in np.flatnonzero(support):
                    encoded_name = str(encoded_names[encoded_index])
                    feature = map_encoded_name(encoded_name, raw_scope_columns)
                    score = float(scores[encoded_index]) if encoded_index < len(scores) else float("nan")
                    records.append({"scope": scope, "selection_method": method, "fold": fold_number,
                                    "selected_encoded_feature": encoded_name,
                                    "feature_name": feature, "modality": modality_of(feature),
                                    "selection_score": score, "selected_encoded_count_this_fold": int(support.sum())})
            logger.write(f"Feature selection completed: {scope} / {method}; {SELECTOR_FOLDS} training folds")
    detail = pd.DataFrame(records)
    detail.to_csv(OUTPUT / "selected_features_by_fold.csv", index=False)
    metrics = pd.DataFrame(fold_metrics)
    metrics.to_csv(OUTPUT / "feature_selection_fold_metrics.csv", index=False)
    if detail.empty:
        return detail, metrics
    summary = detail.groupby(["scope", "selection_method", "feature_name", "modality"], as_index=False).agg(
        folds_selected=("fold", "nunique"), mean_selection_score=("selection_score", "mean"),
        mean_selected_encoded_count=("selected_encoded_count_this_fold", "mean"),
    )
    summary["selection_frequency"] = summary.folds_selected / SELECTOR_FOLDS
    summary["stable_in_at_least_4_of_5_folds"] = summary.folds_selected >= 4
    group_keys = ["scope", "selection_method"]
    summary["unique_source_features_across_folds"] = summary.groupby(group_keys)["feature_name"].transform("nunique")
    summary.to_csv(OUTPUT / "selected_features.csv", index=False)
    selection_counts = detail.groupby(group_keys).agg(
        unique_source_features_across_folds=("feature_name", "nunique"),
    ).reset_index()
    per_fold_counts = detail.groupby(group_keys + ["fold"])["selected_encoded_feature"].nunique().groupby(level=[0, 1]).mean()
    selection_counts["mean_encoded_features_selected_per_fold"] = [float(per_fold_counts.loc[(r.scope, r.selection_method)]) for r in selection_counts.itertuples()]
    selection_counts.to_csv(OUTPUT / "feature_selection_summary.csv", index=False)
    logger.write("\nSelected features by scope/selector (one row per source feature; scores are selector-specific):")
    logger.write("Selected feature counts:")
    logger.write(selection_counts.to_string(index=False, float_format=lambda value: f"{value:.1f}"))
    logger.write(summary.sort_values(["scope", "selection_method", "selection_frequency", "mean_selection_score"],
                                    ascending=[True, True, False, False]).to_string(
                                        index=False, float_format=lambda value: f"{value:.6g}"))
    return summary, metrics


def save_confusion_outputs(scope: str, model_name: str, y_true: Any,
                           y_pred: Any) -> None:
    cm = confusion_matrix(y_true, y_pred, labels=[0, 1])
    stem = f"{scope.replace(' ', '_').lower()}__{model_name}"
    pd.DataFrame(cm, index=["Actual Alive (0)", "Actual Deceased (1)"],
                 columns=["Predicted Alive (0)", "Predicted Deceased (1)"]).to_csv(
        OUTPUT / "confusion_matrices" / f"{stem}.csv", index=True)
    fig, ax = plt.subplots(figsize=(5.2, 4.2))
    image = ax.imshow(cm, cmap="Blues")
    ax.set_xticks([0, 1], ["Alive (0)", "Deceased (1)"])
    ax.set_yticks([0, 1], ["Alive (0)", "Deceased (1)"])
    ax.set_xlabel("Predicted class")
    ax.set_ylabel("True class")
    ax.set_title(f"{scope} / {MODEL_DISPLAY[model_name]}\nHoldout confusion matrix")
    for (row, col), value in np.ndenumerate(cm):
        ax.text(col, row, str(value), ha="center", va="center", color="black", fontsize=12)
    fig.colorbar(image, ax=ax, label="Number of samples")
    fig.tight_layout()
    fig.savefig(OUTPUT / "confusion_matrices" / f"{stem}.png", dpi=300)
    plt.close(fig)


def make_evaluation_plots(test_predictions: dict[str, dict[str, Any]],
                          comparison: pd.DataFrame, logger: RunLog) -> None:
    fig, ax = plt.subplots(figsize=(7, 6))
    for scope, result in test_predictions.items():
        fpr, tpr, _ = roc_curve(result["y_true"], result["probability"])
        auc = roc_auc_score(result["y_true"], result["probability"])
        ax.plot(fpr, tpr, label=f"{scope} (AUC={auc:.3f})")
    ax.plot([0, 1], [0, 1], "k--", linewidth=1)
    ax.set(xlabel="False-positive rate", ylabel="True-positive rate",
           title="ROC curves — held-out test set (descriptive)", xlim=(0, 1), ylim=(0, 1))
    ax.legend(fontsize=8)
    fig.tight_layout()
    fig.savefig(OUTPUT / "evaluation_plots" / "roc_curves.png", dpi=300)
    fig.savefig(OUTPUT / "evaluation_plots" / "roc_curves.pdf")
    plt.close(fig)

    fig, ax = plt.subplots(figsize=(7, 6))
    for scope, result in test_predictions.items():
        precision, recall, _ = precision_recall_curve(result["y_true"], result["probability"])
        ap = average_precision_score(result["y_true"], result["probability"])
        ax.plot(recall, precision, label=f"{scope} (AP={ap:.3f})")
    prevalence = float(np.mean(next(iter(test_predictions.values()))["y_true"]))
    ax.axhline(prevalence, linestyle="--", color="black", label=f"Test prevalence={prevalence:.3f}")
    ax.set(xlabel="Recall — deceased", ylabel="Precision — deceased",
           title="Precision-recall curves — held-out test set (descriptive)", xlim=(0, 1), ylim=(0, 1))
    ax.legend(fontsize=7)
    fig.tight_layout()
    fig.savefig(OUTPUT / "evaluation_plots" / "precision_recall_curves.png", dpi=300)
    fig.savefig(OUTPUT / "evaluation_plots" / "precision_recall_curves.pdf")
    plt.close(fig)

    labels = comparison["scope"].tolist()
    x_pos = np.arange(len(labels))
    width = 0.36
    fig, ax = plt.subplots(figsize=(9, 5))
    ax.bar(x_pos - width / 2, comparison["cv_roc_auc_mean"], width,
           yerr=comparison["cv_roc_auc_std"], capsize=3, label="Nested outer CV ROC-AUC")
    ax.bar(x_pos + width / 2, comparison["test_roc_auc"], width, label="Holdout ROC-AUC (descriptive)")
    ax.set_xticks(x_pos, labels, rotation=15, ha="right")
    ax.set_ylabel("ROC-AUC")
    ax.set_ylim(0, 1)
    ax.set_title("Selected model by scope: training-CV selection and holdout evaluation")
    ax.legend(fontsize=8)
    fig.tight_layout()
    fig.savefig(OUTPUT / "evaluation_plots" / "model_comparison.png", dpi=300)
    fig.savefig(OUTPUT / "evaluation_plots" / "model_comparison.pdf")
    plt.close(fig)
    logger.write("Saved ROC, precision-recall, and model comparison figures.")


def main() -> int:
    global OUTPUT
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--data", type=Path, default=DEFAULT_DATA,
                        help=f"Exact Dataset C CSV to use (default: {DEFAULT_DATA})")
    parser.add_argument("--output", type=Path, default=OUTPUT,
                        help=f"Results folder (default: {OUTPUT})")
    args = parser.parse_args()
    source = args.data.expanduser().resolve()
    output_path = args.output.expanduser().resolve()
    OUTPUT = output_path
    for folder in [OUTPUT, OUTPUT / "confusion_matrices", OUTPUT / "evaluation_plots", OUTPUT / "trained_models"]:
        folder.mkdir(parents=True, exist_ok=True)
    logger = RunLog(OUTPUT / "execution_log.txt")
    start = datetime.now().astimezone()
    try:
        logger.write("Dataset C reproducibility run")
        logger.write(f"Start: {start.isoformat(timespec='seconds')}")
        logger.write(f"Python={sys.version.split()[0]} platform={platform.platform()}")
        logger.write(f"numpy={np.__version__}; pandas={pd.__version__}; scikit-learn={sklearn.__version__}; matplotlib={matplotlib.__version__}; joblib={joblib.__version__}")
        logger.write(f"Seed={SEED}; split=test_size={TEST_SIZE}; nested CV={OUTER_FOLDS} outer x {INNER_FOLDS} inner; selector CV={SELECTOR_FOLDS}")
        if not source.is_file():
            raise FileNotFoundError(f"The explicitly selected Dataset C file does not exist: {source}. No alternate dataset will be substituted.")
        logger.write(f"Reading exact input: {source}")
        df = pd.read_csv(source)
        x, y = audit_schema(df)
        qc_dataset(df, x, y, source, logger)

        # Use and verify the exact stratified split from the final optimization.
        train_positions, test_positions = train_test_split(
            np.arange(len(df)), test_size=TEST_SIZE, random_state=SEED, stratify=y
        )
        X_train = x.iloc[train_positions].reset_index(drop=True)
        y_train = y.iloc[train_positions].reset_index(drop=True)
        X_test = x.iloc[test_positions].reset_index(drop=True)
        y_test = y.iloc[test_positions].reset_index(drop=True)
        logger.write(f"Split: training={len(train_positions)} (alive={(y_train == 0).sum()}, deceased={(y_train == 1).sum()}); test={len(test_positions)} (alive={(y_test == 0).sum()}, deceased={(y_test == 1).sum()}).")
        split_frame = pd.DataFrame({"original_row_position": np.arange(len(df)),
                                    "split": np.where(np.isin(np.arange(len(df)), train_positions), "train", "test"),
                                    TARGET: y.to_numpy()})
        split_frame.to_csv(OUTPUT / "train_test_indices.csv", index=False)

        outer_df, inner_df, candidate_df, oof_df, scope_winners = nested_evaluation(
            X_train, y_train, list(x.columns), logger
        )
        logger.write("\nCandidate inner-CV scores (training only; fold SD is descriptive across overlapping outer training partitions):")
        candidate_display = candidate_df[["scope", "model", "inner_roc_auc_mean", "inner_roc_auc_sd_fold_scores",
                                          "inner_pr_auc_mean", "inner_balanced_accuracy_mean",
                                          "inner_macro_f1_mean", "inner_deceased_recall_mean",
                                          "mean_composite_rank", "scope_selected"]]
        logger.write(candidate_display.to_string(index=False, float_format=lambda v: f"{v:.4f}"))
        dummy_result = run_dummy_baseline(X_train, y_train, X_test, y_test, test_positions, logger)

        # Scope algorithm winners and the global winner use training inner CV
        # only. The final test scores below never enter ranking or selection.
        eligible = candidate_df[candidate_df.scope_selected.astype(bool)].copy()
        if eligible.empty:
            raise RuntimeError("No scope has a complete model candidate from training-only nested CV.")
        global_candidates: dict[str, dict[str, float]] = {}
        for _, row in eligible.iterrows():
            key = f"{row['scope']}__{row['model']}"
            global_candidates[key] = {m: float(row[f"inner_{m}_mean"]) for m in PRIORITY_METRICS}
        global_ranks = average_rank(global_candidates)
        global_winner = sorted(global_candidates,
                               key=lambda name: (global_ranks[name], -global_candidates[name]["roc_auc"],
                                                 -global_candidates[name]["pr_auc"]))[0]
        best_scope, best_model = global_winner.rsplit("__", 1)
        logger.write(f"Global winner selected by training-only composite validation rank: {best_scope} / {MODEL_DISPLAY[best_model]}")

        test_predictions: dict[str, dict[str, Any]] = {}
        comparison_rows: list[dict[str, Any]] = []
        classification_rows: list[dict[str, Any]] = []
        all_test_predictions: list[dict[str, Any]] = []
        outer_summary = outer_df.groupby("scope").agg({
            m: ["mean", "std"] for m in ["roc_auc", "pr_auc_average_precision", "balanced_accuracy",
                                          "macro_f1", "deceased_recall", "specificity_alive"]
        })
        outer_summary.columns = ["_".join(col) for col in outer_summary.columns]
        outer_summary = outer_summary.reset_index()
        for scope, model_name in scope_winners.items():
            try:
                final_pipeline = build_pipeline(scope, model_name, list(x.columns), seed=SEED)
                fit_pipeline(final_pipeline, X_train, y_train)
                predicted = final_pipeline.predict(X_test)
                probability = final_pipeline.predict_proba(X_test)[:, 1]
            except Exception as exc:
                message = f"{type(exc).__name__}: {exc}"
                logger.write(f"ERROR final fit for {scope}/{model_name}: {message}")
                logger.write(traceback.format_exc().rstrip())
                raise RuntimeError(f"Could not fit final scope model {scope}/{model_name}; see execution_log.txt") from exc
            metrics = metric_values(y_test, predicted, probability)
            metrics_row = outer_summary[outer_summary.scope == scope].iloc[0].to_dict()
            model_key = f"{scope}__{model_name}"
            cv_info = candidate_df[(candidate_df.scope == scope) & (candidate_df.model == model_name)].iloc[0]
            comparison_rows.append({
                "scope": scope, "selected_model": MODEL_DISPLAY[model_name],
                "selected_by_training_CV": True, "global_best_by_training_CV": model_key == global_winner,
                "training_composite_rank": float(global_ranks[global_winner] if model_key == global_winner else global_ranks[model_key]),
                "candidate_inner_roc_auc_mean": float(cv_info.inner_roc_auc_mean),
                "candidate_inner_pr_auc_mean": float(cv_info.inner_pr_auc_mean),
                "cv_roc_auc_mean": float(metrics_row.get("roc_auc_mean", np.nan)),
                "cv_roc_auc_std": float(metrics_row.get("roc_auc_std", np.nan)),
                "cv_pr_auc_mean": float(metrics_row.get("pr_auc_average_precision_mean", np.nan)),
                "cv_balanced_accuracy_mean": float(metrics_row.get("balanced_accuracy_mean", np.nan)),
                "cv_macro_f1_mean": float(metrics_row.get("macro_f1_mean", np.nan)),
                "cv_deceased_recall_mean": float(metrics_row.get("deceased_recall_mean", np.nan)),
                **{f"test_{k}": v for k, v in metrics.items()},
            })
            test_predictions[scope] = {"y_true": y_test.to_numpy(), "prediction": predicted,
                                       "probability": probability, "model": model_name}
            for original_position, true_value, pred_value, prob_value in zip(test_positions, y_test, predicted, probability):
                all_test_predictions.append({"scope": scope, "model": model_name,
                                             "original_row_position": int(original_position),
                                             "true_vital_status": int(true_value),
                                             "predicted_vital_status": int(pred_value),
                                             "probability_deceased": float(prob_value)})
            report = classification_report(y_test, predicted, labels=[0, 1],
                                           target_names=["Alive (0)", "Deceased (1)"],
                                           output_dict=True, zero_division=0)
            logger.write(f"\nClassification report — {scope} / {MODEL_DISPLAY[model_name]} (Deceased=positive):")
            logger.write(classification_report(y_test, predicted, labels=[0, 1],
                                               target_names=["Alive (0)", "Deceased (1)"],
                                               zero_division=0, digits=4))
            for class_name, values in report.items():
                if isinstance(values, dict):
                    classification_rows.append({"scope": scope, "model": MODEL_DISPLAY[model_name],
                                                "class_or_average": class_name,
                                                "precision": values.get("precision"),
                                                "recall": values.get("recall"),
                                                "f1_score": values.get("f1-score"),
                                                "support": values.get("support")})
            save_confusion_outputs(scope, model_name, y_test, predicted)
            model_path = OUTPUT / "trained_models" / f"{model_key}.joblib"
            joblib.dump(final_pipeline, model_path, compress=3)
            config = {
                "scope": scope, "model": model_name, "display_name": MODEL_DISPLAY[model_name],
                "seed": SEED, "training_rows": int(len(train_positions)), "test_rows": int(len(test_positions)),
                "selection": "Nested stratified 3x3 training-only CV; mean rank on ROC-AUC, PR-AUC, balanced accuracy, macro-F1, deceased recall",
                "predictor_columns": selected_columns(scope, list(x.columns)),
                "target_excluded": TARGET not in selected_columns(scope, list(x.columns)),
                "target_mapping": {str(k): v for k, v in TARGET_MAPPING.items()},
                "duplicate_CNV_aliases": "retained; no CNV deduplication in primary model",
                "preprocessing": {
                    "CNV": "most-frequent imputation; one-hot, unknown categories ignored",
                    "RNA": "median imputation and standard scaling",
                    "Mutation": "median imputation and training-fitted constant-feature variance filter",
                    "Protein": "median imputation and standard scaling",
                    "fit_boundary": "inside each training fold and final training fit only",
                },
                "estimator_parameters": final_pipeline.named_steps["model"].get_params(),
                "original_train_positions": [int(v) for v in train_positions],
                "original_test_positions": [int(v) for v in test_positions],
            }
            (OUTPUT / "trained_models" / f"{model_key}.json").write_text(
                json.dumps(config, indent=2, default=str), encoding="utf-8"
            )
            logger.write(f"Final descriptive holdout evaluation: {scope} / {MODEL_DISPLAY[model_name]}: ROC-AUC={metrics['roc_auc']:.4f}; PR-AUC={metrics['pr_auc_average_precision']:.4f}; deceased recall={metrics['deceased_recall']:.4f}; specificity={metrics['specificity_alive']:.4f}.")

        comparison = pd.DataFrame(comparison_rows)
        comparison.to_csv(OUTPUT / "model_comparison.csv", index=False)
        pd.DataFrame(classification_rows).to_csv(OUTPUT / "classification_reports.csv", index=False)
        pd.DataFrame(all_test_predictions).to_csv(OUTPUT / "heldout_predictions.csv", index=False)
        outer_summary.to_csv(OUTPUT / "nested_cv_scope_summary.csv", index=False)
        pd.DataFrame([{"scope": k, "selected_model": v, "display_name": MODEL_DISPLAY[v],
                       "global_best": f"{k}__{v}" == global_winner,
                       "training_composite_rank": float(global_ranks[f"{k}__{v}"])}
                      for k, v in scope_winners.items()]).to_csv(OUTPUT / "selected_scope_models.csv", index=False)

        selected_feature_summary, selection_metrics = run_feature_selection(X_train, y_train, list(x.columns), logger)
        make_evaluation_plots(test_predictions, comparison, logger)

        logger.write("\n=== Final model comparison (holdout metrics are descriptive) ===")
        printable = comparison[["scope", "selected_model", "cv_roc_auc_mean", "cv_roc_auc_std",
                                "cv_pr_auc_mean", "cv_balanced_accuracy_mean", "test_accuracy",
                                "test_precision_deceased", "test_deceased_recall", "test_specificity_alive",
                                "test_f1_deceased", "test_macro_f1", "test_roc_auc",
                                "test_pr_auc_average_precision", "test_tn", "test_fp", "test_fn", "test_tp"]]
        logger.write(printable.to_string(index=False, float_format=lambda v: f"{v:.4f}"))
        logger.write(f"\nBest by prespecified training-CV composite rank: {best_scope} / {MODEL_DISPLAY[best_model]}; rank={global_ranks[global_winner]:.3f}.")
        logger.write("Held-out test results were not used to choose the best model.")
        logger.write("Positive class: Deceased (1). Confusion matrix order: [Alive (0), Deceased (1)] = TN, FP, FN, TP.")

        final_summary = {
            "run_started": start.isoformat(timespec="seconds"),
            "run_finished": datetime.now().astimezone().isoformat(timespec="seconds"),
            "input_file": str(source), "input_sha256": file_sha256(source),
            "dataset_shape_including_target": list(df.shape), "predictor_count": int(x.shape[1]),
            "modality_counts": {m: int(sum(c.startswith(p) for c in x.columns)) for m, p in MODALITIES.items()},
            "target": TARGET, "target_mapping": TARGET_MAPPING,
            "target_counts": {str(k): int(v) for k, v in y.value_counts().sort_index().items()},
            "train_samples": int(len(train_positions)), "test_samples": int(len(test_positions)),
            "seed": SEED, "nested_cv": {"outer_folds": OUTER_FOLDS, "inner_folds": INNER_FOLDS},
            "candidate_models": MODELS, "scope_winners": scope_winners,
            "global_training_cv_winner": {"scope": best_scope, "model": best_model,
                                           "composite_rank": float(global_ranks[global_winner])},
            "candidate_count": int(len(candidate_df)),
            "model_fit_errors": int(len(pd.read_csv(OUTPUT / "model_fit_errors.csv"))),
            "test_results_are_descriptive_only": True,
            "patient_level_independence_verified": False,
            "survival_time_prediction": False,
            "clinically_validated": False,
            "limitations": [
                "Recorded vital status classification only; no survival duration or censoring.",
                "No patient IDs/barcodes; row-level patient independence cannot be verified.",
                "Only 94 deceased cases overall and 19 in the fixed held-out test set.",
                "No independent external validation; exploratory and not clinically validated.",
            ],
        }
        (OUTPUT / "run_metadata.json").write_text(json.dumps(final_summary, indent=2, default=str), encoding="utf-8")
        # Plain-text handover includes measured, not copied, results.
        lines = [
            "Dataset C reproducibility run complete",
            f"Input: {source}", f"Input SHA-256: {final_summary['input_sha256']}",
            f"Shape: {df.shape}; samples={len(df)}; predictors={x.shape[1]}",
            f"Target counts: Alive={int((y == 0).sum())}; Deceased={int((y == 1).sum())}",
            f"Nested CV: {OUTER_FOLDS} outer x {INNER_FOLDS} inner folds; seed={SEED}",
            f"Best by training CV: {best_scope} / {MODEL_DISPLAY[best_model]}",
            "",
            printable.to_string(index=False, float_format=lambda v: f"{v:.4f}"),
            "",
            "Scientific interpretation: exploratory classification of recorded vital.status, not survival-time prediction or clinical validation.",
            "Patient-level independence cannot be verified because patient identifiers are unavailable.",
        ]
        (OUTPUT / "final_summary.txt").write_text("\n".join(lines) + "\n", encoding="utf-8")
        logger.write("\n".join(lines))
        logger.write(f"\nOutputs written to: {OUTPUT}")
        logger.write("Completed without substituting another dataset or using test data for model selection.")
        return 0
    except Exception as exc:
        logger.write(f"FATAL {type(exc).__name__}: {exc}")
        logger.write(traceback.format_exc().rstrip())
        return 1
    finally:
        logger.close()


if __name__ == "__main__":
    raise SystemExit(main())
