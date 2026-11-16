# Dataset overview

## Archive inventory

The inspected ZIP contains six files: two CSV datasets, one Jupyter notebook, two `.model` QLattice artifacts, and one README. It contains no Excel workbook, independent clinical metadata table, separate layer-specific CSVs, or standalone Python script. The archive was inspected without running notebook code or models.

## Dataset tables

| File | Rows / samples | Columns | Numeric columns | Categorical columns | Duplicate rows | Exact duplicate columns (by values) | Missing cells |
|---|---:|---:|---:|---:|---:|---:|---:|
| `brca_data_w_subtypes.csv` | 705 | 1941 | 1937 | 4 | 0 | 62 groups / 99 extra columns | 389 |
| `data.csv` | 705 | 1937 | 1937 | 0 | 0 | 62 groups / 99 extra columns | 0 |

`brca_data_w_subtypes.csv` contains 1,937 feature/target columns shared with `data.csv`, plus four columns: `PR.Status`, `ER.Status`, `HER2.Final.Status`, and `histological.type`. It has 705 rows and 1,941 columns. `data.csv` has 705 rows and 1,937 columns.

## Feature blocks observed

The shared columns contain 604 `rs_` variables, 860 `cn_` variables, 249 `mu_` variables, 223 `pp_` variables, and `vital.status`. These labels are inferred from column-name prefixes; the files provide no data dictionary to confirm assay details or units. All shared feature/target columns are numeric. The four added columns in the subtype-enriched table are categorical.

## Target and candidate tasks

Both tables contain `vital.status` with values `0` (611/705, 86.67%) and `1` (94/705, 13.33%). The table itself does not define what these codes mean, event timing, censoring, or follow-up. The notebook code explicitly uses `vital.status` as the output for a QLattice binary classification run; that establishes the intended notebook task, not the clinical interpretation or validity of the target. A provisional binary prediction task is possible if the target definition, time point, and data leakage are clarified.

The additional columns in `brca_data_w_subtypes.csv` are receptor-status markers and histological type, not a single breast-cancer molecular subtype field. Receptor status categories include missing/indeterminate/not-performed states. There is no explicit survival time, recurrence, prognosis endpoint, or named subtype label. `vital.status` alone cannot support survival analysis.

## IDs and provenance

No patient/sample identifier column was found in either CSV. The two tables have exactly matching rows in exactly matching order for all 1,937 shared columns; their four additional clinical columns therefore align by row position. This is reproducible table alignment, but it is not a patient-ID join and does not establish patient-level linkage to an external source. The README gives a Kaggle dataset title/link and says the source is ?Kaggle / Computational Oncology (likely derived from TCGA or similar)?; the parenthetical is explicitly uncertain. No original cohort citation, TCGA project code, accession, codebook, or record-level provenance is present in the archive. TCGA provenance is unverified.

## Other project artifacts

The notebook imports `feyn`, loads `data.csv`, splits on `vital.status`, and invokes `QLattice.auto_run`. Notebook narrative reports prior scores and model interpretations; they have not been independently reproduced. Both `.model` files are JSON-encoded QLattice artifacts with `vital.status` outputs and feature names including `mu_TNXB`, `rs_KRT23`, `rs_APOB`, `cn_ANKRD30B`, and `rs_VWDE`. The artifacts support that QLattice models were saved, not that their performance is valid. See `DATASET_ISSUES_AND_RECOMMENDATIONS.md`.
