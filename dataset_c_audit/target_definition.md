# Target definition: `vital.status`

## Verified meaning and coding

The peer-reviewed QLattice study explicitly defines the 705-row breast-cancer target as vital status: code **0 = survived/alive** (611 patients) and code **1 = fatal/deceased** (94 patients). The provided CSV reproduces these counts exactly. This is a binary all-cause vital-status outcome as described by the source; it is not documented as breast-cancer-specific mortality. TCGA/GDC defines vital status as living versus deceased; the event timing is separate.

| Code | Meaning | Count | Share |
|---:|---|---:|---:|
| 0 | Survived/alive at recorded status | 611 | 86.67% |
| 1 | Fatal/deceased | 94 | 13.33% |

## Is it classification or survival analysis?

The notebook deliberately calls QLattice with `output_name="vital.status"` and `kind="classification"`, so a binary classification task is defined and technically possible. However, a status-only endpoint discards time-to-event information. The supplied CSVs contain no days-to-death, days-to-last-follow-up, censoring time, diagnosis date, or status-update date. Survival analysis cannot be run from these two CSVs alone. If an authoritative ID mapping and event/follow-up durations can be recovered from the TCGA clinical source, a time-to-event endpoint would generally preserve more information and handle censoring; otherwise, status classification must be described as status-at-recorded-follow-up, with no fixed prediction horizon.

## Leakage assessment

- `vital.status` is included in `data.csv` as a column. It must never be included among predictors when it is the output. The notebook passes `output_name="vital.status"` to Feyn, which is intended to reserve that column as the outcome. A downstream custom pipeline should explicitly remove it from X.
- The notebook splits rows before QLattice training and stratifies by the target. It contains no patient IDs to support a patient-group split, but the published source says rows represent different patients.
- The model artifacts inspected have `vital.status` as output and their listed input features do not include `vital.status`.
- The four receptor/histology outcomes are present only in `brca_data_w_subtypes.csv`; the notebook loads `data.csv`, which omits those columns. If any of those outcomes is later predicted, do not use the other outcome columns as predictors without a clearly justified prediction time; this risk is not relevant to the notebook's `data.csv` input as written.
- Measurement times and outcome assessment time are absent, so temporal leakage from follow-up or treatment variables cannot be fully assessed. The model inputs shown are molecular feature columns; no age, stage, treatment, or survival-duration covariates are included in `data.csv`.

## Decision

The label mapping is verified and the current task is binary classification. It is not a survival-time analysis dataset in its current exported form. A more scientifically informative survival analysis requires re-linkable TCGA clinical follow-up times and censoring indicators, which are not available in this ZIP.

## References and source URLs

1. Christensen NJ et al. *Identifying interactions in omics data for clinical biomarker discovery using symbolic regression.* Bioinformatics 38(15), 3749-3758 (2022). DOI: https://doi.org/10.1093/bioinformatics/btac405
2. Ciriello G et al.; TCGA Research Network. *Comprehensive Molecular Portraits of Invasive Lobular Breast Cancer.* Cell 163(2), 506-519 (2015). DOI: https://doi.org/10.1016/j.cell.2015.09.033
3. Original paper code/data repository: https://github.com/abzu-ai/QLattice-clinical-omics (README states the repository holds the analysis and data for the Bioinformatics paper; `data/` includes `brca_data_w_subtypes.csv` and `brca_data_w_meta.csv`; notebook `notebooks/bc_case_multi.ipynb`).
4. Kaggle dataset identified by the paper: https://www.kaggle.com/datasets/samdemharter/brca-multiomics-tcga. The Dataset C README links to a second Kaggle upload: https://www.kaggle.com/datasets/mannekuntanagendra/multi-omics-breast-cancer-prediction-qlattice.
5. Christensen et al. arXiv version, *Understanding Breast Cancer Survival: Using Causality and Language Models on Multi-omics Data*: https://arxiv.org/abs/2305.18410. Used as a secondary cross-check, not as a substitute for the original data source.
6. `curatedTCGAData` package documentation: https://bioconductor.org/packages/curatedTCGAData/ (TCGA assays and clinical data coordinated using patient/sample identifiers).
7. GDC TCGA metadata definitions: https://docs.cancergenomicscloud.org/docs/tcga-metadata (defines vital status as living/deceased; days-to-death and follow-up are distinct fields).
