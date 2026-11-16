# Dataset issues and recommendations

## Observed issues

1. **No patient/sample IDs.** No record linkage is possible across external data or independently verified across feature blocks. The two supplied tables align exactly by row, but there is no ID-based match.
2. **Unverified provenance.** The README links to a Kaggle dataset and labels TCGA origin as ?likely?; no primary source citation, cohort/accession, codebook, or metadata confirms it.
3. **Undocumented target.** `vital.status` is 0/1 with 611/94 counts, but the meaning of the encoding, event definition, follow-up horizon, and censoring are not given. The notebook treats it as binary classification.
4. **Subtype claims exceed the columns.** There is no explicit molecular subtype variable. `PR.Status`, `ER.Status`, `HER2.Final.Status`, and `histological.type` are not interchangeable with a validated subtype label. Receptor fields include absent and non-result categories.
5. **Prefix-based assay interpretation only.** There are 604 `rs_`, 860 `cn_`, 249 `mu_`, and 223 `pp_` columns. Prefixes suggest expression, copy-number, mutation, and protein/phosphoprotein blocks, but no schema document confirms them. In particular, `rs_` does not by itself establish SNP data.
6. **Substantial label missingness in subtype-enriched table.** `PR.Status` and `ER.Status` each have 122 blank cells; `HER2.Final.Status` has 145. There are no blank cells in `data.csv`. Other non-result categories are retained values, not blanks.
7. **Duplicate-valued columns.** Each table has 62 groups of columns with identical value vectors (99 extra columns beyond one per group). These often occur in `cn_` variables. This is an observed redundancy; do not drop them without investigating provenance and feature semantics.
8. **No duplicate rows.** Each table has 705 distinct rows.
9. **Small, imbalanced target.** The positive-coded class is 94/705 (13.33%). Accuracy alone would be misleading; any future evaluation should use stratified validation and report class-sensitive metrics with uncertainty.
10. **Saved experiment evidence is not self-validating.** The notebook includes narrative scores, but no reproducible environment/locked dependencies or saved split predictions. The notebook uses an older single train/test split, so its narrative performance claims should not be treated as confirmed.

## Recommendations before modeling

1. Obtain the original source citation, data dictionary, target codebook, event/follow-up definition, preprocessing details, and identifier mapping from the dataset publisher. Confirm whether each row is one patient, one tumor, or another unit.
2. Confirm the clinical meaning of `vital.status` and whether prediction is intended at diagnosis or another time point. If the target is survival status at variable follow-up, formulate a proper time-to-event analysis only when time and censoring data are available.
3. Confirm the four feature blocks and their assay platforms/units. Obtain patient-level IDs or an authoritative sample manifest before claiming multi-omics integration or joining to new data.
4. For an eventual classification task, define a leakage-safe prediction time and split by patient (and by cohort/site where relevant). Use nested or repeated stratified validation appropriate to the small positive class; report ROC-AUC and PR-AUC, sensitivity, specificity, calibration, and uncertainty.
5. Treat `brca_data_w_subtypes.csv` as an extension of `data.csv`, not an independently matched cohort. Decide whether receptor/histology columns are targets or predictors only after defining the clinical question; they may leak outcome-related information.
6. Resolve duplicate-valued feature columns only after checking the source schema. Preserve an audit trail for any later preprocessing.
7. Do not combine this project with the separate image-classification project at this stage.

## Suggested next modeling step

The next step is **not model fitting yet**: establish provenance and target semantics, obtain identifiers/data dictionary, and verify that each feature block is patient-aligned. Once resolved, write a prespecified analysis plan for one outcome (the notebook?s apparent candidate is binary `vital.status`) and an independent validation strategy. Then build a reproducible baseline on the verified cohort before any QLattice comparison.
