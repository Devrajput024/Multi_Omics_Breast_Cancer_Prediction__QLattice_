# Patient ID and cross-file matching

## Finding

Neither CSV contains a column that can be identified from its name as a patient, case, sample, barcode, or record ID. No identifiers are available to report, and no patient-level join can be verified. No sample values are reproduced here.

## Exact comparison of the two tables

The tables each contain 705 rows. The first 1,937 headers match in the same order; every row matches exactly across these shared columns. `brca_data_w_subtypes.csv` then adds `PR.Status`, `ER.Status`, `HER2.Final.Status`, and `histological.type`. Therefore, those four labels correspond to rows in `data.csv` by position in the supplied files. This establishes a direct row-aligned duplicate/extension in this archive, not matching through patient IDs.

## Can omics layers be combined?

The feature matrix has four prefix groups (`rs_`, `cn_`, `mu_`, `pp_`) in the same row table, so those columns can technically be used together as a single aligned feature matrix. However, the archive contains no IDs or layer-specific sample manifests to prove that each block came from the same patients, and no data dictionary defining the prefixes, assays, transformations, or missingness conventions. Cross-layer patient matching cannot be independently verified from the available files. Do not join these data to external omics by row order.

## Privacy

Reports contain aggregate counts and column names only. They do not include patient/sample values or row-level identifiers.
