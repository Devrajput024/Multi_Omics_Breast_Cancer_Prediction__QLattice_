# Molecular/biomarker table audit (Dataset 3)

## File and structure

File: `dataset/dataset3/molecular_biomarker_dataset.csv` (83,235 bytes in the ZIP), read from the archive without extraction. It is a comma-separated UTF-8 CSV with 780 rows and 10 columns. No missing/blank cells or duplicate rows were observed. `Patient ID` is unique for all 780 rows. All non-ID columns are categorical text annotations; there are no raw continuous expression, sequencing, copy-number, or other assay-value columns.

| Column | Type | Distinct values | Observed values/counts |
|---|---|---:|---|
| Patient ID | identifier/text | 780 | one per row |
| Cancer Type | categorical | 1 | Breast Cancer 780 |
| Cancer Type Detailed | categorical | 5 | Invasive Ductal Carcinoma 614; Mixed Ductal and Lobular 94; Invasive Lobular 59; Invasive Mixed Mucinous 10; Breast 3 |
| Cellularity | categorical | 3 | High 394; Moderate 299; Low 87 |
| Chemotherapy | categorical | 2 | No 589; Yes 191 |
| Pam50 + Claudin-low subtype | categorical | 7 | LumA 324; LumB 194; claudin-low 79; Her2 74; Basal 62; Normal 46; NC 1 |
| Integrative Cluster | categorical | 11 | 1, 2, 3, 4ER+, 4ER-, 5, 6, 7, 8, 9, 10 |
| Oncotree Code | categorical | 5 | IDC 614; MDLC 94; ILC 59; IMMC 10; BREAST 3 |
| 3-Gene classifier subtype | categorical | 4 | ER+/HER2- Low Prolif 292; ER+/HER2- High Prolif 282; ER-/HER2- 121; HER2+ 85 |
| PR Status | categorical | 2 | Positive 405; Negative 375 |

## Interpretation

Despite the filename and dataset-level description, the available table contains precomputed categorical subtype/clinical annotations, not raw molecular measurements. PAM50/Claudin-low, integrative cluster, and 3-gene classifier values are labels; the underlying assay values, algorithms, versions, thresholds, and provenance are absent. Cancer type and subtype annotations align in counts with the clinical histology fields, but that consistency alone does not establish independent molecular assay data or provenance. See [patient_id_matching.md](patient_id_matching.md) for exact linkage to Dataset 2 and images.
