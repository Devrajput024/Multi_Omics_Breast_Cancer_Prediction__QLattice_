# Clinical table audit (Dataset 2)

## File and structure

File: `dataset/dataset2/patient_history_dataset.csv` (142,053 bytes in the ZIP), read directly without extraction. It is a comma-separated UTF-8 CSV with 780 data rows and 25 columns. There are 780 unique Patient IDs, no duplicate rows, and no missing/blank cells in the 25 columns. `Patient ID` is unique on every row. Column type below is inferred from observed values (numeric means parseable as numbers; categorical means text labels):

| Column | Inferred type | Distinct values |
|---|---|---:|
| Patient ID | identifier/text | 780 |
| Age at Diagnosis | numeric | 720 |
| Type of Breast Surgery | categorical | 2 |
| Cohort | numeric-coded | 3 |
| ER status measured by IHC | categorical | 2 |
| ER Status | categorical | 2 |
| Neoplasm Histologic Grade | numeric-coded | 3 |
| HER2 status measured by SNP6 | categorical | 3 |
| HER2 Status | categorical | 2 |
| Tumor Other Histologic Subtype | categorical | 7 |
| Hormone Therapy | categorical | 2 |
| Inferred Menopausal State | categorical | 2 |
| Primary Tumor Laterality | categorical | 2 |
| Lymph nodes examined positive | numeric | 25 |
| Tumor Size | numeric | 61 |
| Tumor Stage | numeric-coded | 4 |
| Mutation Count | numeric | 26 |
| Nottingham prognostic index | numeric | 178 |
| Overall Survival (Months) | numeric | 739 |
| Overall Survival Status | categorical | 2 |
| Relapse Free Status (Months) | numeric | 743 |
| Relapse Free Status | categorical | 2 |
| Sex | categorical | 1 |
| Patient's Vital Status | categorical | 3 |
| class | categorical | 3 |

The numerical fields span: age 21.93–96.29; cohort 1–3; histologic grade 1–3; positive nodes 0–41; tumor size 1–180; tumor stage 1–4; mutation count 1–46; Nottingham index 2.002–6.36; overall-survival time 1.233333333–337.0333333; relapse-free time 1.22–298.88. The file does not declare units for age, tumor size, mutation count, or other numeric fields except the two survival column names, which say “Months.”

## Categories and internal checks

- `class`: benign 437, malignant 210, normal 133.
- `Overall Survival Status`: Deceased 407, Living 373.
- `Patient's Vital Status`: Living 373, Died of Disease 256, Died of Other Causes 151. The `Living` survival-status count corresponds exactly to `Living` vital-status values; all Deceased records split into the two deceased vital-status categories.
- ER Status: Positive 613, Negative 167. The IHC field has a source spelling anomaly: `Positve` (614), `Negative` (166); it does not exactly agree with the separate ER Status field for two records.
- HER2 SNP6 field: Neutral 566, Gain 173, Loss 41; HER2 Status: Negative 685, Positive 95.
- Histologic subtype: Ductal/NST 594, Mixed 94, Lobular 59, Tubular/ cribriform 12, Mucinous 10, Medullary 8, Other 3.
- Surgery: Mastectomy 483, Breast Conserving 297. Hormone Therapy: Yes 486, No 294. Menopausal state: Post 589, Pre 191. Laterality: Left 411, Right 369. Sex is Female for all 780 rows.
- Recurrence status: Not Recurred 461, Recurred 319. See [survival_variable_audit.md](survival_variable_audit.md) for cross-tabulation, endpoint limits, and units evidence.

`Cohort`, grade, stage, and mutation count are stored as numbers; their exact coding/unit definitions are undocumented in the ZIP. Names and value labels should not be treated as a complete codebook.
