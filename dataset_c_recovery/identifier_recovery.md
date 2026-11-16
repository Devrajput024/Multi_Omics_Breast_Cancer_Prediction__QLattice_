# Identifier recovery and row comparison

Repository snapshot: `main` at commit `fc29704243010890969b524ebe3e7c2be1f7ef53`; retrieved 2026-11-10. Files were saved under this folder. Raw GitHub URLs at that pinned revision are used below.

## Identifier search result

No patient/sample identifier field was found in either recovered CSV header or the investigated notebook. Specifically, no TCGA barcode, case UUID, sample UUID, or stable patient key is present. The histology field `patient.clinical.cqcf.histological.type` is a phenotype label, not an identifier. No identifiers were reconstructed or guessed.

## Reproducible content comparison

Compared row signatures using SHA-256 over the ordered values in the 1,937 columns shared by each local file and repository `brca_data_w_meta.csv`. This is a content equality check, not a patient-ID match and not a row-position join.

| Local file | Rows in local/repository | Unique signatures | Matched signatures | Unmatched | Ambiguous |
|---|---:|---:|---:|---:|---:|
| Dataset C `data.csv` | 705 / 705 | 705 / 705 | 705 | 0 | 0 |
| Dataset C `brca_data_w_subtypes.csv` | 705 / 705 | 705 / 705 | 705 | 0 | 0 |

The 1,937 shared columns comprise the molecular feature columns plus `vital.status`. The repository metadata CSV adds exactly one histology column. The repository `brca_data_w_subtypes.csv` matches local `data.csv` byte-for-byte. The local subtype-enriched file adds `ER.Status`, `PR.Status`, `HER2.Final.Status`, and `histological.type` compared with `data.csv`.

There are zero verified identifier matches because the source files provide no identifiers. The 705 one-to-one content matches establish that the same rows/profiles occur in the files, but cannot map any row to a named TCGA case. There were no duplicate molecular-row signatures in the compared sets. Duplicate identifier counts are not applicable because identifiers are absent.

## Exact reconstruction status

The feature/target row content can be matched exactly to the repository's metadata version by content signature. The local histology label disagrees for one of those rows (see `histology_discrepancy.md`). The cohort cannot be reconstructed as an identifier-verified TCGA patient list, nor can a patient-level clinical join be performed from these sources.

## Sources and retrieval

Repository snapshot and source links: [repository inspection](repository_inspection.md). Local and downloaded files remain unmodified. Retrieval date: 2026-11-10.
