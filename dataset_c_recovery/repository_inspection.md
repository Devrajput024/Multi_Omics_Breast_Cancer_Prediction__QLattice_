# Repository inspection

Repository snapshot: `main` at commit `fc29704243010890969b524ebe3e7c2be1f7ef53`; retrieved 2026-11-10. Files were saved under this folder. Raw GitHub URLs at that pinned revision are used below.

## Scope and recovered files

Inspected the original repository README, `data/brca_data_w_meta.csv`, `data/brca_data_w_subtypes.csv`, and `notebooks/bc_case_multi.ipynb`. Downloaded copies and hashes:

| File | Bytes | SHA-256 |
|---|---:|---|
| `README.md` | 5,257 | `06a69f22e996c47d9564d4fd70d43f9bbad8dc7798ef97ac84e3e7eea8533a63` |
| `data__brca_data_w_meta.csv` | 10,809,133 | `bd5f42e878359afc0a436fe2bf3aac83f405d89d7277eea4072fce92227c91ad` |
| `data__brca_data_w_subtypes.csv` | 10,787,813 | `da120be71091411cda4185fc45f3c204c4d244a4060d8c7b35f42bbd4c3f69d5` |
| `notebooks__bc_case_multi.ipynb` | 440,140 | `0807b9ea23899ac4bdb2cb6cc71abd983f9298dd04dbb86919480cca24887a8f` |

## Findings

- `brca_data_w_meta.csv` has 705 rows and 1,938 columns: 1,937 columns shared with Dataset C `data.csv`, plus `patient.clinical.cqcf.histological.type`.
- The downloaded repository `brca_data_w_subtypes.csv` has 705 rows and 1,937 columns and is byte-for-byte identical to Dataset C `data.csv`. Despite its filename, it has no subtype or histology fields.
- The metadata CSV has histology and `vital.status`; neither it nor the other CSV has an explicit patient ID, sample ID, TCGA barcode, or UUID column. The word “patient” in the histology column name does not make that column an identifier.
- The notebook loads the prepared metadata CSV (`pd.read_csv('../data/brca_data_w_meta.csv')`) and sets `target = vital.status`. Its relevant code does not construct the table, retrieve raw assays, map identifiers, or merge assays and clinical metadata.
- The README identifies the publication as Christensen et al., *Identifying interactions in omics data for clinical biomarker discovery using symbolic regression*, *Bioinformatics* (2022), DOI `10.1093/bioinformatics/btac405`. It says the paper's results were generated with Feyn 2.1.2 and a now-deprecated remote Abzu server. The README's separate AD proteomics tutorial is not Dataset C documentation.
- The 2022 paper reports TCGA breast cancer data accessed through `curatedTCGAData`; this establishes the stated source cohort for the research data. It does not restore sample-level identifiers in these derivative CSVs.

## Sources

- Repository: https://github.com/abzu-ai/QLattice-clinical-omics
- Pinned README: https://github.com/abzu-ai/QLattice-clinical-omics/blob/fc29704243010890969b524ebe3e7c2be1f7ef53/README.md
- Pinned metadata CSV: https://github.com/abzu-ai/QLattice-clinical-omics/blob/fc29704243010890969b524ebe3e7c2be1f7ef53/data/brca_data_w_meta.csv
- Pinned subtype-named CSV: https://github.com/abzu-ai/QLattice-clinical-omics/blob/fc29704243010890969b524ebe3e7c2be1f7ef53/data/brca_data_w_subtypes.csv
- Pinned notebook: https://github.com/abzu-ai/QLattice-clinical-omics/blob/fc29704243010890969b524ebe3e7c2be1f7ef53/notebooks/bc_case_multi.ipynb
- Christensen et al. (2022): https://doi.org/10.1093/bioinformatics/btac405
- `curatedTCGAData`: https://bioconductor.org/packages/curatedTCGAData
