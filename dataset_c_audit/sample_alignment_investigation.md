# Sample identity and row-alignment investigation

## Search findings

A search across both CSV header sets, all notebook cells and outputs, README, and both serialized model files found no patient ID, sample ID, TCGA barcode, TCGA submitter ID, assay barcode, or row-key column. The notebook has no `merge`, `join`, concat, identifier mapping, barcode parsing, or table-construction code. It simply reads `/content/data.csv`; the files are already prepared before the notebook runs.

The source paper states that the matrix contains 705 tumor samples from different patients, and says it came from TCGA via `curatedTCGAData`, a framework designed to coordinate assays and patient identifiers. This is source-level evidence that the original study assembled the assays at patient level. It does not provide the mapping from each current CSV row to a patient barcode, nor independently demonstrate the cross-layer sample correspondence in the distributed CSV.

## Direct table comparison

- `data.csv`: 705 rows, 1,937 columns.
- `brca_data_w_subtypes.csv`: 705 rows, 1,941 columns.
- The 1,937 shared columns have identical names and order.
- All 705 rows match exactly at the same positions across all shared columns.
- The enriched table appends `PR.Status`, `ER.Status`, `HER2.Final.Status`, and `histological.type`.

Thus the enriched file is an exact row-aligned extension of `data.csv`. This is direct positional evidence for these two files only. It is not an ID-based join and does not prove that each modality column belongs to a particular TCGA patient. The notebook appears to assume the prepared row matrix it reads is ready to model; it does not create or validate the integration.

## Integration decision

The four feature families are in one wide table, so they are computationally co-located. The original research paper describes these as four TCGA-BRCA modalities assembled through `curatedTCGAData`, supporting an intended patient-level integration. However, because IDs/manifests were discarded from this export, a reviewer cannot reproduce the mapping, identify unmatched samples, detect repeated patient aliquots, or verify barcode-level sample selection. Do not join this export to other TCGA data by row order. Request the original ID-bearing `brca_data_w_meta.csv`/manifest or recreate from TCGA with explicit barcode linkage before claiming independently verified integration.

## References and source URLs

1. Christensen NJ et al. *Identifying interactions in omics data for clinical biomarker discovery using symbolic regression.* Bioinformatics 38(15), 3749-3758 (2022). DOI: https://doi.org/10.1093/bioinformatics/btac405
2. Ciriello G et al.; TCGA Research Network. *Comprehensive Molecular Portraits of Invasive Lobular Breast Cancer.* Cell 163(2), 506-519 (2015). DOI: https://doi.org/10.1016/j.cell.2015.09.033
3. Original paper code/data repository: https://github.com/abzu-ai/QLattice-clinical-omics (README states the repository holds the analysis and data for the Bioinformatics paper; `data/` includes `brca_data_w_subtypes.csv` and `brca_data_w_meta.csv`; notebook `notebooks/bc_case_multi.ipynb`).
4. Kaggle dataset identified by the paper: https://www.kaggle.com/datasets/samdemharter/brca-multiomics-tcga. The Dataset C README links to a second Kaggle upload: https://www.kaggle.com/datasets/mannekuntanagendra/multi-omics-breast-cancer-prediction-qlattice.
5. Christensen et al. arXiv version, *Understanding Breast Cancer Survival: Using Causality and Language Models on Multi-omics Data*: https://arxiv.org/abs/2305.18410. Used as a secondary cross-check, not as a substitute for the original data source.
6. `curatedTCGAData` package documentation: https://bioconductor.org/packages/curatedTCGAData/ (TCGA assays and clinical data coordinated using patient/sample identifiers).
7. GDC TCGA metadata definitions: https://docs.cancergenomicscloud.org/docs/tcga-metadata (defines vital status as living/deceased; days-to-death and follow-up are distinct fields).
