# Second-stage conclusions and readiness decision

## Decision: B. Ready only after specific corrections

**Why B rather than A:** scientific provenance and label coding are now supported by the peer-reviewed QLattice paper and the original source repository: this is TCGA-BRCA-derived data, the four prefix families have documented modality meanings, and `vital.status` is coded 0=alive/survived and 1=deceased/fatal. Those core unknowns are resolved. However, the exported files remain insufficient for a defensible model-development run until key integrity and design details are corrected or recovered.

### Required before modeling

1. Recover the original ID-bearing metadata/manifest or recreate the subset from TCGA through `curatedTCGAData`. The current two tables have no IDs, barcodes, release/version information, or follow-up times. They align exactly to each other by row position, but patient-level links and aliquot choices cannot be audited from this export.
2. Resolve the target horizon and obtain days-to-death/last-follow-up and censoring fields if the scientific question is prognosis. If only binary alive/dead classification is pursued, describe it as recorded vital status without a fixed horizon or cancer-specific interpretation.
3. Resolve the histology discrepancy: the local CSV has 574 ductal and 131 lobular values; a related paper describing the 705-case cohort reports 490 IDC, 127 ILC, and 88 mixed. Obtain the exact data version/label transformation. Do not silently recode one to fit the other.
4. Confirm the exact assay preprocessing and original feature dictionary, including GISTIC version/cutoffs, mutation-caller encoding, RNA-seq normalization, protein assay/normalization, and how repeated identical CN vectors map to genomic segments.
5. Treat the two model files as unverified historical artifacts. The publication repository documents Feyn 2.1.2 with a deprecated remote service, while the supplied notebook switches to Feyn 3.4.1. Reproduce and validate after source data IDs, target, preprocessing, version compatibility and evaluation splits have been established; do not rely on unlinked notebook narrative metrics.

### Positive integrity evidence

- 705 rows and 1,936 documented molecular features plus target are consistent with the QLattice article's Dataset C description.
- The two CSVs share 1,937 columns and match exactly across all 705 row positions; the enriched table only adds four clinical/pathology fields.
- 62 exact duplicate-vector groups (99 extra columns) are reproducible, all within GISTIC-style `cn_` columns; this points to a feature-level redundancy that can be handled transparently after verifying segment mappings.
- The saved QLattice models have binary `vital.status` outputs, use actual available features and do not include the target as an input.

This decision does not authorize model training in this stage; none was performed.

## References and source URLs

1. Christensen NJ et al. *Identifying interactions in omics data for clinical biomarker discovery using symbolic regression.* Bioinformatics 38(15), 3749-3758 (2022). DOI: https://doi.org/10.1093/bioinformatics/btac405
2. Ciriello G et al.; TCGA Research Network. *Comprehensive Molecular Portraits of Invasive Lobular Breast Cancer.* Cell 163(2), 506-519 (2015). DOI: https://doi.org/10.1016/j.cell.2015.09.033
3. Original paper code/data repository: https://github.com/abzu-ai/QLattice-clinical-omics (README states the repository holds the analysis and data for the Bioinformatics paper; `data/` includes `brca_data_w_subtypes.csv` and `brca_data_w_meta.csv`; notebook `notebooks/bc_case_multi.ipynb`).
4. Kaggle dataset identified by the paper: https://www.kaggle.com/datasets/samdemharter/brca-multiomics-tcga. The Dataset C README links to a second Kaggle upload: https://www.kaggle.com/datasets/mannekuntanagendra/multi-omics-breast-cancer-prediction-qlattice.
5. Christensen et al. arXiv version, *Understanding Breast Cancer Survival: Using Causality and Language Models on Multi-omics Data*: https://arxiv.org/abs/2305.18410. Used as a secondary cross-check, not as a substitute for the original data source.
6. `curatedTCGAData` package documentation: https://bioconductor.org/packages/curatedTCGAData/ (TCGA assays and clinical data coordinated using patient/sample identifiers).
7. GDC TCGA metadata definitions: https://docs.cancergenomicscloud.org/docs/tcga-metadata (defines vital status as living/deceased; days-to-death and follow-up are distinct fields).
