# Dataset C provenance investigation

## Conclusion

Dataset C is a derived, reduced feature table from **TCGA Breast Cancer (TCGA-BRCA)**, not an independently documented cohort. This conclusion is supported by the attached README's Kaggle link/title, the matching feature schema, and, most decisively, the peer-reviewed QLattice paper whose multi-omics case states the data were obtained from Ciriello et al. (2015) and extracted from TCGA through `curatedTCGAData`. The paper describes 705 breast-tumor samples from different patients.

The likely publication behind the data is Ciriello et al. (2015), *Comprehensive Molecular Portraits of Invasive Lobular Breast Cancer*. The data in this ZIP is a 705-row selected/preprocessed subset, not the full TCGA BRCA release or the entire 817-tumor collection in the Ciriello publication. The QLattice paper reports that raw input was variance-thresholded per data type. It does not identify the exact TCGA aliquot/barcode or patient IDs in this exported CSV.

## Direct evidence inspected

- The Dataset C README links to the Kaggle title "Multi-Omics Breast Cancer Prediction QLattice" and says TCGA derivation is "likely"; that README itself is not sufficient evidence of provenance.
- The Kaggle dataset title used by the original paper is *BRCA Multi-Omics (TCGA)*. The published QLattice paper names Ciriello et al. (2015), TCGA, and the `curatedTCGAData` package as the source and names the four data types. The paper's public companion repository identifies itself as the code/data for that article and includes a `brca_data_w_subtypes.csv` file.
- This Dataset C CSV has 705 rows and the four feature families/counts reported by the paper: 604 expression, 249 mutation, 860 copy-number, and 223 protein columns (1,936 features total), plus `vital.status`.
- The source paper's description is not perfectly consistent with this CSV's subtype labels: a related paper citing the same Kaggle data reports 490 IDC, 127 ILC, and 88 mixed cases, whereas the supplied `brca_data_w_subtypes.csv` has only two observed `histological.type` values: 574 infiltrating ductal and 131 infiltrating lobular. The paper counts should not be substituted for the observed row counts. Whether this difference reflects a different data version or label recoding is unresolved.
- No patient demographics, clinical identifiers, TCGA barcodes, accession/version manifest, or follow-up time is present in the two supplied CSVs. The population can be identified as a TCGA-BRCA tumor subset; age, race, treatment, site, and patient-level cohort membership cannot be reconstructed from this ZIP.

## Provenance limits

The scientific source establishes TCGA-BRCA derivation and a 705-patient tumor cohort, but not the exact data release, barcode level, patient-specific row mapping, preprocessing script/version, or the reason the supplied histology counts differ from the related paper's breakdown. The README's generic references to expression/mutation/clinical files are inaccurate for the archive contents: the ZIP actually has two wide CSV matrices, one notebook, two models, and the README.

## References and source URLs

1. Christensen NJ et al. *Identifying interactions in omics data for clinical biomarker discovery using symbolic regression.* Bioinformatics 38(15), 3749-3758 (2022). DOI: https://doi.org/10.1093/bioinformatics/btac405
2. Ciriello G et al.; TCGA Research Network. *Comprehensive Molecular Portraits of Invasive Lobular Breast Cancer.* Cell 163(2), 506-519 (2015). DOI: https://doi.org/10.1016/j.cell.2015.09.033
3. Original paper code/data repository: https://github.com/abzu-ai/QLattice-clinical-omics (README states the repository holds the analysis and data for the Bioinformatics paper; `data/` includes `brca_data_w_subtypes.csv` and `brca_data_w_meta.csv`; notebook `notebooks/bc_case_multi.ipynb`).
4. Kaggle dataset identified by the paper: https://www.kaggle.com/datasets/samdemharter/brca-multiomics-tcga. The Dataset C README links to a second Kaggle upload: https://www.kaggle.com/datasets/mannekuntanagendra/multi-omics-breast-cancer-prediction-qlattice.
5. Christensen et al. arXiv version, *Understanding Breast Cancer Survival: Using Causality and Language Models on Multi-omics Data*: https://arxiv.org/abs/2305.18410. Used as a secondary cross-check, not as a substitute for the original data source.
6. `curatedTCGAData` package documentation: https://bioconductor.org/packages/curatedTCGAData/ (TCGA assays and clinical data coordinated using patient/sample identifiers).
7. GDC TCGA metadata definitions: https://docs.cancergenomicscloud.org/docs/tcga-metadata (defines vital status as living/deceased; days-to-death and follow-up are distinct fields).
