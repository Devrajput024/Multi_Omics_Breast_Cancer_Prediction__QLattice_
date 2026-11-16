# Dataset C feature dictionary and biological interpretation

The prefix semantics below are explicitly defined in the peer-reviewed QLattice article and match the uploaded schema and observed value types. Each feature is a gene-level variable measured per tumor/patient row. Column values are not raw patient IDs; they are molecular feature values. The exact assay pipeline, units, normalization and aliquot provenance are not recorded in the supplied files.

| Prefix | Features | Supported modality | What the variable represents | Observed values in supplied data | Evidence limits |
|---|---:|---|---|---|---|
| `cn_` | 860 | Copy-number variation | Per-gene copy-number alteration categories derived from GISTIC scores | Integer-like categories among -2, -1, 0, 1, 2; each column has 4 or 5 observed states | GISTIC score categories are documented; exact version/cutoffs and gene-to-segment mapping are not included |
| `mu_` | 249 | Somatic mutation status | Binary indicator of whether the named gene has a mutation | All 249 columns have two observed values, 0/1 | Paper documents mutated/not-mutated categorical encoding; variant type, allele, pathogenicity, and mutation caller are absent |
| `rs_` | 604 | Gene expression / transcriptomics | RNA-sequencing expression for gene symbols | Continuous numeric values | Article documents RNA-seq; units, normalization method, gene-symbol version, and exact preprocessing are not specified in archive |
| `pp_` | 223 | Protein expression / proteomics | Protein-level abundance/measurement for gene/protein symbols | Continuous numeric values | Protein expression is documented; platform, units, normalization and assay batch metadata are not available in ZIP |

The source paper says genes are represented by HUGO symbols and reports variance-threshold preprocessing per modality. `rs_` therefore means expression in this documented dataset, not SNP genotype; `mu_` is somatic mutation status. Names such as `pp_p70S6K.pT389` indicate a protein/phospho-protein-like feature label, but the archive lacks assay documentation to establish its exact analyte or measurement protocol.

## Clinical/phenotype columns

`vital.status` is the verified binary vital-status outcome (0 alive/survived, 1 deceased/fatal). `brca_data_w_subtypes.csv` additionally contains `PR.Status`, `ER.Status`, `HER2.Final.Status`, and `histological.type`. These are clinical/pathology labels; the exact coding conventions and missing/non-result states should be retained as present rather than recoded without the original codebook. No age, race, sex, stage, treatment, diagnosis date, patient identifier, or survival duration columns are present in either supplied CSV.

## Patient-level versus gene-level

Rows represent tumor samples from different breast-cancer patients according to the source paper; feature columns are gene-level summaries/measurements for each row. The files have no independent sample IDs, so this row unit is literature-verified at cohort level but cannot be linked to a named TCGA patient or aliquot from this export.

## References and source URLs

1. Christensen NJ et al. *Identifying interactions in omics data for clinical biomarker discovery using symbolic regression.* Bioinformatics 38(15), 3749-3758 (2022). DOI: https://doi.org/10.1093/bioinformatics/btac405
2. Ciriello G et al.; TCGA Research Network. *Comprehensive Molecular Portraits of Invasive Lobular Breast Cancer.* Cell 163(2), 506-519 (2015). DOI: https://doi.org/10.1016/j.cell.2015.09.033
3. Original paper code/data repository: https://github.com/abzu-ai/QLattice-clinical-omics (README states the repository holds the analysis and data for the Bioinformatics paper; `data/` includes `brca_data_w_subtypes.csv` and `brca_data_w_meta.csv`; notebook `notebooks/bc_case_multi.ipynb`).
4. Kaggle dataset identified by the paper: https://www.kaggle.com/datasets/samdemharter/brca-multiomics-tcga. The Dataset C README links to a second Kaggle upload: https://www.kaggle.com/datasets/mannekuntanagendra/multi-omics-breast-cancer-prediction-qlattice.
5. Christensen et al. arXiv version, *Understanding Breast Cancer Survival: Using Causality and Language Models on Multi-omics Data*: https://arxiv.org/abs/2305.18410. Used as a secondary cross-check, not as a substitute for the original data source.
6. `curatedTCGAData` package documentation: https://bioconductor.org/packages/curatedTCGAData/ (TCGA assays and clinical data coordinated using patient/sample identifiers).
7. GDC TCGA metadata definitions: https://docs.cancergenomicscloud.org/docs/tcga-metadata (defines vital status as living/deceased; days-to-death and follow-up are distinct fields).
