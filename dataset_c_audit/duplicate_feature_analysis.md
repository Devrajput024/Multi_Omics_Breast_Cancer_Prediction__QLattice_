# Exact duplicate feature analysis

## Reproduction

I compared each column's complete 705-value vector in `data.csv` (exact string values, row order preserved) and repeated the check in `brca_data_w_subtypes.csv`. Both yield **62 groups of identical columns**, totaling **99 additional columns** beyond the first representative in each group. All 62 groups are `cn_` variables; no `rs_`, `mu_`, `pp_`, or outcome column participates. The same groups occur unchanged in both CSVs. There are no duplicate rows.

The CN values are categorical GISTIC-style states. These exact-vector duplicates likely occur because multiple genes share the same copy-number segment/state pattern across this cohort. This explanation is biologically and structurally plausible, but the export contains no genomic coordinates or preprocessing script to prove the precise cause. It is not a repeated import of the same named column: the column names differ. The notebook does not perform the feature construction or merge that created them.

## Groups and candidate redundant columns

The first column listed in each row is the representative retained by original header order; remaining names are the 99 redundant columns if the aim is one numeric vector per unique pattern. Preserve the aliases because they may represent distinct genes covered by a common genomic event.

| Representative | Identical-valued aliases (candidate redundant columns) | Group size |
|---|---|---:|
| `cn_ISG15` | `cn_PLCH2`, `cn_SAMD11`, `cn_TNFRSF18` | 4 |
| `cn_KCNQ4` | `cn_CITED4` | 2 |
| `cn_GSTM1` | `cn_GSTM2` | 2 |
| `cn_FMO5` | `cn_PDZK1P1`, `cn_PDZK1` | 3 |
| `cn_CXCL5` | `cn_CXCL2` | 2 |
| `cn_CXCL9` | `cn_ART3`, `cn_CXCL10`, `cn_CXCL11` | 4 |
| `cn_LRRC14B` | `cn_PLEKHG4B` | 2 |
| `cn_C5orf38` | `cn_IRX2` | 2 |
| `cn_CDC20B` | `cn_CCNO` | 2 |
| `cn_PCDHA1` | `cn_PCDHA3`, `cn_PCDHA4`, `cn_PCDHA10`, `cn_PCDHA11`, `cn_PCDHA5`, `cn_PCDHA6` | 7 |
| `cn_PCDHA12` | `cn_PCDHA13` | 2 |
| `cn_DACT2` | `cn_FRMD1` | 2 |
| `cn_ARC` | `cn_LY6D` | 2 |
| `cn_PSCA` | `cn_SLURP1` | 2 |
| `cn_IL33` | `cn_GLDC` | 2 |
| `cn_DNAI1` | `cn_ENHO` | 2 |
| `cn_CCL19` | `cn_CCL21` | 2 |
| `cn_AKR1E2` | `cn_AKR1C1`, `cn_AKR1C2` | 3 |
| `cn_PLAC9` | `cn_SFTPA2` | 2 |
| `cn_CTSW` | `cn_SPDYC`, `cn_ZNHIT2` | 3 |
| `cn_ACY3` | `cn_ALDH3B2` | 2 |
| `cn_MMP10` | `cn_MMP1` | 2 |
| `cn_MMP3` | `cn_MMP12`, `cn_MMP13` | 3 |
| `cn_C11orf92` | `cn_POU2AF1` | 2 |
| `cn_TMPRSS4` | `cn_SCN2B` | 2 |
| `cn_CD3D` | `cn_CD3G` | 2 |
| `cn_HEPACAM` | `cn_HEPN1` | 2 |
| `cn_B3GAT1` | `cn_GLB1L2` | 2 |
| `cn_ZIC5` | `cn_ZIC2` | 2 |
| `cn_TCL1B` | `cn_TCL1A` | 2 |
| `cn_CRIP1` | `cn_ADAM6`, `cn_KIAA0125` | 3 |
| `cn_ARHGDIG` | `cn_HBA1`, `cn_HBA2`, `cn_RGS11` | 4 |
| `cn_METRN` | `cn_CCDC78`, `cn_HAGHL`, `cn_GNG13`, `cn_MSLN` | 5 |
| `cn_CACNA1H` | `cn_TPSAB1`, `cn_TPSB2`, `cn_TPSD1`, `cn_TPSG1` | 5 |
| `cn_PRSS33` | `cn_PRSS21`, `cn_ZG16B` | 3 |
| `cn_ABCC6` | `cn_MPV17L`, `cn_MYH11`, `cn_PLA2G10` | 4 |
| `cn_ASPHD1` | `cn_DOC2A`, `cn_FAM57B`, `cn_SEZ6L2` | 4 |
| `cn_MT1M` | `cn_MT1A` | 2 |
| `cn_SDR42E1` | `cn_HSD17B2` | 2 |
| `cn_CPNE7` | `cn_DPEP1`, `cn_TUBB3` | 3 |
| `cn_CCL11` | `cn_CCL13`, `cn_CCL7`, `cn_CCL8` | 4 |
| `cn_CCL18` | `cn_TBC1D3C`, `cn_TBC1D3G` | 3 |
| `cn_KRT13` | `cn_KRT15` | 2 |
| `cn_KRT16` | `cn_KRT17` | 2 |
| `cn_GFAP` | `cn_KIF18B` | 2 |
| `cn_CRHR1` | `cn_MAPT` | 2 |
| `cn_HOXB5` | `cn_HOXB6` | 2 |
| `cn_CBLN2` | `cn_NETO1` | 2 |
| `cn_CFD` | `cn_KISS1R` | 2 |
| `cn_LRG1` | `cn_PLIN4`, `cn_PLIN5` | 3 |
| `cn_FUT3` | `cn_FUT6` | 2 |
| `cn_DMKN` | `cn_SBSN` | 2 |
| `cn_CYP2A6` | `cn_CYP2A7` | 2 |
| `cn_CEACAM7` | `cn_CEACAM5` | 2 |
| `cn_KLK7` | `cn_KLK11`, `cn_KLK12` | 3 |
| `cn_KLK8` | `cn_KLK10` | 2 |
| `cn_TNNT1` | `cn_TNNI3`, `cn_SYT5` | 3 |
| `cn_CDC45` | `cn_PRODH`, `cn_SLC7A4`, `cn_TBX1` | 4 |
| `cn_VPREB3` | `cn_MMP11` | 2 |
| `cn_INPP5J` | `cn_PLA2G3` | 2 |
| `cn_SULT4A1` | `cn_PNPLA3` | 2 |
| `cn_KLHDC7B` | `cn_MAPK8IP2`, `cn_MLC1` | 3 |

## Safe future handling

Do not alter either source CSV. Before downstream modeling, retain an alias map from each duplicate vector to all gene names and verify the values against the original GISTIC/genomic-coordinate source. For prediction, a deterministic one-representative-per-identical-vector reduction is mathematically lossless for the matrix but changes gene-level interpretation; it should be documented and applied consistently. For biological interpretation, report the whole alias group or map to the shared segment rather than attributing a signal to one arbitrarily selected gene. Never treat the 99 aliases as independent evidence or importance signals.
## References and source URLs

1. Christensen NJ et al. *Identifying interactions in omics data for clinical biomarker discovery using symbolic regression.* Bioinformatics 38(15), 3749-3758 (2022). DOI: https://doi.org/10.1093/bioinformatics/btac405
2. Ciriello G et al.; TCGA Research Network. *Comprehensive Molecular Portraits of Invasive Lobular Breast Cancer.* Cell 163(2), 506-519 (2015). DOI: https://doi.org/10.1016/j.cell.2015.09.033
3. Original paper code/data repository: https://github.com/abzu-ai/QLattice-clinical-omics (README states the repository holds the analysis and data for the Bioinformatics paper; `data/` includes `brca_data_w_subtypes.csv` and `brca_data_w_meta.csv`; notebook `notebooks/bc_case_multi.ipynb`).
4. Kaggle dataset identified by the paper: https://www.kaggle.com/datasets/samdemharter/brca-multiomics-tcga. The Dataset C README links to a second Kaggle upload: https://www.kaggle.com/datasets/mannekuntanagendra/multi-omics-breast-cancer-prediction-qlattice.
5. Christensen et al. arXiv version, *Understanding Breast Cancer Survival: Using Causality and Language Models on Multi-omics Data*: https://arxiv.org/abs/2305.18410. Used as a secondary cross-check, not as a substitute for the original data source.
6. `curatedTCGAData` package documentation: https://bioconductor.org/packages/curatedTCGAData/ (TCGA assays and clinical data coordinated using patient/sample identifiers).
7. GDC TCGA metadata definitions: https://docs.cancergenomicscloud.org/docs/tcga-metadata (defines vital status as living/deceased; days-to-death and follow-up are distinct fields).
