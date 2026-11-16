# Assay and preprocessing provenance

Repository snapshot: `main` at commit `fc29704243010890969b524ebe3e7c2be1f7ef53`; retrieved 2026-11-10. Files were saved under this folder. Raw GitHub URLs at that pinned revision are used below.

## What is supported by the research publication

The Christensen et al. 2022 paper describes a TCGA-BRCA multi-omics dataset obtained using `curatedTCGAData` and identifies four modalities: RNA expression (`rs_`), copy-number alterations (`cn_`), somatic mutations (`mu_`), and protein expression (`pp_`). It reports variance-based feature filtering before modeling. Feature names and paper descriptions support these broad modality interpretations; the processed CSV values alone do not establish raw assay provenance.

The related Ciriello et al. paper describes TCGA invasive breast tumor assays including RNA sequencing, exome sequencing, SNP arrays, and RPPA; it notes RPPA coverage for 633 of 817 tumors. This is useful cohort context but does not prove the exact platform/version/normalization of every `pp_` feature in this 705-row derivative table.

## Not specified in the available preparation artifacts

The repository contains prepared CSVs and an analysis notebook, not the pipeline that retrieved or transformed the assays. The notebook reads the finished CSV; it does not show raw-data preparation. The available files do not verify:

- GISTIC release, thresholds, or segment-to-gene mapping used for `cn_` values;
- RNA-seq assay release, normalization method, units, or any transformation for `rs_` values;
- mutation caller, variant filtering, or exact binary/other encoding for `mu_` values;
- exact protein assay platform/version, antibody panel, normalization, or units for `pp_` values;
- the exact variance threshold or feature-selection calculation used in the published preparation.

The source paper's modality descriptions and variance filtering are not enough to fill those implementation details. The user's earlier duplicate-feature audit found 62 groups of identical CNV columns (99 extra columns); the published/repository preparation artifacts do not map those vectors to genomic segments. Shared segment-level GISTIC calls are a plausible explanation for gene columns being identical, but this remains a hypothesis without coordinates and transformation code.

## Safe future handling

Preserve the source columns and retain a machine-readable mapping from each original feature to its assay annotation. Investigate redundant CNV vectors against verified genomic coordinates and GISTIC segment definitions before dropping or consolidating any column. Recover assay-specific processing metadata and reference versions before interpreting values biologically. Do not infer that duplicate columns are import artifacts or remove them solely because their values are identical.

## Sources

- Christensen et al. (2022): https://doi.org/10.1093/bioinformatics/btac405
- Ciriello et al. (2015): https://pmc.ncbi.nlm.nih.gov/articles/PMC4603750/
- `curatedTCGAData`: https://bioconductor.org/packages/curatedTCGAData
- Repository README/notebook and pinned sources: [repository inspection](repository_inspection.md)
