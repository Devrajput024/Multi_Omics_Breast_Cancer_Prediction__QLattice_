# Recovery feasibility conclusion

Repository snapshot: `main` at commit `fc29704243010890969b524ebe3e7c2be1f7ef53`; retrieved 2026-11-10. Files were saved under this folder. Raw GitHub URLs at that pinned revision are used below.

## Decision: C — Not ready because essential information is unavailable

The source repository confirms a 705-row processed TCGA multi-omics cohort and provides a binary vital-status target, a histology field, and evidence sufficient to match all rows across files by exact shared content. It does **not** provide identifiers or the preparation scripts needed to verify patient identity, restore clinical follow-up, resolve the histology recoding, or reproduce assay processing. These gaps block a validated patient-level clinical recovery and a fully provenance-secure modeling dataset.

## Recovery status by item

| Item | Status | Evidence / limit |
|---|---|---|
| Repository provenance and published source | **Recovered and verified** | Repository README cites Christensen et al. (2022); paper identifies TCGA data via `curatedTCGAData`. |
| Cohort size and common row content | **Recovered and verified** | 705 rows; 705 unique content-signature matches to repo metadata, with no unmatched/ambiguous signatures. This is not an ID match. |
| Patient/sample identifiers or TCGA barcodes | **Not recoverable from available sources** | No identifier fields in the supplied CSVs or investigated notebook. |
| Binary vital-status semantics | **Recovered and verified** | Paper maps 0 to alive/survived and 1 to fatal/deceased; dataset counts 611/94. |
| Days to death / last follow-up / diagnosis timing | **Not recoverable from available sources** | Not in the derivatives; no IDs for a safe GDC lookup. No clinical join was attempted. |
| Histology labels in repository versus local | **Recovered and verified** | 574/131 local versus 575/130 repo; one matched row differs from lobular to ductal. |
| Explanation for published 490/127/88 versus binary labels | **Partially recovered** | 490+127+88=705 and mixed-label collapse is plausible, but no IDs or recoding script verify the actual rule. |
| Broad assay modality names | **Recovered and verified** | Paper documents expression, CNV, mutations, and protein features; prefixes are consistent with its feature groups. |
| Assay versions, thresholds, normalization, caller, segment maps | **Not recoverable from available sources** | No derivation pipeline or sufficiently detailed assay metadata in the repo artifacts. |
| QLattice notebook target/workflow context | **Recovered and verified** | Notebook reads prepared CSV and targets `vital.status`; README states published results used Feyn 2.1.2 remote service. |

## Next recovery requirements

Obtain an authoritative sample-ID-to-row mapping or original cohort table and the original assay preparation scripts/metadata. With verified identifiers, retrieve current GDC clinical records and derive event time and censoring according to documented rules; do not join by row order. Resolve the histology discrepancy against pathology-reviewed case-level labels while retaining mixed/other categories. Then document feature time origin and exclusions before deciding between binary classification and survival analysis.

No original Dataset C files were modified. No clinical records were downloaded or joined, and no models were trained.

## References

- Repository pinned at `fc29704243010890969b524ebe3e7c2be1f7ef53`: https://github.com/abzu-ai/QLattice-clinical-omics/tree/fc29704243010890969b524ebe3e7c2be1f7ef53
- Christensen et al. (2022): https://doi.org/10.1093/bioinformatics/btac405
- Ciriello et al. (2015): https://doi.org/10.1016/j.cell.2015.09.033
- GDC clinical fields: https://docs.gdc.cancer.gov/Data_Portal/Users_Guide/clinical_data_analysis/
- `curatedTCGAData`: https://bioconductor.org/packages/curatedTCGAData
