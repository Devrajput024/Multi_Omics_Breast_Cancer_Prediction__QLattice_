# Clinical follow-up investigation

Repository snapshot: `main` at commit `fc29704243010890969b524ebe3e7c2be1f7ef53`; retrieved 2026-11-10. Files were saved under this folder. Raw GitHub URLs at that pinned revision are used below.

## What is present

The repository metadata CSV includes `vital.status` and a histology label, but no days-to-death, days-to-last-follow-up, time-to-event, diagnosis date, or follow-up date field. The notebook uses `vital.status` as its prediction target and contains no clinical data retrieval or time-to-event processing. Content-signature comparison confirms that its `vital.status` values agree with Dataset C across all 705 rows.

The 2022 source paper defines the binary coding as 0 = alive/survived and 1 = fatal/deceased and reports 611 and 94 observations respectively. In the provided table, `vital.status` alone has no time horizon, censoring time, or event duration. It is an observed vital-status outcome, not a complete survival endpoint. A binary classification target is possible as defined in the paper, but its interpretation is limited by unknown follow-up and censoring; time-to-event analysis would generally be more appropriate if valid event/censoring times can be recovered.

## Can TCGA follow-up be recovered?

GDC clinical records expose vital status and time fields such as `demographic.days_to_death`, `diagnoses.days_to_last_follow_up`, and in current data structures `follow_ups.days_to_follow_up`. Bioconductor's `curatedTCGAData` packages TCGA assays with clinical/sample metadata, and TCGA utilities document barcode mapping. These sources show that relevant clinical information exists in the authoritative TCGA ecosystem.

However, Dataset C's derivative CSVs and the investigated repository files contain no sample barcode, case ID, or UUID. Therefore the cohort rows cannot be linked to individual TCGA clinical records with verified identifiers. No external clinical data were downloaded or joined. Diagnosis date / follow-up fields are **not recovered for this cohort**.

## Leakage assessment

The table contains only a binary vital-status target, with no survival duration columns. From the inspected files alone, direct survival-time leakage is not evidenced. The provenance and temporal relationship of every feature to the outcome are not documented well enough to rule out indirect leakage or guarantee that all measurements precede the outcome. Before modeling, obtain cohort IDs and define a prediction time origin; exclude post-outcome or outcome-derived variables if present.

## Sources

- Paper target definition: https://doi.org/10.1093/bioinformatics/btac405
- GDC clinical analysis fields: https://docs.gdc.cancer.gov/Data_Portal/Users_Guide/clinical_data_analysis/
- GDC available clinical fields: https://docs.gdc.cancer.gov/API/Users_Guide/Appendix_A_Available_Fields/
- `curatedTCGAData`: https://bioconductor.org/packages/curatedTCGAData
- `TCGAutils` barcode and sample mapping: https://bioconductor.org/packages/TCGAutils
- Repository and pinned file links: [repository inspection](repository_inspection.md)
