# Histology discrepancy

Repository snapshot: `main` at commit `fc29704243010890969b524ebe3e7c2be1f7ef53`; retrieved 2026-11-10. Files were saved under this folder. Raw GitHub URLs at that pinned revision are used below.

## Observed labels

The published Ciriello et al. study reports a pathology-reviewed TCGA cohort of 817 tumors: 490 IDC, 127 ILC, 88 mixed IDC/ILC, and 112 other histologies. It also says that only 633 cases had RPPA data. Dataset C has 705 rows and only two histology labels: 574 infiltrating ductal and 131 infiltrating lobular. The recovered repository metadata CSV has 575 ductal and 130 lobular.

Using shared molecular/target content signatures, all 705 rows matched uniquely without using row order. The cross-tab is:

| Dataset C histology | Repository metadata histology | Count |
|---|---|---:|
| Ductal | Ductal | 574 |
| Lobular | Lobular | 130 |
| Lobular | Ductal | 1 |

Thus the one-label difference between local and repository histology is verified. No label has been recoded.

## What likely explains 705 versus the publication's counts?

The published subtotal 490 + 127 + 88 equals 705 exactly. This is consistent with a cohort formed from IDC, ILC, and mixed tumors while excluding the 112 other histologies. That arithmetic is suggestive, not proof that these CSV rows represent exactly those pathology-reviewed cases, because the derivative CSVs lack IDs and the repo has no cohort-construction script.

A binary 574/131 split could arise if mixed cases were reassigned to ductal/lobular categories. Relative to the publication counts, it would imply allocating approximately 84 of 88 mixed cases to ductal and 4 to lobular; repository metadata's 575/130 split would imply 85 and 3. These are arithmetic implications only. The available files do not establish that this is what happened, the exact recoding rule, or which case was labeled differently.

A one-row label change between repo metadata and local export demonstrates at least a label/version discrepancy. It does not explain the aggregate difference from the publication. Different cohort/version, histology filtering, and collapsing mixed labels remain plausible but unverified explanations.

## Needed to resolve

Recover the original sample IDs and the actual cohort extraction / label-preparation script or authoritative supplemental sample table. Compare case-level pathology labels, retaining mixed and other categories; do not recode until that mapping is evidenced.

## Sources

- Ciriello et al. (2015), *Comprehensive Molecular Portraits of Invasive Lobular Breast Cancer*: https://pmc.ncbi.nlm.nih.gov/articles/PMC4603750/ (DOI: https://doi.org/10.1016/j.cell.2015.09.033)
- Christensen et al. (2022): https://doi.org/10.1093/bioinformatics/btac405
- Source repository and pinned CSVs: [repository inspection](repository_inspection.md)
