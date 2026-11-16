# Patient ID matching audit

## Exact matching results

Compared the literal `Patient ID` values (no trimming, normalization, or fuzzy matching) in the clinical and molecular CSVs. Also compared each PNG basename after removing `.png` against those IDs.

| Comparison | Unique IDs left | Unique IDs right | Exact intersection | Left only | Right only | Duplicate ID occurrences |
|---|---:|---:|---:|---:|---:|---:|
| Dataset 2 clinical vs Dataset 3 molecular | 780 | 780 | 780 | 0 | 0 | 0 in either table |
| Tables vs Dataset 1 image basenames | 780 | 780 | 780 | 0 | 0 | 0 image basename duplicates |
| Tables vs Dataset 1 mask basenames | 780 | 780 | 780 | 0 | 0 | 0 mask basename duplicates |

The two table ID sequences are also identical in row order. Both tables have one row per ID. For every image/mask basename, the corresponding clinical `class` exactly matches its folder category (benign/malignant/normal): 780/780 category matches.

## What the match does and does not establish

The archive provides a deterministic join key across its three directory/table components: all 780 image IDs and all 780 mask IDs have one record in each table, with no unmatched or repeated IDs. This supports a file-level join within this archive. It does not establish that these are independently collected patient modalities, that the PNG represents the patient’s diagnostic slide or tissue, or that patient-level clinical/molecular labels are valid targets for each image. There is no source slide identifier, specimen identifier, image metadata, collection date, site, image protocol, assay provenance, or linkage documentation. The `MB-####` naming scheme and exact class agreement should be reported as observed metadata, not used to infer TCGA/METABRIC or identity overlap with Dataset C.
