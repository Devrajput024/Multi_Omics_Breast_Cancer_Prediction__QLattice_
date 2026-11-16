# Archive inventory

## Scope and method

Audited `C:\Users\asus\Downloads\archive (10).zip` in place on 2026-11-10. The ZIP was enumerated and CSV/PNG content was read directly from the archive; it was not extracted. No source file was changed. The downloaded ZIP is 204,419,676 bytes (195.0 MiB by decimal/binary conversion: 194.95 MiB). It contains 1,562 files, no explicit directory entries, and no README, codebook, or documentation files.

## Structure and counts

| Archive path | Files | Format | Uncompressed bytes |
|---|---:|---|---:|
| `dataset/dataset1/benign/images/` | 437 | PNG | 148,088,597 |
| `dataset/dataset1/benign/masks/` | 437 | PNG | 239,231 |
| `dataset/dataset1/malignant/images/` | 210 | PNG | 67,882,887 |
| `dataset/dataset1/malignant/masks/` | 210 | PNG | 187,367 |
| `dataset/dataset1/normal/images/` | 133 | PNG | 49,623,545 |
| `dataset/dataset1/normal/masks/` | 133 | PNG | 151,062 |
| `dataset/dataset2/patient_history_dataset.csv` | 1 | CSV | 142,053 |
| `dataset/dataset3/molecular_biomarker_dataset.csv` | 1 | CSV | 83,235 |

Total: 1,560 PNG files and 2 CSV files. All archive data paths use the three dataset folders shown; Dataset 1 has benign, malignant, and normal categories, each with `images` and `masks` subfolders. Each PNG has a same-basename image/mask partner in its category. See [image_audit.md](image_audit.md) for validation findings.

The ZIP entry timestamps are 2026-03-24. Sizes above are uncompressed member sizes; archive total size is the downloaded file size noted above. File extensions are exclusively `.png` and `.csv`.

## Source and evidence boundaries

The source identified by the user is [Kaggle: Multi-Modal Breast Cancer Dataset](https://www.kaggle.com/datasets/ajithdari/multi-modal-breast-cancer-dataset). The page could be opened but provided no readable description or documentation to this audit. The archive itself contains no source citation, image acquisition details, mask-generation protocol, data dictionary, or publication references. Do not infer that this archive is TCGA, METABRIC, or a cohort matching Dataset C from its contents alone.
