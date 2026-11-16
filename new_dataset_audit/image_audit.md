# Histopathology image and mask audit

## Inventory by category

| Category | Images | Masks | Same-name pairs | Image bytes | Mask bytes |
|---|---:|---:|---:|---:|---:|
| benign | 437 | 437 | 437 | 148,088,597 | 239,231 |
| malignant | 210 | 210 | 210 | 67,882,887 | 187,367 |
| normal | 133 | 133 | 133 | 49,623,545 | 151,062 |
| **Total** | **780** | **780** | **780** | **265,595,029** | **577,660** |

All 1,560 image/mask files have `.png` extensions and valid PNG signatures/headers. No corrupt image files were detected while parsing headers. Image dimensions vary (not a single fixed resolution); examples/range depend on individual files. Images span roughly hundreds of pixels per side; mask files are far smaller than corresponding source images in aggregate, consistent with compact masks but not proving their intended semantics. No mask annotation guide or class/pixel encoding is included.

For each category, image and mask basename sets match exactly: 437/437 benign, 210/210 malignant, 133/133 normal. There are no unmatched image/mask names. The basename is the same `MB-####` Patient ID appearing in both tables; all 780 image IDs and mask IDs map exactly once to both tables. Clinical `class` also agrees with directory category for every image ID.

## Duplicate check and limits

SHA-256 comparison found three exact duplicate-content pairs among all PNGs: `benign/images/MB-2990.png` is byte-identical to `malignant/images/MB-3298.png`; `normal/masks/MB-4959.png` equals `normal/masks/MB-4974.png`; and `normal/masks/MB-4986.png` equals `normal/masks/MB-5011.png`. No other duplicate file contents were found. This means one exact duplicated image occurs across benign and malignant labels, a label-quality issue to review. The two mask duplicates may be legitimate identical empty/simple masks or duplication errors; contents need semantic inspection with documented mask encoding to adjudicate.

Filenames include no slide/specimen accession or patient metadata beyond `MB-####`; the ZIP contains no sidecar metadata. The deterministic ID basename match supports an archive-internal join only. Without documentation, one cannot verify that each PNG is a patient-level histopathology image, the sampling/patch process, tissue source, stain/scanner, or the masks’ meaning. The requested category image/mask counts and joins are directly verifiable; clinical validity is not.

## Histological subtype labels

Yes, the tables include subtype-like labels: Dataset 2 has `Tumor Other Histologic Subtype` including Ductal/NST, Lobular, and Mixed; Dataset 3 has detailed cancer type and Oncotree Code values corresponding to ductal/IDC, lobular/ILC, mixed ductal-lobular/MDLC, mixed mucinous/IMMC, plus generic Breast/BREAST. Counts by Dataset 3 annotation are IDC 614, MDLC 94, ILC 59, IMMC 10, BREAST 3. IDC/ILC/mixed can be identified from these table labels, but they are absent from filenames and are not image-level pathology annotations. No explicit `IDC`, `ILC`, or mixed label is attached to the image records except by joining Patient ID to these categorical table labels.
