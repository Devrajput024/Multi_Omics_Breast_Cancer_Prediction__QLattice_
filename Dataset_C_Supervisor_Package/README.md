# Dataset C experiment handover

This folder contains the script, its required Dataset C input, pinned Python dependencies, and a small set of **reference outputs** from a previous run. The script calculates metrics from model execution each time it runs; it does not read the reference outputs to construct results.

## Contents

- `dataset_C_reproducibility.py` — complete experiment and output-generation script.
- `requirements.txt` — pinned third-party dependencies.
- `dataset_c_extracted/data.csv` — the exact primary Dataset C matrix used by the project (705 samples, 1,936 predictors, `vital.status`).
- `reference_outputs/` — an example final summary, dataset summary and comparison from a previous run. These are provided for comparison and are not outputs from the supervisor's run.

The experiment does not use the subtype-augmented CSV, QLattice notebooks or model files, previous trained pipelines, or any files outside this package. It uses no VS Code/Codex features or laptop-specific absolute paths. By default, input and output paths are resolved relative to the script directory.

## Requirements

- Python 3.11 (recommended and used for the local reproduction).
- A terminal with Python available as `python` (or use `py -3.11` on Windows).
- Internet access for the initial pip dependency installation.

## Install and run (Windows PowerShell)

Open PowerShell in this folder and run:

```powershell
py -3.11 -m venv .venv
.\.venv\Scripts\Activate.ps1
python -m pip install --upgrade pip
python -m pip install -r requirements.txt
python dataset_C_reproducibility.py
```

If PowerShell blocks activation, skip activation and call the virtual environment's interpreter directly:

```powershell
py -3.11 -m venv .venv
.\.venv\Scripts\python.exe -m pip install -r requirements.txt
.\.venv\Scripts\python.exe dataset_C_reproducibility.py
```

On macOS/Linux, create and activate the environment with `python3.11 -m venv .venv` and `source .venv/bin/activate`, then install requirements and run `python dataset_C_reproducibility.py`.

## Dataset placement

The included dataset is already in the expected location: `dataset_c_extracted/data.csv`. Keep that relative path and run the command from anywhere; the script resolves the default input relative to its own file. If the CSV is missing, the script stops with an explicit error and does not search for or substitute another dataset. An alternate copy of the same verified dataset can be given with `--data PATH`, but its schema and labels must match the validated Dataset C matrix.

## Results and reproducibility

The default run creates `dataset_C_reproducibility_results/` beside the script. It includes the dataset/schema audit, modality feature names, candidate and nested-CV scores, model-selection details, held-out metrics and predictions, classification reports, selected features and scores, confusion matrices, ROC and precision-recall plots, model-comparison plots, trained pipelines, `execution_log.txt`, and `final_summary.txt`.

The verified method uses a stratified 80/20 row split (seed `20261110`), then nested stratified 3-fold outer × 3-fold inner CV on training rows. Model choice ranks ROC-AUC, PR-AUC, balanced accuracy, macro-F1 and deceased recall. The holdout is evaluated after selection and does not select features or models. Preprocessing and selectors are fit within training folds. The five candidates are Logistic Regression, Random Forest, Extra Trees, RBF SVM and HistGradientBoosting. Duplicate CNV vectors are retained, as in the final project methodology.

This is exploratory classification of recorded vital status (`0=Alive`, `1=Deceased`), not survival-time prediction or clinical validation. There are no patient identifiers to establish patient-level independence.

## Troubleshooting

- **Missing package:** confirm the environment is active and run `python -m pip install -r requirements.txt` using that environment's Python.
- **Dataset missing:** restore `dataset_c_extracted/data.csv` inside this folder. The script deliberately refuses substitute files.
- **Schema/target mismatch:** use the verified primary matrix; the subtype-augmented file is not the experiment input.
- **SVC deprecation warning:** scikit-learn 1.9 may warn about the preserved `SVC(probability=True)` option. It is a warning, not a fit error; it is retained to match the archived final experiment.
- **Fit error on another setup:** check `execution_log.txt` and `model_fit_errors.csv`. The run logs exceptions and does not silently replace failed models.
- **Different numerical results:** confirm Python/package versions, input SHA-256 in `run_metadata.json`, and that the exact same data file was used. Small platform/library differences can affect floating-point results.

The reference outputs are comparison aids only. Running the script writes newly computed outputs into the results folder and does not use reference metrics as inputs.
