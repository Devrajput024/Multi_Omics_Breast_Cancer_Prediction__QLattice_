# Existing QLattice model inspection

No model was loaded into Feyn, no prediction/inference was executed, and no performance is claimed. Both files are parsed as serialized graph metadata only.

## `model1_200its_mc7.model`

- Output/target: `vital.status` (`out:lr`).
- Input features: `mu_TNXB`, `rs_KRT23`, `rs_APOB`, `cn_ANKRD30B`.
- Graph operators: `multiply`, `gaussian2`, `add`.
- Internal serialization version: `2021-08-30`.
- Saved graph metadata: `{'pid': 10961630358094799780, 'ppid': 5768885376330529529, 'generation': 8, 'action': 'c', 'ix': 4}`.
- The output node is a logistic/binary-classification output in the Feyn/QLattice serialization. Input nodes include learned scale, weight, and bias parameters, so some scaling is embedded in the graph. No separate preprocessing pipeline, imputer, feature-selection object, training split, predictions, or sample IDs are stored.

## `model_200its_md2_bic_0.model`

- Output/target: `vital.status` (`out:lr(i)->b`).
- Input features: `rs_APOB`, `cn_ANKRD30B`, `rs_VWDE`, `mu_TNXB`.
- Graph operators: `cell:add(i,i)->i`, `cell:add(i,i)->i`, `cell:multiply(i,i)->i`.
- Internal serialization version: `2020-02-07`.
- Saved graph metadata: `{'node_count': 8, 'link_count': 7}`.
- The output node is a logistic/binary-classification output in the Feyn/QLattice serialization. Input nodes include learned scale, weight, and bias parameters, so some scaling is embedded in the graph. No separate preprocessing pipeline, imputer, feature-selection object, training split, predictions, or sample IDs are stored.

## Relationship to the supplied dataset and notebook

Both graphs use `vital.status` as output. Every listed input feature exists in `data.csv` and belongs to a documented feature family. The notebook reads `data.csv`, generates a separate QLattice search on a stratified split, and later explicitly loads `model1_200its_mc7.model`; it then saves it back to the same path. The second model file is not referenced in notebook code. The notebook does not show how either artifact was originally trained or which source rows/version were used. Neither model contains a dataset hash, patient IDs, feature-matrix hash, train/test membership, date of training, or validation record. Thus training against exactly these CSV bytes cannot be established, even though file names, output labels and input features are compatible.

The companion repository states that the published paper results were generated with Feyn 2.1.2 using a remote Abzu server, which is now deprecated. The supplied notebook instead records a failed install attempt for Feyn 2.1.4 followed by Feyn 3.4.1. The artifacts carry older internal serialization version strings (2021-08-30 and 2020-02-07). This version history means exact compatibility/reproduction needs to be checked in a controlled environment; it does not by itself prove the models cannot be loaded. The model file names include strings such as `200its` and `bic`, but those are filenames, not reliable evidence of search epochs or criterion. The notebook's Markdown contains performance claims, but they are not tied to these serialized models by artifact IDs or saved prediction tables and are not independently reproducible from the available files. Do not cite those scores as verified performance.

## Leakage and reuse notes

- Exclude `vital.status` from predictors explicitly, even when the QLattice API receives it separately as `output_name`. The inspected graphs themselves do not list it as an input.
- Confirm input scaling and categorical encoding against the original Feyn version before using the graph. The artifact parameters encode input transforms, but the exact training data type declarations/preprocessing are not saved.
- These models are historical hypotheses, not validated deployable models. Refit/re-evaluate only after source/sample alignment, label time horizon, folds and data version are fixed.
## References and source URLs

1. Christensen NJ et al. *Identifying interactions in omics data for clinical biomarker discovery using symbolic regression.* Bioinformatics 38(15), 3749-3758 (2022). DOI: https://doi.org/10.1093/bioinformatics/btac405
2. Ciriello G et al.; TCGA Research Network. *Comprehensive Molecular Portraits of Invasive Lobular Breast Cancer.* Cell 163(2), 506-519 (2015). DOI: https://doi.org/10.1016/j.cell.2015.09.033
3. Original paper code/data repository: https://github.com/abzu-ai/QLattice-clinical-omics (README states the repository holds the analysis and data for the Bioinformatics paper; `data/` includes `brca_data_w_subtypes.csv` and `brca_data_w_meta.csv`; notebook `notebooks/bc_case_multi.ipynb`).
4. Kaggle dataset identified by the paper: https://www.kaggle.com/datasets/samdemharter/brca-multiomics-tcga. The Dataset C README links to a second Kaggle upload: https://www.kaggle.com/datasets/mannekuntanagendra/multi-omics-breast-cancer-prediction-qlattice.
5. Christensen et al. arXiv version, *Understanding Breast Cancer Survival: Using Causality and Language Models on Multi-omics Data*: https://arxiv.org/abs/2305.18410. Used as a secondary cross-check, not as a substitute for the original data source.
6. `curatedTCGAData` package documentation: https://bioconductor.org/packages/curatedTCGAData/ (TCGA assays and clinical data coordinated using patient/sample identifiers).
7. GDC TCGA metadata definitions: https://docs.cancergenomicscloud.org/docs/tcga-metadata (defines vital status as living/deceased; days-to-death and follow-up are distinct fields).
