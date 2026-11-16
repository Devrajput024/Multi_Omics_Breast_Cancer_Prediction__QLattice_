# Dataset C feature-quality report

The primary matrix has 705 rows and 1936 predictors. No source values were modified. The complete per-feature table is `feature_quality_summary.csv`; high absolute correlations (|Pearson r| ≥ 0.95) are in `correlation_analysis.csv`.

## Missingness and constants

Missing cells in the primary modelling matrix: 0. Constant features: 0; features with one value in at least 95% of rows: 230. Feature-level dominance, zero fractions, unique counts and distribution quantiles are in the summary CSV.

## Duplicate and redundant features

The previous audit found 62 groups of exactly identical vectors and 99 extra CNV aliases. These are identical observed data vectors, not biologically identical genes. The exhaustive numeric correlation screen reports 1116 pairs at |r|≥0.95, including exact CNV aliases; see CSV for identifiers and correlations.

## Modality distributions

- CNV observed categories and per-feature counts are in `tables/cnv_category_distribution.csv`; the modelling pipeline treats them categorically.
- Mutation matrix zero fraction is 0.9688; the mean nonzero mutation-feature count per sample is 7.76. Observed values outside binary 0/1 are listed in `reports/mutation_sparsity.json` if any.
- RNA and protein feature distributions, quantiles and IQR-screened outlier counts are in `feature_quality_summary.csv`. IQR flags are descriptive, not automatic errors; no rows or values were removed.

## Correlation and outlier caveats

Pearson correlations were computed across all numeric columns, including categorical GISTIC-style CNV codes treated numerically for this screening only. This highlights redundancy but is not a biological association analysis. IQR outlier counts were computed only for RNA/protein continuous features; CNV states and sparse mutation indicators were not treated as continuous outlier candidates.
