# Final figure captions

1. **Modality feature distribution.** Counts of the 1,936 prefixed feature columns in the primary matrix: CNV 860, RNA 604, mutation 249, protein 223.
2. **Recorded vital-status distribution.** Alive code 0: 611; deceased code 1: 94. The target is status only, with no follow-up time.
3. **Individual and integrated modality comparison.** Repeated 5×5 training CV ROC-AUC means and standard deviations across five repeat means for the prespecified architecture per scope. Error bars describe resampling variability, not independent-sample CIs.
4. **Integrated versus strongest individual modality.** Repeated-CV ROC-AUC comparison for RNA Random Forest and integrated Extra Trees, selected by prior training-CV results.
5. **ROC curves.** Frozen predictions on the common 141-row exploratory holdout for Dummy and the CV-selected model architecture in each modality. Test curves are descriptive.
6. **Precision–recall curves.** Frozen holdout precision/recall curves; the dashed horizontal line marks deceased prevalence in the test set.
7. **Confusion matrices.** Held-out confusion counts at threshold 0.5; rows/columns use class order alive (0), deceased (1).
8. **Feature importance.** Top repeated-fold tree importance summaries for RNA and integrated models; importance is associative and model-dependent.
9. **Feature stability.** Frequency with which variables ranked in the top 30 across 25 training folds; the dashed threshold marks the report’s descriptive 60% stability cutoff.
10. **Calibration curves.** Uncalibrated probabilities from the saved holdout predictions and probabilities from sigmoid calibrators fitted using training-only cross-validation; curves use the same held-out cases (n=141, 19 deceased).

Each figure is exported as a 300-dpi PNG and a vector PDF.
