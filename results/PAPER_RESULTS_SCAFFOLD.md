# HIV-ShiftGuard: Mutation-Novelty-Aware Uncertainty-Calibrated Prediction of HIV-1 Antiretroviral Susceptibility

> **AUTO-GENERATED METHODS AND RESULTS SCAFFOLD. NOT A SUBMISSION-READY PAPER.**
> Verify experiment quality, compare literature and write scholarly narrative before submission.

## Abstract — to be written after validation

Background: Existing HIV-1 genotype-to-phenotype prediction can be overconfident on unfamiliar viral mutation patterns.
Methods: Drug-specific mutation-based models were benchmarked against baselines using four grouped partitions and three independent seeds, with mutation-novelty-aware conformal inference and selective review.
Results: 32313 observed genotype–drug pairs, 23 evaluated drugs, ensemble mean per-drug-seed RMSE 1.017 log2 fold; observed pooled novelty-scaled interval coverage 90.0%. These are in-vitro outcomes only.
Conclusion: Report whether improvement over all baselines is actually supported by the per-drug tables and discuss failures. **Never imply clinical treatment selection or HIV cure discovery.**

## Introduction
Describe established ART and the limits of existing genotypic resistance prediction. Position this work against recent protein language model baselines, including 2026 HIV-ESM-2 research; do not claim all component methods are novel.

## Data
Stanford HIVDB genotype-phenotype research datasets: PI, NRTI, NNRTI and INI (direct links in README). Fold susceptibility is laboratory IC50 fold-change, not patient response. For the retrieved versions, 32313 valid drug-isolate measurements were retained and 7148 gene-genotype isolates identified. Source files are SHA256-tracked.

## Methods
Mutation substitution tokens encode each observed position (ignore missing '.' and ambiguous multi-amino-acid calls). For each drug independently, a linear ridge baseline, full-substitution ExtraTrees regressor, and position-only ablation were fitted. A fusion weight was chosen using validation error, not calibration or test data. Exact-genotype / sequence-ID linked components ensure no exact overlap between train, validation, calibration, and test groups. The four sets are nominally 60/10/10/20 percent but may differ in row-level sizes because the splits are group-based. Fold-change outcomes were transformed as log2(FC) and evaluated with RMSE, MAE, R², Pearson r and Spearman rho. The unscaled split-conformal residual quantile and novelty-scaled variant were calibrated on separate held-out calibration samples. Novelty is the nearest-training-mutation Jaccard distance. The novelty gate is an experimental abstention analysis and not a clinical decision rule.

## Results: ablation (mean ± SD across drug-seed runs)

| Model | RMSE mean | RMSE SD | MAE mean | Runs |
|---|---:|---:|---:|---:|
| median_baseline | 2.334 | 0.779 | 1.872 | 67 |
| linear_mutations | 1.113 | 0.365 | 0.824 | 67 |
| position_only_ablation | 1.267 | 0.377 | 0.862 | 67 |
| mutation_forest | 1.111 | 0.307 | 0.759 | 67 |
| shiftguard_ensemble | 1.017 | 0.283 | 0.719 | 67 |

Observed pooled standard interval coverage: 90.1%; novelty-scaled coverage: 90.0% (nominal 90%). Do not claim novelty is superior unless both error and interval efficiency support the claim, with uncertainty intervals and appropriate comparisons.

## Figures
- Figure 1 — system architecture (<code>figures/Figure_01_Research_Architecture.svg</code>)
- Figure 2 — model/ablation comparison
- Figure 3 — drug-specific generalization
- Figure 4 — predicted vs observed fold change
- Figure 5 — marginal interval coverage
- Figure 6 — novelty and test error
- Figure 7 — novelty subgroup sensitivity

## Limitations
- Data are phenotype measurements from virus isolates, not clinical trial outcomes or cure efficacy.
- Patient-unique identifiers, longitudinal response, treatment adherence, adverse effects, viral reservoirs and drug–drug interactions are not modeled.
- Mutation novelty does not in itself guarantee reliable intervals or robust generalization.
- Several older drugs are included because they have laboratory phenotype data; rankings are not treatment recommendations.
- Group splits guard against direct exact-genotype overlap, not all patient-level, study-level or temporal dependencies.
- This is a proposed combination/benchmark, not proven unique in the published literature.

## References to investigate and cite accurately
1. Stanford HIV Drug Resistance Database: https://hivdb.stanford.edu/pages/genopheno.dataset.html
2. NIH HIV drug resistance guidance: https://hivinfo.nih.gov/understanding-hiv/fact-sheets/drug-resistance
3. WHO 2026 HIV recommendations: https://www.who.int/news/item/07-01-2026-who-releases-updated-recommendations-on-hiv-clinical-management
4. Farquhar, 2026, Protein language model embeddings improve HIV drug resistance prediction: https://pubmed.ncbi.nlm.nih.gov/42105210/
5. Original split-conformal prediction literature, and clinically validated HIVdb reference papers.

Do not submit without author validation, literature review and correct venue formatting.
