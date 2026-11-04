# HIV-ShiftGuard — preregisterable implementation study protocol

## Venue and provisional paper

Target venue: ICIRSET 2026, AI and Data Science / Emerging Computing track, 20–21 November 2026. Verify current deadline, fees, virtual author presentation and proceedings indexing directly with the organizing committee.

**Working manuscript title:** *HIV-ShiftGuard: Mutation-Novelty-Gated Conformal Prediction of Antiretroviral Susceptibility under Genotypic Distribution Shift*

## Purpose and clinical boundary

The study is about the response of HIV-1 isolates to established antiretroviral drugs in laboratory phenotypic assays. It is not a study of patient-level virologic suppression, drug safety, treatment regimen optimization, reservoir eradication, a functional cure or a sterilizing cure. No human treatment recommendations will be generated.

## Testable hypotheses

- **H1 (point prediction):** Substitution-aware nonlinear genotype models improve *in-vitro* fold-susceptibility prediction over a training-median baseline under exact-genotype-heldout evaluation. If the linear baseline wins, report that instead.
- **H2 (uncertainty):** Novelty-scaled conformal prediction intervals alter observed coverage and mean interval width relative to ordinary split conformal intervals. Improvement is not assumed in advance.
- **H3 (selective prediction):** A gate combining genotypic novelty and interval width can reduce error on accepted (non-flagged) research predictions, with a measurable tradeoff in retention. Gate rules are selected without access to test labels.
- **H4 (ablation):** Amino acid substitution identity provides measurable value beyond mutation position alone. The position-only ablation assesses this directly.

## Data management / leakage risks

Four Stanford high-quality filtered genotype–phenotype tables, accessed automatically from official Stanford URLs. Each genomic sequence genotype is grouped by its mutation pattern (within gene) and Stanford SeqID (across genes) so that exact repeats cannot cross splits. Lack of patient and source-study identifiers prevents claims of patient- or study-independent validation. Duplicate gene/SeqID/drug measurements are collapsed by the median and counted. Missing drugs, ambiguous mutations and other source limitations remain explicitly recorded.

## Model design

Training-only mutation vocabulary; one independent supervised regressor per drug. Fit a ridge baseline, nonlinear substitution ExtraTrees, and position-only tree ablation. Obtain a convex ridge/forest fusion coefficient from validation RMSE only. In separate calibration data, fit (i) absolute residual conformal intervals, (ii) fold-aware novelty-scaled intervals with predeclared gamma=1.0. Mutation novelty uses minimum exact Jaccard distance to a training genotype. No test label is observed during fitting, tuning or calibration.

## Repetitions and metrics

Three fixed seeds: 11, 23, 37. Train/validation/calibration/test target ratios: 60/10/10/20 percent on groups; actual sample fractions are logged. Quantitative outputs: log2-fold RMSE, MAE, R2, Pearson r, Spearman rho; per-drug performance; 90% interval coverage and mean interval width; selective-prediction retention/accepted RMSE; 300 group-cluster bootstrap resamples per drug-seed for approximate 95% test RMSE and baseline-improvement intervals. Multiple models evaluated against the exact same test split. Missing/insufficient drug sets are transparently skipped.

## Expected artifacts (without fabricating experiments)

`results/latest_report.html`, `results/summary.json`, 10 CSV tables, 7 publication-useful PNG/SVG figures, per-seed checkpoints and test prediction files, and `PAPER_RESULTS_SCAFFOLD.md`. All result numbers come from downloads and model outputs at run time; the ZIP ships with no prefilled experimental results.

## Novelties to evaluate critically

The proposed combination of held-out genotype novelty, calibrated uncertainty, and research abstention is *a proposed integrated study design*. Existing HIVDR ML, conformal prediction, and protein-language-model techniques need proper comparisons and attribution. Treat a mixed or negative result as publishable evidence, not an excuse to claim universal superiority.

## Publication and reproducibility checklist

[ ] Confirm venue still accepting papers, hybrid author presentation, publication agreement, fees and Scopus status independently.
[ ] Obtain original Stanford HIVDB data and preserve exact hashes.
[ ] Run all three seeds, check errors/skips, review every result and double-check leakage checks.
[ ] Validate whether reported fold measurement units and assay methodology support all claims.
[ ] Compare with recent 2026 HIV protein language model research and HIVdb algorithms; do not claim novelty without literature review.
[ ] Human author verifies design, interpretation, references and reproducibility; seek expert review from HIV pharmacology specialists if available.
[ ] Insert measured figures and tables in actual submission template; do not use placeholder results or claim therapeutic efficacy.
[ ] Document no patient efficacy/cure prediction and data/source restrictions.

## Sources

Stanford HIVDB: https://hivdb.stanford.edu/pages/genopheno.dataset.html  
NIH clinical HIV drug resistance guidance: https://hivinfo.nih.gov/understanding-hiv/fact-sheets/drug-resistance  
Farquhar, Bioinformatics (2026), protein language model benchmarking: https://pubmed.ncbi.nlm.nih.gov/42105210/  
WHO updated treatment guidance (2026): https://www.who.int/news/item/07-01-2026-who-releases-updated-recommendations-on-hiv-clinical-management
