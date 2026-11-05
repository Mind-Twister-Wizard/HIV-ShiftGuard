# HIV-ShiftGuard — Study design and reproducibility

**Current manuscript:** HIV-ShiftGuard: Mutation-Aware Ensemble Regression and an Uncertainty Audit of HIV-1 Antiretroviral Susceptibility.

**Prepared format:** IEEE ICAAIC 2027 author-named manuscript (six pages). This designation does not indicate conference acceptance. See [paper metadata](../paper/README.md).

## Research scope

Retrospective prediction of HIV-1 antiretroviral susceptibility measured in vitro as log2 fold-change. The model is not validated for patient outcomes, drug selection, treatment safety, or HIV cure prediction.

## Dataset and partitions

Four Stanford HIVDB genotype–phenotype files (PI, NRTI, NNRTI, INI) supplied **32,313 valid drug–isolate measurements**. The study reports **7,148** gene-genotype records and **5,043** linked groups. Exact sequence IDs and shared within-gene mutation signatures define non-overlapping split groups. GroupShuffleSplit uses seeds **11, 23, 37** and target proportions of **60/10/10/20** for training, validation, calibration, and test respectively. These splits do not establish patient-level or temporal independence.

## Models and evaluation

- Baselines: training-target median, mutation ridge, mutation-based extremely randomized trees, and a position-only variant.
- Ensemble: validation-tuned convex combination of ridge and mutation forest; selection uses validation error, not test labels.
- Prediction intervals: 90% nominal split conformal based on calibration residuals, compared with a Jaccard genotype-novelty-scaled version.
- Research review gate: uses calibration interval width and validation novelty thresholds; it is not a clinical recommendation mechanism.
- Evaluation: 67 completed drug–seed runs across 23 drugs; test RMSE, MAE, and additional archived metrics; 300 group-cluster bootstrap draws per run.

The position-only comparator also changes model capacity. It is an **exploratory comparison**, not a single-variable-controlled ablation.

## Recorded findings

| Model | Mean RMSE | Mean MAE |
|---|---:|---:|
| Median baseline | 2.334 | 1.872 |
| Mutation ridge | 1.113 | 0.824 |
| Position-only trees | 1.267 | 0.862 |
| Mutation forest | 1.111 | 0.759 |
| HIV-ShiftGuard ensemble | **1.017** | **0.719** |

At a nominal 90% target, pooled interval coverage was **90.12%** unscaled and **90.00%** novelty-scaled; the novelty-scaled intervals were slightly wider. The latter did **not** improve pooled coverage–width efficiency.

## Reproducible artifacts

[Verification script](../scripts/verify_results.py) · [Archived results](../results/) · [Run provenance](../results/PROVENANCE.json) · [One-click guide](REPRODUCING.md).

The recorded fingerprint is `42374487aabb5e1a`. Original data and per-observation records are not redistributed. No unverified new experiments or clinical claims are introduced here.
