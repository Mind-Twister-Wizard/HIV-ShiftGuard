# HIV-ShiftGuard — IEEE Manuscript

**Title:** HIV-ShiftGuard: Mutation-Aware Ensemble Regression and an Uncertainty Audit of HIV-1 Antiretroviral Susceptibility

**Version:** Author-named IEEE-format manuscript prepared for ICAAIC 2027 (six pages). **No conference acceptance or publication is claimed.**

## Authors (manuscript order)

1. **Jasmine Kaur** — jaledarjasmine@gmail.com
2. **Hossein Tabasi** — Hosseintabasi@shooliniuniversity.com
3. **Alok Raj** — corresponding author, alokraj2027@gmail.com · [ORCID 0009-0003-7870-7065](https://orcid.org/0009-0003-7870-7065)
4. **Pankaj Vaidya** — pankaj.vaidya@shooliniuniversity.com · [ORCID 0000-0003-1304-630X](https://orcid.org/0000-0003-1304-630X)
5. **Anurag Rana** — anuragrana.anu@gmail.com · [ORCID 0000-0003-0247-8908](https://orcid.org/0000-0003-0247-8908)

**Shared affiliation:** Yogananda School of AI, Computers and Data Sciences, Shoolini University of Biotechnology and Management Sciences, Solan, Himachal Pradesh, India.

## Abstract

Antiretroviral resistance complicates HIV treatment, but in vitro susceptibility does not predict clinical outcomes. This study evaluates HIV-ShiftGuard, a mutation-aware regression ensemble with genotype-novelty-dependent intervals and a research-only review gate. Four Stanford HIV Drug Resistance Database genotype–phenotype files provided 32,313 valid drug–isolate observations. Drug-specific ridge and extremely randomized tree models were fused using validation-selected weights and compared with median, single-model, and position-only alternatives. Across 67 drug–seed evaluations spanning 23 drugs, the ensemble achieved a mean test RMSE of 1.017 log2 fold-change units, compared with 1.113 for mutation ridge, 1.111 for mutation forest, and 2.334 for the median baseline. Unscaled conformal intervals achieved 90.12% pooled coverage at a nominal 90% level, compared with 90.00% for novelty-scaled intervals; the latter were slightly wider. These findings support mutation-aware ensemble benchmarking but do not establish clinical efficacy, improved interval efficiency, or a cure for HIV.

**Keywords:** HIV-1 drug resistance; genotype–phenotype prediction; ensemble learning; conformal prediction; uncertainty quantification.

## Study artifacts

- [Verified results](../results/) · [Methods](../docs/RESEARCH_PROTOCOL.md) · [Reproducibility](../docs/REPRODUCING.md) · [Citation metadata](../CITATION.cff).
- The existing `HIV_ShiftGuard_ICIRSET2026_Manuscript.docx` and `.pdf` files are **historical manuscript snapshots**, not the revised IEEE submission. Their older author lists are part of that archive.
- The revised IEEE Word and PDF documents are prepared separately and are **not yet uploaded to this directory**. When published here, use filenames `HIV_ShiftGuard_ICAAIC_2027_IEEE_Submission.docx` and `HIV_ShiftGuard_ICAAIC_2027_IEEE_Submission.pdf`.

The repository records experimental software and retrospective *in vitro* results; it is not a clinical decision system.
