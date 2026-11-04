# HIV-ShiftGuard  
### Mutation-novelty-aware uncertainty calibration for HIV-1 *in vitro* antiretroviral susceptibility

> **IMPORTANT SCIENTIFIC LIMITATION:** This is research software for predicting laboratory susceptibility measurements. It **does not** discover a cure, select an antiretroviral therapy, rank treatments for an individual, or replace qualified clinical care. All modeled fold-changes come from HIV-1 isolates, not patient treatment outcomes.

## Start in one click (Windows 11)

1. Extract this entire ZIP anywhere, for example `C:\Users\User\Desktop\Development\HIV_ShiftGuard_OneClick`.
2. Double-click **`RUN_ALL.bat`**. It checks for Python 3.10+, creates a private `.venv`, installs libraries, runs unit tests, downloads verified Stanford HIVDB data, and performs all research experiments.
3. Open **`results\latest_report.html`** when finished. Inspect `results\figures`, `results\tables`, and `results\PAPER_RESULTS_SCAFFOLD.md`.
4. If your laptop restarts, the battery dies or the Internet disconnects, **double-click `RUN_ALL.bat` again**. Completed drug/seed checkpoints are reused when source data, code, and settings are unchanged.

**Python auto-install:** If no usable Python is present, the launcher attempts `winget` Python 3.12 installation. Close and rerun the BAT after Windows updates its PATH. Python download, dependency installation, and Stanford data downloads require Internet access. It does not silently switch to toy data.

**No CUDA required:** The default experiments use scikit-learn on CPU, compatible with a laptop carrying an NVIDIA RTX 3050 Ti 4GB GPU. No 16GB VRAM language model is downloaded. Config has 3 seeds and capped ExtraTrees models for practical execution. The full runtime depends on your processor/Internet; no fixed duration is promised.

## Novel research *hypothesis* — not an already proven novelty claim

**Proposed framing:** Most in-vitro genotype–phenotype ML pipelines report a point estimate; this project asks whether predictions should additionally expose (i) distance from observed training mutation combinations, (ii) two conformal interval methods, and (iii) a controlled research-level review/abstention gate under unobserved genotype patterns. The study couples these elements with a strict exact-genotype/SeqID-connected group split and ablation analysis. Each component is established; **the research proposal is a particular integrated evaluation, not a claim that Jaccard novelty, ensembles, or conformal prediction were invented here.** Compare against recent 2026 HIV PLM benchmarks before claiming novelty or superiority.

## Verified input sources (Stanford HIV Drug Resistance Database)

The launcher downloads these **unmodified public research files** at runtime. Data are NOT included in the ZIP and may be subject to Stanford HIVDB usage terms and citation requirements. Keep source attribution.

| Gene / drug class | Direct source | Measurement |
|---|---|---|
| HIV-1 protease, PI | https://hivdb.stanford.edu/download/GenoPhenoDatasets/PI_DataSet.txt | Viral inhibitor fold susceptibility |
| HIV-1 RT, NRTI | https://hivdb.stanford.edu/download/GenoPhenoDatasets/NRTI_DataSet.txt | Viral inhibitor fold susceptibility |
| HIV-1 RT, NNRTI | https://hivdb.stanford.edu/download/GenoPhenoDatasets/NNRTI_DataSet.txt | Viral inhibitor fold susceptibility |
| HIV-1 integrase, INI | https://hivdb.stanford.edu/download/GenoPhenoDatasets/INI_DataSet.txt | Viral inhibitor fold susceptibility |

Download listing and documentation: https://hivdb.stanford.edu/pages/genopheno.dataset.html

**Please note**: Some drug measurements are missing; rows are not all assayed against all compounds. Source files include a mutation column per gene position, observed fold change, SeqID and `CompMutList`. Target is `log2(fold-change)`; high values can indicate decreased *in vitro* susceptibility, not automatic clinical failure. Older drugs in this research database are not current first-line treatment recommendations.

## Automatic research protocol

- **Ingestion:** HTTP streaming with retries and header checks; SHA256 and source timestamps recorded in `results/PROVENANCE.json`.
- **Quality control:** Nonpositive/non-numeric susceptibility removed; exact duplicate genotype/drug isolate records collapsed by median; missing residue calls and mixed variants excluded from mutation features.
- **Leakage auditing:** Connected component grouping of identical mutation signatures within genes and any shared Stanford SeqID across genes. Splits are train/validation/calibration/test with target 60%/10%/10%/20% *group-based*. Repeated genotype patterns never cross these four partitions.
- **Preprocessing:** Build mutation tokens (`P<position><mutant residue>`) using training-only feature vocabulary. No pretrained model or proprietary API needed.
- **Models:** Baseline per-drug training median; mutation ridge regression; position-only mutation ablation; nonlinear substitution-aware ExtraTrees regression; main validation-weighted ensemble (ridge + ExtraTrees). Each drug is evaluated separately; no patient-target model is built.
- **Evaluation:** RMSE, MAE, R2, Pearson, Spearman, per-drug 300-resample genotype-cluster bootstrap confidence intervals, test scatter, drug-specific results and three independent seeds. Tests are untouched by model selection.
- **Uncertainty:** Per-drug held-out split conformal intervals, nominal 90% target, plus scaled intervals based on nearest mutation Jaccard distance from training patterns. Quantiles fit **only on calibration data**. Under genotype shift, empirical coverage may fall below 90%; this is explicitly tested, not promised.
- **Review gate:** A research-only selective-prediction threshold from validation novelty and calibration interval width, with held-out retention and accepted-sample RMSE. It does not make clinical decisions.
- **Ablations:** Median/no genotype; position-only; linear mutation; nonlinear mutation; validation-weighted ensemble. Each method shares the exact test partition for its drug/seed.
- **Reporting:** CSV tables, 7 PNG + SVG figures including architecture, an HTML report, and a methods/results Markdown scaffold. Failed tasks are marked rather than replaced with imaginary results.

### Outputs

```
results/
  latest_report.html             <- browser report (automatic)
  summary.json
  PROVENANCE.json
  PAPER_RESULTS_SCAFFOLD.md    <- evidence-based draft scaffold; NOT a finished paper
  tables/                     <- 10 publication-useful CSVs
  figures/                    <- 7 scientific figures PNG and SVG
  runs/<fingerprint>/
    checkpoints/              <- fitted scikit-learn model bundles (.joblib)
    seeds/                    <- per-drug + seed CSV and JSON, fully resumable
    splits/                   <- every observation's grouping and partition
    tables/ figures/          <- immutable run-specific outputs
```

The run fingerprint includes every raw-file SHA256, Python source file and config JSON; changes intentionally begin a new run instead of silently mixing incompatible checkpoints. No dataset modification is needed between runs.

### Manual terminal operation

In the extracted project root, after environment installation:

```bat
set PYTHONPATH=%CD%\src
.venv\Scripts\python.exe -m unittest discover -s tests -v
.venv\Scripts\python.exe -m hivshiftguard.pipeline --config config.json
```

### When downloads fail

If your university Internet blocks Stanford, download the four original `.txt` files yourself into `data\raw\` using the URLs above, keeping the filenames. Then rerun the BAT. **No synthetic research dataset is substituted.** The unit tests have a tiny isolated temporary fixture; no temporary test measurements enter the results pipeline.

### Interruption and repairs

- Run the same BAT. It will skip completed models with valid per-drug JSON, CSV and joblib checkpoints.
- Failed drugs are listed in `results\tables\Table_06_Skipped_Or_Failed_Drugs.csv`.
- If you edit `config.json` or source files, a new run ID is created. Your older run stays intact in `results\runs\`.
- Do not rename the Stanford text files; use documented exact filenames.

### Scientific limitations and responsible use

1. Exact-genotype/SeqID groups may not completely eliminate patient/study overlap since no patient IDs or study batches are supplied.
2. Stanford filtered phenotypes are not clinical outcome labels. A resistance phenotype is neither a probability of cure nor a safe recommended medication.
3. The limited drug-specific training size, subtype mix, and missing phenotypes may restrict the findings. No claim of external clinical validation is made.
4. Conformal validity relies on assumptions, including exchangeability; it can fail under shifted test genotypes. Analyze observed intervals rather than assuming guaranteed coverage.
5. Quantitative novelty/superiority claims require a systematic review and careful comparison with existing drug-specific state-of-the-art models; an attractive architecture figure does not establish novelty.
6. Any therapeutic translation requires wet-lab research, broader cohorts, expert interpretation, and eventually appropriately regulated clinical validation. The HIV latent reservoir is not modeled; thus **this project cannot discover or establish an AIDS cure**.

### Citations / background

- Stanford HIV Drug Resistance Database: https://hivdb.stanford.edu/
- Stanford research phenotype data: https://hivdb.stanford.edu/pages/genopheno.dataset.html
- NIH HIV drug resistance fact sheet: https://hivinfo.nih.gov/understanding-hiv/fact-sheets/drug-resistance
- WHO 2026 ART update: https://www.who.int/news/item/07-01-2026-who-releases-updated-recommendations-on-hiv-clinical-management
- Recent benchmark (Farquhar 2026): https://pubmed.ncbi.nlm.nih.gov/42105210/
- NIH 2026 ART guidelines: https://clinicalinfo.hiv.gov/en/guidelines/hiv-clinical-guidelines-adult-and-adolescent-arv/what-start-initial-combination-regimens

**License:** Research scripts are made available under the MIT license in `LICENSE`. Stanford data are not bundled and retain their own applicable terms.
