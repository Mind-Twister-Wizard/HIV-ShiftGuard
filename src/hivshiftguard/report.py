import html
import json
from pathlib import Path

import numpy as np
import pandas as pd

from .figures import build_figures

VARIANTS = ["median_baseline", "linear_mutations", "position_only_ablation", "mutation_forest", "shiftguard_ensemble"]


def reduce_runs(run_dir, successes, preds, provenance, cfg, failures):
    run_dir = Path(run_dir)
    tables = run_dir / "tables"; figures = run_dir / "figures"
    tables.mkdir(exist_ok=True, parents=True)
    figure_dir = figures
    raw = pd.DataFrame(successes)
    pred = pd.concat(preds, ignore_index=True) if preds else pd.DataFrame()
    metric_records, interval_records, novelty_records = [], [], []
    for r in successes:
        for model, m in r["models"].items():
            metric_records.append({"seed": r["seed"], "gene": r["gene"], "drug": r["drug"],
                                   "model": model, **m, **{f"n_{key}": num for key, num in r["split_counts"].items()}})
        interval_records.append({"seed": r["seed"], "gene": r["gene"], "drug": r["drug"], **r["intervals"]})
        novelty_records.append({"seed": r["seed"], "gene": r["gene"], "drug": r["drug"], **r["novelty"]})
    bootstrap_records = [{"seed": r["seed"], "gene": r["gene"], "drug": r["drug"],
                          **(r["group_bootstrap_95ci"] or {})} for r in successes
                         if r.get("group_bootstrap_95ci")]
    pd.DataFrame(bootstrap_records).to_csv(tables / "Table_10_Group_Cluster_Bootstrap_CI.csv", index=False)
    metrics_df = pd.DataFrame(metric_records)
    intervals_df = pd.DataFrame(interval_records)
    novelty_df = pd.DataFrame(novelty_records)
    metrics_df.to_csv(tables / "Table_01_Evaluation_By_Drug_And_Seed.csv", index=False)
    intervals_df.to_csv(tables / "Table_02_Conformal_Diagnostics.csv", index=False)
    novelty_df.to_csv(tables / "Table_03_Novelty_Review_Gate.csv", index=False)
    pd.DataFrame(provenance["datasets"]).to_csv(tables / "Table_04_Dataset_Inventory.csv", index=False)
    if successes:
        pd.DataFrame([{"seed": r["seed"], "gene": r["gene"], "drug": r["drug"],
                       "validation_optimal_ridge_weight": r["ridge_weight_selected_on_validation"],
                       "train_samples": r["split_counts"]["train"],
                       "mutation_features": r["mutation_feature_count"]} for r in successes]).to_csv(
            tables / "Table_05_Model_Selections.csv", index=False)
    pd.DataFrame(failures if failures else [{"status": "No failures"}]).to_csv(tables / "Table_06_Skipped_Or_Failed_Drugs.csv", index=False)
    abl = (metrics_df.groupby("model", as_index=False).agg(
           mean_rmse=("rmse", "mean"), sd_rmse=("rmse", "std"), mean_mae=("mae", "mean"),
           experiments=("rmse", "size")))
    abl["model"] = pd.Categorical(abl["model"], VARIANTS, ordered=True)
    abl = abl.sort_values("model")
    abl.to_csv(tables / "Table_07_Ablation_Summary.csv", index=False)
    group_summary = []
    for (seed, gene, drug), subset in pred.groupby(["seed", "gene", "drug"]):
        threshold = subset["novelty_jaccard"].median()
        for subgroup, p in (("lower_novelty", subset[subset.novelty_jaccard <= threshold]),
                            ("higher_novelty", subset[subset.novelty_jaccard > threshold])):
            if len(p):
                group_summary.append({"seed": seed, "gene": gene, "drug": drug, "group": subgroup,
                                      "n": len(p), "rmse": float(np.sqrt(np.mean((p.observed_log2fc - p.predicted_log2fc)**2))),
                                      "mean_novelty": float(p.novelty_jaccard.mean())})
    pd.DataFrame(group_summary).to_csv(tables / "Table_08_Novelty_Subgroup_Sensitivity.csv", index=False)
    pred.to_csv(tables / "Table_09_All_Heldout_Predictions.csv", index=False)
    build_figures(figure_dir, metrics_df, pred)
    primary = abl[abl["model"] == "shiftguard_ensemble"].iloc[0]
    baseline = abl[abl["model"] == "median_baseline"].iloc[0]
    coverage = np.mean((pred.observed_log2fc >= pred.interval_lower_novelty) &
                       (pred.observed_log2fc <= pred.interval_upper_novelty))
    standard = np.mean((pred.observed_log2fc >= pred.interval_lower_unscaled) &
                       (pred.observed_log2fc <= pred.interval_upper_unscaled))
    details = {"completed_experiments": int(len(successes)), "distinct_drugs": int(len(set(r["drug"] for r in successes))),
               "seeds": sorted(set(int(r["seed"]) for r in successes)), "observations": provenance["valid_pairs"],
               "genotypes": provenance["unique_gene_genotypes"], "groups": provenance["connected_isolate_groups"],
               "mean_drug_seed_rmse": float(primary["mean_rmse"]),
               "mean_drug_seed_median_baseline_rmse": float(baseline["mean_rmse"]),
               "pooled_nominal_target": cfg["coverage_target"], "pooled_standard_coverage": float(standard),
               "pooled_novelty_coverage": float(coverage), "test_predictions": len(pred),
               "number_of_failures": len(failures)}
    (run_dir / "summary.json").write_text(json.dumps(details, indent=2), encoding="utf-8")
    write_html(run_dir, details, abl, metrics_df, intervals_df, novelty_df, failures)
    write_manuscript_scaffold(run_dir, details, abl, cfg)
    return details


def frame_to_html(df, limit=25):
    return df.head(limit).to_html(index=False, classes="paper-table", border=0,
                                float_format=lambda x: f"{x:.4f}", escape=True)


def write_html(run_dir, result, abl, metrics, interval, novelty, failures):
    intro = "<strong>Research only:</strong> susceptibility measurements are in vitro. These predictions do not identify an HIV cure, select a clinical regimen, or imply safety or effectiveness for an individual."
    figures = ["Figure_01_Research_Architecture", "Figure_02_Model_Ablation_RMSE",
               "Figure_03_Drugwise_Test_RMSE", "Figure_04_Predicted_vs_Observed",
               "Figure_05_Interval_Coverage", "Figure_06_Error_vs_Novelty", "Figure_07_Novelty_Quintiles"]
    gallery = "\n".join(f'<figure><img src="figures/{html.escape(name)}.png" alt="{html.escape(name)}"><figcaption>{html.escape(name.replace("_", " "))}</figcaption></figure>' for name in figures)
    script = f"""<!doctype html><html lang='en'><head><meta charset='utf-8'>
<meta name='viewport' content='width=device-width, initial-scale=1'><title>HIV-ShiftGuard | Research audit</title>
<style>
:root{{--ink:#172c45;--accent:#177b90;--bg:#f3f7fa;--muted:#506477}}
*{{box-sizing:border-box}}body{{margin:0;background:var(--bg);color:var(--ink);font:15px/1.55 Segoe UI,Arial,sans-serif}}
main{{max-width:1200px;margin:auto;padding:0 25px 50px}}
header{{background:linear-gradient(110deg,#152943,#14607c);color:white;padding:47px 28px;margin-bottom:26px}}
header .inner{{max-width:1160px;margin:auto}}h1{{font-size:34px;margin:3px 0}}h2{{font-size:21px;margin-top:10px}}h3{{font-size:16px}}
.subtitle{{color:#b8d8e5}}.notice{{border-left:4px solid #c47d34;background:#fff5e9;padding:14px 20px;border-radius:5px}}
.card{{background:white;border:1px solid #e3eaf0;border-radius:12px;padding:24px;margin:22px 0;box-shadow:0 3px 14px #233d5210}}
.stats{{display:grid;grid-template-columns:repeat(auto-fit,minmax(190px,1fr));gap:14px}}
.stat{{background:white;border-radius:9px;padding:17px;border:1px solid #d9e3ec}}.stat strong{{display:block;font-size:25px;color:#176c8b}}
.stat span{{font-size:12px;color:#607487}}figure{{margin:20px auto;max-width:950px;text-align:center}}
figure img{{width:100%;height:auto;border:1px solid #e5ecf2;border-radius:6px}}figcaption{{color:#657789;font-size:12px}}
.paper-table{{border-collapse:collapse;width:100%;font-size:13px}}.paper-table th{{background:#e9f1f6;text-align:left}}
.paper-table th,.paper-table td{{padding:8px;border-bottom:1px solid #e3eaf0}}.scroll{{overflow-x:auto}}
a{{color:#176d95}}.foot{{font-size:12px;color:#637588}}code{{background:#eaf0f4;padding:2px 5px}}
</style></head><body><header><div class='inner'><div class='subtitle'>REPRODUCIBLE RESEARCH • AUTOMATICALLY COMPUTED FROM REAL DATA</div>
<h1>HIV-ShiftGuard</h1><p>Mutation-novelty-gated uncertainty analysis of HIV-1 antiretroviral drug susceptibility</p></div></header>
<main><div class='notice'>{intro}</div><section class='stats'>
<div class='stat'><strong>{result['observations']:,}</strong><span>Real observed genotype–drug pairs</span></div>
<div class='stat'><strong>{result['distinct_drugs']}</strong><span>Evaluated antiretroviral drugs</span></div>
<div class='stat'><strong>{result['mean_drug_seed_rmse']:.3f}</strong><span>ShiftGuard mean per-experiment RMSE (log2 fold)</span></div>
<div class='stat'><strong>{result['pooled_novelty_coverage']:.1%}</strong><span>Observed coverage, novelty-scaled intervals</span></div>
<div class='stat'><strong>{result['completed_experiments']}</strong><span>Completed drug-by-seed experiments</span></div></section>
<section class='card'><h2>Research protocol</h2><p>Stanford HIVDB filtered in-vitro phenotype datasets; drug-specific models; four disjoint partitions
(train, validation, calibration, test), grouped by shared SeqID and exact within-gene mutation signatures. The
unseen test partition is never used to select fusion weights, intervals or novelty thresholds. Prediction target:
log2 of reported fold-change in susceptibility versus wild-type reference.</p>
<p>Comparators: training-set median; mutation ridge; position-only ExtraTrees ablation; full-mutation ExtraTrees;
validation-weighted ensemble. Two conformal approaches are compared: constant-width and genotype-novelty-scaled.
No exchangeability guarantee is asserted under distribution shift.</p>
<p>Connected groups guard against exact-genotype leakage but cannot guarantee patient-disjoint splits because
patient identifiers are not present. Reported performance is isolate-level in vitro only.</p></section>
<section class='card'><h2>Table: model / ablation results</h2><div class='scroll'>{frame_to_html(abl)}</div></section>
<section class='card'><h2>Table: test RMSE by drug (first 40)</h2><div class='scroll'>{frame_to_html(metrics[metrics.model=='shiftguard_ensemble'].sort_values(['gene','drug','seed']),40)}</div></section>
<section class='card'><h2>Uncertainty summary</h2><p>Target marginal coverage: {result['pooled_nominal_target']:.0%}; observed standard intervals: {result['pooled_standard_coverage']:.1%}; observed novelty-scaled intervals: {result['pooled_novelty_coverage']:.1%}. A difference in width and coverage must both be assessed.</p><div class='scroll'>{frame_to_html(interval,12)}</div></section>
<section class='card'><h2>Figures (PNG + editable vector SVG)</h2>{gallery}</section>
<section class='card'><h2>Artifacts and checks</h2><p>All numerical results are in <code>tables/</code>; split assignments and sample-level held-out predictions
are preserved. Each fitted drug-level model has a checkpoint in <code>checkpoints/</code>. Use the checkpoint directory with care: it is research software, not a clinical device.</p>
<p>Skipped or failed experiments: {len(failures)}. Investigate any skips before making manuscript claims. Models should be compared within the same drug/seed; averaging across drugs gives equal weight to each experiment and differs from pooled sample-level RMSE.</p>
<p class='foot'>Sources: <a href='https://hivdb.stanford.edu/pages/genopheno.dataset.html'>Stanford HIVDB</a> •
<a href='https://hivinfo.nih.gov/understanding-hiv/fact-sheets/drug-resistance'>NIH HIV drug resistance</a> •
<a href='https://clinicalinfo.hiv.gov/en/guidelines/hiv-clinical-guidelines-adult-and-adolescent-arv/what-start-initial-combination-regimens'>NIH ART guidelines</a>.</p></section></main></body></html>"""
    (run_dir / "FINAL_REPORT.html").write_text(script, encoding="utf-8")


def write_manuscript_scaffold(run_dir, result, abl, cfg):
    model_rows = "\n".join(f"| {str(r.model)} | {r.mean_rmse:.3f} | {r.sd_rmse:.3f} | {r.mean_mae:.3f} | {int(r.experiments)} |" for _, r in abl.iterrows())
    manuscript = f"""# HIV-ShiftGuard: Mutation-Novelty-Aware Uncertainty-Calibrated Prediction of HIV-1 Antiretroviral Susceptibility

> **AUTO-GENERATED METHODS AND RESULTS SCAFFOLD. NOT A SUBMISSION-READY PAPER.**
> Verify experiment quality, compare literature and write scholarly narrative before submission.

## Abstract — to be written after validation

Background: Existing HIV-1 genotype-to-phenotype prediction can be overconfident on unfamiliar viral mutation patterns.
Methods: Drug-specific mutation-based models were benchmarked against baselines using four grouped partitions and three independent seeds, with mutation-novelty-aware conformal inference and selective review.
Results: {result['observations']} observed genotype–drug pairs, {result['distinct_drugs']} evaluated drugs, ensemble mean per-drug-seed RMSE {result['mean_drug_seed_rmse']:.3f} log2 fold; observed pooled novelty-scaled interval coverage {result['pooled_novelty_coverage']:.1%}. These are in-vitro outcomes only.
Conclusion: Report whether improvement over all baselines is actually supported by the per-drug tables and discuss failures. **Never imply clinical treatment selection or HIV cure discovery.**

## Introduction
Describe established ART and the limits of existing genotypic resistance prediction. Position this work against recent protein language model baselines, including 2026 HIV-ESM-2 research; do not claim all component methods are novel.

## Data
Stanford HIVDB genotype-phenotype research datasets: PI, NRTI, NNRTI and INI (direct links in README). Fold susceptibility is laboratory IC50 fold-change, not patient response. For the retrieved versions, {result['observations']} valid drug-isolate measurements were retained and {result['genotypes']} gene-genotype isolates identified. Source files are SHA256-tracked.

## Methods
Mutation substitution tokens encode each observed position (ignore missing '.' and ambiguous multi-amino-acid calls). For each drug independently, a linear ridge baseline, full-substitution ExtraTrees regressor, and position-only ablation were fitted. A fusion weight was chosen using validation error, not calibration or test data. Exact-genotype / sequence-ID linked components ensure no exact overlap between train, validation, calibration, and test groups. The four sets are nominally 60/10/10/20 percent but may differ in row-level sizes because the splits are group-based. Fold-change outcomes were transformed as log2(FC) and evaluated with RMSE, MAE, R², Pearson r and Spearman rho. The unscaled split-conformal residual quantile and novelty-scaled variant were calibrated on separate held-out calibration samples. Novelty is the nearest-training-mutation Jaccard distance. The novelty gate is an experimental abstention analysis and not a clinical decision rule.

## Results: ablation (mean ± SD across drug-seed runs)

| Model | RMSE mean | RMSE SD | MAE mean | Runs |
|---|---:|---:|---:|---:|
{model_rows}

Observed pooled standard interval coverage: {result['pooled_standard_coverage']:.1%}; novelty-scaled coverage: {result['pooled_novelty_coverage']:.1%} (nominal {cfg['coverage_target']:.0%}). Do not claim novelty is superior unless both error and interval efficiency support the claim, with uncertainty intervals and appropriate comparisons.

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
"""
    (run_dir / "PAPER_RESULTS_SCAFFOLD.md").write_text(manuscript, encoding="utf-8")
