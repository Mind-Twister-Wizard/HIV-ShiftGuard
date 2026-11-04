"""Robust restartable one-click research pipeline. No generated / synthetic research findings."""
import argparse
import hashlib
import json
import shutil
import sys
import time
import traceback
from datetime import datetime, timezone
from pathlib import Path

import pandas as pd

from . import __version__
from .data import load_stanford_data, save_prepared_data, split_grouped
from .io import atomic_json, ensure_data, sha256_file
from .ml import train_drug
from .report import reduce_runs


def run(root, cfg):
    t0 = time.monotonic()
    print("\nHIV-SHIFTGUARD | Verifiable in-vitro research, NOT a clinical or cure prediction tool\n", flush=True)
    records = ensure_data(root, cfg)
    obs, provenance = load_stanford_data(root)
    save_prepared_data(root, obs, provenance)
    fingerprint = hashlib.sha256()
    for record in records:
        fingerprint.update(record["sha256"].encode())
    for script in sorted((root / "src" / "hivshiftguard").glob("*.py")):
        fingerprint.update(sha256_file(script).encode())
    fingerprint.update(json.dumps(cfg, sort_keys=True).encode())
    signature = fingerprint.hexdigest()[:16]
    run_dir = root / "results" / "runs" / signature
    run_dir.mkdir(parents=True, exist_ok=True)
    atomic_json(run_dir / "PROVENANCE.json", {"package_version": __version__,
                                                 "software_and_data_fingerprint": signature,
                                                 "input_sources": records, "data_audit": provenance,
                                                 "configuration": cfg,
                                                 "started_or_resumed_utc": datetime.now(timezone.utc).isoformat()})
    print(f"Valid observed pairs: {len(obs):,}; sequence groups: {provenance['connected_isolate_groups']:,}")
    print(f"Reproducible run ID: {signature} (will reuse valid prior checkpoints)", flush=True)
    completed, predictions, failures = [], [], []
    count_total = 0
    for seed in cfg["seeds"]:
        assignment, split_stats = split_grouped(obs, seed, (cfg["test_fraction"], cfg["validation_fraction"], cfg["calibration_fraction"]))
        splits = obs[["observation_id", "gene", "drug", "seqid", "sample_id", "split_group"]].copy()
        splits["partition"] = assignment
        target = run_dir / "splits" / f"seed_{seed}_split_assignments.csv"
        target.parent.mkdir(parents=True, exist_ok=True)
        if not target.exists():
            splits.to_csv(target, index=False)
        print(f"\n[SEED {seed}] {split_stats}", flush=True)
        for (gene, drug), frame in obs.groupby(["gene", "drug"], sort=True):
            if len(frame) < cfg["minimum_samples_per_drug"]:
                failures.append({"seed": seed, "gene": gene, "drug": drug,
                                 "status": "SKIPPED", "reason": f"Only {len(frame)} valid samples"})
                continue
            count_total += 1
            base = run_dir / "seeds" / str(seed)
            base.mkdir(parents=True, exist_ok=True)
            tag = f"{gene}_{drug}".replace("/", "_")
            metric_path, pred_path = base / (tag + ".json"), base / (tag + "_predictions.csv")
            checkpoint = run_dir / "checkpoints" / f"seed_{seed}" / (tag + ".joblib")
            if metric_path.exists() and pred_path.exists() and checkpoint.exists():
                try:
                    obj = json.loads(metric_path.read_text(encoding="utf-8"))
                    pred = pd.read_csv(pred_path)
                    if len(pred) == obj["split_counts"]["test"] and len(pred):
                        completed.append(obj)
                        predictions.append(pred)
                        print(f"[RESUME SKIP] {seed} {tag} valid outputs exist", flush=True)
                        continue
                except Exception:
                    print(f"[RESUME FIX] {seed} {tag} damaged output will be rebuilt", flush=True)
            print(f"[TRAIN] seed {seed}  {tag}  n={len(frame)}", flush=True)
            try:
                local_labels = assignment[frame.index.to_numpy()]
                outcome, preds, reason = train_drug(frame, local_labels, cfg, seed, checkpoint)
                if reason:
                    failures.append({"seed": seed, "gene": gene, "drug": drug,
                                     "status": "SKIPPED", "reason": reason})
                    print(f"[SKIP] {tag}: {reason}", flush=True)
                    continue
                tmp = pred_path.with_suffix(".csv.tmp")
                preds.to_csv(tmp, index=False)
                tmp.replace(pred_path)
                atomic_json(metric_path, outcome)  # write JSON last to signal completion
                completed.append(outcome)
                predictions.append(preds)
                print(f"[DONE] {tag}  heldout RMSE={outcome['models']['shiftguard_ensemble']['rmse']:.4f}"
                      f"  novelty interval coverage={outcome['intervals']['novelty_scaled_coverage']:.1%}", flush=True)
            except Exception as exc:
                error = {"seed": seed, "gene": gene, "drug": drug, "status": "FAILED",
                         "reason": repr(exc), "traceback": traceback.format_exc()}
                failures.append(error)
                print(f"[ERROR] {tag}: {exc}", file=sys.stderr, flush=True)
                (base / (tag + ".error.txt")).write_text(error["traceback"], encoding="utf-8")
    if not completed:
        raise RuntimeError("No experiment succeeded. Check data, split sizes, and errors; no figures will be fabricated.")
    details = reduce_runs(run_dir, completed, predictions, provenance, cfg, failures)
    details["software_and_data_fingerprint"] = signature
    details["elapsed_seconds"] = round(time.monotonic() - t0, 2)
    atomic_json(run_dir / "summary.json", details)
    latest = root / "results"
    # Preserve full run folder and copy report/graphics/tables to accessible top-level paths.
    for name in ("figures", "tables"):
        dst = latest / name
        if dst.exists():
            shutil.rmtree(dst)
        shutil.copytree(run_dir / name, dst)
    for name in ("FINAL_REPORT.html", "PAPER_RESULTS_SCAFFOLD.md", "PROVENANCE.json", "summary.json"):
        shutil.copy2(run_dir / name, latest / ("latest_report.html" if name == "FINAL_REPORT.html" else name))
    atomic_json(latest / "latest_run.json", {"run_id": signature,
                                               "full_run_directory": str(run_dir), "summary": details})
    print("\n" + "=" * 65)
    print(f"SUCCESS: {len(completed)}/{count_total} eligible drug/seed experiments produced results.")
    if failures:
        print(f"WARNING: {len(failures)} skipped/failed; inspect tables/Table_06_Skipped_Or_Failed_Drugs.csv")
    print(f"Report: {latest / 'latest_report.html'}")
    print(f"Tables: {latest / 'tables'}")
    print(f"Figures: {latest / 'figures'}")
    print(f"Manuscript scaffold: {latest / 'PAPER_RESULTS_SCAFFOLD.md'}")
    print(f"Resumable model checkpoints: {run_dir / 'checkpoints'}")
    print("No clinical recommendations or HIV-cure claims are generated.", flush=True)
    return details


class TeeStream:
    def __init__(self, console, logfile):
        self.console = console
        self.logfile = logfile

    def write(self, text):
        self.console.write(text)
        self.logfile.write(text)
        self.logfile.flush()

    def flush(self):
        self.console.flush()
        self.logfile.flush()

    def isatty(self):
        return self.console.isatty()


def main():
    ap = argparse.ArgumentParser(description="HIV-ShiftGuard reproducible and resumable research pipeline")
    ap.add_argument("--config", default="config.json")
    args = ap.parse_args()
    root = Path(__file__).resolve().parents[2]
    cfg = json.loads((root / args.config).read_text(encoding="utf-8"))
    if cfg.get("mode") != "research_only":
        raise ValueError("Only research_only mode is supported. This is not a clinical tool.")
    log_path = root / "results" / "run.log"
    log_path.parent.mkdir(parents=True, exist_ok=True)
    with log_path.open("a", encoding="utf-8") as log:
        old_out, old_err = sys.stdout, sys.stderr
        sys.stdout, sys.stderr = TeeStream(old_out, log), TeeStream(old_err, log)
        try:
            print("\n" + "-" * 70 + "\nSession started: " + datetime.now(timezone.utc).isoformat(), flush=True)
            run(root, cfg)
        except Exception:
            print(traceback.format_exc(), file=sys.stderr, flush=True)
            raise
        finally:
            sys.stdout, sys.stderr = old_out, old_err


if __name__ == "__main__":
    main()
