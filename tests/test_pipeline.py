import math
import tempfile
import unittest
from pathlib import Path

import numpy as np
import pandas as pd

from hivshiftguard.data import mutation_tokens, load_stanford_data, split_grouped
from hivshiftguard.io import SOURCES, validate_header
from hivshiftguard.ml import conformal_q, metrics, train_drug
from hivshiftguard.report import reduce_runs


class ResearchPipelineTests(unittest.TestCase):
    def test_mutation_columns_ignore_missing(self):
        mutation, missing = mutation_tokens({"P1":"-", "P2":"K", "P3":".", "P4":"KV"}, ["P1", "P2", "P3", "P4"])
        self.assertEqual(mutation, ("2K",))
        self.assertEqual(missing, 2)

    def test_conformal_quantile(self):
        self.assertEqual(conformal_q([1,2,3,4,5], .2), 5.)
        self.assertEqual(conformal_q([1,2,3,4], .50), 3.)

    def test_same_genotype_group_does_not_leak(self):
        df = pd.DataFrame({"split_group": [k for k in range(100) for _ in range(3)]})
        split, _ = split_grouped(df, 42)
        groups = [set(df.split_group[split == name]) for name in ("train", "validation", "calibration", "test")]
        for i in range(4):
            for j in range(i + 1, 4):
                self.assertFalse(groups[i] & groups[j])
        self.assertTrue(all(len(s) > 0 for s in groups))

    def test_full_pipeline_with_fixture_is_not_released_as_real_data(self):
        # Synthetic values created in TEMPORARY folder for self-test ONLY; never written
        # into the research data/ or results/ directories.
        with tempfile.TemporaryDirectory() as tmp:
            root = Path(tmp)
            raw = root / "data" / "raw"
            raw.mkdir(parents=True)
            genes = {"PI": ["ATV","DRV"], "NRTI": ["3TC","TDF"],
                     "NNRTI": ["EFV","ETR"], "INI": ["DTG","RAL"]}
            for gene, drugs in genes.items():
                rng = np.random.default_rng(12 + list(genes).index(gene))
                entries = []
                for i in range(140):
                    a, b, c = ("A" if int(i & 1) else "-", "V" if int(i & 2) else "-",
                               "T" if int(i & 4) else "-")
                    # Make each genotype distinct so group holdouts have sufficient groups.
                    d = str(["K","R","Q","M"][i % 4])
                    e = str(["D","E","I","N","S"][i % 5])
                    f = str(["F","H","L","P","Y","W","G"][i % 7])
                    x = 0.4 * (a != "-") + 0.7 * (b != "-") + .2 * (c != "-") + rng.normal(0,.2)
                    r = {"SeqID": str(i + 1000), **{drug: f"{max(.1, 2**(x + j * .1)):.3f}" for j, drug in enumerate(drugs)},
                         "P1": a, "P2": b, "P3": c, "P4": d, "P5": e, "P6": f, **{f"P{k}": "-" for k in range(7, 21)}, "CompMutList": ""}
                    entries.append(r)
                pd.DataFrame(entries).to_csv(raw / f"{gene}_DataSet.txt", sep="\t", index=False)
                validate_header(raw / f"{gene}_DataSet.txt")
            obs, meta = load_stanford_data(root)
            self.assertEqual(meta["valid_pairs"], 1120)
            self.assertEqual(len(obs), 1120)
            assignment, stats = split_grouped(obs, 11)
            frame = obs[(obs.gene == "PI") & (obs.drug == "DRV")]
            cfg = {"minimum_train_samples": 25, "minimum_calibration_samples": 5,
                   "minimum_test_samples": 5, "forest_trees": 16, "forest_max_depth": 8,
                   "min_samples_leaf": 2, "coverage_target": .9, "num_threads": 1, "novelty_scale_gamma": 1.0}
            result, preds, error = train_drug(frame, assignment[frame.index], cfg, 11, root / "checkpoint.joblib")
            self.assertIsNone(error)
            self.assertTrue((root / "checkpoint.joblib").exists())
            self.assertEqual(len(preds), result["split_counts"]["test"])
            self.assertTrue(np.isfinite(preds["predicted_log2fc"]).all())
            self.assertTrue(result["models"]["shiftguard_ensemble"]["rmse"] >= 0)
            config = {"coverage_target": .9}
            report = reduce_runs(root / "reports", [result], [preds], meta, config, [])
            self.assertTrue((root / "reports" / "FINAL_REPORT.html").exists())
            self.assertTrue((root / "reports" / "figures" / "Figure_01_Research_Architecture.svg").exists())
            self.assertEqual(report["completed_experiments"], 1)


if __name__ == "__main__":
    unittest.main()
