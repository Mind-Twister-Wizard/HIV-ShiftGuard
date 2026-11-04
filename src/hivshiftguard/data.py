"""Load observed *in-vitro* susceptibility, not clinical treatment outcomes."""
import hashlib
import json
import re
from collections import defaultdict
from pathlib import Path

import numpy as np
import pandas as pd
from sklearn.model_selection import GroupShuffleSplit

from .io import SOURCES

MUTATION_COLUMN = re.compile(r"^P([1-9][0-9]*)$")
AMINO_ACID = re.compile(r"^[ACDEFGHIKLMNPQRSTVWY]$")


class DSU:
    def __init__(self, n):
        self.parent = list(range(n))
        self.rank = [0] * n

    def find(self, x):
        while self.parent[x] != x:
            self.parent[x] = self.parent[self.parent[x]]
            x = self.parent[x]
        return x

    def union(self, a, b):
        a, b = self.find(a), self.find(b)
        if a == b:
            return
        if self.rank[a] < self.rank[b]:
            a, b = b, a
        self.parent[b] = a
        if self.rank[a] == self.rank[b]:
            self.rank[a] += 1


def mutation_tokens(row, positions):
    """Treat '-' as reference and '.' / mixtures as unavailable, not mutations."""
    found = []
    unavailable = 0
    for column in positions:
        amino = str(row[column]).strip().upper()
        if amino == "-":
            continue
        if AMINO_ACID.fullmatch(amino):
            found.append(f"{column[1:]}{amino}")
        else:
            unavailable += 1
    return tuple(sorted(found)), unavailable


def load_stanford_data(root):
    """Parse Stanford high-quality filtered GenoPhenoDatasets. No synthetic fallback."""
    observations = []
    dataset_info = []
    duplicates = 0
    for gene in SOURCES:
        path = Path(root) / "data" / "raw" / f"{gene}_DataSet.txt"
        df = pd.read_csv(path, sep="\t", dtype=str, keep_default_na=False, encoding="utf-8-sig")
        if "SeqID" not in df.columns:
            raise ValueError(f"Missing SeqID in {path.name}")
        positions = [c for c in df.columns if MUTATION_COLUMN.fullmatch(c)]
        drugs = list(df.columns[1:df.columns.get_loc(positions[0])])
        if not positions or not drugs or "CompMutList" not in df.columns:
            raise ValueError(f"Unrecognized Stanford file layout for {gene}: {path.name}")
        before = len(observations)
        for idx, row in df.iterrows():
            seqid = str(row["SeqID"]).strip()
            if not seqid:
                continue
            muts, missing = mutation_tokens(row, positions)
            mutation_key = ",".join(muts)
            # Identical genotypes must remain in the same split, even if SeqID differs.
            sample_id = f"{gene}:{seqid}:{hashlib.sha1(mutation_key.encode()).hexdigest()[:12]}"
            for drug in drugs:
                fold = pd.to_numeric(row[drug], errors="coerce")
                if pd.isna(fold) or not np.isfinite(fold) or fold <= 0:
                    continue
                observations.append({"gene": gene, "drug": drug, "seqid": seqid,
                                     "sample_id": sample_id, "signature": mutation_key,
                                     "mutations": muts, "num_mutations": len(muts),
                                     "unavailable_positions": missing, "fold": float(fold),
                                     "log2fc": float(np.log2(fold))})
        dataset_info.append({"gene": gene, "source_rows": len(df), "drug_columns": drugs,
                             "observed_drug_isolate_pairs": len(observations) - before})
    if not observations:
        raise ValueError("No valid observed genotype/phenotype pairs. Source data unavailable or malformed.")
    obs = pd.DataFrame(observations)
    # One genotype-isolate/drug label: median of exact duplicates; do not leak duplicates across splits.
    pre = len(obs)
    obs = (obs.groupby(["gene", "drug", "seqid", "sample_id", "signature"], as_index=False)
              .agg(fold=("fold", "median"), num_replicates=("fold", "size"),
                   unavailable_positions=("unavailable_positions", "max")))
    duplicates = pre - len(obs)
    obs["log2fc"] = np.log2(obs["fold"].astype(float))
    obs["mutations"] = obs["signature"].map(lambda s: tuple(s.split(",")) if s else tuple())
    obs["num_mutations"] = obs["mutations"].map(len)
    # Connected components join any reused SeqID (even across genes) and any identical genotype
    # within a gene. This is stricter than random sample splitting.
    samples = obs[["sample_id", "gene", "seqid", "signature"]].drop_duplicates("sample_id").reset_index(drop=True)
    uf = DSU(len(samples))
    seen_sequence, seen_genotype = {}, {}
    for i, r in samples.iterrows():
        k1 = r["seqid"]
        k2 = (r["gene"], r["signature"])
        if k1 in seen_sequence:
            uf.union(i, seen_sequence[k1])
        else:
            seen_sequence[k1] = i
        if k2 in seen_genotype:
            uf.union(i, seen_genotype[k2])
        else:
            seen_genotype[k2] = i
    samples["split_group"] = [uf.find(i) for i in range(len(samples))]
    obs = obs.merge(samples[["sample_id", "split_group"]], on="sample_id", validate="many_to_one")
    obs = obs.sort_values(["gene", "drug", "sample_id"]).reset_index(drop=True)
    obs["observation_id"] = np.arange(len(obs), dtype=int)
    provenance = {"datasets": dataset_info, "valid_pairs": int(len(obs)),
                  "unique_gene_genotypes": int(len(samples)),
                  "connected_isolate_groups": int(samples["split_group"].nunique()),
                  "exact_duplicate_pair_measurements_collapsed": int(duplicates),
                  "outcome": "log2 of observed phenotypic fold-change (in vitro)",
                  "not_clinical_outcome": True,
                  "patient_identifiers_available": False}
    return obs, provenance


def split_grouped(obs, seed, fractions=(.2, .1, .1)):
    """Perform 60/10/10/20 train/validation/calibration/test group-disjoint split."""
    test, val, cal = fractions
    groups = obs["split_group"].values
    dummy = np.zeros(len(obs))
    splitter = GroupShuffleSplit(n_splits=1, test_size=test, random_state=seed)
    initial, test_idx = next(splitter.split(dummy, groups=groups))
    remainder_f = (val + cal) / (1 - test)
    internal = GroupShuffleSplit(n_splits=1, test_size=remainder_f, random_state=seed + 100003)
    train_local, middle_local = next(internal.split(dummy[initial], groups=groups[initial]))
    train_idx, middle_idx = initial[train_local], initial[middle_local]
    cal_fraction = cal / (val + cal)
    final = GroupShuffleSplit(n_splits=1, test_size=cal_fraction, random_state=seed + 200003)
    val_local, cal_local = next(final.split(dummy[middle_idx], groups=groups[middle_idx]))
    indexes = {"train": train_idx, "validation": middle_idx[val_local],
               "calibration": middle_idx[cal_local], "test": test_idx}
    # Ensure identical isolate (or genotype pattern) cannot appear in two splits.
    split_groups = {name: set(groups[ix]) for name, ix in indexes.items()}
    for a in split_groups:
        for b in split_groups:
            if a != b:
                assert not split_groups[a].intersection(split_groups[b]), "Group leakage detected"
    assignment = np.full(len(obs), "", dtype=object)
    for name, idx in indexes.items():
        assignment[idx] = name
    assert (assignment != "").all()
    return assignment, {k: {"pairs": len(v), "groups": len(split_groups[k])} for k, v in indexes.items()}


def save_prepared_data(root, obs, provenance):
    processed = Path(root) / "data" / "processed"
    processed.mkdir(parents=True, exist_ok=True)
    export = obs.drop(columns=["mutations"]).copy()
    export.to_csv(processed / "phenotypes_clean.csv", index=False)
    (processed / "provenance.json").write_text(json.dumps(provenance, indent=2), encoding="utf-8")
