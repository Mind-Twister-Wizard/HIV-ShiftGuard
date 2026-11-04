"""Mutation-level benchmark, ablations, uncertainty intervals and novelty gating."""
import math
from pathlib import Path

import joblib
import numpy as np
import pandas as pd
from scipy import sparse
from scipy.stats import pearsonr, spearmanr
from sklearn.ensemble import ExtraTreesRegressor
from sklearn.feature_extraction import DictVectorizer
from sklearn.linear_model import Ridge
from sklearn.metrics import mean_absolute_error, mean_squared_error, r2_score


def features(mutations, mode="full"):
    record = {"BIAS": 1.0}
    if mode == "count":
        record["mutation_count"] = float(len(mutations))
    else:
        for item in mutations:
            token = item if mode == "full" else item[:-1]  # position-only ablation
            record[f"mut:{token}"] = 1.0
    return record


def novelty_jaccard(x_eval, x_train, names, query_mutations):
    """Exact minimum mutation-set Jaccard distance to training genotype patterns.

    Crucially, mutation variants *unseen* during training remain counted in the
    query's full set cardinality, even though DictVectorizer ignores their columns.
    """
    names = np.asarray(names)
    keep = np.where(np.char.startswith(names.astype(str), "mut:"))[0]
    true_query_size = np.asarray([len(set(tokens)) for tokens in query_mutations], dtype=float)
    if len(keep) == 0:
        return np.where(true_query_size == 0, 0.0, 1.0)
    train = x_train[:, keep].tocsr().astype(np.float32)
    query = x_eval[:, keep].tocsr().astype(np.float32)
    train_size = np.asarray(train.sum(axis=1)).ravel()
    out = np.empty(query.shape[0], dtype=float)
    for start in range(0, query.shape[0], 160):
        chunk = query[start:start + 160]
        inter = (chunk @ train.T).toarray()
        union = true_query_size[start:start + chunk.shape[0], None] + train_size[None, :] - inter
        similarities = np.divide(inter, union, out=np.zeros_like(inter), where=union > 0)
        similarities[union == 0] = 1.0  # both mutation sets empty: identical
        out[start:start + chunk.shape[0]] = 1.0 - similarities.max(axis=1)
    return np.clip(out, 0.0, 1.0)


def finite_or_none(x):
    return float(x) if np.isfinite(x) else None


def metrics(y, pred):
    y = np.asarray(y, dtype=float)
    pred = np.asarray(pred, dtype=float)
    if len(y) == 0:
        return {"n": 0}
    mae = float(mean_absolute_error(y, pred))
    rmse = float(np.sqrt(mean_squared_error(y, pred)))
    spearman = float(spearmanr(y, pred).statistic) if len(y) >= 3 and np.std(y) > 1e-10 and np.std(pred) > 1e-10 else float("nan")
    pearson = float(pearsonr(y, pred).statistic) if len(y) >= 3 and np.std(y) > 1e-10 and np.std(pred) > 1e-10 else float("nan")
    return {"n": int(len(y)), "mae": mae, "rmse": rmse,
            "r2": finite_or_none(r2_score(y, pred)) if len(y) > 1 else None,
            "pearson": finite_or_none(pearson), "spearman": finite_or_none(spearman)}


def conformal_q(residual_scores, alpha):
    vals = np.sort(np.asarray(residual_scores, dtype=float))
    if len(vals) == 0:
        raise ValueError("Cannot calibrate on an empty set")
    k = max(1, min(len(vals), math.ceil((len(vals) + 1) * (1.0 - alpha))))
    return float(vals[k - 1])


def train_drug(drug_obs, partition_labels, config, seed, root_checkpoint):
    """Fit all models for one antiretroviral and produce never-before-seen test scores."""
    pieces = {name: drug_obs.loc[partition_labels == name].reset_index(drop=True)
              for name in ("train", "validation", "calibration", "test")}
    n = {key: len(p) for key, p in pieces.items()}
    if n["train"] < config["minimum_train_samples"] or n["calibration"] < config["minimum_calibration_samples"] or n["test"] < config["minimum_test_samples"] or n["validation"] < 10:
        return None, None, f"Insufficient split sizes: {n}"
    train = pieces["train"]
    y = {k: p["log2fc"].to_numpy(dtype=float) for k, p in pieces.items()}
    tokens = {k: [tuple(m) for m in p["mutations"]] for k, p in pieces.items()}
    vec = DictVectorizer(sparse=True)
    X = {"train": vec.fit_transform([features(m) for m in tokens["train"]])}
    for k in ("validation", "calibration", "test"):
        X[k] = vec.transform([features(m) for m in tokens[k]])
    ridge = Ridge(alpha=12.0, solver="lsqr")
    ridge.fit(X["train"], y["train"])
    forest = ExtraTreesRegressor(n_estimators=config["forest_trees"], max_depth=config["forest_max_depth"],
                                 min_samples_leaf=config["min_samples_leaf"], max_features=0.75,
                                 random_state=seed, n_jobs=config["num_threads"])
    forest.fit(X["train"], y["train"])
    y_val = y["validation"]
    predicted = {name: {"ridge": ridge.predict(x), "forest": forest.predict(x)} for name, x in X.items() if name != "train"}
    candidates = np.linspace(0.0, 1.0, 11)
    scores = [np.mean((weight * predicted["validation"]["ridge"] + (1 - weight) * predicted["validation"]["forest"] - y_val) ** 2)
              for weight in candidates]
    alpha_ridge = float(candidates[int(np.argmin(scores))])
    fused = {k: alpha_ridge * v["ridge"] + (1 - alpha_ridge) * v["forest"] for k, v in predicted.items()}
    position_vec = DictVectorizer(sparse=True)
    position_X = position_vec.fit_transform([features(m, "position") for m in tokens["train"]])
    position_test = position_vec.transform([features(m, "position") for m in tokens["test"]])
    pos_forest = ExtraTreesRegressor(n_estimators=max(24, config["forest_trees"] // 2),
                                      max_depth=config["forest_max_depth"], min_samples_leaf=config["min_samples_leaf"],
                                      random_state=seed + 907, n_jobs=config["num_threads"])
    pos_forest.fit(position_X, y["train"])
    median = float(np.median(y["train"]))
    predictions = {"median_baseline": np.full(n["test"], median),
                   "linear_mutations": predicted["test"]["ridge"],
                   "position_only_ablation": pos_forest.predict(position_test),
                   "mutation_forest": predicted["test"]["forest"],
                   "shiftguard_ensemble": fused["test"]}
    nov = {}
    for k in ("validation", "calibration", "test"):
        nov[k] = novelty_jaccard(X[k], X["train"], vec.get_feature_names_out(), tokens[k])
    conf_alpha = 1.0 - config["coverage_target"]
    unscaled_q = conformal_q(np.abs(y["calibration"] - fused["calibration"]), conf_alpha)
    gamma = float(config["novelty_scale_gamma"])
    calibration_scale = 1 + gamma * nov["calibration"]
    scaled_q = conformal_q(np.abs(y["calibration"] - fused["calibration"]) / calibration_scale, conf_alpha)
    lower_u, upper_u = fused["test"] - unscaled_q, fused["test"] + unscaled_q
    half_scaled = scaled_q * (1 + gamma * nov["test"])
    lower_s, upper_s = fused["test"] - half_scaled, fused["test"] + half_scaled
    cover_unscaled = (y["test"] >= lower_u) & (y["test"] <= upper_u)
    cover_scaled = (y["test"] >= lower_s) & (y["test"] <= upper_s)
    gate_novelty_threshold = float(np.quantile(nov["validation"], .90))
    gate_width_threshold = float(np.quantile(scaled_q * (1 + gamma * nov["calibration"]), .75))
    accept = (nov["test"] <= gate_novelty_threshold) & (half_scaled <= gate_width_threshold)
    result = {"gene": str(train["gene"].iloc[0]), "drug": str(train["drug"].iloc[0]),
              "seed": int(seed), "split_counts": n, "ridge_weight_selected_on_validation": alpha_ridge,
              "models": {name: metrics(y["test"], estimate) for name, estimate in predictions.items()},
              "intervals": {"target_coverage": float(config["coverage_target"]),
                            "unscaled_coverage": float(cover_unscaled.mean()), "unscaled_mean_width": float(2 * unscaled_q),
                            "novelty_scaled_coverage": float(cover_scaled.mean()), "novelty_scaled_mean_width": float(np.mean(2 * half_scaled)),
                            "n_calibration": n["calibration"]},
              "novelty": {"mean_test": float(np.mean(nov["test"])), "median_test": float(np.median(nov["test"])),
                          "gate_novelty_threshold": gate_novelty_threshold, "gate_halfwidth_threshold": gate_width_threshold,
                          "gate_retention": float(np.mean(accept)),
                          "gate_accepted_rmse": float(np.sqrt(np.mean((y["test"][accept] - fused["test"][accept]) ** 2))) if accept.any() else None,
                          "gate_rejected_count": int((~accept).sum())},
              "mutation_feature_count": int(X["train"].shape[1] - 1),
              "train_median_log2_fold": median}
    obs_ids = pieces["test"]["observation_id"].astype(int).to_numpy()
    pred_df = pd.DataFrame({"observation_id": obs_ids, "gene": result["gene"], "drug": result["drug"],
                            "seed": seed, "split_group": pieces["test"]["split_group"].to_numpy(),
                            "observed_log2fc": y["test"], "predicted_log2fc": fused["test"],
                            "drug_median_pred": predictions["median_baseline"],
                            "ridge_pred": predictions["linear_mutations"],
                            "position_only_pred": predictions["position_only_ablation"],
                            "forest_pred": predictions["mutation_forest"],
                            "novelty_jaccard": nov["test"], "interval_lower_unscaled": lower_u,
                            "interval_upper_unscaled": upper_u, "interval_lower_novelty": lower_s,
                            "interval_upper_novelty": upper_s, "gate_accept_research": accept})
    # Group-cluster bootstrap: preserves within-genotype dependence and avoids
    # treating all repeated measurements as independent for interval estimation.
    group_values = pred_df["split_group"].to_numpy()
    unique_groups = np.unique(group_values)
    if len(unique_groups) >= 3:
        rng = np.random.default_rng(seed * 991 + sum(ord(c) for c in result["drug"]))
        group_indices = [np.flatnonzero(group_values == group) for group in unique_groups]
        bootstrap_rmse, bootstrap_delta = [], []
        for _ in range(300):
            drawn = rng.integers(0, len(unique_groups), len(unique_groups))
            chosen = np.concatenate([group_indices[g] for g in drawn])
            full_error = y["test"][chosen] - fused["test"][chosen]
            naive_error = y["test"][chosen] - predictions["median_baseline"][chosen]
            rms_full = np.sqrt(np.mean(full_error ** 2))
            rms_naive = np.sqrt(np.mean(naive_error ** 2))
            bootstrap_rmse.append(rms_full)
            bootstrap_delta.append(rms_naive - rms_full)
        result["group_bootstrap_95ci"] = {
            "method": "300 group-cluster bootstrap resamples of untouched test set",
            "ensemble_rmse_lower": float(np.quantile(bootstrap_rmse, .025)),
            "ensemble_rmse_upper": float(np.quantile(bootstrap_rmse, .975)),
            "improvement_over_median_lower": float(np.quantile(bootstrap_delta, .025)),
            "improvement_over_median_upper": float(np.quantile(bootstrap_delta, .975)),
            "clusters": int(len(unique_groups))}
    else:
        result["group_bootstrap_95ci"] = None
    ckpt = {"gene": result["gene"], "drug": result["drug"], "vectorizer": vec, "ridge": ridge,
            "forest": forest, "ridge_weight": alpha_ridge, "q_unscaled": unscaled_q,
            "q_scaled": scaled_q, "gamma": gamma, "train_mutation_tokens": tokens["train"],
            "gate_novelty_threshold": gate_novelty_threshold, "gate_halfwidth_threshold": gate_width_threshold,
            "research_only": True, "outcome": "in-vitro log2 fold susceptibility"}
    root_checkpoint = Path(root_checkpoint)
    root_checkpoint.parent.mkdir(parents=True, exist_ok=True)
    temp = root_checkpoint.with_name(root_checkpoint.name + ".tmp")
    joblib.dump(ckpt, temp, compress=2)
    temp.replace(root_checkpoint)
    return result, pred_df, None
