"""Independently verify public, aggregate archived-run results (stdlib only)."""
from pathlib import Path
import csv, json, math

ROOT = Path(__file__).resolve().parents[1]
TAB = ROOT / 'results' / 'tables'
summary = json.loads((ROOT / 'results' / 'summary.json').read_text(encoding='utf-8'))
with (TAB / 'Table_01_Evaluation_By_Drug_And_Seed.csv').open(newline='', encoding='utf-8') as fp:
    rows = list(csv.DictReader(fp))
models = {name: [] for name in (
    'median_baseline', 'linear_mutations', 'position_only_ablation',
    'mutation_forest', 'shiftguard_ensemble')}
for row in rows:
    if row['model'] in models:
        models[row['model']].append(float(row['rmse']))

assert len(models['shiftguard_ensemble']) == summary['completed_experiments'] == 67
assert len({(r['gene'],r['drug']) for r in rows if r['model']=='shiftguard_ensemble'}) == summary['distinct_drugs'] == 23
expected = {
    'median_baseline': 2.3341239923243005,
    'linear_mutations': 1.1133184956460418,
    'position_only_ablation': 1.2666077309877437,
    'mutation_forest': 1.1110467204952688,
    'shiftguard_ensemble': 1.0174300110644539,
}
for name, values in models.items():
    mean = sum(values) / len(values)
    assert math.isclose(mean, expected[name], rel_tol=1e-10), (name, mean)
    print(f'{name:24s}: {mean:.3f} RMSE (n={len(values)})')
assert math.isclose(summary['pooled_standard_coverage'], 0.9012028374627326)
assert math.isclose(summary['pooled_novelty_coverage'], 0.9000205613241493)
assert summary['software_and_data_fingerprint'] == '42374487aabb5e1a'
print('PASS: Archived metrics agree with the manuscript source run.')
