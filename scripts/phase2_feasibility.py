"""Phase-2 feasibility: predicting the EXPERIMENTAL optical gap of ChalMolDB polymer records.

Questions
  Q1  How well can the experimental optical gap be predicted, and from what?
      inputs: (a) nothing (mean), (b) oligomer DFT gap only (linear calibration),
              (c) structure descriptors only, (d) structure + DFT values,
              (e) delta-learning: linear DFT calibration + model of the residual.
  Q2  Does a model predict the S->Se change of the experimental gap better than raw DFT?

Validation
  5-fold GroupKFold, groups = connected clusters of records whose oligomer structures have
  Morgan-fingerprint Tanimoto >= 0.7 (so near-identical backbones never span train/test);
  the InChIKey grouping is reported as the looser reference.
"""
import argparse, json, sys
from pathlib import Path
import numpy as np, pandas as pd
from rdkit import Chem, DataStructs, RDLogger
from rdkit.Chem import AllChem, Descriptors, rdMolDescriptors
from sklearn.ensemble import RandomForestRegressor
from sklearn.linear_model import LinearRegression, RidgeCV
from sklearn.metrics import mean_absolute_error, r2_score
from sklearn.model_selection import GroupKFold
from sklearn.pipeline import make_pipeline
from sklearn.preprocessing import StandardScaler
import xgboost as xgb

RDLogger.DisableLog("rdApp.*")
ap = argparse.ArgumentParser(); ap.add_argument("--root", required=True); ap.add_argument("--out", required=True)
a = ap.parse_args()
ROOT, OUT = Path(a.root), Path(a.out); OUT.mkdir(parents=True, exist_ok=True)
sys.path.insert(0, str(ROOT))
from descriptor_engine import calculate_chalcogen_aware_descriptors  # ChalMolDB S/Se/Te layer

df = pd.read_csv(ROOT / "data/chalcogen_database_v1_master.csv.gz", low_memory=False)
pol = df[df.Record_ID.str.startswith("POL_")].copy()
for c in ["Experimental_Eg_eV", "Eg_eV", "HOMO_eV", "LUMO_eV", "TDDFT_S1_eV", "Experimental_HOMO_eV", "Experimental_LUMO_eV"]:
    pol[c] = pd.to_numeric(pol[c], errors="coerce")
pol = pol[pol.Experimental_Eg_eV.notna() & pol.Canonical_SMILES.notna() & pol.Eg_eV.notna()].reset_index(drop=True)
log = {"n_records": int(len(pol))}

# ---------------- features
GEN = {"MolWt": Descriptors.MolWt, "TPSA": Descriptors.TPSA, "LogP": Descriptors.MolLogP,
       "RotB": Descriptors.NumRotatableBonds, "Rings": rdMolDescriptors.CalcNumRings,
       "AroRings": rdMolDescriptors.CalcNumAromaticRings, "Hetero": rdMolDescriptors.CalcNumHeteroatoms,
       "FracSP3": rdMolDescriptors.CalcFractionCSP3, "BertzCT": Descriptors.BertzCT,
       "HBA": rdMolDescriptors.CalcNumHBA, "HeavyAtoms": lambda m: m.GetNumHeavyAtoms()}
rows, fps, mols = [], [], []
for smi in pol.Canonical_SMILES:
    m = Chem.MolFromSmiles(smi); mols.append(m)
    r = {k: f(m) for k, f in GEN.items()}
    for el, z in [("N", 7), ("O", 8), ("F", 9), ("S", 16), ("Se", 34)]:
        r[f"n_{el}"] = sum(1 for at in m.GetAtoms() if at.GetAtomicNum() == z)
    try:
        r.update({f"ch_{k}": v for k, v in calculate_chalcogen_aware_descriptors(m).items() if isinstance(v, (int, float))})
    except Exception as e:
        log.setdefault("chalcogen_errors", 0); log["chalcogen_errors"] += 1
    rows.append(r)
    fps.append(AllChem.GetMorganFingerprintAsBitVect(m, 2, 1024))
X_desc = pd.DataFrame(rows).fillna(0.0)
X_fp = np.array([list(fp) for fp in fps], dtype=np.uint8)
dft_cols = ["Eg_eV", "HOMO_eV", "LUMO_eV"]
X_dft = pol[dft_cols].to_numpy()
y = pol.Experimental_Eg_eV.to_numpy()
log["n_descriptors"] = int(X_desc.shape[1])

# ---------------- groups: Tanimoto >= 0.7 clusters (union-find), and InChIKey
n = len(fps); parent = list(range(n))
def find(i):
    while parent[i] != i:
        parent[i] = parent[parent[i]]; i = parent[i]
    return i
for i in range(n):
    sims = DataStructs.BulkTanimotoSimilarity(fps[i], fps[i + 1:])
    for j, s in enumerate(sims):
        if s >= 0.7:
            ri, rj = find(i), find(i + 1 + j)
            if ri != rj: parent[ri] = rj
g_sim = np.array([find(i) for i in range(n)])
g_ik = pd.factorize(pol.InChIKey)[0]
log["groups_tanimoto07"] = int(len(set(g_sim))); log["largest_group"] = int(pd.Series(g_sim).value_counts().max())
log["groups_inchikey"] = int(len(set(g_ik)))

def xgbm(): return xgb.XGBRegressor(n_estimators=600, max_depth=5, learning_rate=0.03, subsample=0.8, colsample_bytree=0.6, random_state=0, n_jobs=-1)
def rfm(): return RandomForestRegressor(n_estimators=500, min_samples_leaf=2, random_state=0, n_jobs=-1)

def cv(feature_fn, model_fn, groups, delta=False):
    pred = np.zeros_like(y)
    for tr, te in GroupKFold(5).split(y, y, groups):
        if delta:
            cal = LinearRegression().fit(X_dft[tr, :1], y[tr])
            base_tr, base_te = cal.predict(X_dft[tr, :1]), cal.predict(X_dft[te, :1])
            m = model_fn().fit(feature_fn(tr), y[tr] - base_tr)
            pred[te] = base_te + m.predict(feature_fn(te))
        else:
            m = model_fn().fit(feature_fn(tr), y[tr]); pred[te] = m.predict(feature_fn(te))
    return pred

Xd = X_desc.to_numpy(); Xfull = np.hstack([Xd, X_fp]); Xfull_dft = np.hstack([Xfull, X_dft])
F = {
    "mean": None,
    "dft_gap_linear": lambda idx: X_dft[idx, :1],
    "dft_hlg_linear": lambda idx: X_dft[idx],
    "struct_xgb": lambda idx: Xfull[idx],
    "struct+dft_xgb": lambda idx: Xfull_dft[idx],
    "struct+dft_rf": lambda idx: Xfull_dft[idx],
    "delta_struct+dft_xgb": lambda idx: Xfull_dft[idx],
}
results, preds = {}, {}
for gname, groups in [("tanimoto07", g_sim), ("inchikey", g_ik)]:
    for name, fn in F.items():
        if name == "mean":
            p = np.zeros_like(y)
            for tr, te in GroupKFold(5).split(y, y, groups): p[te] = y[tr].mean()
        elif name.endswith("_linear"):
            p = cv(fn, LinearRegression, groups)
        elif name.startswith("delta_"):
            p = cv(fn, xgbm, groups, delta=True)
        elif name.endswith("_rf"):
            p = cv(fn, rfm, groups)
        else:
            p = cv(fn, xgbm, groups)
        results[f"{gname}|{name}"] = {"MAE": round(mean_absolute_error(y, p), 4), "R2": round(r2_score(y, p), 3)}
        preds[f"{gname}|{name}"] = p
        print(gname, name, results[f"{gname}|{name}"], flush=True)

# ---------------- Q2: S->Se pairs, out-of-fold predictions under the strict grouping
an = json.load(open(ROOT / "data/build/polymer_chalcogen_analysis.json"))["se_s_swap_pairs"]["pairs"]
idx = {r: i for i, r in enumerate(pol.Record_ID)}
q2 = {}
for key in ["tanimoto07|struct+dft_xgb", "tanimoto07|delta_struct+dft_xgb", "tanimoto07|dft_gap_linear"]:
    p = preds[key]; d_true, d_pred = [], []
    for pr in an:
        if pr["Se_record"] in idx and pr["S_record"] in idx:
            i, j = idx[pr["Se_record"]], idx[pr["S_record"]]
            d_true.append(y[i] - y[j]); d_pred.append(p[i] - p[j])
    d_true, d_pred = np.array(d_true), np.array(d_pred)
    q2[key] = {"pairs": int(len(d_true)), "MAE_delta": round(float(np.abs(d_true - d_pred).mean()), 4),
               "sign_agreement": round(float(np.mean(np.sign(d_true) == np.sign(d_pred))), 3),
               "pearson_delta": round(float(np.corrcoef(d_true, d_pred)[0, 1]), 3)}
    print("Q2", key, q2[key], flush=True)
# raw DFT difference as reference
d_true, d_raw = [], []
for pr in an:
    if pr["Se_record"] in idx and pr["S_record"] in idx:
        i, j = idx[pr["Se_record"]], idx[pr["S_record"]]
        d_true.append(y[i] - y[j]); d_raw.append(X_dft[i, 0] - X_dft[j, 0])
d_true, d_raw = np.array(d_true), np.array(d_raw)
q2["raw_dft_difference"] = {"MAE_delta": round(float(np.abs(d_true - d_raw).mean()), 4), "sign_agreement": round(float(np.mean(np.sign(d_true) == np.sign(d_raw))), 3),
                            "pearson_delta": round(float(np.corrcoef(d_true, d_raw)[0, 1]), 3)}
print("Q2 raw", q2["raw_dft_difference"])
# per-chalcogen error for the best model
best = "tanimoto07|delta_struct+dft_xgb"; p = preds[best]
has_se = pol.Canonical_SMILES.str.contains(r"\[se\]|\[Se\]")
log["per_chalcogen_MAE"] = {"S-only": round(mean_absolute_error(y[~has_se], p[~has_se]), 4), "Se": round(mean_absolute_error(y[has_se], p[has_se]), 4), "n_Se": int(has_se.sum())}
json.dump({"log": log, "Q1": results, "Q2": q2}, open(OUT / "phase2_feasibility.json", "w"), indent=1)
print(json.dumps(log, indent=1))
