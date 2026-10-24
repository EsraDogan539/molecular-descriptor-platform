"""Phase-2: pair models for the S->Se change of the experimental optical gap.

Target: d_exp = Eg_exp(Se polymer) - Eg_exp(S partner) for the 139 ChalMolDB swap pairs.
Baselines (no training):  zero; constant mean shift (fitted in CV); raw DFT difference;
                          calibrated DFT difference (slope fitted in CV).
Pair models (all predict d_exp directly from pair features):
  P1  ridge  on the 139 S->Se pairs only
  P2  xgb    on the 139 S->Se pairs only
  P3  xgb    pre-trained on many matched polymer pairs (any two polymers with Tanimoto >= 0.5),
             then the S->Se pairs added; i.e. it learns "how does the measured gap change between
             two similar polymers" from thousands of pairs, not only 139.
  P4  delta: calibrated DFT difference + P3-style model of the residual.
Pair features: x(B) - x(A) for descriptors and DFT values, x(A), and the position of the
introduced Se (selenadiazole N[Se]N acceptor ring vs selenophene C[Se]C ring).
Validation: 5-fold GroupKFold over connected components of the pair graph (a record never
appears in both a training and a test pair); general pairs touching any test record are dropped.
"""
import argparse, json, sys, itertools
from pathlib import Path
import numpy as np, pandas as pd
from rdkit import Chem, DataStructs, RDLogger
from rdkit.Chem import AllChem, Descriptors, rdMolDescriptors
from sklearn.linear_model import RidgeCV, LinearRegression
from sklearn.model_selection import GroupKFold
from sklearn.pipeline import make_pipeline
from sklearn.preprocessing import StandardScaler
import xgboost as xgb

RDLogger.DisableLog("rdApp.*")
ap = argparse.ArgumentParser(); ap.add_argument("--root", required=True); ap.add_argument("--out", required=True); a = ap.parse_args()
ROOT, OUT = Path(a.root), Path(a.out); OUT.mkdir(parents=True, exist_ok=True)

df = pd.read_csv(ROOT / "data/chalcogen_database_v1_master.csv.gz", low_memory=False)
pol = df[df.Record_ID.str.startswith("POL_")].copy()
for c in ["Experimental_Eg_eV", "Eg_eV", "HOMO_eV", "LUMO_eV"]:
    pol[c] = pd.to_numeric(pol[c], errors="coerce")
pol = pol[pol.Experimental_Eg_eV.notna() & pol.Canonical_SMILES.notna() & pol.Eg_eV.notna()].reset_index(drop=True)
rid = {r: i for i, r in enumerate(pol.Record_ID)}
y = pol.Experimental_Eg_eV.to_numpy()

PAT_SEDIAZ = Chem.MolFromSmarts("[#7][Se,se][#7]")
PAT_SEPHEN = Chem.MolFromSmarts("[#6][Se,se][#6]")
PAT_SDIAZ = Chem.MolFromSmarts("[#7][S,s][#7]")
PAT_SPHEN = Chem.MolFromSmarts("[#6;a][s][#6;a]")
GEN = {"MolWt": Descriptors.MolWt, "TPSA": Descriptors.TPSA, "LogP": Descriptors.MolLogP, "RotB": Descriptors.NumRotatableBonds,
       "Rings": rdMolDescriptors.CalcNumRings, "AroRings": rdMolDescriptors.CalcNumAromaticRings, "Hetero": rdMolDescriptors.CalcNumHeteroatoms,
       "FracSP3": rdMolDescriptors.CalcFractionCSP3, "BertzCT": Descriptors.BertzCT, "HBA": rdMolDescriptors.CalcNumHBA,
       "HeavyAtoms": lambda m: m.GetNumHeavyAtoms()}
rows, fps = [], []
for smi in pol.Canonical_SMILES:
    m = Chem.MolFromSmiles(smi)
    r = {k: f(m) for k, f in GEN.items()}
    for el, z in [("N", 7), ("O", 8), ("F", 9), ("S", 16), ("Se", 34)]:
        r[f"n_{el}"] = sum(1 for at in m.GetAtoms() if at.GetAtomicNum() == z)
    r["n_Se_diazole"] = len(m.GetSubstructMatches(PAT_SEDIAZ)); r["n_Se_phene"] = len(m.GetSubstructMatches(PAT_SEPHEN))
    r["n_S_diazole"] = len(m.GetSubstructMatches(PAT_SDIAZ)); r["n_S_phene"] = len(m.GetSubstructMatches(PAT_SPHEN))
    rows.append(r); fps.append(AllChem.GetMorganFingerprintAsBitVect(m, 2, 1024))
D = pd.DataFrame(rows).fillna(0.0)
D["dft_gap"] = pol.Eg_eV.values; D["dft_homo"] = pol.HOMO_eV.values; D["dft_lumo"] = pol.LUMO_eV.values
X = D.to_numpy(float); cols = list(D.columns)
i_gap = cols.index("dft_gap")

def pair_feat(i, j):  # j = new (Se) polymer, i = reference (S) polymer
    return np.concatenate([X[j] - X[i], X[i], [DataStructs.TanimotoSimilarity(fps[i], fps[j])]])

an = json.load(open(ROOT / "data/build/polymer_chalcogen_analysis.json"))["se_s_swap_pairs"]["pairs"]
SP = [(rid[p["S_record"]], rid[p["Se_record"]]) for p in an if p["S_record"] in rid and p["Se_record"] in rid]
SP = np.array(SP); n_sp = len(SP)
Xs = np.array([pair_feat(i, j) for i, j in SP]); ds = np.array([y[j] - y[i] for i, j in SP]); dd = np.array([X[j, i_gap] - X[i, i_gap] for i, j in SP])

# connected components of the S->Se pair graph -> groups
par = {}
def f(x):
    par.setdefault(x, x)
    while par[x] != x: par[x] = par[par[x]]; x = par[x]
    return x
for i, j in SP: par[f(i)] = f(j)
grp = np.array([f(i) for i, _ in SP])

# general matched pairs (Tanimoto >= 0.5, both directions), built once
gen_i, gen_j = [], []
for i in range(len(fps)):
    s = DataStructs.BulkTanimotoSimilarity(fps[i], fps[i + 1:])
    for k, v in enumerate(s):
        if v >= 0.5:
            gen_i.append(i); gen_j.append(i + 1 + k)
gen_i, gen_j = np.array(gen_i), np.array(gen_j)
print("general pairs", len(gen_i), flush=True)
Xg = np.vstack([np.array([pair_feat(i, j) for i, j in zip(gen_i, gen_j)]), np.array([pair_feat(j, i) for i, j in zip(gen_i, gen_j)])])
dg = np.concatenate([y[gen_j] - y[gen_i], y[gen_i] - y[gen_j]])
gi_all = np.concatenate([gen_i, gen_j]); gj_all = np.concatenate([gen_j, gen_i])

def xgbm(n=500): return xgb.XGBRegressor(n_estimators=n, max_depth=4, learning_rate=0.03, subsample=0.8, colsample_bytree=0.6, min_child_weight=3, random_state=0, n_jobs=-1)

P = {k: np.zeros(n_sp) for k in ["zero", "mean_shift", "raw_dft", "cal_dft", "P1_ridge_pairs", "P2_xgb_pairs", "P3_xgb_general+pairs", "P4_delta_general+pairs"]}
for tr, te in GroupKFold(5).split(Xs, ds, grp):
    test_records = set(SP[te].ravel())
    P["zero"][te] = 0.0
    P["mean_shift"][te] = ds[tr].mean()
    P["raw_dft"][te] = dd[te]
    cal = LinearRegression().fit(dd[tr, None], ds[tr]); P["cal_dft"][te] = cal.predict(dd[te, None])
    r = make_pipeline(StandardScaler(), RidgeCV(alphas=np.logspace(-2, 3, 20))).fit(Xs[tr], ds[tr]); P["P1_ridge_pairs"][te] = r.predict(Xs[te])
    m = xgbm(300).fit(Xs[tr], ds[tr]); P["P2_xgb_pairs"][te] = m.predict(Xs[te])
    keep = np.array([(a not in test_records) and (b not in test_records) for a, b in zip(gi_all, gj_all)])
    Xtr = np.vstack([Xg[keep], np.repeat(Xs[tr], 5, axis=0)]); ytr = np.concatenate([dg[keep], np.repeat(ds[tr], 5)])
    m3 = xgbm(600).fit(Xtr, ytr); P["P3_xgb_general+pairs"][te] = m3.predict(Xs[te])
    # delta: calibrated DFT difference + residual model (general pairs calibrated with their own DFT diffs)
    ddg = Xg[:, i_gap]  # first block = X(B)-X(A), dft_gap difference
    calg = LinearRegression().fit(np.concatenate([ddg[keep], np.repeat(dd[tr], 5)])[:, None], ytr)
    res = ytr - calg.predict(np.concatenate([ddg[keep], np.repeat(dd[tr], 5)])[:, None])
    m4 = xgbm(600).fit(Xtr, res); P["P4_delta_general+pairs"][te] = calg.predict(dd[te, None]) + m4.predict(Xs[te])

res = {}
for k, p in P.items():
    res[k] = {"MAE": round(float(np.abs(ds - p).mean()), 4),
              "sign_agreement": round(float(np.mean(np.sign(ds) == np.sign(p))), 3) if np.any(p != 0) else None,
              "pearson": round(float(np.corrcoef(ds, p)[0, 1]), 3) if np.std(p) > 0 else None}
    print(k, res[k], flush=True)
# where is the Se? per-position error of best model
pos = np.where(Xs[:, cols.index("n_Se_diazole")] > 0, "acceptor_selenadiazole", np.where(Xs[:, cols.index("n_Se_phene")] > 0, "donor_selenophene", "other"))
best = min([k for k in P if k.startswith("P")], key=lambda k: res[k]["MAE"])
res["by_position"] = {p: {"n": int((pos == p).sum()), "MAE_best": round(float(np.abs(ds - P[best])[pos == p].mean()), 4),
                          "MAE_raw_dft": round(float(np.abs(ds - dd)[pos == p].mean()), 4), "mean_d_exp": round(float(ds[pos == p].mean()), 4)}
                      for p in set(pos)}
res["best"] = best; res["n_pairs"] = int(n_sp); res["n_general_pairs"] = int(len(gen_i)); res["n_groups"] = int(len(set(grp)))
print(json.dumps(res["by_position"], indent=1))
json.dump(res, open(OUT / "phase2_pairs.json", "w"), indent=1)
