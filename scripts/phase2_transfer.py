"""Phase-2: does the large computed development collection help predict EXPERIMENTAL polymer gaps?

Source task : DEV collection (3,088 D-A molecules, DFT HOMO-LUMO gap, B3LYP/6-311+G(d)).
Target task : POL collection (1,323 polymers, experimental optical gap).
Same feature space for both (RDKit descriptors + chalcogen-position counts + Morgan r2/1024).

Settings, each with and without the polymer's own oligomer DFT values:
  base        XGB trained on POL only
  feat        XGB on POL with an extra feature = prediction of a DEV-trained model ("DEV-style gap")
  finetune    XGB first trained on DEV gaps, then boosting continued on POL experimental gaps
  stack       same as 'feat' but DEV model trained on DEV + EXT (more chalcogen diversity)
Validation: 5-fold GroupKFold over Butina clusters (Tanimoto >= 0.6) of the polymer oligomers,
repeated with 3 seeds of cluster-to-fold assignment; mean and SD of MAE reported.
Also a learning curve: MAE of base vs feat when only 10/25/50/100 % of the POL training folds are used.
"""
import argparse, json
from pathlib import Path
import numpy as np, pandas as pd
from rdkit import Chem, DataStructs, RDLogger
from rdkit.Chem import AllChem, Descriptors, rdMolDescriptors
from rdkit.ML.Cluster import Butina
from sklearn.linear_model import LinearRegression
from sklearn.metrics import mean_absolute_error, r2_score
import xgboost as xgb

RDLogger.DisableLog("rdApp.*")
ap = argparse.ArgumentParser(); ap.add_argument("--root", required=True); ap.add_argument("--out", required=True); a = ap.parse_args()
ROOT, OUT = Path(a.root), Path(a.out); OUT.mkdir(parents=True, exist_ok=True)
df = pd.read_csv(ROOT / "data/chalcogen_database_v1_master.csv.gz", low_memory=False)
for c in ["Experimental_Eg_eV", "Eg_eV", "HOMO_eV", "LUMO_eV"]:
    df[c] = pd.to_numeric(df[c], errors="coerce")

PATS = {k: Chem.MolFromSmarts(v) for k, v in {"Se_diaz": "[#7][Se,se][#7]", "Se_phen": "[#6][Se,se][#6]", "S_diaz": "[#7][S,s][#7]",
        "S_phen": "[#6;a][s][#6;a]", "Te_any": "[Te,te]", "O_furan": "[#6;a][o][#6;a]"}.items()}
GEN = {"MolWt": Descriptors.MolWt, "TPSA": Descriptors.TPSA, "LogP": Descriptors.MolLogP, "RotB": Descriptors.NumRotatableBonds,
       "Rings": rdMolDescriptors.CalcNumRings, "AroRings": rdMolDescriptors.CalcNumAromaticRings, "Hetero": rdMolDescriptors.CalcNumHeteroatoms,
       "FracSP3": rdMolDescriptors.CalcFractionCSP3, "BertzCT": Descriptors.BertzCT, "HBA": rdMolDescriptors.CalcNumHBA,
       "HeavyAtoms": lambda m: m.GetNumHeavyAtoms()}

def featurize(smiles):
    rows, fps = [], []
    for smi in smiles:
        m = Chem.MolFromSmiles(smi)
        r = {k: f(m) for k, f in GEN.items()}
        for el, z in [("N", 7), ("O", 8), ("F", 9), ("S", 16), ("Se", 34), ("Te", 52)]:
            r[f"n_{el}"] = sum(1 for at in m.GetAtoms() if at.GetAtomicNum() == z)
        for k, p in PATS.items(): r[f"n_{k}"] = len(m.GetSubstructMatches(p))
        r["conj_frac"] = sum(b.GetIsConjugated() for b in m.GetBonds()) / max(1, m.GetNumBonds())
        rows.append(r); fps.append(AllChem.GetMorganFingerprintAsBitVect(m, 2, 1024))
    D = pd.DataFrame(rows).fillna(0.0)
    return np.hstack([D.to_numpy(float), np.array([list(f) for f in fps], float)]), fps

dev = df[df.Record_ID.str.startswith(("DEV_", "EROL_")) & df.Canonical_SMILES.notna() & df.Eg_eV.notna()]
ext = df[df.Record_ID.str.startswith(("EXT_", "HAKAN_")) & df.Canonical_SMILES.notna() & df.Eg_eV.notna() & (df.Unit_Type.astype(str) != "Polymer")]
pol = df[df.Record_ID.str.startswith("POL_") & df.Canonical_SMILES.notna() & df.Eg_eV.notna() & df.Experimental_Eg_eV.notna()].reset_index(drop=True)
Xdev, _ = featurize(dev.Canonical_SMILES); ydev = dev.Eg_eV.to_numpy()
Xext, _ = featurize(ext.Canonical_SMILES); yext = ext.Eg_eV.to_numpy()
Xpol, fps = featurize(pol.Canonical_SMILES); y = pol.Experimental_Eg_eV.to_numpy()
DFT = pol[["Eg_eV", "HOMO_eV", "LUMO_eV"]].to_numpy()
print("sizes dev/ext/pol", len(ydev), len(yext), len(y), flush=True)

# Butina clusters on polymer oligomers
n = len(fps); dists = []
for i in range(1, n):
    s = DataStructs.BulkTanimotoSimilarity(fps[i], fps[:i]); dists.extend([1 - x for x in s])
clusters = Butina.ClusterData(dists, n, 0.4, isDistData=True)
cl = np.zeros(n, int)
for k, c in enumerate(clusters):
    for i in c: cl[i] = k
sizes = pd.Series(cl).value_counts()
print("butina clusters", len(clusters), "largest", int(sizes.max()), flush=True)

def folds(seed):
    rng = np.random.default_rng(seed); ids = rng.permutation(len(clusters))
    # greedy balanced assignment of clusters to 5 folds
    load = np.zeros(5); assign = {}
    for k in sorted(ids, key=lambda k: -len(clusters[k])):
        f = int(np.argmin(load)); assign[k] = f; load[f] += len(clusters[k])
    fo = np.array([assign[c] for c in cl])
    return [(np.where(fo != f)[0], np.where(fo == f)[0]) for f in range(5)]

P = dict(n_estimators=700, max_depth=5, learning_rate=0.03, subsample=0.8, colsample_bytree=0.5, min_child_weight=3, n_jobs=-1)
dev_model = xgb.XGBRegressor(random_state=0, **P).fit(Xdev, ydev)
devext_model = xgb.XGBRegressor(random_state=0, **P).fit(np.vstack([Xdev, Xext]), np.concatenate([ydev, yext]))
dev_feat = dev_model.predict(Xpol)[:, None]; devext_feat = devext_model.predict(Xpol)[:, None]
print("corr(DEV-style gap, exp gap) =", round(float(np.corrcoef(dev_feat[:, 0], y)[0, 1]), 3),
      " corr(oligomer DFT gap, exp gap) =", round(float(np.corrcoef(DFT[:, 0], y)[0, 1]), 3), flush=True)

def run(setting, use_dft, frac=1.0, seed=0):
    maes = []; pred = np.zeros_like(y)
    for tr, te in folds(seed):
        if frac < 1.0:
            rng = np.random.default_rng(seed + 7); tr = rng.choice(tr, max(30, int(frac * len(tr))), replace=False)
        Xa = Xpol
        if setting == "feat": Xa = np.hstack([Xa, dev_feat])
        if setting == "stack": Xa = np.hstack([Xa, devext_feat])
        if use_dft: Xa = np.hstack([Xa, DFT])
        if setting == "finetune":
            Xs = np.hstack([Xdev, np.zeros((len(Xdev), 3))]) if use_dft else Xdev
            # source model then continued boosting on the target
            src = xgb.XGBRegressor(random_state=seed, **{**P, "n_estimators": 400}).fit(Xs, ydev)
            m = xgb.XGBRegressor(random_state=seed, **{**P, "n_estimators": 400})
            m.fit(Xa[tr], y[tr], xgb_model=src.get_booster())
        else:
            m = xgb.XGBRegressor(random_state=seed, **P).fit(Xa[tr], y[tr])
        p = m.predict(Xa[te]); pred[te] = p; maes.append(mean_absolute_error(y[te], p))
    return float(np.mean(maes)), pred

res = {"n_clusters": len(clusters), "largest_cluster": int(sizes.max())}
for use_dft in (False, True):
    for setting in ("base", "feat", "stack", "finetune"):
        ms, r2s = [], []
        for seed in (0, 1, 2):
            m, p = run(setting, use_dft, seed=seed); ms.append(mean_absolute_error(y, p)); r2s.append(r2_score(y, p))
        key = f"{'struct+DFT' if use_dft else 'struct_only'}|{setting}"
        res[key] = {"MAE_mean": round(float(np.mean(ms)), 4), "MAE_sd": round(float(np.std(ms)), 4), "R2": round(float(np.mean(r2s)), 3)}
        print(key, res[key], flush=True)
# learning curve (struct only, where transfer should matter most)
lc = {}
for frac in (0.1, 0.25, 0.5, 1.0):
    lc[str(frac)] = {s: round(float(np.mean([mean_absolute_error(y, run(s, False, frac, seed)[1]) for seed in (0, 1, 2)])), 4) for s in ("base", "feat")}
    print("learning curve", frac, lc[str(frac)], flush=True)
res["learning_curve_struct_only"] = lc
json.dump(res, open(OUT / "phase2_transfer.json", "w"), indent=1)
