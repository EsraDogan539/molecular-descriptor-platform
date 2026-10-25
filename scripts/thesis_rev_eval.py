"""Revision run 2: (1) D-A molecules, five splits x five seeds (fold assignments), XGBoost on
    3D descriptors from DFT geometries / from MMFF geometries / RDKit 2D descriptors / Morgan counts, plus one-hot identity;
    mean and SD over seeds. (2) Polymers with oligomer 3D descriptors: delta-learning with and without the 67 descriptors on
    the 1,064 polymers with geometries (cluster-grouped CV, 3 seeds). (3) Out-of-fold SHAP consistency for section 1.
"""
import argparse, json, sys, time
from pathlib import Path
import numpy as np, pandas as pd
from rdkit import Chem, DataStructs, RDLogger
from rdkit.Chem import AllChem, Descriptors, rdMolDescriptors
from rdkit.ML.Cluster import Butina
from sklearn.linear_model import LinearRegression, Ridge
from sklearn.metrics import mean_absolute_error, r2_score
from sklearn.model_selection import KFold
from sklearn.preprocessing import OneHotEncoder
import xgboost as xgb, shap
RDLogger.DisableLog("rdApp.*")
ap = argparse.ArgumentParser(); ap.add_argument("--root", required=True); ap.add_argument("--out", required=True); a = ap.parse_args()
ROOT, OUT = Path(a.root), Path(a.out)
sys.path.insert(0, str(ROOT)); sys.path.insert(0, str(ROOT / "benchmarks"))
import grouped_cv_donor_acceptor as G
t0 = time.time(); res = {}
dev, _ = G.load_development_collection(ROOT / "data/chalcogen_database_v1_master.csv.gz")
dg, ag, _ = G.fragment_groups(dev); y = dev["Eg_eV"].to_numpy(float)
def load(fn):
    d = pd.read_csv(OUT / fn); d = d.drop(columns=[c for c in ("Record_ID", "Eg_eV") if c in d]); d = d.loc[:, d.notna().any()]
    d = d.apply(lambda c: c.fillna(c.median())); return d.loc[:, d.std() > 0].to_numpy(float)
X = {"3D (DFT geometry)": load("s1_descriptors.csv.gz"), "3D (MMFF geometry)": load("rev_desc_mmff.csv.gz"),
     "RDKit 2D descriptors": load("rev_desc_rdkit2d.csv.gz"), "Morgan counts (r=2, 2048)": load("rev_morgan_counts.csv.gz")}
X["3D (MMFF) + RDKit 2D"] = np.hstack([X["3D (MMFF geometry)"], X["RDKit 2D descriptors"]])
print({k: v.shape for k, v in X.items()}, flush=True)
frag = pd.DataFrame({"d": dg.values, "a": ag.values})
P = dict(n_estimators=600, learning_rate=0.05, max_depth=6, subsample=0.8, colsample_bytree=0.6, n_jobs=-1, tree_method="hist")
SCH = ["random", "donor", "acceptor", "donor+acceptor"]; SEEDS = [0, 1, 2]
def splits(s, seed): return list(KFold(5, shuffle=True, random_state=seed).split(y)) if s == "random" else list(G.split_indices(dev, s, dg, ag, seed=seed))
res["s1_seeds"] = {}
for s in SCH:
    out = {}
    for name in ["one-hot identity"] + list(X):
        maes, r2s = [], []
        for seed in SEEDS:
            yt, yp = [], []
            for tr, te in splits(s, seed):
                if name == "one-hot identity":
                    enc = OneHotEncoder(handle_unknown="ignore"); p = Ridge(alpha=1.0).fit(enc.fit_transform(frag.iloc[tr]), y[tr]).predict(enc.transform(frag.iloc[te]))
                else:
                    p = xgb.XGBRegressor(random_state=seed, **P).fit(X[name][tr], y[tr]).predict(X[name][te])
                yt.append(y[te]); yp.append(p)
            yt, yp = np.concatenate(yt), np.concatenate(yp); maes.append(mean_absolute_error(yt, yp)); r2s.append(r2_score(yt, yp))
        out[name] = {"MAE_mean": round(float(np.mean(maes)), 4), "MAE_sd": round(float(np.std(maes)), 4), "R2_mean": round(float(np.mean(r2s)), 3)}
    res["s1_seeds"][s] = out; print(s, json.dumps(out), round(time.time() - t0), "s", flush=True)
    json.dump(res, open(OUT / "rev_eval.json", "w"), indent=1)

# (3) out-of-fold SHAP (random and acceptor splits, seed 0) vs in-sample ranking
Xd = X["3D (DFT geometry)"]; cols = [c for c in pd.read_csv(OUT / "s1_descriptors.csv.gz", nrows=5).columns if c not in ("Record_ID", "Eg_eV")]
Fd = pd.read_csv(OUT / "s1_descriptors.csv.gz").drop(columns=["Record_ID", "Eg_eV"]); Fd = Fd.loc[:, Fd.std() > 0]; cols = list(Fd.columns)
def shap_rank(pairs):
    sv = np.zeros_like(Xd)
    for tr, te in pairs:
        m = xgb.XGBRegressor(random_state=0, **P).fit(Xd[tr], y[tr]); sv[te] = shap.TreeExplainer(m).shap_values(Xd[te])
    return pd.Series(np.abs(sv).mean(0), index=cols).sort_values(ascending=False)
ins = shap_rank([(np.arange(len(y)), np.arange(len(y)))])
res["shap_oof"] = {}
for s in ["random", "acceptor"]:
    r = shap_rank(splits(s, 0)); k = [c for c in r.index if c in ins.index]
    from scipy.stats import spearmanr
    res["shap_oof"][s] = {"top10": list(r.index[:10]), "top10_overlap_with_insample": len(set(r.index[:10]) & set(ins.index[:10])),
                          "spearman_all": round(float(spearmanr(r[k], ins[k]).correlation), 3)}
print("shap oof", json.dumps(res["shap_oof"]), flush=True)
json.dump(res, open(OUT / "rev_eval.json", "w"), indent=1)

# (2) polymers with oligomer 3D descriptors
df = pd.read_csv(ROOT / "data/chalcogen_database_v1_master.csv.gz", low_memory=False)
for c in ["Experimental_Eg_eV", "Eg_eV", "HOMO_eV", "LUMO_eV"]: df[c] = pd.to_numeric(df[c], errors="coerce")
P3 = pd.read_csv(OUT / "rev_desc_polymer3d.csv.gz")
pol = df[df.Record_ID.str.startswith("POL_") & df.Canonical_SMILES.notna() & df.Eg_eV.notna() & df.Experimental_Eg_eV.notna()]
pol = pol.merge(P3, on="Record_ID", how="inner").reset_index(drop=True)
d3cols = [c for c in P3.columns if c != "Record_ID" and c in pol]
D3 = pol[d3cols].apply(lambda c: c.fillna(c.median()) if c.notna().any() else c.fillna(0)); D3 = D3.loc[:, D3.std() > 0].to_numpy(float)
PATS = {k: Chem.MolFromSmarts(v) for k, v in {"Se_diaz": "[#7][Se,se][#7]", "Se_phen": "[#6][Se,se][#6]", "S_diaz": "[#7][S,s][#7]", "S_phen": "[#6;a][s][#6;a]", "O_furan": "[#6;a][o][#6;a]"}.items()}
GEN = {"MolWt": Descriptors.MolWt, "TPSA": Descriptors.TPSA, "LogP": Descriptors.MolLogP, "RotB": Descriptors.NumRotatableBonds, "Rings": rdMolDescriptors.CalcNumRings,
       "AroRings": rdMolDescriptors.CalcNumAromaticRings, "Hetero": rdMolDescriptors.CalcNumHeteroatoms, "FracSP3": rdMolDescriptors.CalcFractionCSP3,
       "BertzCT": Descriptors.BertzCT, "HBA": rdMolDescriptors.CalcNumHBA, "HeavyAtoms": lambda m: m.GetNumHeavyAtoms()}
rows, fps = [], []
for smi in pol.Canonical_SMILES:
    m = Chem.MolFromSmiles(smi); r = {k: f(m) for k, f in GEN.items()}
    for el, z in [("N", 7), ("O", 8), ("F", 9), ("S", 16), ("Se", 34)]: r[f"n_{el}"] = sum(at.GetAtomicNum() == z for at in m.GetAtoms())
    for k, p in PATS.items(): r[f"n_{k}"] = len(m.GetSubstructMatches(p))
    r["conj_frac"] = sum(b.GetIsConjugated() for b in m.GetBonds()) / max(1, m.GetNumBonds()); rows.append(r); fps.append(AllChem.GetMorganFingerprintAsBitVect(m, 2, 1024))
Xs = np.hstack([pd.DataFrame(rows).fillna(0).to_numpy(float), np.array([list(f) for f in fps], float)])
yp_ = pol.Experimental_Eg_eV.to_numpy(); DFT = pol[["Eg_eV", "HOMO_eV", "LUMO_eV"]].to_numpy(); n = len(yp_)
dists = []
for i in range(1, n): dists.extend([1 - s for s in DataStructs.BulkTanimotoSimilarity(fps[i], fps[:i])])
clusters = Butina.ClusterData(dists, n, 0.4, isDistData=True); cl = np.zeros(n, int)
for k, c in enumerate(clusters):
    for i in c: cl[i] = k
def folds(seed):
    rng = np.random.default_rng(seed); ids = rng.permutation(len(clusters)); load_ = np.zeros(5); asg = {}
    for k in sorted(ids, key=lambda k: -len(clusters[k])): f = int(np.argmin(load_)); asg[k] = f; load_[f] += len(clusters[k])
    fo = np.array([asg[c] for c in cl]); return [(np.where(fo != f)[0], np.where(fo == f)[0]) for f in range(5)]
PP = dict(n_estimators=700, max_depth=5, learning_rate=0.03, subsample=0.8, colsample_bytree=0.5, min_child_weight=3, n_jobs=-1)
res["polymer3d"] = {"n": int(n), "n_3d_desc": int(D3.shape[1])}
for name, XX in [("2D+Morgan+DFT", np.hstack([Xs, DFT])), ("2D+Morgan+DFT+3D", np.hstack([Xs, DFT, D3])), ("3D+DFT", np.hstack([D3, DFT]))]:
    maes, r2s = [], []
    for seed in (0, 1, 2):
        p = np.zeros(n)
        for tr, te in folds(seed):
            cal = LinearRegression().fit(DFT[tr, :1], yp_[tr]); m = xgb.XGBRegressor(random_state=seed, **PP).fit(XX[tr], yp_[tr] - cal.predict(DFT[tr, :1]))
            p[te] = cal.predict(DFT[te, :1]) + m.predict(XX[te])
        maes.append(mean_absolute_error(yp_, p)); r2s.append(r2_score(yp_, p))
    res["polymer3d"][f"delta|{name}"] = {"MAE_mean": round(float(np.mean(maes)), 4), "MAE_sd": round(float(np.std(maes)), 4), "R2": round(float(np.mean(r2s)), 3)}
    print("polymer3d", name, res["polymer3d"][f"delta|{name}"], flush=True)
json.dump(res, open(OUT / "rev_eval.json", "w"), indent=1)
print("done", round(time.time() - t0), "s")
