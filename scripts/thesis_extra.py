"""Thesis paper, supporting runs.
(1) Polymer optical gap: random 5-fold vs cluster-grouped CV for the same models, on the 1,096 Nagasawa/Liu subset
    and on all 1,323 records, to place the grouped result next to the published random-split value.
(2) Section 1 out-of-fold predictions (XGB on all descriptor families, one-hot baseline) for the parity figures.
"""
import argparse, json, sys
from pathlib import Path
import numpy as np, pandas as pd
from rdkit import Chem, DataStructs, RDLogger
from rdkit.Chem import AllChem, Descriptors, rdMolDescriptors
from rdkit.ML.Cluster import Butina
from sklearn.linear_model import LinearRegression, Ridge
from sklearn.metrics import mean_absolute_error, r2_score
from sklearn.model_selection import KFold
from sklearn.preprocessing import OneHotEncoder
import xgboost as xgb

RDLogger.DisableLog("rdApp.*")
ap = argparse.ArgumentParser(); ap.add_argument("--root", required=True); ap.add_argument("--out", required=True); a = ap.parse_args()
ROOT, OUT = Path(a.root), Path(a.out)
sys.path.insert(0, str(ROOT)); sys.path.insert(0, str(ROOT / "benchmarks")); sys.path.insert(0, str(Path(__file__).parent))
df = pd.read_csv(ROOT / "data/chalcogen_database_v1_master.csv.gz", low_memory=False)
for c in ["Experimental_Eg_eV", "Eg_eV", "HOMO_eV", "LUMO_eV"]:
    df[c] = pd.to_numeric(df[c], errors="coerce")
GEN = {"MolWt": Descriptors.MolWt, "TPSA": Descriptors.TPSA, "LogP": Descriptors.MolLogP, "RotB": Descriptors.NumRotatableBonds,
       "Rings": rdMolDescriptors.CalcNumRings, "AroRings": rdMolDescriptors.CalcNumAromaticRings, "Hetero": rdMolDescriptors.CalcNumHeteroatoms,
       "FracSP3": rdMolDescriptors.CalcFractionCSP3, "BertzCT": Descriptors.BertzCT, "HBA": rdMolDescriptors.CalcNumHBA, "HeavyAtoms": lambda m: m.GetNumHeavyAtoms()}
PATS = {k: Chem.MolFromSmarts(v) for k, v in {"Se_diaz": "[#7][Se,se][#7]", "Se_phen": "[#6][Se,se][#6]", "S_diaz": "[#7][S,s][#7]",
        "S_phen": "[#6;a][s][#6;a]", "Te_any": "[Te,te]", "Te_diaz": "[#7][Te,te][#7]", "Te_phen": "[#6][Te,te][#6]", "O_furan": "[#6;a][o][#6;a]"}.items()}
def featurize(smiles):
    rows, fps = [], []
    for smi in smiles:
        m = Chem.MolFromSmiles(smi); r = {k: f(m) for k, f in GEN.items()}
        for el, z in [("N", 7), ("O", 8), ("F", 9), ("S", 16), ("Se", 34), ("Te", 52)]: r[f"n_{el}"] = sum(1 for at in m.GetAtoms() if at.GetAtomicNum() == z)
        for k, p in PATS.items(): r[f"n_{k}"] = len(m.GetSubstructMatches(p))
        r["conj_frac"] = sum(b.GetIsConjugated() for b in m.GetBonds()) / max(1, m.GetNumBonds())
        rows.append(r); fps.append(AllChem.GetMorganFingerprintAsBitVect(m, 2, 1024))
    return np.hstack([pd.DataFrame(rows).fillna(0).to_numpy(float), np.array([list(f) for f in fps], float)]), fps
P = dict(n_estimators=700, max_depth=5, learning_rate=0.03, subsample=0.8, colsample_bytree=0.5, min_child_weight=3, n_jobs=-1)
res = {}
polall = df[df.Record_ID.str.startswith("POL_") & df.Canonical_SMILES.notna() & df.Eg_eV.notna() & df.Experimental_Eg_eV.notna()].reset_index(drop=True)
for name, pol in [("all_1323", polall), ("nagasawa_1096", polall[polall.Record_ID.str[4:].astype(int) <= 1096].reset_index(drop=True))]:
    X, fps = featurize(pol.Canonical_SMILES); y = pol.Experimental_Eg_eV.to_numpy(); DFT = pol[["Eg_eV", "HOMO_eV", "LUMO_eV"]].to_numpy(); Xd = np.hstack([X, DFT])
    n = len(y); dists = []
    for i in range(1, n): dists.extend([1 - s for s in DataStructs.BulkTanimotoSimilarity(fps[i], fps[:i])])
    clusters = Butina.ClusterData(dists, n, 0.4, isDistData=True); cl = np.zeros(n, int)
    for k, c in enumerate(clusters):
        for i in c: cl[i] = k
    def grouped(seed):
        rng = np.random.default_rng(seed); ids = rng.permutation(len(clusters)); load = np.zeros(5); asg = {}
        for k in sorted(ids, key=lambda k: -len(clusters[k])): f = int(np.argmin(load)); asg[k] = f; load[f] += len(clusters[k])
        fo = np.array([asg[c] for c in cl]); return [(np.where(fo != f)[0], np.where(fo == f)[0]) for f in range(5)]
    def rand(seed): return list(KFold(5, shuffle=True, random_state=seed).split(y))
    res[name] = {"n": int(n)}
    for sname, sp in [("random", rand), ("cluster", grouped)]:
        for kind in ["dft_linear", "struct+dft_xgb", "delta"]:
            maes, r2s = [], []
            for seed in (0, 1, 2):
                p = np.zeros(n)
                for tr, te in sp(seed):
                    if kind == "dft_linear": p[te] = LinearRegression().fit(DFT[tr, :1], y[tr]).predict(DFT[te, :1]); continue
                    if kind == "delta":
                        cal = LinearRegression().fit(DFT[tr, :1], y[tr]); m = xgb.XGBRegressor(random_state=seed, **P).fit(Xd[tr], y[tr] - cal.predict(DFT[tr, :1]))
                        p[te] = cal.predict(DFT[te, :1]) + m.predict(Xd[te])
                    else: p[te] = xgb.XGBRegressor(random_state=seed, **P).fit(Xd[tr], y[tr]).predict(Xd[te])
                maes.append(mean_absolute_error(y, p)); r2s.append(r2_score(y, p))
            res[name][f"{sname}|{kind}"] = {"MAE": round(float(np.mean(maes)), 4), "R2": round(float(np.mean(r2s)), 3)}
            print(name, sname, kind, res[name][f"{sname}|{kind}"], flush=True)

# (2) S1 OOF predictions
import grouped_cv_donor_acceptor as G
dev, _ = G.load_development_collection(ROOT / "data/chalcogen_database_v1_master.csv.gz")
dg, ag, _ = G.fragment_groups(dev)
F = pd.read_csv(OUT / "s1_descriptors.csv.gz"); y = F.pop("Eg_eV").to_numpy(); F.pop("Record_ID"); X = F.to_numpy(float)
frag = pd.DataFrame({"d": dg.values, "a": ag.values}); oof = {"y": y}
for scheme in ["random", "acceptor", "donor+acceptor"]:
    sp = list(KFold(5, shuffle=True, random_state=0).split(X)) if scheme == "random" else list(G.split_indices(dev, scheme, dg, ag))
    p1, p2 = np.zeros(len(y)), np.zeros(len(y))
    for tr, te in sp:
        p1[te] = xgb.XGBRegressor(n_estimators=600, learning_rate=0.05, max_depth=6, subsample=0.8, colsample_bytree=0.6, random_state=0, n_jobs=-1).fit(X[tr], y[tr]).predict(X[te])
        enc = OneHotEncoder(handle_unknown="ignore"); p2[te] = Ridge(alpha=1.0).fit(enc.fit_transform(frag.iloc[tr]), y[tr]).predict(enc.transform(frag.iloc[te]))
    oof[f"xgb_{scheme}"] = p1; oof[f"onehot_{scheme}"] = p2
pd.DataFrame(oof).to_csv(OUT / "s1_oof.csv", index=False)
json.dump(res, open(OUT / "extra_results.json", "w"), indent=1)
