"""Thesis paper, section 1 add-on: learned-representation benchmark for the D-A molecules.
SOAP (DScribe, averaged over atoms, unit-normalized) compared with the interpretable 3D descriptor families under the
same five splits and against the same one-hot identity baseline.
Kernel ridge regression with precomputed Laplacian / Gaussian kernels; alpha and kernel width chosen by 3-fold inner CV
on the training folds only. XGBoost on fold-wise PCA(SOAP). Results are written after every split.
"""
import argparse, gzip, json, sys, time
from pathlib import Path
import numpy as np, pandas as pd
from ase import Atoms
from dscribe.descriptors import SOAP
from rdkit import Chem, RDLogger
from scipy.spatial.distance import cdist
from sklearn.decomposition import PCA
from sklearn.metrics import mean_absolute_error, r2_score
from sklearn.model_selection import KFold
from sklearn.preprocessing import StandardScaler
import xgboost as xgb

RDLogger.DisableLog("rdApp.*")
ap = argparse.ArgumentParser(); ap.add_argument("--root", required=True); ap.add_argument("--out", required=True); a = ap.parse_args()
ROOT, OUT = Path(a.root), Path(a.out)
sys.path.insert(0, str(ROOT)); sys.path.insert(0, str(ROOT / "benchmarks"))
import grouped_cv_donor_acceptor as G
dev, _ = G.load_development_collection(ROOT / "data/chalcogen_database_v1_master.csv.gz")
dg, ag, _ = G.fragment_groups(dev); y = dev["Eg_eV"].to_numpy(float); n = len(y)
pub = dev["Record_ID"].astype(str).str.replace("EROL_", "DEV_", regex=False)
mols = {}
with gzip.open(ROOT / "data/structures_3d_v1.sdf.gz", "rb") as fh:
    for m in Chem.ForwardSDMolSupplier(fh, removeHs=False, sanitize=False):
        if m is None: continue
        mols[m.GetProp("ChalMolDB_Record_ID") if m.HasProp("ChalMolDB_Record_ID") else m.GetProp("_Name")] = m
atoms = [Atoms([at.GetSymbol() for at in mols[r].GetAtoms()], positions=mols[r].GetConformer().GetPositions()) for r in pub]
species = sorted({s for at in atoms for s in at.get_chemical_symbols()})
t0 = time.time()
Xs = SOAP(species=species, r_cut=5.0, n_max=6, l_max=4, average="inner", periodic=False, sparse=False).create(atoms, n_jobs=4)
Xs = Xs[:, Xs.std(0) > 0]; Xs = Xs / np.linalg.norm(Xs, axis=1, keepdims=True)
print("species", species, "SOAP dim", Xs.shape, round(time.time() - t0, 1), "s", flush=True)
F = pd.read_csv(OUT / "s1_descriptors.csv.gz"); F.pop("Eg_eV"); F.pop("Record_ID"); Xd = F.to_numpy(float)
Xd_std = StandardScaler().fit_transform(Xd)  # label-free; only used for the descriptor-KRR comparison

# precomputed distance matrices (label-free)
D = {"SOAP|L1": cdist(Xs, Xs, "cityblock"), "SOAP|L2sq": cdist(Xs, Xs, "sqeuclidean"),
     "3D|L1": cdist(Xd_std, Xd_std, "cityblock")}
print("distances", round(time.time() - t0, 1), "s", flush=True)

def krr_fit_predict(Dm, kind, tr, te, sigmas, alphas=(1e-4, 1e-3, 1e-2, 1e-1)):
    scale = np.median(Dm[np.ix_(tr, tr)])
    def solve(itr, ite, s, al):
        K = np.exp(-Dm[np.ix_(itr, itr)] / (s * scale)); mu = y[itr].mean()
        c = np.linalg.solve(K + al * np.eye(len(itr)), y[itr] - mu)
        return np.exp(-Dm[np.ix_(ite, itr)] / (s * scale)) @ c + mu
    best, bs = None, 1e9
    inner = list(KFold(3, shuffle=True, random_state=1).split(tr))
    for s in sigmas:
        for al in alphas:
            e = np.mean([mean_absolute_error(y[tr[b]], solve(tr[a_], tr[b], s, al)) for a_, b in inner])
            if e < bs: bs, best = e, (s, al)
    return solve(tr, te, *best)

def splits(scheme):
    return list(KFold(5, shuffle=True, random_state=0).split(Xd)) if scheme == "random" else list(G.split_indices(dev, scheme, dg, ag))

def evaluate(fn, scheme):
    yt, yp = [], []
    for tr, te in splits(scheme):
        yt.append(y[te]); yp.append(fn(tr, te))
    yt, yp = np.concatenate(yt), np.concatenate(yp)
    return {"MAE": round(mean_absolute_error(yt, yp), 4), "R2": round(r2_score(yt, yp), 3)}

def xgb_pca(tr, te, extra=False):
    p = PCA(100, random_state=0).fit(Xs[tr]); A, B = p.transform(Xs[tr]), p.transform(Xs[te])
    if extra: A, B = np.hstack([A, Xd[tr]]), np.hstack([B, Xd[te]])
    return xgb.XGBRegressor(n_estimators=600, learning_rate=0.05, max_depth=6, subsample=0.8, colsample_bytree=0.6, random_state=0, n_jobs=-1).fit(A, y[tr]).predict(B)

res = {"soap_dim": int(Xs.shape[1]), "species": species, "soap_params": {"r_cut": 5.0, "n_max": 6, "l_max": 4, "average": "inner", "normalized": True}}
for scheme in ["random", "inchikey", "donor", "acceptor", "donor+acceptor"]:
    r = {"SOAP|KRR-laplacian": evaluate(lambda tr, te: krr_fit_predict(D["SOAP|L1"], "lap", tr, te, (0.05, 0.1, 0.3, 1.0)), scheme),
         "SOAP|KRR-gaussian": evaluate(lambda tr, te: krr_fit_predict(D["SOAP|L2sq"], "gau", tr, te, (0.01, 0.03, 0.1, 0.3)), scheme),
         "SOAP|XGB(PCA100)": evaluate(lambda tr, te: xgb_pca(tr, te), scheme),
         "3D descriptors|KRR-laplacian": evaluate(lambda tr, te: krr_fit_predict(D["3D|L1"], "lap", tr, te, (0.05, 0.1, 0.3, 1.0)), scheme),
         "SOAP(PCA100)+3D|XGB": evaluate(lambda tr, te: xgb_pca(tr, te, extra=True), scheme)}
    res[scheme] = r; print(scheme, json.dumps(r), round(time.time() - t0), "s", flush=True)
    json.dump(res, open(OUT / "s1_soap_results.json", "w"), indent=1)
