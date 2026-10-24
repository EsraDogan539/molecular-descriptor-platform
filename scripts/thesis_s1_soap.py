"""Thesis paper, section 1 add-on: learned-representation benchmark for the D-A molecules.
SOAP (DScribe, averaged over atoms) and MBTR-free Coulomb-type baseline compared with the interpretable 3D descriptor
families under the same five splits and against the same one-hot identity baseline.
Models: kernel ridge regression (RBF / Laplacian; alpha and gamma tuned by 3-fold inner CV on the training folds) and XGBoost.
"""
import argparse, gzip, json, sys, time
from pathlib import Path
import numpy as np, pandas as pd
from ase import Atoms
from dscribe.descriptors import SOAP
from rdkit import Chem, RDLogger
from sklearn.decomposition import PCA
from sklearn.kernel_ridge import KernelRidge
from sklearn.metrics import mean_absolute_error, r2_score
from sklearn.model_selection import GridSearchCV, KFold
from sklearn.pipeline import make_pipeline
from sklearn.preprocessing import StandardScaler
import xgboost as xgb

RDLogger.DisableLog("rdApp.*")
ap = argparse.ArgumentParser(); ap.add_argument("--root", required=True); ap.add_argument("--out", required=True); a = ap.parse_args()
ROOT, OUT = Path(a.root), Path(a.out)
sys.path.insert(0, str(ROOT)); sys.path.insert(0, str(ROOT / "benchmarks"))
import grouped_cv_donor_acceptor as G
dev, _ = G.load_development_collection(ROOT / "data/chalcogen_database_v1_master.csv.gz")
dg, ag, _ = G.fragment_groups(dev); y = dev["Eg_eV"].to_numpy(float)
pub = dev["Record_ID"].astype(str).str.replace("EROL_", "DEV_", regex=False)
mols = {}
with gzip.open(ROOT / "data/structures_3d_v1.sdf.gz", "rb") as fh:
    for m in Chem.ForwardSDMolSupplier(fh, removeHs=False, sanitize=False):
        if m is None: continue
        mols[m.GetProp("ChalMolDB_Record_ID") if m.HasProp("ChalMolDB_Record_ID") else m.GetProp("_Name")] = m
atoms = []
for r in pub:
    m = mols[r]; atoms.append(Atoms([at.GetSymbol() for at in m.GetAtoms()], positions=m.GetConformer().GetPositions()))
species = sorted({s for at in atoms for s in at.get_chemical_symbols()})
t0 = time.time()
soap = SOAP(species=species, r_cut=5.0, n_max=6, l_max=4, average="inner", periodic=False, sparse=False)
Xs = soap.create(atoms, n_jobs=4); Xs = Xs[:, Xs.std(0) > 0]
print("species", species, "SOAP dim", Xs.shape, round(time.time() - t0, 1), "s", flush=True)
F = pd.read_csv(OUT / "s1_descriptors.csv.gz"); F.pop("Eg_eV"); F.pop("Record_ID"); Xd = F.to_numpy(float)

def splits(scheme):
    return list(KFold(5, shuffle=True, random_state=0).split(Xd)) if scheme == "random" else list(G.split_indices(dev, scheme, dg, ag))

def run(X, model, scheme, pca=None, extra=None):
    yt, yp = [], []
    for tr, te in splits(scheme):
        Xtr, Xte = X[tr], X[te]
        if pca:
            sc = StandardScaler().fit(Xtr); p = PCA(pca, random_state=0).fit(sc.transform(Xtr)); Xtr, Xte = p.transform(sc.transform(Xtr)), p.transform(sc.transform(Xte))
        if extra is not None: Xtr, Xte = np.hstack([Xtr, extra[tr]]), np.hstack([Xte, extra[te]])
        if model == "xgb":
            m = xgb.XGBRegressor(n_estimators=600, learning_rate=0.05, max_depth=6, subsample=0.8, colsample_bytree=0.6, random_state=0, n_jobs=-1).fit(Xtr, y[tr])
        else:
            d = Xtr.shape[1]
            pipe = make_pipeline(StandardScaler(), KernelRidge(kernel=model))
            gam = [g / d for g in ([0.1, 0.3, 1.0, 3.0] if model == "rbf" else [0.03, 0.1, 0.3, 1.0])]
            m = GridSearchCV(pipe, {"kernelridge__alpha": [1e-3, 1e-2, 1e-1], "kernelridge__gamma": gam}, cv=3, scoring="neg_mean_absolute_error", n_jobs=-1).fit(Xtr, y[tr])
        yt.append(y[te]); yp.append(m.predict(Xte))
    yt, yp = np.concatenate(yt), np.concatenate(yp)
    return {"MAE": round(mean_absolute_error(yt, yp), 4), "R2": round(r2_score(yt, yp), 3)}

res = {"soap_dim": int(Xs.shape[1]), "species": species, "soap_params": {"r_cut": 5.0, "n_max": 6, "l_max": 4, "average": "inner"}}
for scheme in ["random", "inchikey", "donor", "acceptor", "donor+acceptor"]:
    r = {}
    r["SOAP|KRR-laplacian"] = run(Xs, "laplacian", scheme)
    r["SOAP|KRR-rbf"] = run(Xs, "rbf", scheme)
    r["SOAP|XGB(PCA200)"] = run(Xs, "xgb", scheme, pca=200)
    r["3D descriptors|KRR-laplacian"] = run(Xd, "laplacian", scheme)
    r["SOAP+3D|XGB"] = run(Xs, "xgb", scheme, pca=100, extra=Xd)
    res[scheme] = r; print(scheme, json.dumps(r), round(time.time() - t0), "s", flush=True)
json.dump(res, open(OUT / "s1_soap_results.json", "w"), indent=1)
