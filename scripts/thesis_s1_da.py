"""Thesis paper, section 1: interpretable 3D-structure descriptors for the DFT HOMO-LUMO gap of the
3,088 donor-acceptor molecules (ChalMolDB development collection), re-evaluated under grouped splits.

Descriptor families (computed from the published DFT geometries, B3LYP/6-311+G(d)):
  base      composition, ring/size/topology, global XYZ geometry (Rg, interatomic distances, vdW size)
  bondtype  bond counts by element pair and bond order, and their ratios
  bondlen   mean and SD of bond lengths per element pair (C-C, C-N, C-S, C-Se, ...), bond-length alternation
  planar    out-of-plane metrics: RMS deviation from the best plane, z-range, smallest principal moment ratio
Splits: random 5-fold (as in the original draft), InChIKey, leave-donor-out, leave-acceptor-out, donor+acceptor.
Baselines: training mean; ridge on one-hot donor+acceptor identity ("fragment" memorisation baseline).
Models: SVR (RBF), Extra Trees, XGBoost. Ablation and SHAP (XGBoost) on the full and grouped settings.
"""
import argparse, gzip, json, re, sys
from pathlib import Path
import numpy as np, pandas as pd
from rdkit import Chem, RDLogger
from rdkit.Chem import Descriptors, rdMolDescriptors
from sklearn.dummy import DummyRegressor
from sklearn.ensemble import ExtraTreesRegressor
from sklearn.linear_model import Ridge
from sklearn.metrics import mean_absolute_error, mean_squared_error, r2_score
from sklearn.model_selection import KFold
from sklearn.pipeline import make_pipeline
from sklearn.preprocessing import OneHotEncoder, StandardScaler
from sklearn.svm import SVR
import xgboost as xgb

RDLogger.DisableLog("rdApp.*")
ap = argparse.ArgumentParser(); ap.add_argument("--root", required=True); ap.add_argument("--out", required=True); a = ap.parse_args()
ROOT, OUT = Path(a.root), Path(a.out); OUT.mkdir(parents=True, exist_ok=True)
sys.path.insert(0, str(ROOT)); sys.path.insert(0, str(ROOT / "benchmarks"))
import grouped_cv_donor_acceptor as G

dev, excluded = G.load_development_collection(ROOT / "data/chalcogen_database_v1_master.csv.gz")
donor_group, acceptor_group, _ = G.fragment_groups(dev)
y = dev["Eg_eV"].to_numpy(float)

# ---------------- 3D structures by record id
mols = {}
with gzip.open(ROOT / "data/structures_3d_v1.sdf.gz", "rb") as fh:
    for m in Chem.ForwardSDMolSupplier(fh, removeHs=False, sanitize=False):
        if m is None: continue
        rid = m.GetProp("ChalMolDB_Record_ID") if m.HasProp("ChalMolDB_Record_ID") else m.GetProp("_Name")
        mols[rid] = m
pub = dev["Record_ID"].astype(str).str.replace("EROL_", "DEV_", regex=False)
have = [r in mols for r in pub]
print("records", len(dev), "with 3D", sum(have), flush=True)

EL = ["C", "H", "N", "O", "S", "Se", "Te", "F", "Cl", "Si"]
PAIRS = ["C-C", "C-H", "C-N", "C-O", "C-S", "C-Se", "C-Te", "N-S", "N-Se", "N-Te", "C-F", "N-N", "N-O", "C-Si"]

def feats(m, smi):
    f = {}
    syms = [at.GetSymbol() for at in m.GetAtoms()]
    pos = m.GetConformer().GetPositions()
    heavy = [i for i, s in enumerate(syms) if s != "H"]
    # base
    for e in EL: f[f"base_n_{e}"] = syms.count(e)
    mm = Chem.MolFromSmiles(smi)
    f["base_heavy"] = len(heavy); f["base_rings"] = rdMolDescriptors.CalcNumRings(mm); f["base_arom_rings"] = rdMolDescriptors.CalcNumAromaticRings(mm)
    f["base_rotb"] = rdMolDescriptors.CalcNumRotatableBonds(mm); f["base_molwt"] = Descriptors.MolWt(mm); f["base_tpsa"] = Descriptors.TPSA(mm)
    f["base_fsp3"] = rdMolDescriptors.CalcFractionCSP3(mm)
    P = pos[heavy]; c = P - P.mean(0)
    f["base_rg"] = float(np.sqrt((c ** 2).sum(1).mean()))
    d = np.sqrt(((P[:, None] - P[None]) ** 2).sum(-1)); iu = np.triu_indices(len(P), 1)
    f["base_dmax"] = float(d[iu].max()); f["base_dmean"] = float(d[iu].mean()); f["base_dstd"] = float(d[iu].std())
    # planarity
    w, v = np.linalg.eigh(c.T @ c / len(c))
    f["planar_rms"] = float(np.sqrt(w[0])); f["planar_ratio"] = float(w[0] / max(w[2], 1e-9))
    f["planar_zrange"] = float(np.ptp(c @ v[:, 0]))
    f["planar_mid_ratio"] = float(w[1] / max(w[2], 1e-9))
    # bonds: in the SDF bond orders come from the source file
    lens = {p: [] for p in PAIRS}; order_cnt = {"single": 0, "double": 0, "triple": 0, "aromatic": 0}
    cc_single, cc_double = [], []
    for b in m.GetBonds():
        s1, s2 = sorted([syms[b.GetBeginAtomIdx()], syms[b.GetEndAtomIdx()]], key=lambda s: (s != "C", s))
        key = f"{s1}-{s2}" if f"{s1}-{s2}" in lens else f"{s2}-{s1}"
        L = float(np.linalg.norm(pos[b.GetBeginAtomIdx()] - pos[b.GetEndAtomIdx()]))
        bt = b.GetBondTypeAsDouble()
        order_cnt[{1.0: "single", 2.0: "double", 3.0: "triple", 1.5: "aromatic"}.get(bt, "single")] += 1
        if key in lens: lens[key].append(L)
        if key == "C-C":
            (cc_double if bt >= 1.5 else cc_single).append(L)
    nb = max(1, m.GetNumBonds())
    for p in PAIRS:
        f[f"bondtype_n_{p}"] = len(lens[p]); f[f"bondtype_frac_{p}"] = len(lens[p]) / nb
    for k, v_ in order_cnt.items(): f[f"bondtype_{k}_frac"] = v_ / nb
    for p in PAIRS:
        if p == "C-H": continue
        f[f"bondlen_mean_{p}"] = float(np.mean(lens[p])) if lens[p] else np.nan
        f[f"bondlen_sd_{p}"] = float(np.std(lens[p])) if len(lens[p]) > 1 else np.nan
    f["bondlen_BLA_CC"] = (float(np.mean(cc_single)) - float(np.mean(cc_double))) if cc_single and cc_double else np.nan
    return f

rows = []
for r, smi in zip(pub, dev["Canonical_SMILES"]):
    rows.append(feats(mols[r], smi) if r in mols else {})
F = pd.DataFrame(rows)
# missing bond-length columns: impute with column median, keep the count columns as the indicator
F = F.apply(lambda col: col.fillna(col.median()) if col.notna().any() else col.fillna(0.0))
F = F.loc[:, F.std() > 0]
fam = {k: [c for c in F.columns if c.startswith(k + "_")] for k in ["base", "bondtype", "bondlen", "planar"]}
print({k: len(v) for k, v in fam.items()}, flush=True)

def models():
    return {"SVR": lambda: make_pipeline(StandardScaler(), SVR(C=10.0, epsilon=0.05, gamma="scale")),
            "ExtraTrees": lambda: ExtraTreesRegressor(n_estimators=500, n_jobs=-1, random_state=0),
            "XGB": lambda: xgb.XGBRegressor(n_estimators=600, learning_rate=0.05, max_depth=6, subsample=0.8, colsample_bytree=0.6, random_state=0, n_jobs=-1)}

def splits(scheme):
    if scheme == "random":
        return list(KFold(5, shuffle=True, random_state=0).split(F))
    return list(G.split_indices(dev, scheme, donor_group, acceptor_group))

SCHEMES = ["random", "inchikey", "donor", "acceptor", "donor+acceptor"]
frag = pd.DataFrame({"d": donor_group.values, "a": acceptor_group.values})

def evaluate(cols, model_fn, scheme):
    X = F[cols].to_numpy(float); yt, yp = [], []
    for tr, te in splits(scheme):
        m = model_fn().fit(X[tr], y[tr]); yt.append(y[te]); yp.append(m.predict(X[te]))
    yt, yp = np.concatenate(yt), np.concatenate(yp)
    return {"MAE": round(mean_absolute_error(yt, yp), 4), "RMSE": round(float(np.sqrt(mean_squared_error(yt, yp))), 4), "R2": round(r2_score(yt, yp), 3)}

def baseline(kind, scheme):
    yt, yp = [], []
    for tr, te in splits(scheme):
        if kind == "mean": p = np.full(len(te), y[tr].mean())
        else:
            enc = OneHotEncoder(handle_unknown="ignore"); Xtr = enc.fit_transform(frag.iloc[tr]); Xte = enc.transform(frag.iloc[te])
            p = Ridge(alpha=1.0).fit(Xtr, y[tr]).predict(Xte)
        yt.append(y[te]); yp.append(p)
    yt, yp = np.concatenate(yt), np.concatenate(yp)
    return {"MAE": round(mean_absolute_error(yt, yp), 4), "RMSE": round(float(np.sqrt(mean_squared_error(yt, yp))), 4), "R2": round(r2_score(yt, yp), 3)}

res = {"n": int(len(y)), "excluded": len(excluded), "families": {k: len(v) for k, v in fam.items()}, "main": {}, "ablation": {}}
ALL = sum(fam.values(), [])
for scheme in SCHEMES:
    res["main"][scheme] = {"mean": baseline("mean", scheme), "fragment_onehot": baseline("frag", scheme)}
    for name, fn in models().items():
        res["main"][scheme][name] = evaluate(ALL, fn, scheme)
    print(scheme, json.dumps(res["main"][scheme]), flush=True)

ABL = {"base": fam["base"], "base+planar": fam["base"] + fam["planar"], "base+bondtype": fam["base"] + fam["bondtype"],
       "base+bondlen": fam["base"] + fam["bondlen"], "all": ALL}
for scheme in ["random", "donor", "acceptor"]:
    res["ablation"][scheme] = {k: evaluate(v, models()["XGB"], scheme) for k, v in ABL.items()}
    print("ablation", scheme, json.dumps(res["ablation"][scheme]), flush=True)

# SHAP (XGB, all data) - top 15 by mean |SHAP|
try:
    import shap
    m = models()["XGB"]().fit(F[ALL].to_numpy(float), y)
    sv = shap.TreeExplainer(m).shap_values(F[ALL].to_numpy(float))
    imp = pd.Series(np.abs(sv).mean(0), index=ALL).sort_values(ascending=False)
    res["shap_top15"] = {k: round(float(v), 4) for k, v in imp.head(15).items()}
    res["shap_by_family"] = {k: round(float(imp[[c for c in v]].sum() / imp.sum()), 3) for k, v in fam.items()}
    print("shap", json.dumps(res["shap_top15"]), json.dumps(res["shap_by_family"]), flush=True)
except Exception as e:
    res["shap_error"] = repr(e)
F.assign(Record_ID=pub.values, Eg_eV=y).to_csv(OUT / "s1_descriptors.csv.gz", index=False)
json.dump(res, open(OUT / "s1_results.json", "w"), indent=1)
