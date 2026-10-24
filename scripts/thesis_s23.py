"""Thesis paper, sections 2 and 3.

S2  Experimental optical gap of the ChalMolDB polymer records (POL_): final model.
    Butina clusters (Tanimoto distance 0.4) -> balanced 5-fold group CV, 3 seeds.
    Models: mean; linear DFT calibration; struct XGB; struct+DFT XGB; delta (linear DFT calibration + XGB residual).
    Uncertainty: split-conformal 90 % intervals (calibrated on out-of-fold residuals of the other folds).
    Applicability domain: error vs. maximum Tanimoto similarity of a test record to its training folds.
    SHAP of the delta residual model (descriptor features individually, Morgan bits summed).
S3  Where the chalcogen sits.
    (a) 139 experimental S->Se swap pairs: donor selenophene vs acceptor selenadiazole; bootstrap 95 % CI and
        Mann-Whitney for d(Eg_exp), d(Eg_DFT), d(HOMO_exp), d(LUMO_exp).
    (b) The same question in the 3,086 DFT donor-acceptor molecules: matched pairs that differ by one S->Se or
        S->Te replacement (identical skeleton after mapping Se/Te -> S), split by ring type.
    (c) Helicene reference set (Kumar et al. 2025): does each DFT functional reproduce the S->Se and S->Te gap
        shift of the CCSD(T) IP-EA gap?
"""
import argparse, json, sys
from pathlib import Path
import numpy as np, pandas as pd
from scipy.stats import mannwhitneyu
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
R = {}
rng0 = np.random.default_rng(0)

PATS = {k: Chem.MolFromSmarts(v) for k, v in {"Se_diaz": "[#7][Se,se][#7]", "Se_phen": "[#6][Se,se][#6]", "S_diaz": "[#7][S,s][#7]",
        "S_phen": "[#6;a][s][#6;a]", "Te_any": "[Te,te]", "Te_diaz": "[#7][Te,te][#7]", "Te_phen": "[#6][Te,te][#6]", "O_furan": "[#6;a][o][#6;a]"}.items()}
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
    return D, np.array([list(f) for f in fps], float), fps

# ======================= S2
pol = df[df.Record_ID.str.startswith("POL_") & df.Canonical_SMILES.notna() & df.Eg_eV.notna() & df.Experimental_Eg_eV.notna()].reset_index(drop=True)
D, FP, fps = featurize(pol.Canonical_SMILES)
y = pol.Experimental_Eg_eV.to_numpy(); DFT = pol[["Eg_eV", "HOMO_eV", "LUMO_eV"]].to_numpy()
desc_cols = list(D.columns); fp_cols = [f"fp{i}" for i in range(FP.shape[1])]
Xs = np.hstack([D.to_numpy(float), FP]); Xsd = np.hstack([Xs, DFT]); names_sd = desc_cols + fp_cols + ["dft_gap", "dft_homo", "dft_lumo"]
n = len(fps); dists = []
for i in range(1, n):
    dists.extend([1 - x for x in DataStructs.BulkTanimotoSimilarity(fps[i], fps[:i])])
clusters = Butina.ClusterData(dists, n, 0.4, isDistData=True)
cl = np.zeros(n, int)
for k, c in enumerate(clusters):
    for i in c: cl[i] = k

def folds(seed):
    rng = np.random.default_rng(seed); ids = rng.permutation(len(clusters)); load = np.zeros(5); asg = {}
    for k in sorted(ids, key=lambda k: -len(clusters[k])):
        f = int(np.argmin(load)); asg[k] = f; load[f] += len(clusters[k])
    fo = np.array([asg[c] for c in cl]); return fo

P = dict(n_estimators=700, max_depth=5, learning_rate=0.03, subsample=0.8, colsample_bytree=0.5, min_child_weight=3, n_jobs=-1)
def oof(kind, seed):
    fo = folds(seed); pred = np.zeros(n)
    for f in range(5):
        tr, te = np.where(fo != f)[0], np.where(fo == f)[0]
        if kind == "mean": pred[te] = y[tr].mean(); continue
        if kind == "dft_linear": pred[te] = LinearRegression().fit(DFT[tr, :1], y[tr]).predict(DFT[te, :1]); continue
        X = Xs if kind == "struct_xgb" else Xsd
        if kind == "delta":
            cal = LinearRegression().fit(DFT[tr, :1], y[tr])
            m = xgb.XGBRegressor(random_state=seed, **P).fit(X[tr], y[tr] - cal.predict(DFT[tr, :1]))
            pred[te] = cal.predict(DFT[te, :1]) + m.predict(X[te])
        else:
            pred[te] = xgb.XGBRegressor(random_state=seed, **P).fit(X[tr], y[tr]).predict(X[te])
    return pred, fo

s2 = {"n": int(n), "n_clusters": len(clusters), "largest_cluster": int(pd.Series(cl).value_counts().max()), "models": {}}
keep_pred = {}
for kind in ["mean", "dft_linear", "struct_xgb", "struct+dft_xgb", "delta"]:
    maes, r2s, rmses = [], [], []
    for seed in (0, 1, 2):
        p, fo = oof(kind, seed); maes.append(mean_absolute_error(y, p)); r2s.append(r2_score(y, p)); rmses.append(float(np.sqrt(np.mean((y - p) ** 2))))
        if seed == 0: keep_pred[kind] = (p, fo)
    s2["models"][kind] = {"MAE": round(float(np.mean(maes)), 4), "MAE_sd": round(float(np.std(maes)), 4), "RMSE": round(float(np.mean(rmses)), 4), "R2": round(float(np.mean(r2s)), 3)}
    print("S2", kind, s2["models"][kind], flush=True)

p, fo = keep_pred["delta"]; res = np.abs(y - p)
# conformal: for fold f use |residuals| of the other folds as calibration
lo, hi = np.zeros(n), np.zeros(n)
for f in range(5):
    q = np.quantile(res[fo != f], 0.9 * (1 + 1 / (fo != f).sum())); te = fo == f; lo[te] = p[te] - q; hi[te] = p[te] + q
s2["conformal90"] = {"coverage": round(float(np.mean((y >= lo) & (y <= hi))), 3), "half_width_eV": round(float(np.mean((hi - lo) / 2)), 3)}
# applicability domain
maxsim = np.zeros(n)
for f in range(5):
    tr = np.where(fo != f)[0]; trf = [fps[i] for i in tr]
    for i in np.where(fo == f)[0]: maxsim[i] = max(DataStructs.BulkTanimotoSimilarity(fps[i], trf))
bins = [0, 0.4, 0.5, 0.6, 0.7, 1.01]; s2["applicability"] = []
for b0, b1 in zip(bins[:-1], bins[1:]):
    m = (maxsim >= b0) & (maxsim < b1)
    if m.sum(): s2["applicability"].append({"maxsim": f"{b0}-{min(b1,1.0)}", "n": int(m.sum()), "MAE": round(float(res[m].mean()), 4),
                                            "coverage90": round(float(np.mean((y[m] >= lo[m]) & (y[m] <= hi[m]))), 3)})
s2["corr_maxsim_abs_error"] = round(float(np.corrcoef(maxsim, res)[0, 1]), 3)
has = {"Se": pol.Canonical_SMILES.str.contains(r"\[se\]|\[Se\]").to_numpy(), "Te": pol.Canonical_SMILES.str.contains(r"\[te\]|\[Te\]").to_numpy()}
s2["by_chalcogen"] = {"S/O only": {"n": int((~has["Se"] & ~has["Te"]).sum()), "MAE": round(float(res[~has["Se"] & ~has["Te"]].mean()), 4)}}
for k, m in has.items():
    if m.sum(): s2["by_chalcogen"][k] = {"n": int(m.sum()), "MAE": round(float(res[m].mean()), 4)}
print("S2 conformal/AD", s2["conformal90"], s2["applicability"], s2["by_chalcogen"], flush=True)
# SHAP for the delta residual model on all data
try:
    import shap
    cal = LinearRegression().fit(DFT[:, :1], y); s2["dft_calibration"] = {"slope": round(float(cal.coef_[0]), 3), "intercept": round(float(cal.intercept_), 3)}
    m = xgb.XGBRegressor(random_state=0, **P).fit(Xsd, y - cal.predict(DFT[:, :1]))
    sv = np.abs(shap.TreeExplainer(m).shap_values(Xsd)).mean(0)
    imp = pd.Series(sv, index=names_sd); fp_sum = float(imp[fp_cols].sum())
    top = imp.drop(fp_cols).sort_values(ascending=False).head(12)
    s2["shap_residual_top"] = {k: round(float(v), 4) for k, v in top.items()}; s2["shap_morgan_total"] = round(fp_sum, 4)
    s2["shap_share"] = {"descriptors": round(float(imp[desc_cols].sum() / imp.sum()), 3), "morgan": round(fp_sum / float(imp.sum()), 3),
                        "dft": round(float(imp[["dft_gap", "dft_homo", "dft_lumo"]].sum() / imp.sum()), 3)}
    print("S2 shap", s2["shap_residual_top"], s2["shap_share"], flush=True)
except Exception as e:
    s2["shap_error"] = repr(e)
pd.DataFrame({"Record_ID": pol.Record_ID, "y_exp": y, "dft_gap": DFT[:, 0], "pred_delta": p, "lo90": lo, "hi90": hi, "maxsim": maxsim, "fold": fo,
              "pred_dft_linear": keep_pred["dft_linear"][0], "pred_struct": keep_pred["struct_xgb"][0]}).to_csv(OUT / "s2_oof.csv", index=False)
R["S2"] = s2

# ======================= S3a experimental pairs by position
def boot_ci(x, B=5000):
    x = np.asarray(x); bs = rng0.choice(x, (B, len(x)), replace=True).mean(1)
    return [round(float(np.mean(x)), 4), round(float(np.percentile(bs, 2.5)), 4), round(float(np.percentile(bs, 97.5)), 4)]
pairs = json.load(open(ROOT / "data/build/polymer_chalcogen_analysis.json"))["se_s_swap_pairs"]["pairs"]
smi = dict(zip(df.Record_ID, df.Canonical_SMILES))
rows = []
for pr in pairs:
    mse, ms = Chem.MolFromSmiles(smi[pr["Se_record"]]), Chem.MolFromSmiles(smi[pr["S_record"]])
    dz = len(mse.GetSubstructMatches(PATS["Se_diaz"])) - len(ms.GetSubstructMatches(PATS["Se_diaz"]))
    dp = len(mse.GetSubstructMatches(PATS["Se_phen"])) - len(ms.GetSubstructMatches(PATS["Se_phen"]))
    posn = "acceptor (selenadiazole)" if dz > 0 and dp - dz <= 0 else ("donor (selenophene)" if dp > 0 and dz <= 0 else "both/other")
    rows.append({**pr, "position": posn, "n_Se_added": int(max(dz, 0) + max(dp - max(dz, 0), 0))})
PP = pd.DataFrame(rows); PP.to_csv(OUT / "s3_exp_pairs.csv", index=False)
s3a = {}
for posn, g in PP.groupby("position"):
    s3a[posn] = {"n": int(len(g))}
    for k in ["dEg_exp", "dEg_dft", "dHOMO_exp", "dLUMO_exp"]:
        v = pd.to_numeric(g[k], errors="coerce").dropna()
        if len(v) > 2: s3a[posn][k] = {"mean_CI95": boot_ci(v), "n": int(len(v)), "share_negative": round(float((v < 0).mean()), 3)}
A, Dn = PP[PP.position.str.startswith("acceptor")], PP[PP.position.str.startswith("donor")]
s3a["mann_whitney"] = {}
for k in ["dEg_exp", "dEg_dft", "dHOMO_exp", "dLUMO_exp"]:
    a_, d_ = pd.to_numeric(A[k], errors="coerce").dropna(), pd.to_numeric(Dn[k], errors="coerce").dropna()
    if len(a_) > 2 and len(d_) > 2:
        s3a["mann_whitney"][k] = {"U": float(mannwhitneyu(a_, d_).statistic), "p": float(mannwhitneyu(a_, d_).pvalue), "ratio_means_acc_over_don": round(float(a_.mean() / d_.mean()), 2)}
# per added Se atom (donor swaps can add more than one Se)
for posn, g in PP.groupby("position"):
    nn = g.n_Se_added.clip(lower=1); s3a[posn]["dEg_exp_per_Se_atom"] = boot_ci(g.dEg_exp / nn)
print("S3a", json.dumps(s3a), flush=True)
R["S3a"] = s3a

# ======================= S3b DFT molecules: matched single-atom swaps
dev = df[df.Record_ID.str.startswith(("DEV_", "EROL_", "EXT_", "HAKAN_")) & df.Canonical_SMILES.notna() & df.Eg_eV.notna() & (df.Unit_Type.astype(str) != "Polymer")].reset_index(drop=True)
def skeleton(m):
    rw = Chem.RWMol(m)
    for at in rw.GetAtoms():
        if at.GetAtomicNum() in (34, 52): at.SetAtomicNum(16)
    try: return Chem.MolToSmiles(rw)
    except Exception: return None
info = []
for r in dev.itertuples():
    m = Chem.MolFromSmiles(r.Canonical_SMILES)
    if m is None: continue
    info.append({"rid": r.Record_ID, "coll": "DEV" if r.Record_ID.startswith(("DEV_", "EROL_")) else "EXT", "Eg": r.Eg_eV, "HOMO": r.HOMO_eV, "LUMO": r.LUMO_eV,
                 "skel": skeleton(m), "nSe": sum(at.GetAtomicNum() == 34 for at in m.GetAtoms()), "nTe": sum(at.GetAtomicNum() == 52 for at in m.GetAtoms()),
                 "SeD": len(m.GetSubstructMatches(PATS["Se_diaz"])), "TeD": len(m.GetSubstructMatches(PATS["Te_diaz"])),
                 "SeP": len(m.GetSubstructMatches(PATS["Se_phen"])), "TeP": len(m.GetSubstructMatches(PATS["Te_phen"]))})
I = pd.DataFrame(info).dropna(subset=["skel"])
mp = []
for sk, g in I.groupby("skel"):
    if len(g) < 2: continue
    base = g[(g.nSe == 0) & (g.nTe == 0)]
    for b in base.itertuples():
        for o in g[(g.nSe + g.nTe) == 1].itertuples():
            el = "Se" if o.nSe else "Te"
            ring = "diazole (acceptor)" if (o.SeD + o.TeD) > 0 else ("chalcogenophene (donor-type)" if (o.SeP + o.TeP) > 0 else "other")
            mp.append({"coll": o.coll, "element": el, "ring": ring, "dEg": o.Eg - b.Eg, "dHOMO": o.HOMO - b.HOMO, "dLUMO": o.LUMO - b.LUMO, "S_rid": b.rid, "X_rid": o.rid})
MP = pd.DataFrame(mp); MP.to_csv(OUT / "s3_dft_matched_pairs.csv", index=False)
s3b = {"n_pairs": int(len(MP))}
if len(MP):
    for (el, ring), g in MP.groupby(["element", "ring"]):
        if len(g) >= 3:
            s3b[f"{el}|{ring}"] = {"n": int(len(g)), "dEg": boot_ci(g.dEg), "dHOMO": boot_ci(g.dHOMO), "dLUMO": boot_ci(g.dLUMO)}
    for el in ["Se", "Te"]:
        a_, d_ = MP[(MP.element == el) & MP.ring.str.startswith("diazole")].dEg, MP[(MP.element == el) & MP.ring.str.startswith("chalcogenophene")].dEg
        if len(a_) > 2 and len(d_) > 2:
            s3b[f"mann_whitney_{el}"] = {"p": float(mannwhitneyu(a_, d_).pvalue), "ratio_means_acc_over_don": round(float(a_.mean() / d_.mean()), 2)}
print("S3b", json.dumps(s3b), flush=True)
R["S3b"] = s3b

# ======================= S3c helicenes: functional vs CCSD(T) chalcogen shift
H = ROOT / "data/sources/kumar_2025_helicenes"
orb = pd.read_csv(H / "helicene_orbital_energies_13_functionals.csv"); orb["gap"] = orb.LUMO_eV - orb.HOMO_eV
cc = pd.read_csv(H / "helicene_ccsdt_tables_s5_s8_s11.csv").set_index(["chalcogen", "n"]).IP_minus_EA_eV
s3c = {}
for (fn, st), g in orb.groupby(["functional", "set"]):
    gg = g.set_index(["chalcogen", "n"]).gap; out = {}
    for el in ["Se", "Te"]:
        e = []
        for nn in range(1, 5):
            if (el, nn) in gg.index and ("S", nn) in gg.index and (el, nn) in cc.index and ("S", nn) in cc.index:
                e.append((gg[(el, nn)] - gg[("S", nn)]) - (cc[(el, nn)] - cc[("S", nn)]))
        if e: out[f"shift_err_{el}_mean"] = round(float(np.mean(e)), 3); out[f"shift_err_{el}_MAE"] = round(float(np.mean(np.abs(e))), 3)
    s3c[f"{fn}|{st}"] = out
ccs = {el: [round(float(cc[(el, nn)] - cc[("S", nn)]), 3) for nn in range(1, 5) if (el, nn) in cc.index] for el in ["Se", "Te"]}
s3c["_ccsdt_shift_vs_S_n1to4"] = ccs
print("S3c", json.dumps(s3c), flush=True)
R["S3c"] = s3c
json.dump(R, open(OUT / "s23_results.json", "w"), indent=1, default=float)
