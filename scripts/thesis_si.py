"""Thesis paper: SHAP figure data and Supporting Information tables for section 1.
(1) SHAP values of the XGBoost model on all 67 descriptors (all data) -> top features for a beeswarm plot.
(2) Eight-model comparison of the original draft (Linear, Ridge, Lasso, SVR, RF, Extra Trees, Gradient Boosting, MLP)
    under the five validation schemes.
(3) Low-gap region (Eg < 2.0 and < 1.5 eV): MAE and mean signed error of out-of-fold predictions per split.
"""
import argparse, json, sys, time
from pathlib import Path
import numpy as np, pandas as pd
from sklearn.ensemble import ExtraTreesRegressor, GradientBoostingRegressor, RandomForestRegressor
from sklearn.linear_model import Lasso, LinearRegression, Ridge
from sklearn.metrics import mean_absolute_error, r2_score
from sklearn.model_selection import KFold
from sklearn.neural_network import MLPRegressor
from sklearn.pipeline import make_pipeline
from sklearn.preprocessing import StandardScaler
from sklearn.svm import SVR
import xgboost as xgb, shap

ap = argparse.ArgumentParser(); ap.add_argument("--root", required=True); ap.add_argument("--out", required=True); a = ap.parse_args()
ROOT, OUT = Path(a.root), Path(a.out)
sys.path.insert(0, str(ROOT)); sys.path.insert(0, str(ROOT / "benchmarks"))
import grouped_cv_donor_acceptor as G
dev, _ = G.load_development_collection(ROOT / "data/chalcogen_database_v1_master.csv.gz")
dg, ag, _ = G.fragment_groups(dev)
F = pd.read_csv(OUT / "s1_descriptors.csv.gz"); y = F.pop("Eg_eV").to_numpy(); F.pop("Record_ID"); cols = list(F.columns); X = F.to_numpy(float)
t0 = time.time(); res = {}

# (1) SHAP
XGBP = dict(n_estimators=600, learning_rate=0.05, max_depth=6, subsample=0.8, colsample_bytree=0.6, random_state=0, n_jobs=-1)
m = xgb.XGBRegressor(**XGBP).fit(X, y)
sv = shap.TreeExplainer(m).shap_values(X)
imp = pd.Series(np.abs(sv).mean(0), index=cols).sort_values(ascending=False)
top = list(imp.index[:15])
pd.DataFrame(sv[:, [cols.index(c) for c in top]], columns=top).to_csv(OUT / "si_shap_values_top15.csv.gz", index=False)
F[top].to_csv(OUT / "si_shap_features_top15.csv.gz", index=False)
res["shap_importance_all"] = {k: round(float(v), 5) for k, v in imp.items()}
res["shap_direction_top15"] = {c: round(float(np.corrcoef(F[c], sv[:, cols.index(c)])[0, 1]), 3) for c in top}
print("shap done", round(time.time() - t0), "s", flush=True)

# (2) eight models
def sc(est): return make_pipeline(StandardScaler(), est)
MODELS = {"Linear Regression": lambda: sc(LinearRegression()), "Ridge": lambda: sc(Ridge(alpha=1.0)), "Lasso": lambda: sc(Lasso(alpha=0.001, max_iter=20000)),
          "SVR": lambda: sc(SVR(C=10.0, epsilon=0.05, gamma="scale")), "Random Forest": lambda: RandomForestRegressor(n_estimators=500, n_jobs=-1, random_state=0),
          "Extra Trees": lambda: ExtraTreesRegressor(n_estimators=500, n_jobs=-1, random_state=0), "Gradient Boosting": lambda: GradientBoostingRegressor(random_state=0),
          "MLP": lambda: sc(MLPRegressor(hidden_layer_sizes=(128, 64), alpha=1e-3, early_stopping=True, max_iter=1000, random_state=0)),
          "XGBoost": lambda: xgb.XGBRegressor(**XGBP)}
SCH = ["random", "inchikey", "donor", "acceptor", "donor+acceptor"]
def splits(s): return list(KFold(5, shuffle=True, random_state=0).split(X)) if s == "random" else list(G.split_indices(dev, s, dg, ag))
oof = {}
res["eight_models"] = {}
for name, fn in MODELS.items():
    res["eight_models"][name] = {}
    for s in SCH:
        p = np.full(len(y), np.nan); tr_r2 = []
        for tr, te in splits(s):
            mm = fn().fit(X[tr], y[tr]); p[te] = mm.predict(X[te]); tr_r2.append(r2_score(y[tr], mm.predict(X[tr])))
        k = ~np.isnan(p)
        res["eight_models"][name][s] = {"MAE": round(mean_absolute_error(y[k], p[k]), 4), "RMSE": round(float(np.sqrt(np.mean((y[k] - p[k]) ** 2))), 4),
                                        "R2": round(r2_score(y[k], p[k]), 3), "train_R2": round(float(np.mean(tr_r2)), 3)}
        if name in ("SVR", "XGBoost"): oof[f"{name}|{s}"] = p
    print(name, json.dumps(res["eight_models"][name]), round(time.time() - t0), "s", flush=True)
    json.dump(res, open(OUT / "si_results.json", "w"), indent=1)

# (3) low-gap region
res["low_gap"] = {}
for key, p in oof.items():
    k = ~np.isnan(p); r = {}
    for lab, mask in [("all", k), ("Eg<2.0", k & (y < 2.0)), ("Eg<1.5", k & (y < 1.5)), ("Eg>=3.5", k & (y >= 3.5))]:
        if mask.sum(): r[lab] = {"n": int(mask.sum()), "MAE": round(float(np.abs(y[mask] - p[mask]).mean()), 4), "mean_signed_error": round(float((p[mask] - y[mask]).mean()), 4)}
    res["low_gap"][key] = r
print("low gap", json.dumps(res["low_gap"]), flush=True)
json.dump(res, open(OUT / "si_results.json", "w"), indent=1)
