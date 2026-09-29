"""Grouped cross-validation of Eg models on the ChalMolDB development collection.

The development collection is a combinatorial donor x acceptor space, so a
split that only keeps identical structures (InChIKey) together still lets every
test donor and acceptor appear in training. This script compares:

  inchikey        5-fold GroupKFold by standardized structure (manuscript setting)
  donor           5-fold leave-donor-group-out
  acceptor        5-fold leave-acceptor-group-out
  donor+acceptor  5-fold "cold-cold": test pairs whose donor AND acceptor are unseen

Donors/acceptors that are chemically identical in the source library (they
produce the same InChIKey) are merged into one group before splitting, so an
alias of a test donor cannot leak into training.

Representations: general descriptors, chalcogen-aware descriptors, combined.
Models: random forest, SVR, XGBoost, plus two reference baselines:
  mean      - training-set mean (no information)
  fragments - ridge on one-hot donor + acceptor identity (memorisation baseline)

Usage:
    python benchmarks/grouped_cv_donor_acceptor.py --out benchmarks/results
"""

import argparse
import re
import sys
import time
from pathlib import Path

import numpy as np
import pandas as pd
from sklearn.dummy import DummyRegressor
from sklearn.ensemble import RandomForestRegressor
from sklearn.linear_model import Ridge
from sklearn.metrics import mean_absolute_error, mean_squared_error, r2_score
from sklearn.model_selection import GroupKFold
from sklearn.pipeline import make_pipeline
from sklearn.preprocessing import OneHotEncoder, StandardScaler
from sklearn.svm import SVR

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT))

DATABASE_PATH = ROOT / "data" / "chalcogen_database_v1_master.csv.gz"
TARGET = "Eg_eV"
N_FOLDS = 5
SEED = 0

GENERAL_COLUMNS = [
    "Molecular Weight", "Exact Molecular Weight", "LogP", "TPSA",
    "H-Bond Donors", "H-Bond Acceptors", "Rotatable Bonds", "Ring Count",
    "Aromatic Ring Count", "Aliphatic Ring Count", "Heavy Atom Count",
    "Heteroatom Count", "Fraction Csp3", "Formal Charge", "Bertz CT",
    "Balaban J", "Carbon Count", "Nitrogen Count", "Oxygen Count",
    "Phosphorus Count", "Halogen Count",
]
CHALCOGEN_COLUMNS = [
    "Sulfur Count", "Selenium Count", "Tellurium Count",
    "Target Chalcogen Count", "Target Chalcogen Fraction",
    "Contains S", "Contains Se", "Contains Te",
    "Aromatic Chalcogen Count", "NonAromatic Chalcogen Count",
    "Mixed Chalcogen Flag", "Ring Incorporated Chalcogen Count",
    "Chalcogen-C Bond Count", "Chalcogen-Heteroatom Bond Count",
    "Aromatic Neighbor Count", "Conjugated Bond Count",
    "Aromatic Bond Fraction", "Conjugated Atom Fraction",
    "Heteroaromatic Ring Count",
]
REPRESENTATIONS = {
    "General": GENERAL_COLUMNS,
    "Chalcogen-aware": CHALCOGEN_COLUMNS,
    "Combined": GENERAL_COLUMNS + CHALCOGEN_COLUMNS,
}

NAME_PATTERN = re.compile(r"(D\d+)\W*(A\d+)")


# --------------------------------------------------------------------------- data

def load_development_collection(path=DATABASE_PATH):
    df = pd.read_csv(path, low_memory=False)
    dev = df[df["Split_Role"] == "Development/Training"].copy()
    dev = dev[dev[TARGET].notna() & dev["Canonical_SMILES"].notna()]

    parsed = dev["Molecule_Name"].astype(str).str.extract(NAME_PATTERN)
    dev["Donor"] = parsed[0]
    dev["Acceptor"] = parsed[1]
    # A few source names lack the acceptor label (e.g. "D21"); they cannot be
    # assigned to an acceptor group, so they are excluded rather than guessed.
    missing = dev["Donor"].isna() | dev["Acceptor"].isna()
    excluded = dev.loc[missing, ["Record_ID", "Molecule_Name"]].reset_index(drop=True)
    return dev[~missing].reset_index(drop=True), excluded


class _UnionFind:
    def __init__(self):
        self.parent = {}

    def find(self, item):
        self.parent.setdefault(item, item)
        while self.parent[item] != item:
            self.parent[item] = self.parent[self.parent[item]]
            item = self.parent[item]
        return item

    def union(self, a, b):
        ra, rb = self.find(a), self.find(b)
        if ra != rb:
            self.parent[max(ra, rb)] = min(ra, rb)


def fragment_groups(dev):
    """Merge donors (acceptors) that yield identical structures with the same partner."""
    donors, acceptors = _UnionFind(), _UnionFind()
    for _, block in dev.groupby("InChIKey"):
        if len(block) < 2:
            continue
        pairs = list(zip(block["Donor"], block["Acceptor"]))
        for d1, a1 in pairs:
            for d2, a2 in pairs:
                if a1 == a2:
                    donors.union(d1, d2)
                if d1 == d2:
                    acceptors.union(a1, a2)
    donor_group = dev["Donor"].map(donors.find)
    acceptor_group = dev["Acceptor"].map(acceptors.find)
    aliases = {
        "donor": sorted({(d, g) for d, g in zip(dev["Donor"], donor_group) if d != g}),
        "acceptor": sorted({(a, g) for a, g in zip(dev["Acceptor"], acceptor_group) if a != g}),
    }
    return donor_group, acceptor_group, aliases


def compute_descriptors(smiles):
    """Descriptor rows from the ChalMolDB descriptor engine (requires RDKit)."""
    from descriptor_engine import calculate_single_molecule_descriptors

    rows = []
    for index, value in enumerate(smiles):
        result, error = calculate_single_molecule_descriptors(f"M{index}", value)
        if result is None:
            raise ValueError(f"Descriptor calculation failed for {value}: {error}")
        rows.append(result)
    return pd.DataFrame(rows)


# --------------------------------------------------------------------------- splits

def _balanced_folds(labels, n_folds, rng):
    unique = np.array(sorted(set(labels)))
    rng.shuffle(unique)
    return {label: i % n_folds for i, label in enumerate(unique)}


def split_indices(dev, scheme, donor_group, acceptor_group, n_folds=N_FOLDS, seed=SEED):
    """Yield (train_idx, test_idx) for a split scheme."""
    n = len(dev)
    if scheme in {"inchikey", "donor", "acceptor"}:
        groups = {
            "inchikey": dev["InChIKey"].fillna(dev["Canonical_SMILES"]),
            "donor": donor_group,
            "acceptor": acceptor_group,
        }[scheme]
        rng = np.random.default_rng(seed)
        fold_of = _balanced_folds(groups, n_folds, rng)
        fold = groups.map(fold_of).to_numpy()
        for k in range(n_folds):
            yield np.flatnonzero(fold != k), np.flatnonzero(fold == k)
        return

    if scheme == "donor+acceptor":
        rng = np.random.default_rng(seed)
        d_fold = donor_group.map(_balanced_folds(donor_group, n_folds, rng)).to_numpy()
        a_fold = acceptor_group.map(_balanced_folds(acceptor_group, n_folds, rng)).to_numpy()
        for k in range(n_folds):
            test = np.flatnonzero((d_fold == k) & (a_fold == k))
            train = np.flatnonzero((d_fold != k) & (a_fold != k))
            if len(test):
                yield train, test
        return

    raise ValueError(f"Unknown split scheme: {scheme}")


# --------------------------------------------------------------------------- models

def make_models():
    models = {
        "RF": lambda: RandomForestRegressor(
            n_estimators=500, min_samples_leaf=1, n_jobs=-1, random_state=SEED
        ),
        "SVR": lambda: make_pipeline(
            StandardScaler(), SVR(C=10.0, epsilon=0.05, gamma="scale")
        ),
    }
    try:
        from xgboost import XGBRegressor

        models["XGB"] = lambda: XGBRegressor(
            n_estimators=600, learning_rate=0.05, max_depth=6, subsample=0.8,
            colsample_bytree=0.8, random_state=SEED, n_jobs=-1,
        )
    except ImportError:
        print("xgboost not installed; XGB skipped", file=sys.stderr)
    return models


def _metrics(y_true, y_pred):
    return {
        "MAE": mean_absolute_error(y_true, y_pred),
        "RMSE": float(np.sqrt(mean_squared_error(y_true, y_pred))),
        "R2": r2_score(y_true, y_pred),
    }


def evaluate(dev, features, donor_group, acceptor_group, schemes):
    y = dev[TARGET].to_numpy(dtype=float)
    fragments = dev[["Donor", "Acceptor"]].astype(str)
    models = make_models()
    records = []

    for scheme in schemes:
        folds = list(split_indices(dev, scheme, donor_group, acceptor_group))
        for fold_id, (train, test) in enumerate(folds):
            base = {"Split": scheme, "Fold": fold_id, "n_train": len(train), "n_test": len(test)}

            dummy = DummyRegressor().fit(np.zeros((len(train), 1)), y[train])
            records.append({**base, "Representation": "-", "Model": "Mean baseline",
                            **_metrics(y[test], dummy.predict(np.zeros((len(test), 1))))})

            frag = make_pipeline(OneHotEncoder(handle_unknown="ignore"), Ridge(alpha=1.0))
            frag.fit(fragments.iloc[train], y[train])
            records.append({**base, "Representation": "Donor/acceptor one-hot",
                            "Model": "Ridge", **_metrics(y[test], frag.predict(fragments.iloc[test]))})

            for rep_name, columns in REPRESENTATIONS.items():
                X = features[columns].to_numpy(dtype=float)
                for model_name, factory in models.items():
                    model = factory().fit(X[train], y[train])
                    records.append({**base, "Representation": rep_name, "Model": model_name,
                                    **_metrics(y[test], model.predict(X[test]))})
            print(f"{scheme} fold {fold_id}: train={len(train)} test={len(test)}", flush=True)

    return pd.DataFrame(records)


def summarise(folds):
    grouped = folds.groupby(["Split", "Representation", "Model"], sort=False)
    summary = grouped[["MAE", "RMSE", "R2"]].agg(["mean", "std"])
    summary.columns = [f"{m}_{s}" for m, s in summary.columns]
    summary["n_folds"] = grouped.size()
    summary["n_test_total"] = grouped["n_test"].sum()
    return summary.reset_index()


def summary_markdown(summary, dev, aliases, elapsed, excluded=None):
    lines = [
        "# Grouped cross-validation: donor/acceptor-aware splits",
        "",
        f"Development records: {len(dev)} · unique InChIKeys: {dev['InChIKey'].nunique()} · "
        f"donors: {dev['Donor'].nunique()} · acceptors: {dev['Acceptor'].nunique()}",
        "",
        "Merged identical fragments (alias -> group): "
        f"donors {aliases['donor'] or 'none'}; acceptors {aliases['acceptor'] or 'none'}",
        "",
        "Excluded (no acceptor label in source name): "
        + (", ".join(f"{r.Record_ID} ({r.Molecule_Name})" for r in excluded.itertuples())
           if excluded is not None and len(excluded) else "none"),
        "",
        "Values are mean ± SD over folds. Target: Eg_eV (calculated HOMO-LUMO gap).",
        "Hyperparameters are fixed (no tuning); see make_models().",
        "",
    ]
    for scheme, block in summary.groupby("Split", sort=False):
        lines += [f"## Split: {scheme}", "",
                  "| Representation | Model | MAE (eV) | RMSE (eV) | R² | test records |",
                  "|---|---|---|---|---|---|"]
        for _, r in block.iterrows():
            lines.append(
                f"| {r['Representation']} | {r['Model']} | "
                f"{r['MAE_mean']:.3f} ± {r['MAE_std']:.3f} | "
                f"{r['RMSE_mean']:.3f} ± {r['RMSE_std']:.3f} | "
                f"{r['R2_mean']:.3f} ± {r['R2_std']:.3f} | {int(r['n_test_total'])} |"
            )
        lines.append("")
    lines.append(f"Runtime: {elapsed:.0f} s")
    return "\n".join(lines) + "\n"


def main(argv=None):
    parser = argparse.ArgumentParser(description=__doc__.splitlines()[0])
    parser.add_argument("--out", default=str(ROOT / "benchmarks" / "results"))
    parser.add_argument("--schemes", nargs="+",
                        default=["inchikey", "donor", "acceptor", "donor+acceptor"])
    args = parser.parse_args(argv)

    start = time.time()
    dev, excluded = load_development_collection()
    donor_group, acceptor_group, aliases = fragment_groups(dev)
    features = compute_descriptors(dev["Canonical_SMILES"].tolist())

    folds = evaluate(dev, features, donor_group, acceptor_group, args.schemes)
    summary = summarise(folds)

    out = Path(args.out)
    out.mkdir(parents=True, exist_ok=True)
    folds.to_csv(out / "grouped_cv_folds.csv", index=False)
    summary.to_csv(out / "grouped_cv_summary.csv", index=False)
    (out / "grouped_cv_summary.md").write_text(
        summary_markdown(summary, dev, aliases, time.time() - start, excluded), encoding="utf-8"
    )
    print((out / "grouped_cv_summary.md").read_text(encoding="utf-8"))


if __name__ == "__main__":
    main()
