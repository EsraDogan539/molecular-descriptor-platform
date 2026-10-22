"""Experimental vs DFT gaps in the polymer collection, and S -> Se swap pairs.

Writes data/build/polymer_chalcogen_analysis.json. Requires RDKit.
A swap pair is two polymer records whose oligomer structures become identical when every Se atom of one
record is replaced by S (same backbone, side chains truncated as in the calculation).
"""

import json
from pathlib import Path

import numpy as np
import pandas as pd
from rdkit import Chem, RDLogger

RDLogger.DisableLog("rdApp.*")
ROOT = Path(__file__).resolve().parents[1]
DB = ROOT / "data" / "chalcogen_database_v1_master.csv.gz"
OUT = ROOT / "data" / "build" / "polymer_chalcogen_analysis.json"
NUM = ["Eg_eV", "HOMO_eV", "LUMO_eV", "TDDFT_S1_eV", "Experimental_Eg_eV", "Experimental_HOMO_eV",
       "Experimental_LUMO_eV"]


def fit(x, y):
    ok = x.notna() & y.notna()
    x, y = x[ok].to_numpy(float), y[ok].to_numpy(float)
    a, b = np.polyfit(x, y, 1)
    pred = a * x + b
    return {"n": int(ok.sum()), "pearson_r": round(float(np.corrcoef(x, y)[0, 1]), 3),
            "slope": round(float(a), 3), "intercept": round(float(b), 3),
            "mae_after_linear_fit": round(float(np.abs(pred - y).mean()), 3),
            "mean_offset_dft_minus_exp": round(float((x - y).mean()), 3)}


def se_to_s_key(smiles):
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return None
    rw = Chem.RWMol(mol)
    changed = False
    for atom in rw.GetAtoms():
        if atom.GetSymbol() == "Se":
            atom.SetAtomicNum(16)
            changed = True
    if not changed:
        return None
    m = rw.GetMol()
    Chem.SanitizeMol(m)
    return Chem.MolToInchiKey(m)


def main():
    df = pd.read_csv(DB, low_memory=False)
    pol = df[df.Record_ID.str.startswith("POL_")].copy()
    for c in NUM:
        pol[c] = pd.to_numeric(pol[c], errors="coerce")
    pol["has_Se"] = pd.to_numeric(pol.Se_Count, errors="coerce").fillna(0) > 0
    res = {
        "records": int(len(pol)),
        "dft_gap_vs_exp_optical_gap": fit(pol.Eg_eV, pol.Experimental_Eg_eV),
        "tddft_s1_vs_exp_optical_gap": fit(pol.TDDFT_S1_eV, pol.Experimental_Eg_eV),
        "dft_homo_vs_exp_homo": fit(pol.HOMO_eV, pol.Experimental_HOMO_eV),
        "dft_lumo_vs_exp_lumo": fit(pol.LUMO_eV, pol.Experimental_LUMO_eV),
        "by_chalcogen": {},
    }
    for label, mask in (("S only", pol.Chalcogen_Type.eq("S")), ("Se-containing", pol.has_Se)):
        sub = pol[mask]
        res["by_chalcogen"][label] = {
            "n": int(len(sub)),
            "mean_exp_gap": round(float(sub.Experimental_Eg_eV.mean()), 3),
            "mean_dft_gap": round(float(sub.Eg_eV.mean()), 3),
            "dft_gap_vs_exp": fit(sub.Eg_eV, sub.Experimental_Eg_eV),
        }
    by_key = pol.groupby("InChIKey")
    pairs = []
    for _, r in pol[pol.has_Se].iterrows():
        key = se_to_s_key(r.Canonical_SMILES)
        if key is None or key not in by_key.groups:
            continue
        for _, s in pol.loc[by_key.groups[key]].iterrows():
            pairs.append({
                "Se_record": r.Record_ID, "S_record": s.Record_ID,
                "Se_name": r.Molecule_Name, "S_name": s.Molecule_Name,
                "dEg_exp": round(r.Experimental_Eg_eV - s.Experimental_Eg_eV, 3),
                "dEg_dft": round(r.Eg_eV - s.Eg_eV, 3),
                "dHOMO_exp": (None if pd.isna(r.Experimental_HOMO_eV) or pd.isna(s.Experimental_HOMO_eV)
                              else round(r.Experimental_HOMO_eV - s.Experimental_HOMO_eV, 3)),
                "dLUMO_exp": (None if pd.isna(r.Experimental_LUMO_eV) or pd.isna(s.Experimental_LUMO_eV)
                              else round(r.Experimental_LUMO_eV - s.Experimental_LUMO_eV, 3)),
            })
    pdf = pd.DataFrame(pairs)
    res["se_s_swap_pairs"] = {
        "n_pairs": int(len(pdf)),
        "n_se_records_with_partner": int(pdf.Se_record.nunique()) if len(pdf) else 0,
        "mean_dEg_exp": round(float(pdf.dEg_exp.mean()), 3) if len(pdf) else None,
        "median_dEg_exp": round(float(pdf.dEg_exp.median()), 3) if len(pdf) else None,
        "share_exp_gap_smaller_with_Se": round(float((pdf.dEg_exp < 0).mean()), 3) if len(pdf) else None,
        "mean_dEg_dft": round(float(pdf.dEg_dft.mean()), 3) if len(pdf) else None,
        "share_dft_gap_smaller_with_Se": round(float((pdf.dEg_dft < 0).mean()), 3) if len(pdf) else None,
        "mean_dHOMO_exp": round(float(pdf.dHOMO_exp.dropna().mean()), 3) if len(pdf) else None,
        "mean_dLUMO_exp": round(float(pdf.dLUMO_exp.dropna().mean()), 3) if len(pdf) else None,
        "pairs": pairs,
    }
    OUT.parent.mkdir(parents=True, exist_ok=True)
    OUT.write_text(json.dumps(res, indent=1))
    print(json.dumps({k: v for k, v in res.items() if k != "se_s_swap_pairs"}, indent=1))
    print({k: v for k, v in res["se_s_swap_pairs"].items() if k != "pairs"})


if __name__ == "__main__":
    main()
