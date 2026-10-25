"""Revision run 1: (a) the 67 descriptors recomputed from cheap RDKit ETKDG + MMFF94 geometries of the D-A molecules
(no DFT input); (b) 2D baselines (RDKit 2D descriptors, Morgan count fingerprint); (c) the 67 descriptors from the
B3LYP-D3/6-31G* oligomer geometries of the polymer collection; (d) DFT vs MMFF bond-length agreement."""
import argparse, gzip, json, sys, time
from multiprocessing import Pool
from pathlib import Path
import numpy as np, pandas as pd
from rdkit import Chem, RDLogger
from rdkit.Chem import AllChem, Descriptors, rdFingerprintGenerator
RDLogger.DisableLog("rdApp.*")
ap = argparse.ArgumentParser(); ap.add_argument("--root", required=True); ap.add_argument("--out", required=True); a = ap.parse_args()
ROOT, OUT = Path(a.root), Path(a.out)
sys.path.insert(0, str(ROOT)); sys.path.insert(0, str(ROOT / "benchmarks")); sys.path.insert(0, str(Path(__file__).parent))
import grouped_cv_donor_acceptor as G
from thesis_feats import feats
dev, _ = G.load_development_collection(ROOT / "data/chalcogen_database_v1_master.csv.gz")
pub = dev["Record_ID"].astype(str).str.replace("EROL_", "DEV_", regex=False).tolist(); smis = dev["Canonical_SMILES"].tolist()
t0 = time.time()
def mmff(arg):
    rid, smi = arg
    try:
        m = Chem.AddHs(Chem.MolFromSmiles(smi)); ps = AllChem.ETKDGv3(); ps.randomSeed = 7
        if AllChem.EmbedMolecule(m, ps) != 0:
            ps.useRandomCoords = True
            if AllChem.EmbedMolecule(m, ps) != 0: return rid, None, "embed_failed"
        st = AllChem.MMFFOptimizeMolecule(m, maxIters=2000) if AllChem.MMFFHasAllMoleculeParams(m) else AllChem.UFFOptimizeMolecule(m, maxIters=2000)
        ff = "MMFF94" if AllChem.MMFFHasAllMoleculeParams(m) else "UFF"
        Chem.Kekulize(m, clearAromaticFlags=True)
        return rid, feats(m, smi), ff
    except Exception as e:
        return rid, None, repr(e)[:80]
with Pool(4) as pool: res = pool.map(mmff, list(zip(pub, smis)), chunksize=8)
ok = {r: f for r, f, s in res if f is not None}; status = pd.Series([s for _, _, s in res]).value_counts().to_dict()
print("MMFF done", len(ok), status, round(time.time() - t0), "s", flush=True)
F = pd.DataFrame([ok.get(r, {}) for r in pub]); F.insert(0, "Record_ID", pub); F.to_csv(OUT / "rev_desc_mmff.csv.gz", index=False)

# DFT vs MMFF bond-length agreement on the shared columns
Fd = pd.read_csv(OUT / "s1_descriptors.csv.gz")
cmp = {}
for c in ["bondlen_mean_N-S", "bondlen_mean_C-S", "bondlen_mean_C-Se", "bondlen_mean_N-Se", "bondlen_mean_C-N", "bondlen_BLA_CC", "planar_rms"]:
    if c in F and c in Fd:
        x, y = Fd[c].to_numpy(float), F[c].to_numpy(float); k = ~(np.isnan(x) | np.isnan(y))
        if k.sum() > 10: cmp[c] = {"pearson": round(float(np.corrcoef(x[k], y[k])[0, 1]), 3), "mean_abs_diff": round(float(np.abs(x[k] - y[k]).mean()), 4), "n": int(k.sum())}
print("DFT vs MMFF", json.dumps(cmp), flush=True)

# 2D baselines
names = [n for n, _ in Descriptors.descList]
gen = rdFingerprintGenerator.GetMorganGenerator(radius=2, fpSize=2048)
D2, FPc = [], []
for smi in smis:
    m = Chem.MolFromSmiles(smi)
    row = []
    for n, f in Descriptors.descList:
        try: row.append(float(f(m)))
        except Exception: row.append(np.nan)
    D2.append(row); FPc.append(gen.GetCountFingerprintAsNumPy(m))
pd.DataFrame(D2, columns=names).assign(Record_ID=pub).to_csv(OUT / "rev_desc_rdkit2d.csv.gz", index=False)
pd.DataFrame(np.array(FPc, dtype=np.uint16)).assign(Record_ID=pub).to_csv(OUT / "rev_morgan_counts.csv.gz", index=False)
print("2D done", round(time.time() - t0), "s", flush=True)

# polymer oligomer 3D descriptors (B3LYP-D3/6-31G* geometries)
df = pd.read_csv(ROOT / "data/chalcogen_database_v1_master.csv.gz", low_memory=False); smi_of = dict(zip(df.Record_ID, df.Canonical_SMILES))
rows = []
with gzip.open(ROOT / "data/structures_3d_polymer.sdf.gz", "rb") as fh:
    for m in Chem.ForwardSDMolSupplier(fh, removeHs=False, sanitize=False):
        if m is None: continue
        rid = m.GetProp("ChalMolDB_Record_ID")
        try: rows.append({"Record_ID": rid, **feats(m, smi_of[rid])})
        except Exception as e: rows.append({"Record_ID": rid})
pd.DataFrame(rows).to_csv(OUT / "rev_desc_polymer3d.csv.gz", index=False)
json.dump({"mmff_status": status, "dft_vs_mmff": cmp, "n_polymer3d": len(rows)}, open(OUT / "rev_geom.json", "w"), indent=1)
print("polymer 3D", len(rows), round(time.time() - t0), "s", flush=True)
