"""3D coordinates for the demo: helicenes (RDKit ETKDG + UFF) and two polymer oligomers from the release SDF."""
import json, gzip
from pathlib import Path
import pandas as pd
from rdkit import Chem
from rdkit.Chem import AllChem
OUT = Path("docs/demo2"); OUT.mkdir(parents=True, exist_ok=True)
df = pd.read_csv("data/chalcogen_database_v1_master.csv.gz", low_memory=False).set_index("Record_ID")
res = {}
def pack(m):
    c = m.GetConformer(); a = []
    for at in m.GetAtoms():
        p = c.GetAtomPosition(at.GetIdx()); a.append([round(p.x, 3), round(p.y, 3), round(p.z, 3), at.GetSymbol()])
    b = [[x.GetBeginAtomIdx(), x.GetEndAtomIdx(), 2 if x.GetBondType() == Chem.BondType.DOUBLE else (1.5 if x.GetIsAromatic() else 1)] for x in m.GetBonds()]
    return {"a": a, "b": b}
for rid in ["LIT_0206", "LIT_0210", "LIT_0214", "LIT_0205"]:
    try:
        smi = df.loc[rid, "Canonical_SMILES"]; m = Chem.AddHs(Chem.MolFromSmiles(smi))
        AllChem.EmbedMolecule(m, randomSeed=7); AllChem.UFFOptimizeMolecule(m, maxIters=2000)
        res[rid] = pack(m); res[rid]["name"] = df.loc[rid, "Molecule_Name"]; print(rid, "ok", smi)
    except Exception as e:
        print(rid, "fail", e)
sup = Chem.ForwardSDMolSupplier(gzip.open("data/structures_3d_polymer.sdf.gz"), removeHs=False)
for m in sup:
    if m is None: continue
    name = m.GetProp("_Name").strip()
    if name in ("POL_0079", "POL_0080"):
        res[name] = pack(m); res[name]["name"] = df.loc[name, "Molecule_Name"]; print(name, "ok", m.GetNumAtoms())
(OUT / "mols.json").write_text(json.dumps(res))
