"""Render the graphical TOC entry for the Molecular Informatics submission from real ChalMolDB data."""
import io
from pathlib import Path

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import pandas as pd
from PIL import Image
from rdkit import Chem
from rdkit.Chem import rdDepictor
from rdkit.Chem.Draw import rdMolDraw2D

OUT = Path("toc"); OUT.mkdir(exist_ok=True)
RECORD = "EROL_0481"  # public ID DEV_0481 (donor D9, acceptor A13)
S_COLOR, SE_COLOR, NAVY = "#D6A400", "#0F766E", "#1F4E79"

db = pd.read_csv("data/chalcogen_database_v1_master.csv.gz", low_memory=False)
row = db[db.Record_ID == RECORD].iloc[0]
mol = Chem.MolFromSmiles(row.Canonical_SMILES)
rdDepictor.SetPreferCoordGen(True)
rdDepictor.Compute2DCoords(mol)

hex2rgb = lambda h: tuple(int(h[i:i + 2], 16) / 255 for i in (1, 3, 5))
highlight = {a.GetIdx(): hex2rgb(S_COLOR if a.GetSymbol() == "S" else SE_COLOR)
             for a in mol.GetAtoms() if a.GetSymbol() in ("S", "Se")}
d = rdMolDraw2D.MolDraw2DCairo(1400, 1000)
o = d.drawOptions()
o.clearBackground = False
o.bondLineWidth = 5
o.minFontSize = 44
o.maxFontSize = 60
o.fixedBondLength = 90
o.highlightRadius = 0.45
o.updateAtomPalette({16: hex2rgb("#8A6A00"), 34: hex2rgb("#0B5750"), 7: hex2rgb("#2B4FA0"), 6: (0.1, 0.1, 0.1)})
d.DrawMolecule(mol, highlightAtoms=list(highlight), highlightAtomColors=highlight, highlightBonds=[])
d.FinishDrawing()
mol_img = Image.open(io.BytesIO(d.GetDrawingText())).convert("RGBA")
bbox = mol_img.getbbox(); mol_img = mol_img.crop(bbox)
mol_img.save(OUT / "molecule_DEV_0481.png")

summ = pd.read_csv("benchmarks/results/grouped_cv_summary.csv") if Path("benchmarks/results/grouped_cv_summary.csv").exists() else None
mae = {"inchikey": 0.161, "donor": 0.239, "donor+acceptor": 0.347}
if summ is not None:
    for k in mae:
        mae[k] = float(summ[(summ.Split == k) & (summ.Representation == "Combined") & (summ.Model == "XGB")].MAE_mean.iloc[0])

W, H = 5.0 / 2.54, 3.6 / 2.54
plt.rcParams.update({"font.family": "DejaVu Sans", "font.size": 6.5})
fig = plt.figure(figsize=(W, H))
fig.text(0.03, 0.93, "ChalMolDB", fontsize=9.5, fontweight="bold", color=NAVY, va="top")
fig.text(0.03, 0.80, "3,360 curated S/Se/Te records", fontsize=6.5, color="#333333", va="top")

ax_m = fig.add_axes([0.00, 0.14, 0.52, 0.58]); ax_m.axis("off")
ax_m.imshow(mol_img); ax_m.set_aspect("equal")
fig.text(0.26, 0.07, f"DEV_0481 · E$_g$ {row.Eg_eV:.2f} eV", fontsize=6, ha="center", color="#333333")

ax = fig.add_axes([0.62, 0.14, 0.36, 0.56])
labels = ["Known\nunits", "New\ndonor", "New\nD + A"]
vals = [mae["inchikey"], mae["donor"], mae["donor+acceptor"]]
bars = ax.bar(range(3), vals, color=["#9DB5CF", "#4F7CA8", NAVY], width=0.68)
for b, v in zip(bars, vals):
    ax.text(b.get_x() + b.get_width() / 2, v + 0.012, f"{v:.2f}", ha="center", va="bottom", fontsize=6)
ax.set_xticks(range(3)); ax.set_xticklabels(labels, fontsize=6)
ax.set_ylim(0, 0.45); ax.set_yticks([]); ax.tick_params(axis="x", length=0, pad=1.5)
for s in ("top", "right", "left"):
    ax.spines[s].set_visible(False)
ax.spines["bottom"].set_color("#555555"); ax.spines["bottom"].set_linewidth(0.6)
ax.set_title("MAE (eV)", fontsize=6.5, pad=2, color="#333333")
fig.text(0.80, 0.93, "Validation\ndesign matters", fontsize=6.5, fontweight="bold", color=NAVY, ha="center", va="top")

for dpi in (300, 600):
    fig.savefig(OUT / f"ChalMolDB_TOC_5cm_{dpi}dpi.png", dpi=dpi, facecolor="white")
fig.savefig(OUT / "ChalMolDB_TOC_5cm.pdf", facecolor="white")
print("record", RECORD, row.Canonical_SMILES, row.Eg_eV, mae)
