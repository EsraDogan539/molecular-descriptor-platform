"""Shared structure-derived descriptor families (same code as thesis_s1_da.py)."""
import numpy as np
from rdkit import Chem
from rdkit.Chem import Descriptors, rdMolDescriptors
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

