"""Build exact structures for external (D–A–D) records from the Kayı-group geometry files.

Input: Gaussian input files in data/sources/kayi_group_geometries/<Family>-DAD_<Monomer|Hexamer>/,
geometries optimized at B3LYP/LANL2DZ (see data/sources/kayi_group_geometries/README.md).

Bond orders cannot be read reliably from Cartesian coordinates, so each system is built from a
family template (known chemistry). The template's identity (SMILES/InChIKey) is accepted when its
full atom/bond graph, including hydrogens, equals either the graph perceived from the coordinates
or the Gaussian connectivity block written in the input file.

Coordinates are exported to the SDF only when, in addition, (i) the file's formula equals the
template formula and (ii) every chalcogen bond length lies in a plausible range for that element
at this level of theory (BOND_RANGES). Check (ii) catches files whose geometry was made by swapping
atom labels on another system's structure without re-optimization (bond lengths then do not change
with the element).

Outputs (data/build/):
  external_structures.csv   one row per geometry file, with SMILES/InChI/InChIKey and checks
  external_structures.sdf   3D structures (template bond orders + file coordinates)
"""

import re
import sys
from pathlib import Path

import pandas as pd
from rdkit import Chem
from rdkit.Chem import rdDetermineBonds
from rdkit.Chem.rdMolDescriptors import CalcMolFormula

ROOT = Path(__file__).resolve().parents[1]
SRC = ROOT / "data" / "sources" / "kayi_group_geometries"
OUT = ROOT / "data" / "build"

SYM = {"O": "o", "S": "s", "Se": "[se]", "Te": "[te]"}
KEKULE_FALLBACK = {"o": "O", "s": "S", "[se]": "[Se]", "[te]": "[Te]"}

# Repeat units with attachment dummies [*:L] (left) / [*:R] (right) on the donor alpha carbons.
TEMPLATES = {
    # 4,7-bis(chalcogenophen-2-yl)-2,1,3-benzochalcogenadiazole
    "A": "c1cc([*:L]){X}c1-c1ccc(-c2ccc([*:R]){X}2)c2n{Y}nc12",
    # 4,7-bis(3,4-ethylenedioxychalcogenophen-2-yl)-2,1,3-benzochalcogenadiazole
    "B": "c1({X}c([*:L])c2OCCOc12)-c1ccc(-c2{X}c([*:R])c3OCCOc23)c2n{Y}nc12",
    # 4,9-bis(chalcogenophen-2-yl)-[1,2,5]chalcogenadiazolo[3,4-g]quinoxaline
    "C": "c1cc([*:L]){X}c1-c1c2nccnc2c(-c2ccc([*:R]){X}2)c2n{Y}nc12",
}
# Plausible ring X-C / X-N bond lengths (Angstrom) at B3LYP/LANL2DZ: the range observed in the
# element-consistent geometries of this set (C family, A/B hexamers) widened by about 0.05 A.
# Only bonds of chalcogens in five-membered (chalcogenophene / chalcogenadiazole) rings are checked.
BOND_RANGES = {
    ("O", "C"): (1.35, 1.47), ("O", "N"): (1.37, 1.47),
    ("S", "C"): (1.76, 1.88), ("S", "N"): (1.70, 1.83),
    ("Se", "C"): (1.87, 2.01), ("Se", "N"): (1.80, 1.93),
    ("Te", "C"): (2.04, 2.19), ("Te", "N"): (1.95, 2.08),
}
FAMILY_LABEL = {
    "A": "benzochalcogenadiazole",
    "B": "B_system / benzochalcogenadiazole-like DAD",
    "C": "chalcogendiazoloquinoxaline",
}


def parse_gjf(path):
    lines = path.read_text(encoding="utf-8", errors="replace").splitlines()
    route = " ".join(l.strip() for l in lines if l.startswith("#"))
    i = next(k for k, l in enumerate(lines) if re.match(r"^\s*-?\d+\s+\d+\s*$", l))
    atoms, coords = [], []
    k = i + 1
    while k < len(lines) and len(lines[k].split()) >= 4:
        p = lines[k].split()
        atoms.append(re.sub(r"[^A-Za-z]", "", p[0]))
        coords.append(tuple(float(x) for x in p[1:4]))
        k += 1
    bonds = None
    if "connectivity" in route.lower():
        bonds = set()
        for line in lines[k + 1:]:
            p = line.split()
            if not p or not p[0].isdigit():
                if bonds:
                    break
                continue
            a = int(p[0]) - 1
            for j in range(1, len(p) - 1, 2):
                b = int(p[j]) - 1
                bonds.add(tuple(sorted((a, b))))
    return route, atoms, coords, bonds


def system_from_name(name):
    m = re.search(r"(Te|Se|O|S)(Te|Se|O|S)(Te|Se|O|S)(\d)", name.replace("-", ""))
    donor, acceptor, donor2, n = m.groups()
    assert donor == donor2, name
    return donor, acceptor, int(n)


def build_template(family, donor, acceptor, n_units, kekule=False):
    def unit(i):
        s = TEMPLATES[family]
        x, y = SYM[donor], SYM[acceptor]
        if kekule:
            x, y = KEKULE_FALLBACK[x], KEKULE_FALLBACK[y]
        s = s.replace("{X}", x).replace("{Y}", y)
        s = s.replace("[*:L]", f"[*:{100 + i}]" if i > 0 else "[H]")
        s = s.replace("[*:R]", f"[*:{101 + i}]" if i < n_units - 1 else "[H]")
        return s

    frags = Chem.MolFromSmiles(".".join(unit(i) for i in range(n_units)))
    if frags is None:
        return None
    mol = Chem.molzip(frags)
    Chem.SanitizeMol(mol)
    return Chem.AddHs(mol)


def skeleton(mol):
    """Element-labelled graph with all bonds single: canonical SMILES = graph identity."""
    rw = Chem.RWMol(mol)
    for a in rw.GetAtoms():
        a.SetIsAromatic(False); a.SetFormalCharge(0); a.SetNoImplicit(True); a.SetNumExplicitHs(0)
    for b in rw.GetBonds():
        b.SetBondType(Chem.BondType.SINGLE); b.SetIsAromatic(False)
    m = rw.GetMol()
    m.UpdatePropertyCache(strict=False)
    Chem.FastFindRings(m)
    return m


def xyz_mol(atoms, coords):
    block = f"{len(atoms)}\n\n" + "\n".join(f"{a} {x:.6f} {y:.6f} {z:.6f}" for a, (x, y, z) in zip(atoms, coords))
    m = Chem.MolFromXYZBlock(block)
    rdDetermineBonds.DetermineConnectivity(m)
    return m


def file_formula(atoms):
    from collections import Counter
    c = Counter(atoms)
    order = ["C", "H"] + sorted(k for k in c if k not in ("C", "H"))
    return "".join(f"{k}{c[k] if c[k] > 1 else ''}" for k in order if k in c)


def bond_length_check(mol, coords):
    """Return list of implausible chalcogen bonds (element pair, length)."""
    import math
    bad = []
    for b in mol.GetBonds():
        a1, a2 = b.GetBeginAtom(), b.GetEndAtom()
        pair = (a1.GetSymbol(), a2.GetSymbol())
        if pair not in BOND_RANGES:
            pair = pair[::-1]
            a1, a2 = a2, a1
        if pair in BOND_RANGES and a1.IsInRingSize(5) and a2.IsInRingSize(5):
            d = math.dist(coords[a1.GetIdx()], coords[a2.GetIdx()])
            lo, hi = BOND_RANGES[pair]
            if not lo <= d <= hi:
                bad.append(f"{pair[0]}-{pair[1]} {d:.3f}")
    return sorted(set(bad))


def graph_from_bonds(atoms, bonds):
    rw = Chem.RWMol()
    for a in atoms:
        rw.AddAtom(Chem.Atom(a))
    for a, b in bonds:
        rw.AddBond(a, b, Chem.BondType.SINGLE)
    return skeleton(rw.GetMol())


def main():
    OUT.mkdir(parents=True, exist_ok=True)
    rows, writer = [], Chem.SDWriter(str(OUT / "external_structures.sdf"))
    for path in sorted(SRC.rglob("*.gjf")):
        folder = path.parent.name  # e.g. B-DAD_Monomer
        family = folder[0]
        unit_type = folder.split("_")[1]
        donor, acceptor, n_units = system_from_name(path.stem)
        assert (n_units == 1) == (unit_type == "Monomer"), path
        route, atoms, coords, bonds = parse_gjf(path)

        tmpl = build_template(family, donor, acceptor, n_units)
        if tmpl is None:
            tmpl = build_template(family, donor, acceptor, n_units, kekule=True)
        sk_t = skeleton(tmpl)
        tmpl_key = Chem.MolToSmiles(sk_t)
        formula_t = CalcMolFormula(tmpl)
        formula_f = file_formula(atoms)
        geo = skeleton(xyz_mol(atoms, coords))
        geom_ok = Chem.MolToSmiles(geo) == tmpl_key
        gauss = graph_from_bonds(atoms, bonds) if bonds else None
        conn_ok = None if gauss is None else Chem.MolToSmiles(gauss) == tmpl_key
        graph_ok = bool(geom_ok or conn_ok)

        notes, bad = [], []
        if formula_f != formula_t:
            notes.append(f"file formula {formula_f} differs from intended {formula_t}")
        if graph_ok:
            ref = geo if geom_ok else gauss
            match = ref.GetSubstructMatch(sk_t)
            assert len(match) == tmpl.GetNumAtoms() == len(atoms)
            mapped = [coords[g] for g in match]
            bad = bond_length_check(tmpl, mapped)
            if bad:
                notes.append("implausible chalcogen bond lengths: " + "; ".join(bad))
            if not geom_ok:
                notes.append("identity from Gaussian connectivity block (distance perception differs)")
            if conn_ok is False:
                notes.append("Gaussian connectivity block incomplete; identity from coordinates")
        export = graph_ok and not bad and formula_f == formula_t

        heavy = Chem.RemoveHs(tmpl)
        smiles = Chem.MolToSmiles(heavy)
        inchi = Chem.MolToInchi(heavy)
        row = {
            "File": f"{folder}/{path.name}", "Family_Code": family, "Family": FAMILY_LABEL[family],
            "Unit_Type": unit_type, "System_Code": f"{donor}{acceptor}{donor}",
            "Formula": formula_t, "Formula_in_file": formula_f, "Atoms_in_file": len(atoms),
            "Geometry_graph_match": geom_ok, "Gaussian_connectivity_match": conn_ok,
            "Identity_confirmed": graph_ok, "Export_3D": export, "Notes": " | ".join(notes),
            "Canonical_SMILES": smiles, "InChI": inchi, "InChIKey": Chem.InchiToInchiKey(inchi),
            "Route": route,
        }
        rows.append(row)

        if export:
            conf = Chem.Conformer(tmpl.GetNumAtoms())
            for t_idx, xyz in enumerate(mapped):
                conf.SetAtomPosition(t_idx, xyz)
            out = Chem.Mol(tmpl)
            out.RemoveAllConformers(); out.AddConformer(conf, assignId=True)
            out.SetProp("_Name", path.stem)
            for k in ("File", "System_Code", "Unit_Type", "Family", "InChIKey"):
                out.SetProp(k, str(row[k]))
            writer.write(out)
    writer.close()
    df = pd.DataFrame(rows)
    df.to_csv(OUT / "external_structures.csv", index=False)
    pd.set_option("display.width", 250); pd.set_option("display.max_colwidth", 120)
    print(df[["File", "Identity_confirmed", "Export_3D", "Notes"]].to_string())
    print("identity confirmed:", int(df.Identity_confirmed.sum()), "/", len(df))
    print("3D exported:", int(df.Export_3D.sum()), "/", len(df))
    print(df.groupby(["Family_Code", "Unit_Type"])[["Identity_confirmed", "Export_3D"]].sum().to_string())


if __name__ == "__main__":
    main()
