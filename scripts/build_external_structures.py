"""Build exact structures for external (D–A–D) records from the Kayı-group geometry files.

Input: Gaussian input files in data/sources/kayi_group_geometries/<Family>-DAD_<Monomer|Hexamer>/,
geometries optimized at B3LYP/LANL2DZ (see data/sources/kayi_group_geometries/README.md).

Bond orders cannot be read reliably from Cartesian coordinates, so each system is built from a
family template (known chemistry) and accepted only if the template's full atom/bond graph,
including hydrogens, is identical to the graph perceived from the coordinates. Where the input
file also carries a Gaussian connectivity block, that graph is checked as well.

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
        geo = xyz_mol(atoms, coords)
        sk_t, sk_g = skeleton(tmpl), skeleton(geo)
        geom_ok = Chem.MolToSmiles(sk_t) == Chem.MolToSmiles(sk_g)
        conn_ok = None
        if bonds:
            conn_ok = Chem.MolToSmiles(graph_from_bonds(atoms, bonds)) == Chem.MolToSmiles(sk_t)

        heavy = Chem.RemoveHs(tmpl)
        smiles = Chem.MolToSmiles(heavy)
        inchi = Chem.MolToInchi(heavy)
        row = {
            "File": f"{folder}/{path.name}", "Family_Code": family, "Family": FAMILY_LABEL[family],
            "Unit_Type": unit_type, "System_Code": f"{donor}{acceptor}{donor}",
            "Formula": CalcMolFormula(tmpl), "Atoms_in_file": len(atoms),
            "Geometry_graph_match": geom_ok, "Gaussian_connectivity_match": conn_ok,
            "Canonical_SMILES": smiles, "InChI": inchi, "InChIKey": Chem.InchiToInchiKey(inchi),
            "Route": route,
        }
        rows.append(row)

        if geom_ok:
            match = sk_g.GetSubstructMatch(sk_t)
            assert len(match) == tmpl.GetNumAtoms() == len(atoms)
            conf = Chem.Conformer(tmpl.GetNumAtoms())
            for t_idx, g_idx in enumerate(match):
                conf.SetAtomPosition(t_idx, coords[g_idx])
            out = Chem.Mol(tmpl)
            out.RemoveAllConformers(); out.AddConformer(conf, assignId=True)
            out.SetProp("_Name", path.stem)
            for k in ("File", "System_Code", "Unit_Type", "Family", "InChIKey"):
                out.SetProp(k, str(row[k]))
            writer.write(out)
    writer.close()
    df = pd.DataFrame(rows)
    df.to_csv(OUT / "external_structures.csv", index=False)
    print(df[["File", "Formula", "Geometry_graph_match", "Gaussian_connectivity_match"]].to_string())
    print("geometry matches:", int(df.Geometry_graph_match.sum()), "/", len(df))
    print("connectivity matches:", df.Gaussian_connectivity_match.value_counts(dropna=False).to_dict())
    if not df.Geometry_graph_match.all():
        sys.exit(1)


if __name__ == "__main__":
    main()
