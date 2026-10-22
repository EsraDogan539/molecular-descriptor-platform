"""Build ChalMolDB Database v1.2: Database v1.1 plus the polymer collection.

Input:
  data/releases/chalcogen_database_v1_1.csv.gz   Database v1.1 as released (tag db-v1.1)
  data/sources/liu_2025_nanoscale/               source files of Liu, Yan & Liu, Nanoscale 2025 (MIT licence)
Output:
  data/chalcogen_database_v1_master.csv.gz       the table the application reads
  data/structures_3d_polymer.sdf.gz              oligomer geometries linked to polymer records
  data/build/polymer_import_report.json          counts and checks
Requires RDKit. Run from the repository root:  python scripts/build_database_v1_2.py

What is added (see docs/DATABASE_CHANGELOG.md):
  * One record per polymer (POL_0001-POL_1096 from the 1,096-polymer set, POL_1097- from the 227 polymers
    reported 2011-2024), with
      - experimental optical gap (Experimental_Eg_eV) and, for the 1,096 set, experimental HOMO and LUMO
        (cyclic voltammetry; new fields Experimental_HOMO_eV, Experimental_LUMO_eV);
      - DFT HOMO, LUMO and HOMO-LUMO gap of the methyl-truncated oligomer (B3LYP-D3/6-31G*, Gaussian 16) in
        HOMO_eV, LUMO_eV, Eg_eV, and the lowest singlet excitation of the oligomer (TD-B3LYP; new field
        TDDFT_S1_eV) where the source gives it;
      - exact structure of the computed oligomer (Canonical_SMILES, InChI, InChIKey, element and atom counts)
        and the source monomer SMILES with full side chains (new field Source_Monomer_SMILES);
      - 3D coordinates for oligomers whose geometry file has the same InChIKey as the oligomer SMILES.
  * Property values are copied as published; nothing is recomputed or rescaled.
"""

import gzip
import io
import json
import sys
import zipfile
from pathlib import Path

import pandas as pd
from rdkit import Chem, RDLogger
from rdkit.Chem import AllChem, rdDetermineBonds
from rdkit.Chem.rdMolDescriptors import CalcMolFormula

RDLogger.DisableLog("rdApp.*")
ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT))
from structures_3d import geometry_descriptors, molblock_atoms  # noqa: E402

V11 = ROOT / "data" / "releases" / "chalcogen_database_v1_1.csv.gz"
SRC = ROOT / "data" / "sources" / "liu_2025_nanoscale"
MASTER = ROOT / "data" / "chalcogen_database_v1_master.csv.gz"
POLY_3D = ROOT / "data" / "structures_3d_polymer.sdf.gz"
REPORT = ROOT / "data" / "build" / "polymer_import_report.json"

NANOSCALE_DOI = "https://doi.org/10.1039/d4nr03702b"
NANOSCALE_REF = "Liu, Yan & Liu, Nanoscale 2025, 17, 7865"
SAEKI_REF = "Nagasawa, Al-Naamani & Saeki, J. Phys. Chem. Lett. 2018, 9, 2639 (doi:10.1021/acs.jpclett.8b00635)"
COORD_SOURCE = ("Liu, Yan & Liu, Nanoscale 2025, 17, 7865 (doi:10.1039/d4nr03702b): B3LYP-D3/6-31G* oligomer "
                "geometry {file} from the authors' GitHub repository")
COORD_LICENSE = "MIT License (Liu Group, University of Florida); see data/sources/liu_2025_nanoscale/LICENSE"
NEW_FIELDS = ["Experimental_HOMO_eV", "Experimental_LUMO_eV", "TDDFT_S1_eV", "Source_Monomer_SMILES"]
COLLECTION = {
    "Dataset_Owner": "Polymer collection",
    "Dataset_Name": "ChalMolDB_polymer_exp1323",
    "Split_Role": "Experimental Benchmark",
    "Paper_Use": "Experimental optical-gap benchmark (conjugated polymers)",
    "Method": "B3LYP-D3 (Gaussian 16)",
    "Basis_Set": "6-31G*",
    "Solvent_or_Conditions": "not stated in source",
    "Unit_Type": "Oligomer (side chains truncated to methyl)",
    "Reference": NANOSCALE_REF,
    "Source_URL": NANOSCALE_DOI,
}


def structure_fields(mol_h):
    heavy = Chem.RemoveHs(mol_h)
    inchi = Chem.MolToInchi(heavy)
    counts = {e: sum(a.GetSymbol() == e for a in mol_h.GetAtoms()) for e in ("O", "S", "Se", "Te")}
    return {
        "Canonical_SMILES": Chem.MolToSmiles(heavy),
        "InChI": inchi,
        "InChIKey": Chem.InchiToInchiKey(inchi),
        "O_Count": counts["O"], "S_Count": counts["S"], "Se_Count": counts["Se"], "Te_Count": counts["Te"],
        "Atom_Count": mol_h.GetNumAtoms(),
        "Bond_Count": mol_h.GetNumBonds(),
        "Heavy_Atom_Count": heavy.GetNumAtoms(),
        "Heteroatom_Count": sum(a.GetSymbol() not in ("C", "H") for a in mol_h.GetAtoms()),
    }


def chalcogen_type(c):
    present = [e for e in ("S", "Se", "Te") if c[f"{e}_Count"] > 0]
    if not present:
        return "None (control)"
    return present[0] if len(present) == 1 else "Mixed"


def scope_flag(c):
    if c["S_Count"] + c["Se_Count"] + c["Te_Count"] > 0:
        return "Core_SSeTe"
    return "O_Control" if c["O_Count"] > 0 else "Control_NonSSeTe"


def read_sources():
    exp = pd.read_csv(SRC / "conjugated-polymer-1096dp.csv")
    olig = pd.read_csv(SRC / "oligomer-DFT-TDDFT-calculations-1096dp.csv")
    assert len(exp) == len(olig) == 1096 and exp.idx.tolist() == olig.idx.tolist()
    a = exp.merge(olig, on="idx", validate="1:1")
    a = pd.DataFrame({
        "set": "1096", "src_idx": a.idx, "name": a.Nickname,
        "exp_homo": a["HOMO(eV)"], "exp_lumo": a["LUMO(eV)"], "exp_gap": a["bandgap(eV)"],
        "exp_source": [f"{SAEKI_REF}, source ref. {r}" for r in a["Ref.No"]],
        "exp_url": "https://doi.org/10.1021/acs.jpclett.8b00635",
        "homo": a["HOMO-B3LYP(eV)"], "lumo": a["LUMO-B3LYP(eV)"], "gap": a["HOMO-LUMO-B3LYP(eV)"],
        "s1": a["S1-state-B3LYP(eV)"], "mono": a.smiles_monomer, "olig": a.smiles_oligomer,
        "year": pd.NA, "group": pd.NA,
    })
    parts = []
    for f in sorted((SRC / "new_227").glob("new-polymer-*.csv")):
        group = f.name.split("-")[2]
        p = pd.read_csv(f).rename(columns={"reference": "Reference"})
        d = pd.read_csv(next((SRC / "new_227").glob(f"oligomer-DFT-{group}-*.csv")))
        assert len(p) == len(d) and p.idx.tolist() == d.idx.tolist(), f.name
        m = p.merge(d, on="idx", validate="1:1")
        parts.append(pd.DataFrame({
            "set": "227", "src_idx": m.idx, "name": m.Name, "group": group, "year": m.Year,
            "exp_homo": pd.NA, "exp_lumo": pd.NA, "exp_gap": m["bandgap(eV)"],
            "exp_source": [f"{r} (via {NANOSCALE_REF}, {group} group)" for r in m.Reference],
            "exp_url": m.Reference,
            "homo": m["HOMO(eV)"], "lumo": m["LUMO(eV)"], "gap": m["HOMO-LUMO(eV)"], "s1": pd.NA,
            "mono": m.smiles_monomer, "olig": m.smiles_oligomer,
        }))
    b = pd.concat(parts, ignore_index=True)
    assert len(b) == 227, len(b)
    return pd.concat([a, b], ignore_index=True)


def xyz_mol(text):
    lines = text.strip("\n").splitlines()
    n = int(lines[0].split()[0])
    atoms = [ln.split() for ln in lines[2:2 + n]]
    block = f"{n}\n\n" + "\n".join(f"{a[0]} {a[1]} {a[2]} {a[3]}" for a in atoms) + "\n"
    mol = Chem.MolFromXYZBlock(block)
    rdDetermineBonds.DetermineConnectivity(mol)
    return mol


def attach_geometry(template_h, text):
    """Return (mol with template bond orders and file coordinates, reason)."""
    raw = xyz_mol(text)
    if CalcMolFormula(raw) != CalcMolFormula(template_h):
        return None, f"formula {CalcMolFormula(raw)} != {CalcMolFormula(template_h)}"
    try:
        mol = AllChem.AssignBondOrdersFromTemplate(template_h, raw)
        Chem.SanitizeMol(mol)
    except Exception as exc:  # connectivity of the file differs from the SMILES
        return None, f"bond assignment failed ({type(exc).__name__})"
    key_file = Chem.MolToInchiKey(Chem.RemoveHs(mol))
    key_smiles = Chem.MolToInchiKey(Chem.RemoveHs(template_h))
    if key_file != key_smiles:
        return None, "InChIKey of geometry differs from oligomer SMILES"
    return mol, None


# ---------------------------------------------------------------------------------------------
# Literature collection (LIT_): small, fully tabulated chalcogen series from recent articles
# ---------------------------------------------------------------------------------------------
HEL = ROOT / "data" / "sources" / "kumar_2025_helicenes"
BIC = ROOT / "data" / "sources" / "buyukmumcu_2026_se_te_bicyclic"
HEL_REF = "Kumar, Kar & Maity, J. Comput. Chem. 2025, 46, e70175"
HEL_URL = "https://doi.org/10.1002/jcc.70175"
BIC_REF = "Büyükmumcu & Selampinar, ACS Omega 2026, 11, 15694"
BIC_URL = "https://doi.org/10.1021/acsomega.5c04257"
HARTREE_EV = 27.2114
LIT_COLLECTION = {
    "Dataset_Owner": "Literature collection",
    "Dataset_Name": "ChalMolDB_literature_series",
    "Split_Role": "Literature Benchmark",
    "Paper_Use": "S/Se/Te series from recent articles (functional benchmark and chalcogen swap)",
}
HEL_NAMES = {"S": ("thia", "TH"), "Se": ("seleno", "SH"), "Te": ("telluro", "TeH")}
FUNCTIONAL_LABEL = {"wB97XD": "ωB97XD", "LC-wPBE": "LC-ωPBE"}
HEL_STRUCTURE_NOTE = (
    "Unsubstituted [n]helicene: the CCSD(T) total energies of the source (SI Tables S5, S8, S11) correspond to the "
    "bare ring system; the ring-fusion pattern is that of the telluro[n]helicene geometries in SI section 1 with "
    "the Br and TMS end groups replaced by H (S and Se analogues: same pattern). No 3D coordinates attached, "
    "because the source geometries belong to the end-substituted molecules."
)


def helicene_core_smiles():
    """{n: SMILES with Te} for n = 1-10 from the source geometries, end groups removed."""
    text = (HEL / "telluro_n_helicenes_end_substituted_B3LYP-D_DCM.xyz").read_text().splitlines()
    out, i = {}, 0
    while i < len(text):
        nat = int(text[i].split()[0])
        n = int(text[i + 1].split("]")[0].strip("["))
        block = f"{nat}\n\n" + "\n".join(text[i + 2:i + 2 + nat]) + "\n"
        mol = Chem.MolFromXYZBlock(block)
        rdDetermineBonds.DetermineConnectivity(mol)
        Chem.FastFindRings(mol)
        ri = mol.GetRingInfo()
        drop = [a.GetIdx() for a in mol.GetAtoms()
                if a.GetSymbol() in ("H", "Br", "Si") or (a.GetSymbol() == "C" and ri.NumAtomRings(a.GetIdx()) == 0)]
        rw = Chem.RWMol(mol)
        for idx in sorted(drop, reverse=True):
            rw.RemoveAtom(idx)
        core = rw.GetMol()
        for a in core.GetAtoms():
            a.SetIsAromatic(True)
            a.SetNoImplicit(False)
            a.SetNumExplicitHs(0)
        for b in core.GetBonds():
            b.SetBondType(Chem.BondType.AROMATIC)
            b.SetIsAromatic(True)
        smi = Chem.MolToSmiles(core, kekuleSmiles=False)
        check = Chem.MolFromSmiles(smi)
        assert check is not None, (n, smi)
        assert CalcMolFormula(Chem.AddHs(check)) == f"C{2 * n + 2}H4Te{n if n > 1 else ''}", (n, CalcMolFormula(check))
        out[n] = Chem.MolToSmiles(check)
        i += 2 + nat
    assert sorted(out) == list(range(1, 11))
    return out


def swap_chalcogen(smiles, element):
    mol = Chem.RWMol(Chem.MolFromSmiles(smiles))
    for a in mol.GetAtoms():
        if a.GetSymbol() == "Te":
            a.SetAtomicNum({"S": 16, "Se": 34, "Te": 52}[element])
    m = mol.GetMol()
    Chem.SanitizeMol(m)
    return Chem.MolToSmiles(m)


def literature_records(columns):
    core = helicene_core_smiles()
    rows = []

    def new(el, n, method, basis, homo, lunmo, eg, notes, src_file, src_id):
        prefix, code = HEL_NAMES[el]
        smi = swap_chalcogen(core[n], el)
        mol_h = Chem.AddHs(Chem.MolFromSmiles(smi))
        fields = structure_fields(mol_h)
        rec = {c: pd.NA for c in columns}
        rec.update(LIT_COLLECTION)
        rec.update(fields)
        rec.update({
            "Source_Record_ID": src_id, "Source_File": src_file,
            "Molecule_Name": f"{prefix}[{n}]helicene", "System_Code": f"[{n}]{code}",
            "Family": "chalcogenophene [n]helicenes (unsubstituted)", "Unit_Type": "Fused oligomer",
            "Oligomer_n": str(n), "HOMO_eV": homo, "LUMO_eV": lunmo, "Eg_eV": eg,
            "Method": method, "Basis_Set": basis, "Solvent_or_Conditions": "gas phase",
            "Reference": HEL_REF, "Source_URL": HEL_URL, "Notes": notes,
            "Chalcogen_Type": chalcogen_type(fields), "Scope_Flag": scope_flag(fields),
            "Structure_Availability": "Exact structure (ring system defined by the source); no 3D coordinates",
            "Curation_Status": "Ready - values transcribed from the source SI and checked",
            "Curation_Note": HEL_STRUCTURE_NOTE, "Duplicate_Flag": False,
        })
        rows.append(rec)

    orb = pd.read_csv(HEL / "helicene_orbital_energies_13_functionals.csv")
    for _, r in orb.iterrows():
        f = FUNCTIONAL_LABEL.get(r.functional, r.functional)
        geom = ("geometry optimized with the same functional" if r.set == "Set 1"
                else "single point on the B3LYP geometry")
        new(r.chalcogen, int(r.n), f"{f} (Gaussian 16; {r.set}: {geom})", "6-311++G(d,p); LANL2DZ for Te",
            r.HOMO_eV, r.LUMO_eV, round(r.LUMO_eV - r.HOMO_eV, 5),
            f"HOMO and LUMO from SI Tables {r.source_tables}; Eg = LUMO - HOMO computed from these values",
            "helicene_orbital_energies_13_functionals.csv", f"{r.chalcogen}:{r.n}:{r.functional}:{r.set}")
    hse = pd.read_csv(HEL / "helicene_hse06_table_s12.csv")
    for _, r in hse.iterrows():
        new(r.chalcogen, int(r.n), "HSE06 (Gaussian 16)", "6-311++G(d,p); LANL2DZ for Te",
            round(r.HOMO_au * HARTREE_EV, 5), round(r.LUMO_au * HARTREE_EV, 5), r.Eg_eV,
            "SI Table S12; HOMO/LUMO converted from hartree (x 27.2114), Eg as published",
            "helicene_hse06_table_s12.csv", f"{r.chalcogen}:{r.n}:HSE06")
    s13 = pd.read_csv(HEL / "helicene_gaps_table_s13.csv")
    for _, r in s13.iterrows():
        method = "B2PLYP-D3(BJ) (Gaussian 16)" if r.method.startswith("B2PLYP") else "ωB97XD (Gaussian 16)"
        new(r.chalcogen, int(r.n), method, "6-311++G(d,p); def2-TZVPP for Te", pd.NA, pd.NA, r.Eg_eV,
            "SI Table S13 (HOMO-LUMO gap only)", "helicene_gaps_table_s13.csv", f"{r.chalcogen}:{r.n}:{r.method}")
    cc = pd.read_csv(HEL / "helicene_ccsdt_tables_s5_s8_s11.csv")
    for _, r in cc.iterrows():
        new(r.chalcogen, int(r.n), "CCSD(T)//CCSD (fundamental gap, IP - EA)", "6-311++G(d,p)", pd.NA, pd.NA,
            round(r.IP_minus_EA_eV, 5),
            f"SI Tables S5/S8/S11: Eg = IP - EA = {r.IP_eV:.4f} - ({r.EA_eV:.4f}) eV (fundamental gap, "
            "reference values of the source benchmark)",
            "helicene_ccsdt_tables_s5_s8_s11.csv", f"{r.chalcogen}:{r.n}:CCSD(T)")

    bic = pd.read_csv(BIC / "monomer_frontier_orbitals_table2.csv")
    for _, r in bic.iterrows():
        sym = {"Se": "[se]", "Te": "[te]"}
        smi = f"c1cc2c{sym[r.X]}cc2{sym[r.Y]}1"
        mol_h = Chem.AddHs(Chem.MolFromSmiles(smi))
        fields = structure_fields(mol_h)
        rec = {c: pd.NA for c in columns}
        rec.update(LIT_COLLECTION)
        rec.update(fields)
        rec.update({
            "Source_Record_ID": f"Table2:{r.system}", "Source_File": "monomer_frontier_orbitals_table2.csv",
            "Molecule_Name": f"{r.system} fused bicyclic monomer (X = {r.X}, Y = {r.Y})", "System_Code": r.system,
            "Family": "chalcogenopheno[3,4-b]chalcogenophenes", "Unit_Type": "Monomer", "Oligomer_n": "1",
            "HOMO_eV": r.HOMO_eV, "LUMO_eV": r.LUMO_eV, "Eg_eV": r.Eg_eV,
            "Method": "B3PW91 (Gaussian 09/16)", "Basis_Set": "6-31G(d,p); LANL2DZ for Se and Te",
            "Solvent_or_Conditions": "gas phase", "Reference": BIC_REF, "Source_URL": BIC_URL,
            "Notes": "Table 2 of the source (monomer frontier orbitals)",
            "Chalcogen_Type": chalcogen_type(fields), "Scope_Flag": scope_flag(fields),
            "Structure_Availability": "Exact structure (Figure 1 of the source); no 3D coordinates",
            "Curation_Status": "Ready - values transcribed from the source table and checked",
            "Curation_Note": "", "Duplicate_Flag": False,
        })
        rows.append(rec)
    for k, rec in enumerate(rows, start=1):
        rec["Record_ID"] = f"LIT_{k:04d}"
    return rows


def main():
    base = pd.read_csv(V11, dtype=str, keep_default_na=False, na_values=[""])
    assert len(base) == 3488, len(base)
    for col in NEW_FIELDS:
        if col not in base.columns:
            base[col] = pd.NA
    src = read_sources()
    xyz = zipfile.ZipFile(SRC / "oligomer-xyz_opt_1096dp.zip")

    rows, blocks, report = [], [], {"structure_parse_failed": [], "geometry_rejected": {}}
    for i, s in src.iterrows():
        rid = f"POL_{i + 1:04d}"
        rec = {c: pd.NA for c in base.columns}
        rec.update(COLLECTION)
        rec.update({
            "Record_ID": rid, "Source_Record_ID": f"{s.set}:{s.src_idx}",
            "Source_File": ("conjugated-polymer-1096dp.csv" if s.set == "1096"
                            else f"new_227/new-polymer-{s.group}-*.csv"),
            "Molecule_Name": s["name"], "Family": ("Saeki compilation (1,096 polymers)" if s.set == "1096"
                                                   else f"Polymers reported 2011-2024 ({s.group} acceptor group)"),
            "HOMO_eV": s.homo, "LUMO_eV": s.lumo, "Eg_eV": s.gap, "TDDFT_S1_eV": s.s1,
            "Experimental_HOMO_eV": s.exp_homo, "Experimental_LUMO_eV": s.exp_lumo,
            "Experimental_Eg_eV": s.exp_gap, "Experimental_Eg_Min_eV": s.exp_gap, "Experimental_Eg_Max_eV": s.exp_gap,
            "Experimental_Eg_Source": s.exp_source, "Source_Monomer_SMILES": s.mono,
            "Duplicate_Flag": False,
            "Notes": ("Experimental optical gap of the polymer; DFT values of the methyl-truncated oligomer "
                      f"({NANOSCALE_REF})" + (f"; reported {int(s.year)}" if pd.notna(s.year) else "")),
            "Curation_Note": "",
        })
        mol = Chem.MolFromSmiles(s.olig)
        if mol is None:
            report["structure_parse_failed"].append(rid)
            rec["Structure_Availability"] = "Source SMILES could not be parsed; values kept, no structure"
            rec["Curation_Status"] = "Values ready; structure missing"
            rec["Scope_Flag"] = "Control_NonSSeTe"
            rec["Chalcogen_Type"] = "None (control)"
            rows.append(rec)
            continue
        mol_h = Chem.AddHs(mol)
        fields = structure_fields(mol_h)
        rec.update(fields)
        rec["Chalcogen_Type"] = chalcogen_type(fields)
        rec["Scope_Flag"] = scope_flag(fields)
        status = ("Exact structure of the computed oligomer (source SMILES)",
                  "Ready - oligomer structure from source SMILES; no 3D coordinates")
        if s.set == "1096":
            fname = f"donor-{int(s.src_idx)}_opt.xyz"
            geom, reason = attach_geometry(mol_h, xyz.read(fname).decode())
            if geom is None:
                report["geometry_rejected"][rid] = reason
                rec["Curation_Note"] = f"Geometry file {fname} not attached: {reason}."
            else:
                geom.SetProp("_Name", rid)
                for k, v in (("ChalMolDB_Record_ID", rid), ("ChalMolDB_InChIKey", fields["InChIKey"]),
                             ("Coordinate_Source", COORD_SOURCE.format(file=fname)),
                             ("Coordinate_License", COORD_LICENSE)):
                    geom.SetProp(k, v)
                block = Chem.MolToMolBlock(geom) + "".join(
                    f"> <{p}>\n{geom.GetProp(p)}\n\n" for p in geom.GetPropNames()) + "$$$$\n"
                blocks.append(block)
                rec.update(geometry_descriptors(molblock_atoms(block)))
                status = ("Exact structure + 3D coordinates (source geometry file)",
                          "Ready - oligomer structure from source SMILES; 3D coordinates included")
        rec["Structure_Availability"], rec["Curation_Status"] = status
        rows.append(rec)

    poly = pd.DataFrame(rows, columns=base.columns)
    lit = pd.DataFrame(literature_records(base.columns), columns=base.columns)
    report["literature_records"] = int(len(lit))
    df = pd.concat([base, poly, lit], ignore_index=True)
    has_key = df.InChIKey.notna()
    df["Duplicate_Flag"] = has_key & df.InChIKey.duplicated(keep=False)
    df["Curation_Note"] = df["Curation_Note"].fillna("")
    df.to_csv(MASTER, index=False, compression={"method": "gzip", "mtime": 0})
    with gzip.GzipFile(POLY_3D, "wb", mtime=0) as fh:
        fh.write("".join(blocks).encode("utf-8"))

    pol = df.Record_ID.str.startswith("POL_")
    keys = lambda m: set(df.loc[m & has_key, "InChIKey"])  # noqa: E731
    report.update({
        "records_total": int(len(df)), "polymer_records": int(pol.sum()),
        "polymer_3d": len(blocks),
        "polymer_chalcogen_type": df.loc[pol, "Chalcogen_Type"].value_counts().to_dict(),
        "polymer_scope": df.loc[pol, "Scope_Flag"].value_counts().to_dict(),
        "polymer_unique_inchikeys": int(df.loc[pol, "InChIKey"].nunique()),
        "polymer_duplicate_records": int(df.loc[pol, "Duplicate_Flag"].sum()),
        "overlap_with_v1_1_inchikeys": len(keys(pol) & keys(~pol)),
        "experimental_eg_total": int(pd.to_numeric(df.Experimental_Eg_eV, errors="coerce").notna().sum()),
        "unique_inchikeys_total": int(df.InChIKey.nunique()),
        "records_with_S_Se_Te": {e: int((pd.to_numeric(df[f"{e}_Count"], errors="coerce") > 0).sum())
                                 for e in ("S", "Se", "Te")},
    })
    REPORT.parent.mkdir(parents=True, exist_ok=True)
    REPORT.write_text(json.dumps(report, indent=1, default=str))
    print(json.dumps({k: v for k, v in report.items() if k != "geometry_rejected"}, indent=1, default=str))
    print("geometry rejected:", len(report["geometry_rejected"]))


if __name__ == "__main__":
    main()
