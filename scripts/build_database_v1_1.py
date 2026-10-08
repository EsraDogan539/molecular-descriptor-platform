"""Build ChalMolDB Database v1.1 from the archived Database v1 table.

Input: data/releases/chalcogen_database_v1.csv.gz (Database v1, as archived on Zenodo).
Output: data/chalcogen_database_v1_master.csv.gz (the table the application reads).
Requires RDKit. Run from the repository root:

    python scripts/build_external_structures.py   # writes data/build/external_structures.{csv,sdf}
    python scripts/build_database_v1_1.py

Changes (see docs/DATABASE_CHANGELOG.md):
  1. Neutral collection metadata (Dataset_Owner, Dataset_Name, Paper_Use).
  2. Development records: Donor_ID / Acceptor_ID for all records, source DOI, DFT level.
  3. External records: the B-system family (112 records) is removed because its geometry files are not
     consistent with the reported systems; exact structures where they can be confirmed and repeat-unit
     SMILES for every remaining record. Record identifiers are not renumbered.
  4. Duplicate_Flag = the record's InChIKey is shared with at least one other record.
  5. 3D coordinates of external records whose geometry file passes all checks
     (data/structures_3d_external.sdf.gz).
  7. Source annotations from the source articles (scripts/annotate_external_sources.py).
  6. Coordinate-derived fields (Planarity_Proxy_Z_Range, Radius_of_Gyration, interatomic distances) for those
     records, with the same definitions as the development collection (structures_3d.geometry_descriptors).
"""

import gzip
import sys
from pathlib import Path

import pandas as pd
from rdkit import Chem
from rdkit.Chem.rdMolDescriptors import CalcMolFormula

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT / "scripts"))
import build_external_structures as bes  # noqa: E402

sys.path.insert(0, str(ROOT))
from structures_3d import geometry_descriptors, molblock_atoms  # noqa: E402
from annotate_external_sources import annotate  # noqa: E402

V1 = ROOT / "data" / "releases" / "chalcogen_database_v1.csv.gz"
MASTER = ROOT / "data" / "chalcogen_database_v1_master.csv.gz"
EXT_CSV = ROOT / "data" / "build" / "external_structures.csv"
EXT_SDF = ROOT / "data" / "build" / "external_structures.sdf"
EXT_3D = ROOT / "data" / "structures_3d_external.sdf.gz"

DEV_DOI = "https://doi.org/10.1021/acs.jcim.5c00345"
DEV_REFERENCE = (
    "T. Haciefendioglu, E. Yildirim, J. Chem. Inf. Model. 2025, 65, 5360-5369, Supporting Information "
    "(Band Gap and Reorganization Energy Prediction of Conducting Polymers by the Integration of Machine "
    "Learning and Density Functional Theory)"
)
EXT_COORD_SOURCE = (
    "H. Kayı group, Ankara University: B3LYP/LANL2DZ geometry file {file} "
    "(data/sources/kayi_group_geometries), released with ChalMolDB by the authors"
)
EXT_COORD_LICENSE = "CC BY 4.0 (https://creativecommons.org/licenses/by/4.0/)"

FAMILY_CODE = {
    "benzochalcogenadiazole": "A",
    "B_system / benzochalcogenadiazole-like DAD": "B",
    "chalcogendiazoloquinoxaline": "C",
}

STATUS = {
    "file_3d": (
        "Exact structure + 3D coordinates (author geometry file)",
        "Ready - exact structure confirmed against author geometry file; 3D coordinates included",
    ),
    "file_no3d": (
        "Exact structure (author geometry file); 3D coordinates not attached",
        "Ready - exact structure confirmed against author geometry file; no 3D coordinates for this record",
    ),
    "template": (
        "Exact structure from family template (verified on n = 1 and n = 6 geometry files)",
        "Ready - exact structure from verified family template; no 3D coordinates",
    ),
    "polymer": (
        "Repeat unit only (Repeat_Unit_SMILES); extrapolated polymer value",
        "External metadata ready; polymer limit, repeat unit given",
    ),
}


def repeat_unit_smiles(family, donor, acceptor):
    s = bes.TEMPLATES[family].replace("{X}", bes.SYM[donor]).replace("{Y}", bes.SYM[acceptor])
    s = s.replace("[*:L]", "[*]").replace("[*:R]", "[*]")
    return Chem.MolToSmiles(Chem.MolFromSmiles(s))


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
        "_formula": CalcMolFormula(mol_h),
    }


A_TABLE = ROOT / "data" / "sources" / "ozkilinc_kayi_2019_table4_oligomers.csv"
A_LEVELS = (("631Gd", "6-31G(d)", "no"), ("LANL2DZ", "LANL2DZ", "no"), ("LANL2DZ_PCM", "LANL2DZ", "yes"))
UNIT_NAMES = {1: "Monomer", 2: "Dimer", 3: "Trimer", 4: "Tetramer", 5: "Pentamer", 6: "Hexamer"}


def add_a_oligomers(df):
    """Monomer-hexamer records of the benzochalcogenadiazole family from Table 4 of the 2019 article."""
    table = pd.read_csv(A_TABLE, dtype=str, keep_default_na=False)
    a_poly = df[df.Family.eq("benzochalcogenadiazole")]
    next_id = max(int(r.split("_")[1]) for r in df.Record_ID if r.startswith("HAKAN_")) + 1
    assert next_id == 273, next_id
    new = []
    for _, t in table.iterrows():
        for suffix, basis, solvent in A_LEVELS:
            eg = t[f"Eg_{suffix}"]
            if eg == "–":
                continue
            like = a_poly[(a_poly.System_Code == t.System_Code) & (a_poly.Basis_Set == basis)
                          & (a_poly.Solvent_or_Conditions == solvent)]
            assert len(like) == 1, (t.System_Code, basis, solvent)
            rec = like.iloc[0].copy()
            for col in df.columns:
                if col not in ("Dataset_Owner", "Dataset_Name", "Split_Role", "Paper_Use", "Scope_Flag", "Family",
                               "System_Code", "Chalcogen_Type", "Donor_Chalcogen", "Acceptor_Chalcogen", "Method",
                               "Basis_Set", "Solvent_or_Conditions", "Reference", "Source_URL", "Notes"):
                    rec[col] = pd.NA
            n = int(t.n)
            rec["Record_ID"] = f"HAKAN_{next_id:04d}"
            rec["Source_File"] = A_TABLE.name
            rec["Unit_Type"] = UNIT_NAMES[n]
            rec["Oligomer_n"] = str(n)
            rec["HOMO_eV"] = float(t[f"HOMO_{suffix}"])
            rec["LUMO_eV"] = float(t[f"LUMO_{suffix}"])
            rec["Eg_eV"] = float(eg)
            rec["Experimental_Eg_eV"] = t.Experimental_Eg or pd.NA
            rec["Notes"] = "Oligomer value from Table 4 of the source article"
            rec["Duplicate_Flag"] = False
            new.append(rec)
            next_id += 1
    out = pd.concat([df, pd.DataFrame(new)], ignore_index=True)
    return out, len(new)


def split_code(code):
    import re
    donor, acceptor, donor2 = re.fullmatch(r"(Te|Se|O|S)(Te|Se|O|S)(Te|Se|O|S)", code).groups()
    assert donor == donor2
    return donor, acceptor


def main():
    df = pd.read_csv(V1, low_memory=False)
    assert len(df) == 3360
    if "Repeat_Unit_SMILES" not in df.columns:
        df.insert(df.columns.get_loc("Structure_Availability"), "Repeat_Unit_SMILES", pd.NA)
    if "Curation_Note" not in df.columns:
        df["Curation_Note"] = pd.NA
    for col in ("Repeat_Unit_SMILES", "Curation_Note", "Canonical_SMILES", "InChI", "InChIKey",
                "Structure_Availability", "Curation_Status", "Reference", "Source_URL", "Method",
                "Basis_Set", "Notes", "Family", "Donor_ID", "Acceptor_ID"):
        df[col] = df[col].astype("object")
    df["Curation_Note"] = pd.NA

    dev = df.Record_ID.str.startswith("EROL_")
    ext = df.Record_ID.str.startswith("HAKAN_")
    assert dev.sum() == 3088 and ext.sum() == 272

    # 1. Neutral collection metadata ---------------------------------------------------------------
    df.loc[dev, "Dataset_Owner"] = "Development collection"
    df.loc[dev, "Dataset_Name"] = "ChalMolDB_development_DA3088"
    df.loc[dev, "Paper_Use"] = "Curated development set and ML benchmark"
    df.loc[ext, "Dataset_Owner"] = "External collection"
    df.loc[ext, "Dataset_Name"] = "ChalMolDB_external_DAD"
    df.loc[ext, "Paper_Use"] = "External transferability benchmark"

    # 2. Development records -----------------------------------------------------------------------
    names = df.loc[dev, "Molecule_Name"].str.extract(r"^\s*(D\d+)\s*,?\s*(A\d+)?\s*$")
    df.loc[dev, "Donor_ID"] = names[0]
    df.loc[dev, "Acceptor_ID"] = names[1]
    a39 = ["EROL_1093", "EROL_1249"]
    a39_smiles = set(df.loc[df.Record_ID.isin(a39), "Acceptor_SMILES"])
    assert len(a39_smiles) == 1
    # Every other record with this acceptor fragment is A39, and A39 is the only acceptor missing for D21/D24.
    for rid, donor in (("EROL_1093", "D21"), ("EROL_1249", "D24")):
        present = set(df.loc[dev & (df.Donor_ID == donor), "Acceptor_ID"].dropna())
        assert "A39" not in present, (rid, donor)
    same_fragment = dev & ~df.Record_ID.isin(a39) & df.Acceptor_SMILES.isin(a39_smiles)
    assert same_fragment.sum() > 0 and set(df.loc[same_fragment, "Acceptor_ID"]) == {"A39"}
    df.loc[df.Record_ID.isin(a39), "Acceptor_ID"] = "A39"
    df.loc[df.Record_ID.isin(a39), "Curation_Note"] = (
        "Acceptor_ID A39 assigned by curation: the source title gives the donor only; the acceptor fragment "
        "is the one used by all other A39 records, and A39 is the only acceptor otherwise missing for this donor."
    )
    assert df.loc[dev, ["Donor_ID", "Acceptor_ID"]].notna().all().all()
    dup_pair = df.loc[dev].duplicated(["Donor_ID", "Acceptor_ID"], keep=False)
    dup_ids = df.loc[dev][dup_pair].Record_ID.tolist()
    assert dup_ids == ["EROL_0316", "EROL_0317"], dup_ids
    df.loc[df.Record_ID.isin(dup_ids), "Curation_Note"] = (
        "Source duplicate: the source SI contains D6A4 twice (Structures_316/317, identical structure) "
        "and no D6A6; both records are kept as published."
    )
    df.loc[dev, "Source_URL"] = DEV_DOI
    df.loc[dev, "Reference"] = DEV_REFERENCE
    df.loc[dev, "Method"] = "B3LYP (Jaguar)"
    df.loc[dev, "Basis_Set"] = "6-311+G(d)"

    # 3. External records --------------------------------------------------------------------------
    b_rows = df.Record_ID.str.startswith("HAKAN_") & df.Family.map(FAMILY_CODE).eq("B")
    assert b_rows.sum() == 112
    df = df[~b_rows].reset_index(drop=True)
    dev = df.Record_ID.str.startswith("EROL_")
    ext = df.Record_ID.str.startswith("HAKAN_")
    assert dev.sum() == 3088 and ext.sum() == 160
    df, n_new = add_a_oligomers(df)
    assert n_new == 240, n_new
    dev = df.Record_ID.str.startswith("EROL_")
    ext = df.Record_ID.str.startswith("HAKAN_")

    built = pd.read_csv(EXT_CSV)
    built["n"] = built.Unit_Type.map({"Monomer": 1, "Hexamer": 6})
    built = built.set_index(["Family_Code", "System_Code", "n"])
    sdf = {m.GetProp("File"): m for m in Chem.SDMolSupplier(str(EXT_SDF), removeHs=False) if m is not None}

    blocks = []
    for idx in df.index[ext]:
        rec = df.loc[idx]
        fam = FAMILY_CODE[rec.Family]
        donor, acceptor = split_code(rec.System_Code)
        assert (donor, acceptor) == (rec.Donor_Chalcogen, rec.Acceptor_Chalcogen), rec.Record_ID
        df.at[idx, "Repeat_Unit_SMILES"] = repeat_unit_smiles(fam, donor, acceptor)
        n = None if str(rec.Oligomer_n) == "polymer" else int(rec.Oligomer_n)

        key = (fam, rec.System_Code, n)
        note = None
        if n is None:
            status, mol = "polymer", None
        elif key in built.index:
            b = built.loc[key]
            assert b.Identity_confirmed, b.File
            mol = bes.build_template(fam, donor, acceptor, n)
            gas_phase_geometry = rec.Solvent_or_Conditions in ("no", "gas phase")
            if b.Export_3D and gas_phase_geometry:
                status = "file_3d"
            else:
                status = "file_no3d"
                note = (f"Geometry file {b.File}: {b.Notes}." if not b.Export_3D else
                        f"Geometry file {b.File} is the gas-phase B3LYP/LANL2DZ geometry; this record used a "
                        "PCM (acetonitrile) re-optimized geometry, so no coordinates are attached.")
        else:
            status, mol = "template", bes.build_template(fam, donor, acceptor, n)

        df.at[idx, "Structure_Availability"], df.at[idx, "Curation_Status"] = STATUS[status]
        if mol is not None:
            fields = structure_fields(mol)
            formula = fields.pop("_formula")
            for k, v in fields.items():
                df.at[idx, k] = v
            assert fields["O_Count"] + fields["S_Count"] + fields["Se_Count"] + fields["Te_Count"] > 0
            if status == "file_3d":
                b = built.loc[key]
                assert b.Formula == formula and b.InChIKey == fields["InChIKey"], rec.Record_ID
                out = Chem.Mol(sdf[b.File])
                public_id = rec.Record_ID.replace("HAKAN_", "EXT_", 1)
                out.SetProp("ChalMolDB_Record_ID", public_id)
                out.SetProp("ChalMolDB_InChIKey", fields["InChIKey"])
                out.SetProp("Coordinate_Source", EXT_COORD_SOURCE.format(file=b.File))
                out.SetProp("Coordinate_License", EXT_COORD_LICENSE)
                blocks.append(Chem.MolToMolBlock(out) + "".join(
                    f"> <{p}>\n{out.GetProp(p)}\n\n" for p in out.GetPropNames()) + "$$$$\n")
                for k, v in geometry_descriptors(molblock_atoms(blocks[-1])).items():
                    df.at[idx, k] = v
        else:
            for k in ("Canonical_SMILES", "InChI", "InChIKey"):
                df.at[idx, k] = pd.NA
        if note:
            df.at[idx, "Curation_Note"] = note

    # 4. Duplicate flag ----------------------------------------------------------------------------
    has_key = df.InChIKey.notna()
    df["Duplicate_Flag"] = has_key & df.InChIKey.duplicated(keep=False)

    df.loc[df.Curation_Note.isna(), "Curation_Note"] = ""
    df = annotate(df)
    assert len(df) == 3248 + 240
    df.to_csv(MASTER, index=False, compression={"method": "gzip", "mtime": 0})
    with gzip.GzipFile(EXT_3D, "wb", mtime=0) as fh:
        fh.write("".join(blocks).encode("utf-8"))

    summary = df.loc[ext].groupby(["Family", "Structure_Availability"]).size()
    print(summary.to_string())
    print("external records with exact structure:", int(df.loc[ext, "InChIKey"].notna().sum()))
    print("external 3D blocks:", len(blocks))
    print("Duplicate_Flag true:", int(df.Duplicate_Flag.sum()),
          "| cross-collection InChIKey overlaps:",
          len(set(df.loc[dev, "InChIKey"].dropna()) & set(df.loc[ext, "InChIKey"].dropna())))


if __name__ == "__main__":
    main()
