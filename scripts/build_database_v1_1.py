"""Build ChalMolDB Database v1.1 from the v1 master table.

All changes are deterministic assignments, so the script can be re-run on its own output.
Requires RDKit. Run from the repository root:

    python scripts/build_external_structures.py   # writes data/build/external_structures.{csv,sdf}
    python scripts/build_database_v1_1.py

Changes (see docs/DATABASE_CHANGELOG.md):
  1. Neutral collection metadata (Dataset_Owner, Dataset_Name, Paper_Use).
  2. Development records: Donor_ID / Acceptor_ID for all records, source DOI, DFT level.
  3. External records: exact structures where they can be confirmed, repeat-unit SMILES for every record,
     provenance of the B-system records, English source-family notes.
  4. Duplicate_Flag = the record's InChIKey is shared with at least one other record.
  5. 3D coordinates of external records whose geometry file passes all checks
     (data/structures_3d_external.sdf.gz).
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
B_REFERENCE = (
    "Kayı group DFT calculations (Ankara University), released with ChalMolDB; "
    "experimental polymer band gaps from E. Poverenov et al., J. Am. Chem. Soc. 2014, 136, 5138-5149, where available"
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
B_FAMILY_LABEL = "ethylenedioxychalcogenophene-benzochalcogenadiazole DAD (B system)"
FAMILY_CODE[B_FAMILY_LABEL] = "B"

STATUS = {
    "file_3d": (
        "Exact structure + 3D coordinates (author geometry file)",
        "Ready - exact structure confirmed against author geometry file; 3D coordinates included",
    ),
    "file_no3d": (
        "Exact structure (author geometry file); 3D coordinates withheld",
        "Exact structure confirmed; geometry file failed coordinate checks (see Curation_Note)",
    ),
    "template": (
        "Exact structure from family template (verified on n = 1 and n = 6 geometry files)",
        "Ready - exact structure from verified family template; no 3D coordinates",
    ),
    "pending": (
        "System-level metadata; exact structure pending author confirmation",
        "External metadata ready; structure and property provenance pending author confirmation",
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


def split_code(code):
    import re
    donor, acceptor, donor2 = re.fullmatch(r"(Te|Se|O|S)(Te|Se|O|S)(Te|Se|O|S)", code).groups()
    assert donor == donor2
    return donor, acceptor


def main():
    df = pd.read_csv(MASTER, low_memory=False)
    n0 = len(df)
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
    # This acceptor fragment occurs in no other record, and A39 is the only acceptor missing for donors D21/D24.
    for rid, donor in (("EROL_1093", "D21"), ("EROL_1249", "D24")):
        present = set(df.loc[dev & (df.Donor_ID == donor), "Acceptor_ID"].dropna())
        assert "A39" not in present, (rid, donor)
    other = df.loc[dev & ~df.Record_ID.isin(a39), "Acceptor_SMILES"]
    assert not other.isin(a39_smiles).any()
    df.loc[df.Record_ID.isin(a39), "Acceptor_ID"] = "A39"
    df.loc[df.Record_ID.isin(a39), "Curation_Note"] = (
        "Acceptor_ID A39 assigned by curation: source title gives the donor only; the acceptor fragment is "
        "unique to these two records and A39 is the only acceptor otherwise missing for this donor."
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
    b_rows = ext & df.Family.map(FAMILY_CODE).eq("B")
    df.loc[b_rows, "Family"] = B_FAMILY_LABEL
    df.loc[b_rows, "Reference"] = B_REFERENCE
    df.loc[b_rows, "Notes"] = df.loc[b_rows, "Notes"].astype(str).str.replace(
        "Source family: Band Gap Data-B Sistemi", "Source family: band-gap data, B system", regex=False
    )

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
            if fam == "B" and n == 6:
                status, mol = "pending", None
                note = (f"Geometry file {b.File} describes {b.Formula_in_file}, not the intended {b.Formula} "
                        "(two H atoms missing on one benzo ring); structure and property values pending "
                        "author confirmation.")
            elif b.Identity_confirmed:
                mol = bes.build_template(fam, donor, acceptor, n)
                status = "file_3d" if b.Export_3D else "file_no3d"
                if not b.Export_3D:
                    note = f"Geometry file {b.File}: {b.Notes}."
                    if fam == "B":
                        note += " Property values pending author confirmation."
            else:
                status, mol = "pending", None
                note = f"Geometry file {b.File}: {b.Notes}."
        elif fam == "C":
            status, mol = "template", bes.build_template(fam, donor, acceptor, n)
        else:
            status, mol = "pending", None
            if fam == "B":
                note = ("No geometry file; B-system monomer and hexamer files failed checks, so the structure "
                        "and property provenance of this series are pending author confirmation.")

        df.at[idx, "Structure_Availability"], df.at[idx, "Curation_Status"] = STATUS[status]
        if note:
            df.at[idx, "Curation_Note"] = note
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
        else:
            for k in ("Canonical_SMILES", "InChI", "InChIKey"):
                df.at[idx, k] = pd.NA

    # 4. Duplicate flag ----------------------------------------------------------------------------
    has_key = df.InChIKey.notna()
    df["Duplicate_Flag"] = has_key & df.InChIKey.duplicated(keep=False)

    df.loc[df.Curation_Note.isna(), "Curation_Note"] = ""
    assert len(df) == n0
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
