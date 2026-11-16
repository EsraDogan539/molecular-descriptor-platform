"""Build ChalMolDB Database v1.3: Database v1.2 plus the external polymers and the literature model compounds.

Input:
  data/releases/chalcogen_database_v1_2.csv.gz                      Database v1.2 as released (tag db-v1.2)
  data/sources/dogan_kayi_2026_external/external_polymers_2015_2026.csv
  data/sources/dogan_kayi_2026_external/literature_model_compounds.csv
Output:
  data/chalcogen_database_v1_master.csv.gz       the table the application reads
  data/build/v1_3_import_report.json             counts and checks
Requires RDKit. Run from the repository root:  python scripts/build_database_v1_3.py

What is added (see docs/DATABASE_CHANGELOG.md):
  * Polymer collection, POL_1324-: linear D-A polymers published 2015-2026 with the optical gap measured in the
    source article (film, or solution where marked). The record structure is the methyl-truncated oligomer built
    from the repeat unit with the convention of the existing polymer records; the repeat unit with full side chains
    is kept in Source_Monomer_SMILES. No DFT values. Polymers that are already in the collection (same name with
    the same backbone, or the same source article) are not added again.
  * Literature collection, LIT_0219-: B3LYP HOMO-LUMO gaps of model compounds of literature polymers (gas phase,
    basis set and number of repeat units stated in the source article).
  * Curation notes on the 15 polymers of the 227-polymer set whose source article could not be verified.
The values of all v1.2 records are unchanged. Only the Curation_Note of these 15 records changes, and the derived
Duplicate_Flag of the v1.2 records whose structure now also occurs among the added records.
"""

import json
import sys
from pathlib import Path

import pandas as pd
from rdkit import Chem, RDLogger

RDLogger.DisableLog("rdApp.*")
ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT / "scripts"))
from build_database_v1_2 import chalcogen_type, scope_flag, structure_fields  # noqa: E402

V12 = ROOT / "data" / "releases" / "chalcogen_database_v1_2.csv.gz"
SRC = ROOT / "data" / "sources" / "dogan_kayi_2026_external"
MASTER = ROOT / "data" / "chalcogen_database_v1_master.csv.gz"
REPORT = ROOT / "data" / "build" / "v1_3_import_report.json"

UNVERIFIED = ["POL_1116", "POL_1117", "POL_1118", "POL_1119", "POL_1122", "POL_1123", "POL_1124", "POL_1125",
              "POL_1135", "POL_1141", "POL_1142", "POL_1143", "POL_1145", "POL_1147", "POL_1148"]
UNVERIFIED_NOTE = ("The source article given for this polymer in the compilation could not be verified "
                   "(checked in v1.3); the values are kept as published in the compilation.")

EXT_COLLECTION = {
    "Dataset_Owner": "Polymer collection",
    "Dataset_Name": "ChalMolDB_polymer_lit2015_2026",
    "Split_Role": "Experimental Benchmark",
    "Paper_Use": "Experimental optical-gap benchmark (polymers published 2015-2026; external test set)",
    "Unit_Type": "Oligomer (side chains truncated to methyl)",
}
LIT_COLLECTION = {
    "Dataset_Owner": "Literature collection",
    "Dataset_Name": "ChalMolDB_literature_models",
    "Split_Role": "Literature Benchmark",
    "Paper_Use": "DFT HOMO-LUMO gaps of model compounds of literature polymers",
}


def url(doi):
    return f"https://doi.org/{doi}" if isinstance(doi, str) and doi else pd.NA


def structure(smiles):
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        raise ValueError(f"unparsable SMILES {smiles}")
    fields = structure_fields(Chem.AddHs(mol))
    fields["Chalcogen_Type"] = chalcogen_type(fields)
    fields["Scope_Flag"] = scope_flag(fields)
    return fields


def external_records(columns, ext):
    rows = []
    for _, r in ext.iterrows():
        rec = {c: pd.NA for c in columns}
        rec.update(EXT_COLLECTION)
        rec.update(structure(r.model_oligomer_smiles))
        state = "film" if r.gap_state == "film" else "solution"
        rec.update({
            "Source_Record_ID": r.id, "Source_File": "external_polymers_2015_2026.csv",
            "Molecule_Name": r.polymer, "Family": f"Polymers reported {r.period} (external test collection)",
            "Oligomer_n": str(int(r.model_oligomer_n)), "Source_Monomer_SMILES": r.repeat_unit_smiles,
            "Experimental_Eg_eV": r.measured_optical_gap_eV, "Experimental_Eg_Min_eV": r.measured_optical_gap_eV,
            "Experimental_Eg_Max_eV": r.measured_optical_gap_eV,
            "Experimental_Eg_Source": r.reference, "Reference": r.reference, "Source_URL": url(r.doi),
            "Solvent_or_Conditions": f"optical gap measured in {state}",
            "Structure_Availability": ("Exact structure of the methyl-truncated oligomer built from the repeat unit; "
                                       "no 3D coordinates"),
            "Curation_Status": "Ready - optical gap and repeat unit taken from the source article",
            "Notes": (f"Experimental optical gap of the polymer ({state}); no DFT values. The oligomer follows the "
                      "convention of the polymer collection (side chains cut to methyl, chain ends capped with methyl)."),
            "Curation_Note": "", "Duplicate_Flag": False,
        })
        rows.append(rec)
    return rows


def literature_records(columns, lit):
    rows = []
    for _, r in lit.iterrows():
        rec = {c: pd.NA for c in columns}
        rec.update(LIT_COLLECTION)
        rec.update(structure(r.model_smiles))
        n = int(r.n_units)
        rec.update({
            "Source_Record_ID": r.id, "Source_File": "literature_model_compounds.csv",
            "Molecule_Name": f"{r.polymer} model ({n} repeat unit{'s' if n > 1 else ''})",
            "Family": "Model compounds of literature polymers (DFT)",
            "Unit_Type": "Oligomer model of a polymer", "Oligomer_n": str(n),
            "Source_Monomer_SMILES": r.repeat_unit_smiles,
            "HOMO_eV": r.homo_eV, "LUMO_eV": r.lumo_eV, "Eg_eV": r.gap_eV,
            "Method": r.functional, "Basis_Set": r.basis, "Solvent_or_Conditions": "gas phase",
            "Reference": r.reference, "Source_URL": url(r.doi),
            "Structure_Availability": "Structure of the model built from the repeat unit; no 3D coordinates",
            "Curation_Status": "Ready - values taken from the source article",
            "Notes": f"HOMO-LUMO gap of the model compound as reported in the source ({r.location_in_source}); "
                     f"model in the source: {r.model_note}",
            "Curation_Note": ("The structure follows the ChalMolDB oligomer convention (side chains cut to methyl); "
                              "the model of the source article may truncate the side chains differently."),
            "Duplicate_Flag": False,
        })
        rows.append(rec)
    return rows


def main():
    base = pd.read_csv(V12, dtype=str, keep_default_na=False, na_values=[""])
    assert len(base) == 5029, len(base)
    ext = pd.read_csv(SRC / "external_polymers_2015_2026.csv")
    ext = ext[ext.backbone != "in ChalMolDB"].reset_index(drop=True)
    lit = pd.read_csv(SRC / "literature_model_compounds.csv")
    lit = lit[lit.used_in_model].reset_index(drop=True)

    pol = pd.DataFrame(external_records(base.columns, ext), columns=base.columns)
    last_pol = int(base.Record_ID[base.Record_ID.str.startswith("POL_")].str[4:].astype(int).max())
    pol["Record_ID"] = [f"POL_{last_pol + k:04d}" for k in range(1, len(pol) + 1)]
    lt = pd.DataFrame(literature_records(base.columns, lit), columns=base.columns)
    last_lit = int(base.Record_ID[base.Record_ID.str.startswith("LIT_")].str[4:].astype(int).max())
    lt["Record_ID"] = [f"LIT_{last_lit + k:04d}" for k in range(1, len(lt) + 1)]

    # curation note on the unverified polymers (only cell changed in v1.2 records)
    m = base.Record_ID.isin(UNVERIFIED)
    assert m.sum() == len(UNVERIFIED)
    old = base.loc[m, "Curation_Note"].fillna("")
    base.loc[m, "Curation_Note"] = [f"{o} {UNVERIFIED_NOTE}".strip() for o in old]

    new = pd.concat([pol, lt], ignore_index=True)
    dup_old = new.InChIKey.isin(set(base.InChIKey.dropna()))
    df = pd.concat([base, new], ignore_index=True)
    # Duplicate_Flag marks every record whose structure (InChIKey) occurs more than once, as in v1.2
    has_key = df.InChIKey.notna()
    df["Duplicate_Flag"] = has_key & df.InChIKey.duplicated(keep=False)
    df["Curation_Note"] = df["Curation_Note"].fillna("")
    df.to_csv(MASTER, index=False, compression={"method": "gzip", "mtime": 0})

    report = {
        "records_total": int(len(df)), "records_v1_2": int(len(base)),
        "added_polymer_records": int(len(pol)), "added_literature_records": int(len(lt)),
        "added_polymers_film": int(ext.gap_state.eq("film").sum()),
        "added_new_structures_sharing_inchikey_with_v1_2": int(dup_old.sum()),
        "unverified_polymers_noted": len(UNVERIFIED),
        "polymer_records": int(df.Record_ID.str.startswith("POL_").sum()),
        "literature_records": int(df.Record_ID.str.startswith("LIT_").sum()),
        "unique_inchikeys_total": int(df.InChIKey.nunique()),
        "experimental_eg_total": int(pd.to_numeric(df.Experimental_Eg_eV, errors="coerce").notna().sum()),
        "core_SSeTe": int(df.Scope_Flag.eq("Core_SSeTe").sum()),
        "records_with_S_Se_Te": {e: int((pd.to_numeric(df[f"{e}_Count"], errors="coerce") > 0).sum())
                                 for e in ("S", "Se", "Te")},
        "added_chalcogen_type": new.Chalcogen_Type.value_counts().to_dict(),
    }
    REPORT.parent.mkdir(parents=True, exist_ok=True)
    REPORT.write_text(json.dumps(report, indent=1))
    print(json.dumps(report, indent=1))


if __name__ == "__main__":
    main()
