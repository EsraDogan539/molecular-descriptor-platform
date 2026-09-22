import os
import sys

import pandas as pd
from rdkit import Chem, RDLogger

RDLogger.DisableLog("rdApp.error")

DATABASE_PATHS = [
    "data/chalcogen_database_v1_master.csv.gz",
    "data/chalcogen_database_v1_master.csv",
    "data/chalcogen_database_preview_v1.csv",
]

PROPERTY_COLUMNS = [
    "HOMO_eV",
    "LUMO_eV",
    "Eg_eV",
    "Experimental_Eg_eV",
]


def load_database():
    for path in DATABASE_PATHS:
        if os.path.exists(path):
            compression = "gzip" if path.endswith(".gz") else "infer"
            return pd.read_csv(path, compression=compression), path
    raise FileNotFoundError("No ChalMolDB database file was found.")


def meaningful(series):
    text = series.fillna("").astype(str).str.strip()
    return text.ne("") & ~text.str.lower().isin({"nan", "none", "null", "n/a", "na"})


def main():
    df, path = load_database()
    failures = []
    warnings = []

    print("ChalMolDB database release audit")
    print(f"Source: {path}")
    print(f"Records: {len(df):,}")
    print("-" * 52)

    if "Record_ID" not in df.columns:
        failures.append("Record_ID column is missing.")
        missing_ids = duplicate_ids = len(df)
    else:
        ids = df["Record_ID"].fillna("").astype(str).str.strip()
        missing_ids = int((~meaningful(df["Record_ID"])).sum())
        duplicate_ids = int((ids.ne("") & ids.duplicated(keep=False)).sum())
        if missing_ids:
            failures.append(f"{missing_ids} record(s) have missing Record_ID.")
        if duplicate_ids:
            failures.append(f"{duplicate_ids} record(s) participate in duplicate Record_ID values.")

    malformed_smiles = 0
    identity_mismatches = 0
    structure_rows = 0
    if "Canonical_SMILES" in df.columns:
        smiles_present = meaningful(df["Canonical_SMILES"])
        structure_rows = int(smiles_present.sum())
        for idx in df.index[smiles_present]:
            smiles = str(df.at[idx, "Canonical_SMILES"]).strip()
            mol = Chem.MolFromSmiles(smiles)
            if mol is None:
                malformed_smiles += 1
                continue

            if "InChIKey" in df.columns:
                stored = df.at[idx, "InChIKey"]
                if pd.notna(stored) and str(stored).strip():
                    try:
                        generated = Chem.InchiToInchiKey(Chem.MolToInchi(mol))
                    except Exception:
                        generated = ""
                    if generated and generated != str(stored).strip():
                        identity_mismatches += 1

        if malformed_smiles:
            failures.append(
                f"{malformed_smiles} non-empty Canonical_SMILES value(s) could not be parsed."
            )
        if identity_mismatches:
            warnings.append(
                f"{identity_mismatches} structure row(s) have a stored InChIKey that differs "
                "from the key regenerated from Canonical_SMILES."
            )

    non_numeric = {}
    for column in PROPERTY_COLUMNS:
        if column not in df.columns:
            continue
        present = meaningful(df[column])
        converted = pd.to_numeric(df[column], errors="coerce")
        bad = int((present & converted.isna()).sum())
        if bad:
            non_numeric[column] = bad
            warnings.append(f"{bad} non-numeric populated value(s) in {column}.")

    print(f"Missing Record IDs: {missing_ids}")
    print(f"Records with duplicate Record IDs: {duplicate_ids}")
    print(f"Rows with standardized structure: {structure_rows}")
    print(f"Malformed populated Canonical SMILES: {malformed_smiles}")
    print(f"SMILES/InChIKey mismatches: {identity_mismatches}")
    if non_numeric:
        print("Non-numeric property values:")
        for column, count in non_numeric.items():
            print(f"  - {column}: {count}")
    else:
        print("Non-numeric populated property values: 0")

    if warnings:
        print("-" * 52)
        print("Warnings")
        for warning in warnings:
            print(f"WARNING: {warning}")

    if failures:
        print("-" * 52)
        print("Critical failures")
        for failure in failures:
            print(f"ERROR: {failure}")
        sys.exit(1)

    print("-" * 52)
    print("Database release audit: PASS")


if __name__ == "__main__":
    main()
