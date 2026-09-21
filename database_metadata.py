import os
import zipfile
import pandas as pd


METADATA_ALIASES = {
    "Compound_Name": ["Compound_Name", "Compound Name", "Name"],
    "HOMO_eV": ["HOMO_eV", "HOMO", "HOMO (eV)"],
    "LUMO_eV": ["LUMO_eV", "LUMO", "LUMO (eV)"],
    "Eg_eV": ["Eg_eV", "Eg", "Band_Gap", "Band Gap", "Bandgap"],
    "Property_Source": ["Property_Source", "Source", "Data_Source"],
    "DOI_or_Reference": ["DOI_or_Reference", "DOI", "Reference", "Citation"],
    "Experimental_or_Computational": [
        "Experimental_or_Computational",
        "Experimental/Computational",
        "Data_Type",
    ],
    "Calculation_Method": ["Calculation_Method", "Method", "DFT_Method"],
    "Basis_Set": ["Basis_Set", "Basis Set"],
    "Solvent_or_Conditions": [
        "Solvent_or_Conditions",
        "Solvent",
        "Conditions",
    ],
    "Data_Quality_Flag": ["Data_Quality_Flag", "Quality", "Quality_Flag"],
}


def _first_available_value(row, aliases):
    for column in aliases:
        if column in row.index:
            value = row[column]
            if pd.notna(value) and str(value).strip() != "":
                return value
    return pd.NA


def normalize_database_metadata(input_df):
    records = []

    for _, row in input_df.iterrows():
        record = {"Molecule_ID": str(row["Molecule_ID"]).strip()}
        for canonical_name, aliases in METADATA_ALIASES.items():
            record[canonical_name] = _first_available_value(row, aliases)
        records.append(record)

    metadata_df = pd.DataFrame(records)

    for numeric_col in ["HOMO_eV", "LUMO_eV", "Eg_eV"]:
        if numeric_col in metadata_df.columns:
            metadata_df[numeric_col] = pd.to_numeric(
                metadata_df[numeric_col], errors="coerce"
            )

    if "Data_Quality_Flag" in metadata_df.columns:
        metadata_df["Data_Quality_Flag"] = (
            metadata_df["Data_Quality_Flag"].fillna("Review")
        )

    return metadata_df


def build_curated_database(input_df, valid_df):
    metadata_df = normalize_database_metadata(input_df)
    database_df = valid_df.merge(
        metadata_df,
        on="Molecule_ID",
        how="left",
        validate="many_to_one",
    )

    preferred_order = [
        "Molecule_ID",
        "Original SMILES",
        "Canonical SMILES",
        "InChI",
        "InChIKey",
        "Molecular Formula",
        "Molecular Weight",
        "Chalcogen Type",
        "Sulfur Count",
        "Selenium Count",
        "Tellurium Count",
        "HOMO_eV",
        "LUMO_eV",
        "Eg_eV",
        "Property_Source",
        "DOI_or_Reference",
        "Experimental_or_Computational",
        "Calculation_Method",
        "Basis_Set",
        "Solvent_or_Conditions",
        "Data_Quality_Flag",
        "Valid SMILES",
        "Duplicate Flag",
    ]

    ordered = [c for c in preferred_order if c in database_df.columns]
    remaining = [c for c in database_df.columns if c not in ordered]
    return database_df[ordered + remaining]


def add_database_export(input_df, results, project_name):
    database_df = build_curated_database(input_df, results["valid_df"])
    safe_project_name = results.get("project_name", str(project_name).strip() or "chalcogen_project")
    output_dir = results.get("output_dir") or os.path.dirname(results["zip_file"]) or "."
    database_file = os.path.join(output_dir, f"{safe_project_name}_curated_database.csv")
    database_df.to_csv(database_file, index=False)

    zip_file = results["zip_file"]
    if os.path.exists(zip_file):
        with zipfile.ZipFile(
            zip_file,
            mode="a",
            compression=zipfile.ZIP_DEFLATED,
        ) as zip_output:
            zip_output.write(database_file, arcname=os.path.basename(database_file))

    results["database_df"] = database_df
    results["database_file"] = database_file
    return results
