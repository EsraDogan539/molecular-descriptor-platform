"""Provenance and citation-ready export helpers for ChalMolDB."""

import re

import pandas as pd

from database_quality import available_reference_column
from release_metadata import DATABASE_VERSION, SCIENTIFIC_CORE_VERSION


DOI_PATTERN = re.compile(r"10\.\d{4,9}/[-._;()/:A-Z0-9]+", re.IGNORECASE)


def extract_doi(value):
    """Extract and normalize a DOI from a reference-like value when present."""
    if value is None:
        return None
    try:
        if pd.isna(value):
            return None
    except (TypeError, ValueError):
        pass

    text = str(value).strip()
    if not text:
        return None

    lowered = text.lower()
    for prefix in ("https://doi.org/", "http://doi.org/", "doi:"):
        if lowered.startswith(prefix):
            text = text[len(prefix):].strip()
            break

    match = DOI_PATTERN.search(text)
    if not match:
        return None

    return match.group(0).rstrip(".,);")


def doi_url(value):
    doi = extract_doi(value)
    return f"https://doi.org/{doi}" if doi else None


def record_reference_value(row):
    row_df = row.to_frame().T
    column = available_reference_column(row_df)
    if column is None:
        return None, None
    return column, row.get(column)


def citation_ready_record(
    row,
    public_record_id,
):
    """Return a citation-ready single-record table preserving source values."""
    reference_column, reference_value = record_reference_value(row)
    doi = extract_doi(reference_value)

    fields = [
        ("Database record ID", public_record_id),
        ("Database version", DATABASE_VERSION),
        ("Scientific core version", SCIENTIFIC_CORE_VERSION),
        ("Dataset / owner", row.get("Dataset_Owner")),
        ("Collection role", row.get("Split_Role")),
        ("Scope", row.get("Scope_Flag")),
        ("Molecule / system", row.get("Molecule_Name") or row.get("System_Code")),
        ("Canonical SMILES", row.get("Canonical_SMILES")),
        ("InChIKey", row.get("InChIKey")),
        ("Chalcogen type", row.get("Chalcogen_Type")),
        ("S count", row.get("S_Count")),
        ("Se count", row.get("Se_Count")),
        ("Te count", row.get("Te_Count")),
        ("HOMO (eV)", row.get("HOMO_eV")),
        ("LUMO (eV)", row.get("LUMO_eV")),
        ("Eg (eV)", row.get("Eg_eV")),
        ("Experimental Eg (eV)", row.get("Experimental_Eg_eV")),
        ("Method", row.get("Method")),
        ("Basis set", row.get("Basis_Set")),
        ("Reference field", reference_column),
        ("Reference / DOI", reference_value),
        ("DOI", doi),
        ("DOI URL", doi_url(reference_value)),
        ("Curation status", row.get("Curation_Status")),
    ]

    return pd.DataFrame(
        [{"Field": field, "Value": value} for field, value in fields]
    )
