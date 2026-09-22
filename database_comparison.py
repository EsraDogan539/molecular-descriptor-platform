"""Curated record comparison helpers for ChalMolDB."""

import pandas as pd

from release_metadata import DATABASE_VERSION, SCIENTIFIC_CORE_VERSION


COMPARISON_FIELDS = [
    ("HOMO_eV", "HOMO (eV)"),
    ("LUMO_eV", "LUMO (eV)"),
    ("Eg_eV", "Eg (eV)"),
    ("Experimental_Eg_eV", "Experimental Eg (eV)"),
    ("S_Count", "S count"),
    ("Se_Count", "Se count"),
    ("Te_Count", "Te count"),
    ("Chalcogen_Type", "Chalcogen type"),
    ("Method", "Method"),
    ("Basis_Set", "Basis set"),
]


def _is_missing(value):
    if value is None:
        return True
    try:
        if pd.isna(value):
            return True
    except (TypeError, ValueError):
        pass
    return str(value).strip().lower() in {"", "nan", "none"}


def _numeric_delta(candidate, reference):
    try:
        if _is_missing(candidate) or _is_missing(reference):
            return None
        return round(float(candidate) - float(reference), 4)
    except (TypeError, ValueError):
        return None


def comparison_table(reference_row, candidate_row):
    """Build a descriptive pairwise comparison without collapsing source values."""
    rows = []
    for field, label in COMPARISON_FIELDS:
        if field not in reference_row.index and field not in candidate_row.index:
            continue

        reference_value = reference_row.get(field)
        candidate_value = candidate_row.get(field)
        rows.append({
            "Field": label,
            "Reference": reference_value,
            "Candidate": candidate_value,
            "Candidate - Reference": _numeric_delta(
                candidate_value,
                reference_value,
            ),
        })

    return pd.DataFrame(rows)


def comparison_export(
    reference_row,
    candidate_row,
    reference_public_id,
    candidate_public_id,
    reference_column=None,
):
    """Return a citation/provenance-aware long-form export table."""
    table = comparison_table(reference_row, candidate_row)

    reference_source = (
        reference_row.get(reference_column)
        if reference_column and reference_column in reference_row.index
        else None
    )
    candidate_source = (
        candidate_row.get(reference_column)
        if reference_column and reference_column in candidate_row.index
        else None
    )

    metadata = pd.DataFrame([
        {
            "Field": "Database record ID",
            "Reference": reference_public_id,
            "Candidate": candidate_public_id,
            "Candidate - Reference": None,
        },
        {
            "Field": "Canonical SMILES",
            "Reference": reference_row.get("Canonical_SMILES"),
            "Candidate": candidate_row.get("Canonical_SMILES"),
            "Candidate - Reference": None,
        },
        {
            "Field": "InChIKey",
            "Reference": reference_row.get("InChIKey"),
            "Candidate": candidate_row.get("InChIKey"),
            "Candidate - Reference": None,
        },
        {
            "Field": "Collection source",
            "Reference": public_collection_label(reference_row.get("Split_Role"), reference_row.get("Dataset_Owner")),
            "Candidate": public_collection_label(candidate_row.get("Split_Role"), candidate_row.get("Dataset_Owner")),
            "Candidate - Reference": None,
        },
        {
            "Field": "Reference / DOI",
            "Reference": reference_source,
            "Candidate": candidate_source,
            "Candidate - Reference": None,
        },
        {
            "Field": "Database version",
            "Reference": DATABASE_VERSION,
            "Candidate": DATABASE_VERSION,
            "Candidate - Reference": None,
        },
        {
            "Field": "Scientific core version",
            "Reference": SCIENTIFIC_CORE_VERSION,
            "Candidate": SCIENTIFIC_CORE_VERSION,
            "Candidate - Reference": None,
        },
    ])

    return pd.concat([metadata, table], ignore_index=True)
