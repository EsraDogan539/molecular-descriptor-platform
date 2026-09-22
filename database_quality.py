"""Database quality and coverage helpers for ChalMolDB v0.9."""

import pandas as pd


REFERENCE_COLUMNS = (
    "DOI_or_Reference",
    "DOI",
    "Reference",
    "Source_Reference",
    "Citation",
)

PROPERTY_COLUMNS = (
    ("HOMO", "HOMO_eV"),
    ("LUMO", "LUMO_eV"),
    ("Eg", "Eg_eV"),
    ("Experimental Eg", "Experimental_Eg_eV"),
)

PROVENANCE_COLUMNS = (
    ("Dataset / owner", "Dataset_Owner"),
    ("Method", "Method"),
    ("Basis set", "Basis_Set"),
    ("Curation status", "Curation_Status"),
)


def _meaningful_mask(series):
    if series is None:
        return pd.Series(dtype=bool)
    text = series.fillna("").astype(str).str.strip()
    return text.ne("") & ~text.str.lower().isin({"nan", "none", "null", "n/a", "na"})



def classify_numeric_measurement(value):
    """Classify a property cell without collapsing multiple reported measurements."""
    if pd.isna(value):
        return "missing"

    text = str(value).strip()
    if text.lower() in {"", "nan", "none", "null", "n/a", "na"}:
        return "missing"

    try:
        float(text)
        return "scalar"
    except ValueError:
        pass

    if ";" in text:
        parts = [part.strip() for part in text.split(";")]
        if len(parts) > 1 and all(parts):
            try:
                [float(part) for part in parts]
                return "multi-valued"
            except ValueError:
                pass

    return "non-numeric"


def available_reference_column(df):
    for column in REFERENCE_COLUMNS:
        if column in df.columns:
            return column
    return None


def reference_mask(df):
    column = available_reference_column(df)
    if column is None:
        return pd.Series(False, index=df.index)
    return _meaningful_mask(df[column]).reindex(df.index, fill_value=False)


def structure_mask(df):
    if "Canonical_SMILES" not in df.columns:
        return pd.Series(False, index=df.index)
    return _meaningful_mask(df["Canonical_SMILES"]).reindex(df.index, fill_value=False)


def method_mask(df):
    if "Method" not in df.columns:
        return pd.Series(False, index=df.index)
    return _meaningful_mask(df["Method"]).reindex(df.index, fill_value=False)


def provenance_label(row):
    reference_present = any(
        column in row.index
        and pd.notna(row.get(column))
        and str(row.get(column)).strip().lower() not in {"", "nan", "none", "null", "n/a", "na"}
        for column in REFERENCE_COLUMNS
    )
    method_present = (
        "Method" in row.index
        and pd.notna(row.get("Method"))
        and str(row.get("Method")).strip().lower() not in {"", "nan", "none", "null", "n/a", "na"}
    )
    owner_present = (
        "Dataset_Owner" in row.index
        and pd.notna(row.get("Dataset_Owner"))
        and str(row.get("Dataset_Owner")).strip().lower() not in {"", "nan", "none", "null", "n/a", "na"}
    )

    if reference_present:
        return "Reference available"
    if method_present:
        return "Method metadata"
    if owner_present:
        return "Dataset provenance"
    return "Limited metadata"


def coverage_table(df):
    total = len(df)
    rows = []

    for label, column in PROPERTY_COLUMNS:
        if column not in df.columns:
            continue

        if column == "Experimental_Eg_eV":
            classes = df[column].map(classify_numeric_measurement)
            available = int(classes.isin({"scalar", "multi-valued"}).sum())
            rows.append({
                "Field": label,
                "Available": available,
                "Missing": total - available,
                "Coverage (%)": round(100 * available / total, 1) if total else 0.0,
            })
            multi_valued = int(classes.eq("multi-valued").sum())
            if multi_valued:
                rows.append({
                    "Field": "Experimental Eg (multi-valued records)",
                    "Available": multi_valued,
                    "Missing": total - multi_valued,
                    "Coverage (%)": round(100 * multi_valued / total, 1) if total else 0.0,
                })
            continue

        available = int(pd.to_numeric(df[column], errors="coerce").notna().sum())
        rows.append({
            "Field": label,
            "Available": available,
            "Missing": total - available,
            "Coverage (%)": round(100 * available / total, 1) if total else 0.0,
        })

    if "Canonical_SMILES" in df.columns:
        available = int(structure_mask(df).sum())
        rows.append({
            "Field": "Standardized structure",
            "Available": available,
            "Missing": total - available,
            "Coverage (%)": round(100 * available / total, 1) if total else 0.0,
        })

    if "Method" in df.columns:
        available = int(method_mask(df).sum())
        rows.append({
            "Field": "Method metadata",
            "Available": available,
            "Missing": total - available,
            "Coverage (%)": round(100 * available / total, 1) if total else 0.0,
        })

    reference_column = available_reference_column(df)
    if reference_column is not None:
        available = int(reference_mask(df).sum())
        rows.append({
            "Field": "Reference / DOI",
            "Available": available,
            "Missing": total - available,
            "Coverage (%)": round(100 * available / total, 1) if total else 0.0,
        })

    return pd.DataFrame(rows)


def database_quality_summary(df):
    duplicate_record_ids = 0
    if "Record_ID" in df.columns:
        ids = df["Record_ID"].fillna("").astype(str).str.strip()
        duplicate_record_ids = int(
            (ids.ne("") & ids.duplicated(keep=False)).sum()
        )

    missing_record_ids = 0
    if "Record_ID" in df.columns:
        missing_record_ids = int((~_meaningful_mask(df["Record_ID"])).sum())

    repeated_structures = 0
    if "InChIKey" in df.columns:
        identities = df["InChIKey"].fillna("").astype(str).str.strip()
        valid = identities.ne("")
        repeated_structures = int(
            identities[valid].duplicated(keep=False).sum()
        )

    return {
        "Records": len(df),
        "Missing Record IDs": missing_record_ids,
        "Duplicate Record IDs": duplicate_record_ids,
        "Repeated Standardized Structures": repeated_structures,
        "Reference Column": available_reference_column(df) or "Not present",
    }
