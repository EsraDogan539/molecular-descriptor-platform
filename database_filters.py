"""Advanced curated-database filtering helpers for ChalMolDB."""

import pandas as pd


def numeric_bounds(df, column):
    """Return numeric min/max bounds for a column, or (None, None)."""
    if column not in df.columns:
        return None, None
    values = pd.to_numeric(df[column], errors="coerce").dropna()
    if values.empty:
        return None, None
    return float(values.min()), float(values.max())


TARGET_CHALCOGENS = ("S", "Se", "Te")
CHALCOGEN_COUNT_COLUMNS = {"S": "S_Count", "Se": "Se_Count", "Te": "Te_Count"}
CHALCOGEN_METADATA_COLUMNS = ("Donor_Chalcogen", "Acceptor_Chalcogen", "Chalcogen_Type")


def chalcogen_presence_mask(df, element):
    """Return a mask of records that contain ``element`` (S, Se or Te).

    Structure-complete records carry atom counts. Metadata-level records
    (no exact structure) carry the chalcogen identity only in donor/acceptor
    or chalcogen-type fields, so presence is taken from either source.
    No atom count is inferred for metadata-level records.
    """
    if element not in CHALCOGEN_COUNT_COLUMNS:
        raise ValueError(f"Unsupported chalcogen: {element}")

    present = pd.Series(False, index=df.index)

    count_column = CHALCOGEN_COUNT_COLUMNS[element]
    if count_column in df.columns:
        present = present | pd.to_numeric(df[count_column], errors="coerce").fillna(0).gt(0)

    # Exact element token, so "S" does not match "Se" and "Mixed" matches nothing.
    pattern = rf"(?:^|[^A-Za-z]){element}(?:[^A-Za-z]|$)"
    for column in CHALCOGEN_METADATA_COLUMNS:
        if column in df.columns:
            metadata = df[column].fillna("").astype(str)
            present = present | metadata.str.contains(pattern, case=True, regex=True)

    return present


def apply_advanced_filters(
    df,
    *,
    numeric_ranges=None,
    minimum_counts=None,
    categorical_filters=None,
    required_chalcogens=None,
):
    """Apply optional filters without imputing missing values."""
    filtered = df.copy()
    numeric_ranges = numeric_ranges or {}
    minimum_counts = minimum_counts or {}
    categorical_filters = categorical_filters or {}
    required_chalcogens = required_chalcogens or []

    for element in required_chalcogens:
        if element in CHALCOGEN_COUNT_COLUMNS:
            filtered = filtered[chalcogen_presence_mask(filtered, element)]

    for column, bounds in numeric_ranges.items():
        if column not in filtered.columns or bounds is None:
            continue
        lower, upper = bounds
        values = pd.to_numeric(filtered[column], errors="coerce")
        mask = values.notna() & values.ge(float(lower)) & values.le(float(upper))
        filtered = filtered[mask]

    for column, minimum in minimum_counts.items():
        if column not in filtered.columns or minimum is None:
            continue
        values = pd.to_numeric(filtered[column], errors="coerce")
        mask = values.notna() & values.ge(float(minimum))
        filtered = filtered[mask]

    for column, selected in categorical_filters.items():
        if column not in filtered.columns or selected in (None, "", "All"):
            continue
        filtered = filtered[filtered[column].astype(str) == str(selected)]

    return filtered
