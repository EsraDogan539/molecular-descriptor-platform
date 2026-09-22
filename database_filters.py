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


def apply_advanced_filters(
    df,
    *,
    numeric_ranges=None,
    minimum_counts=None,
    categorical_filters=None,
):
    """Apply optional filters without imputing missing values."""
    filtered = df.copy()
    numeric_ranges = numeric_ranges or {}
    minimum_counts = minimum_counts or {}
    categorical_filters = categorical_filters or {}

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
