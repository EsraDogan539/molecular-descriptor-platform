"""Public-facing labels and sanitization for ChalMolDB.

Internal source identifiers remain available in the curated source data, while
public UI/export surfaces use neutral collection labels and public record IDs.
"""

import pandas as pd


def public_record_id(value):
    if value is None:
        return value
    text = str(value)
    if text.startswith("EROL_"):
        return text.replace("EROL_", "DEV_", 1)
    if text.startswith("HAKAN_"):
        return text.replace("HAKAN_", "EXT_", 1)
    return text


def public_collection_label(split_role=None, dataset_owner=None):
    """Return a neutral public label based primarily on collection role."""
    role = "" if split_role is None or pd.isna(split_role) else str(split_role).strip()
    if role == "Development/Training":
        return "Development collection"
    if role == "External Validation":
        return "External literature collection"

    owner = "" if dataset_owner is None or pd.isna(dataset_owner) else str(dataset_owner).strip()
    lowered = owner.lower()
    if "erol" in lowered:
        return "Development collection"
    if "hakan" in lowered:
        return "External literature collection"
    return "Curated collection" if owner else "—"


def sanitize_public_dataframe(df):
    """Return a copy with person-linked source labels removed from public output."""
    public = df.copy()

    if "Record_ID" in public.columns:
        public["Record_ID"] = public["Record_ID"].map(public_record_id)

    if "Dataset_Owner" in public.columns:
        if "Split_Role" in public.columns:
            public["Dataset_Owner"] = [
                public_collection_label(role, owner)
                for role, owner in zip(public["Split_Role"], public["Dataset_Owner"])
            ]
        else:
            public["Dataset_Owner"] = public["Dataset_Owner"].map(
                lambda owner: public_collection_label(None, owner)
            )

    return public
