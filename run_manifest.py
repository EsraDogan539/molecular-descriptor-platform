"""Machine-readable reproducibility manifest for ChalMolDB analysis exports."""

import hashlib
import json
import os
import platform
from datetime import datetime, timezone

import numpy as np
import pandas as pd
import rdkit


MANIFEST_SCHEMA_VERSION = "1.0"
SCIENTIFIC_CORE_VERSION = "0.9-development"


def dataframe_sha256(df):
    """Return a deterministic SHA-256 checksum of dataframe content and column order."""
    payload = df.to_csv(index=False, lineterminator="\n").encode("utf-8")
    return hashlib.sha256(payload).hexdigest()


def build_run_manifest(
    input_df,
    project_name,
    summary,
    descriptor_dictionary_version,
    output_files,
    fingerprint_settings,
):
    return {
        "manifest_schema_version": MANIFEST_SCHEMA_VERSION,
        "scientific_core_version": SCIENTIFIC_CORE_VERSION,
        "created_at_utc": datetime.now(timezone.utc).isoformat(),
        "project_name": str(project_name),
        "input": {
            "rows": int(len(input_df)),
            "columns": [str(column) for column in input_df.columns],
            "sha256": dataframe_sha256(input_df),
        },
        "analysis_summary": {
            "total_records": int(summary.get("Total Records", 0)),
            "valid_molecules": int(summary.get("Valid Molecules", 0)),
            "invalid_molecules": int(summary.get("Invalid Molecules", 0)),
            "duplicate_molecules": int(summary.get("Duplicate Molecules", 0)),
            "success_rate_percent": float(summary.get("Success Rate (%)", 0.0)),
        },
        "descriptor_dictionary_version": str(descriptor_dictionary_version),
        "fingerprints": fingerprint_settings,
        "software": {
            "python": platform.python_version(),
            "rdkit": getattr(rdkit, "__version__", "unknown"),
            "pandas": pd.__version__,
            "numpy": np.__version__,
        },
        "outputs": [
            os.path.basename(path)
            for path in output_files
            if path
        ],
    }


def write_run_manifest(manifest, path):
    with open(path, "w", encoding="utf-8") as handle:
        json.dump(
            manifest,
            handle,
            indent=2,
            ensure_ascii=False,
            sort_keys=True,
        )
        handle.write("\n")
    return path
