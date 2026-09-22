"""Release-level smoke checks for ChalMolDB v0.9."""

import json
import os
import sys
import zipfile

ROOT = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
if ROOT not in sys.path:
    sys.path.insert(0, ROOT)

import pandas as pd

from database_metadata import add_database_export
from descriptor_engine import (
    run_molecular_descriptor_platform,
    validate_input_dataframe,
)
from similarity_search import find_similar_molecules


def assert_true(condition, message):
    if not condition:
        raise AssertionError(message)


def smoke_valid_dataset():
    df = pd.DataFrame({
        "Molecule_ID": ["S_001", "SE_001", "TE_001", "MIX_001"],
        "SMILES": [
            "c1ccsc1",
            "c1cc[se]c1",
            "c1cc[te]c1",
            "C[Te]c1ccsc1",
        ],
        "Eg_eV": [2.40, 2.10, 1.85, 1.72],
        "Property_Source": ["Smoke test"] * 4,
    })

    results = run_molecular_descriptor_platform(df, project_name="v09_smoke")
    results = add_database_export(df, results, project_name="v09_smoke")

    assert_true(len(results["valid_df"]) == 4, "Expected 4 valid molecules.")
    assert_true(len(results["invalid_df"]) == 0, "Expected no invalid molecules.")
    assert_true(
        results["summary_df"].iloc[0]["Descriptor Dictionary Version"] == "0.9.0",
        "Descriptor dictionary version missing from summary.",
    )

    with zipfile.ZipFile(results["zip_file"]) as archive:
        names = archive.namelist()
        manifest_name = next(
            name for name in names
            if name.endswith("_run_manifest.json")
        )
        manifest = json.loads(
            archive.read(manifest_name).decode("utf-8")
        )

    assert_true(
        "v09_smoke_curated_database.csv" in names,
        "Curated database export missing from ZIP.",
    )
    assert_true(
        "v09_smoke_descriptor_dictionary_v0.9.0.csv" in names,
        "Descriptor dictionary missing from ZIP.",
    )
    assert_true(
        manifest["analysis_summary"]["valid_molecules"] == 4,
        "Manifest valid-molecule count is incorrect.",
    )
    assert_true(
        "v09_smoke_curated_database.csv" in manifest["outputs"],
        "Manifest package inventory is stale.",
    )

    similarity = find_similar_molecules(
        results["valid_df"],
        reference_id="S_001",
        top_n=3,
        minimum_similarity=0.0,
    )
    assert_true(
        len(similarity) == 3,
        "Similarity search did not return expected result count.",
    )


def smoke_missing_columns():
    df = pd.DataFrame({
        "Molecule_ID": ["A"],
        "Structure": ["c1ccsc1"],
    })
    try:
        validate_input_dataframe(df)
    except ValueError as error:
        assert_true(
            "Missing required columns: SMILES" in str(error),
            "Missing-column validation message changed unexpectedly.",
        )
        return
    raise AssertionError("Missing SMILES column was not rejected.")


def smoke_mixed_invalid():
    df = pd.DataFrame({
        "Molecule_ID": ["S_OK", "BAD_001", "SE_OK", "EMPTY_001", "TE_OK"],
        "SMILES": [
            "c1ccsc1",
            "not_a_smiles",
            "c1cc[se]c1",
            "",
            "c1cc[te]c1",
        ],
    })

    results = run_molecular_descriptor_platform(
        df,
        project_name="v09_mixed",
    )
    summary = results["summary_df"].iloc[0]
    assert_true(int(summary["Total Records"]) == 5, "Submitted count mismatch.")
    assert_true(int(summary["Valid Molecules"]) == 3, "Valid count mismatch.")
    assert_true(int(summary["Invalid Molecules"]) == 2, "Invalid count mismatch.")

    statuses = set(results["invalid_df"]["Status"].astype(str))
    assert_true("Invalid SMILES" in statuses, "Invalid SMILES was not reported.")
    assert_true(
        "Empty SMILES" in statuses or "Missing SMILES" in statuses,
        "Empty SMILES was not reported.",
    )


def main():
    checks = [
        ("valid dataset + export package", smoke_valid_dataset),
        ("missing required column", smoke_missing_columns),
        ("mixed valid/invalid structures", smoke_mixed_invalid),
    ]

    print("ChalMolDB v0.9 release smoke checks")
    print("-" * 48)
    for label, check in checks:
        check()
        print(f"PASS · {label}")
    print("-" * 48)
    print("Release smoke checks: PASS")


if __name__ == "__main__":
    main()
