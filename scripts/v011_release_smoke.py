"""Release-level smoke checks for ChalMolDB v0.11."""

import os
import sys

ROOT = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
if ROOT not in sys.path:
    sys.path.insert(0, ROOT)

import pandas as pd

from database_filters import apply_advanced_filters
from public_labels import public_collection_label, public_record_id, sanitize_public_dataframe
from query_manifest import build_query_manifest, load_query_manifest, replay_configuration


def assert_true(condition, message):
    if not condition:
        raise AssertionError(message)


def smoke_advanced_filters():
    df = pd.DataFrame({
        "Record_ID": ["EROL_0001", "EROL_0002", "HAKAN_0001"],
        "HOMO_eV": [-5.6, -5.0, -4.8],
        "Eg_eV": [3.0, 2.5, 2.1],
        "S_Count": [2, 0, 1],
        "Se_Count": [0, 2, 1],
        "Method": ["DFT", "DFT", "Experiment"],
    })
    filtered = apply_advanced_filters(
        df,
        numeric_ranges={"Eg_eV": (2.4, 3.1)},
        minimum_counts={"S_Count": 1},
        categorical_filters={"Method": "DFT"},
    )
    assert_true(
        filtered["Record_ID"].tolist() == ["EROL_0001"],
        "Advanced filters returned unexpected records.",
    )


def smoke_public_labels():
    df = pd.DataFrame({
        "Record_ID": ["EROL_0001", "HAKAN_0001"],
        "Dataset_Owner": ["Erol dataset", "Hakan dataset"],
        "Split_Role": ["Development/Training", "External Validation"],
    })
    public = sanitize_public_dataframe(df)
    text = public.to_csv(index=False).lower()
    assert_true("erol" not in text and "hakan" not in text, "Person-linked labels leaked.")
    assert_true(public_record_id("EROL_0001") == "DEV_0001", "DEV public ID mapping failed.")
    assert_true(public_record_id("HAKAN_0001") == "EXT_0001", "EXT public ID mapping failed.")
    assert_true(
        public_collection_label("Development/Training", "Erol dataset") == "Development collection",
        "Development collection label failed.",
    )


def smoke_query_replay():
    manifest = build_query_manifest(
        text_query="thiophene",
        structure_query="c1ccsc1",
        standardized_structure_query="c1ccsc1",
        structure_mode="Similarity",
        minimum_similarity=0.55,
        maximum_results=25,
        fingerprint_method="Morgan",
        fingerprint_radius=2,
        fingerprint_bits=2048,
        similarity_metric="Tanimoto",
        filters={
            "collection": "Development/Training",
            "scope": "Core_SSeTe",
            "chalcogen": "S",
            "eg": "Available",
            "reference_or_doi": "Available",
            "advanced_numeric_ranges": {"Eg_eV": [2.0, 3.0]},
            "advanced_minimum_counts": {"S_Count": 1},
            "advanced_categorical": {"Method": "DFT"},
        },
        result_record_ids=["EROL_0001"],
    )
    loaded = load_query_manifest(manifest)
    replay = replay_configuration(loaded)
    assert_true(replay["text_query"] == "thiophene", "Text query replay failed.")
    assert_true(replay["structure_mode"] == "Similarity", "Structure mode replay failed.")
    assert_true(replay["minimum_similarity"] == 0.55, "Similarity threshold replay failed.")
    assert_true(replay["maximum_results"] == 25, "Maximum-result replay failed.")
    assert_true(
        replay["filters"]["advanced_numeric_ranges"]["Eg_eV"] == [2.0, 3.0],
        "Advanced filter replay failed.",
    )


def main():
    checks = [
        ("advanced filters", smoke_advanced_filters),
        ("neutral public labels", smoke_public_labels),
        ("query replay", smoke_query_replay),
    ]

    print("ChalMolDB v0.11 release smoke checks")
    print("-" * 48)
    for label, check in checks:
        check()
        print(f"PASS · {label}")
    print("-" * 48)
    print("v0.11 release smoke checks: PASS")


if __name__ == "__main__":
    main()
