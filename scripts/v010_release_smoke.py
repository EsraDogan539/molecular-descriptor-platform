"""Release-level smoke checks for ChalMolDB v0.10."""

import os
import sys

ROOT = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
if ROOT not in sys.path:
    sys.path.insert(0, ROOT)

import pandas as pd

from database_comparison import comparison_export, comparison_table
from provenance_export import citation_ready_record, extract_doi
from query_manifest import build_query_manifest
from structure_search import structure_search


def assert_true(condition, message):
    if not condition:
        raise AssertionError(message)


def smoke_structure_search():
    df = pd.DataFrame({
        "Record_ID": ["EROL_0001", "EROL_0002", "EROL_0003"],
        "Canonical_SMILES": ["c1ccsc1", "c1sccc1", "c1cc[se]c1"],
        "InChIKey": ["", "", ""],
        "HOMO_eV": [-5.5, -5.4, -5.2],
        "LUMO_eV": [-2.5, -2.4, -2.3],
        "Eg_eV": [3.0, 3.0, 2.9],
    })

    exact, query = structure_search(df, "c1sccc1", "Exact")
    assert_true(
        set(exact["Record_ID"]) == {"EROL_0001", "EROL_0002"},
        "Exact structure search failed equivalent-SMILES matching.",
    )
    assert_true(
        query["canonical_smiles"] == "c1ccsc1",
        "Query canonicalization changed unexpectedly.",
    )

    substructure, _ = structure_search(df, "c1ccsc1", "Substructure")
    assert_true(
        set(substructure["Record_ID"]) == {"EROL_0001", "EROL_0002"},
        "Substructure search returned unexpected records.",
    )

    similarity, _ = structure_search(
        df,
        "c1ccsc1",
        "Similarity",
        minimum_similarity=0.0,
        top_n=3,
    )
    assert_true(len(similarity) == 3, "Similarity search result count mismatch.")
    assert_true(
        float(similarity.iloc[0]["Structure Similarity"]) == 1.0,
        "Identical structure did not rank at similarity 1.0.",
    )


def smoke_comparison_and_provenance():
    reference = pd.Series({
        "Canonical_SMILES": "c1ccsc1",
        "InChIKey": "REFKEY",
        "Dataset_Owner": "Dataset A",
        "HOMO_eV": -5.5,
        "LUMO_eV": -2.5,
        "Eg_eV": 3.0,
        "Experimental_Eg_eV": "1.80; 1.50",
        "S_Count": 1,
        "Se_Count": 0,
        "Te_Count": 0,
        "Chalcogen_Type": "S",
        "Method": "DFT",
        "Basis_Set": "6-31G*",
        "DOI_or_Reference": "10.1000/example",
        "Split_Role": "Development/Training",
        "Scope_Flag": "Core_SSeTe",
        "Molecule_Name": "Reference",
        "Curation_Status": "Ready",
    })
    candidate = reference.copy()
    candidate["Canonical_SMILES"] = "c1cc[se]c1"
    candidate["Eg_eV"] = 2.8
    candidate["Chalcogen_Type"] = "Se"
    candidate["S_Count"] = 0
    candidate["Se_Count"] = 1

    table = comparison_table(reference, candidate)
    eg = table.loc[table["Field"] == "Eg (eV)"].iloc[0]
    assert_true(
        round(float(eg["Candidate - Reference"]), 4) == -0.2,
        "Comparison delta is incorrect.",
    )

    exported = comparison_export(
        reference,
        candidate,
        "DEV_0001",
        "DEV_0002",
        reference_column="DOI_or_Reference",
    )
    assert_true(
        "Scientific core version" in set(exported["Field"]),
        "Comparison export is missing release metadata.",
    )

    citation = citation_ready_record(reference, "DEV_0001")
    values = citation.set_index("Field")["Value"]
    assert_true(
        values["Database record ID"] == "DEV_0001",
        "Citation-ready export has the wrong public record ID.",
    )
    assert_true(
        extract_doi(values["Reference / DOI"]) == "10.1000/example",
        "DOI extraction failed in citation workflow.",
    )


def smoke_query_manifest():
    manifest = build_query_manifest(
        text_query="",
        structure_query="c1ccsc1",
        standardized_structure_query="c1ccsc1",
        structure_mode="Similarity",
        minimum_similarity=0.4,
        maximum_results=50,
        fingerprint_method="Morgan",
        fingerprint_radius=2,
        fingerprint_bits=2048,
        similarity_metric="Tanimoto",
        filters={
            "collection": "All",
            "scope": "All",
            "chalcogen": "All",
            "eg": "All",
            "reference_or_doi": "All",
        },
        result_record_ids=["EROL_0001", "EROL_0002"],
    )

    assert_true(manifest["result_count"] == 2, "Query manifest result count mismatch.")
    assert_true(
        manifest["query"]["fingerprint"]["method"] == "Morgan",
        "Query manifest fingerprint method is missing.",
    )
    assert_true(
        len(manifest["query_checksum_sha256"]) == 64,
        "Query checksum is not SHA-256 length.",
    )
    assert_true(
        len(manifest["result_record_ids_sha256"]) == 64,
        "Result-set checksum is not SHA-256 length.",
    )


def main():
    checks = [
        ("structure search", smoke_structure_search),
        ("comparison + provenance", smoke_comparison_and_provenance),
        ("reproducible query manifest", smoke_query_manifest),
    ]

    print("ChalMolDB v0.10 release smoke checks")
    print("-" * 48)
    for label, check in checks:
        check()
        print(f"PASS · {label}")
    print("-" * 48)
    print("v0.10 release smoke checks: PASS")


if __name__ == "__main__":
    main()
