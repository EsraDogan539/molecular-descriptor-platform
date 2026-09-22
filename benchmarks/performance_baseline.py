import argparse
import time

import pandas as pd

from database_browser import _apply_search, _public_table
from descriptor_engine import create_fingerprint_dataset, process_molecular_dataset
from similarity_search import find_similar_molecules


DATABASE_PATH = "data/chalcogen_database_v1_master.csv.gz"
BASE_SMILES = [
    "c1ccsc1",
    "c1cc[se]c1",
    "c1cc[te]c1",
    "C[Te]c1ccsc1",
]


def timed(label, func):
    start = time.perf_counter()
    result = func()
    elapsed = time.perf_counter() - start
    print(f"{label}: {elapsed:.4f} s")
    return result, elapsed


def build_analysis_frame(rows):
    records = []
    for index in range(rows):
        records.append({
            "Molecule_ID": f"BENCH_{index:05d}",
            "SMILES": BASE_SMILES[index % len(BASE_SMILES)],
        })
    return pd.DataFrame(records)


def main():
    parser = argparse.ArgumentParser(description="ChalMolDB v0.9 performance baseline")
    parser.add_argument("--rows", type=int, default=200)
    args = parser.parse_args()

    print("ChalMolDB performance baseline")
    print(f"Analysis rows: {args.rows}")
    print("-" * 48)

    database_df, _ = timed(
        "Database load",
        lambda: pd.read_csv(DATABASE_PATH, compression="gzip"),
    )

    _, _ = timed(
        "Database search",
        lambda: _apply_search(database_df, "DEV_"),
    )

    _, _ = timed(
        "Public table formatting",
        lambda: _public_table(database_df),
    )

    analysis_df = build_analysis_frame(args.rows)

    processed, descriptor_seconds = timed(
        "Descriptor processing",
        lambda: process_molecular_dataset(analysis_df),
    )
    valid_df, _, _ = processed

    _, fingerprint_seconds = timed(
        "Fingerprint generation",
        lambda: create_fingerprint_dataset(analysis_df),
    )

    if len(valid_df) > 1:
        reference_id = valid_df.iloc[0]["Molecule_ID"]
        _, similarity_seconds = timed(
            "Similarity search",
            lambda: find_similar_molecules(
                valid_df,
                reference_id=reference_id,
                top_n=min(20, len(valid_df) - 1),
                minimum_similarity=0.0,
            ),
        )
    else:
        similarity_seconds = 0.0
        print("Similarity search: skipped")

    print("-" * 48)
    per_row_descriptor = descriptor_seconds / max(args.rows, 1)
    per_row_fingerprint = fingerprint_seconds / max(args.rows, 1)
    print(f"Descriptor per row: {per_row_descriptor:.6f} s")
    print(f"Fingerprint per row: {per_row_fingerprint:.6f} s")
    print(f"Similarity total: {similarity_seconds:.4f} s")


if __name__ == "__main__":
    main()
