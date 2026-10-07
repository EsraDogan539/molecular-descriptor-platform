"""Write the public, programmatic-access copy of the current database release to data/public/.

Outputs (stable paths; versioned through Git tags and Zenodo):
  data/public/chalmoldb_records.csv       all records, public DEV_/EXT_ identifiers
  data/public/MANIFEST.json               release version, row/column counts, SHA-256 of every distributed file
"""

import hashlib
import json
import sys
from pathlib import Path

import pandas as pd

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT))
from public_labels import sanitize_public_dataframe  # noqa: E402
from release_metadata import DATABASE_VERSION, SCIENTIFIC_CORE_VERSION  # noqa: E402

MASTER = ROOT / "data" / "chalcogen_database_v1_master.csv.gz"
OUT = ROOT / "data" / "public"
DISTRIBUTED = [
    "data/public/chalmoldb_records.csv",
    "data/structures_3d_v1.sdf.gz",
    "data/structures_3d_external.sdf.gz",
    "docs/chalcogen_descriptor_dictionary_v1.csv",
]


def sha256(path):
    h = hashlib.sha256()
    with open(path, "rb") as fh:
        for chunk in iter(lambda: fh.read(1 << 20), b""):
            h.update(chunk)
    return h.hexdigest()


def main():
    OUT.mkdir(parents=True, exist_ok=True)
    public = sanitize_public_dataframe(pd.read_csv(MASTER, low_memory=False))
    public.to_csv(OUT / "chalmoldb_records.csv", index=False, lineterminator="\n")
    manifest = {
        "name": "ChalMolDB",
        "database_version": DATABASE_VERSION,
        "scientific_core_version": SCIENTIFIC_CORE_VERSION,
        "records": int(len(public)),
        "fields": list(public.columns),
        "files": {p: {"sha256": sha256(ROOT / p), "bytes": (ROOT / p).stat().st_size} for p in DISTRIBUTED},
    }
    (OUT / "MANIFEST.json").write_text(json.dumps(manifest, indent=2) + "\n", encoding="utf-8")
    print(f"wrote {len(public)} records")


if __name__ == "__main__":
    main()
