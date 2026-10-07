import hashlib
import json
from pathlib import Path

import pandas as pd

from public_labels import public_record_id
from release_metadata import DATABASE_VERSION

ROOT = Path(__file__).resolve().parents[1]


def test_public_csv_matches_master():
    master = pd.read_csv(ROOT / "data/chalcogen_database_v1_master.csv.gz", low_memory=False)
    public = pd.read_csv(ROOT / "data/public/chalmoldb_records.csv", low_memory=False)
    assert list(public.Record_ID) == [public_record_id(r) for r in master.Record_ID]
    assert list(public.columns) == list(master.columns)
    assert not public.Record_ID.str.startswith(("EROL_", "HAKAN_")).any()
    for col in ("Eg_eV", "InChIKey", "Radius_of_Gyration"):
        assert public[col].astype(str).tolist() == master[col].astype(str).tolist(), col


def test_manifest_checksums_and_version():
    manifest = json.loads((ROOT / "data/public/MANIFEST.json").read_text())
    assert manifest["database_version"] == DATABASE_VERSION
    for path, meta in manifest["files"].items():
        data = (ROOT / path).read_bytes()
        assert hashlib.sha256(data).hexdigest() == meta["sha256"], path
