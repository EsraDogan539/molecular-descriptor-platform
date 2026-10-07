import gzip
from pathlib import Path

import pandas as pd

from structures_3d import (
    RECORD_ID_FIELD,
    STRUCTURES_3D_PATH,
    load_structure_blocks,
    parse_sdf_blocks,
    structure_block_for,
)

SAMPLE = (
    "D58,A1\n  test\n\n  1  0  0  0  0  0            999 V2000\n"
    "    0.0000    0.0000    0.0000 C   0  0  0  0  0  0\nM  END\n"
    f"> <{RECORD_ID_FIELD}>\nDEV_0001\n\n$$$$\n"
)


def test_parse_sdf_blocks_keys_by_public_id():
    blocks = parse_sdf_blocks(SAMPLE + SAMPLE.replace("DEV_0001", "DEV_0002"))
    assert set(blocks) == {"DEV_0001", "DEV_0002"}
    assert blocks["DEV_0001"].rstrip().endswith("$$$$")


def test_structure_block_accepts_internal_and_public_ids():
    blocks = parse_sdf_blocks(SAMPLE)
    assert structure_block_for("EROL_0001", blocks) is not None
    assert structure_block_for("DEV_0001", blocks) is not None
    assert structure_block_for("HAKAN_0001", blocks) is None


def test_missing_bundle_returns_empty_mapping(tmp_path):
    assert load_structure_blocks(tmp_path / "absent.sdf.gz") == {}


def test_release_bundle_covers_every_development_record():
    blocks = load_structure_blocks(STRUCTURES_3D_PATH)
    db = pd.read_csv("data/chalcogen_database_v1_master.csv.gz", low_memory=False)
    dev_ids = db.loc[db.Split_Role == "Development/Training", "Record_ID"]
    assert len(blocks) == len(dev_ids) == 3088
    assert all(structure_block_for(rid, blocks) for rid in dev_ids)
