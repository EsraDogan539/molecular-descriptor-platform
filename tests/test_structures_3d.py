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


def test_external_bundle_matches_structured_external_records():
    from structures_3d import EXTERNAL_STRUCTURES_3D_PATH, coordinate_attribution

    blocks = load_structure_blocks(EXTERNAL_STRUCTURES_3D_PATH)
    db = pd.read_csv("data/chalcogen_database_v1_master.csv.gz", low_memory=False)
    expected = db.loc[
        db.Structure_Availability.eq("Exact structure + 3D coordinates (author geometry file)"), "Record_ID"
    ]
    assert len(blocks) == len(expected) == 59
    assert all(rid.startswith("HAKAN_") for rid in expected)
    for rid in expected:
        block = structure_block_for(rid, blocks)
        assert block is not None
        assert "Kay" in coordinate_attribution(block)


def test_attribution_is_per_collection():
    from structures_3d import COORDINATE_ATTRIBUTION, coordinate_attribution

    dev_block = parse_sdf_blocks(SAMPLE)["DEV_0001"]
    assert coordinate_attribution(dev_block) == COORDINATE_ATTRIBUTION
    ext_block = SAMPLE.replace("DEV_0001", "EXT_0001").replace(
        "$$$$", "> <Coordinate_Source>\nGroup X file\n\n> <Coordinate_License>\nCC BY 4.0\n\n$$$$"
    )
    caption = coordinate_attribution(parse_sdf_blocks(ext_block)["EXT_0001"])
    assert caption == "3D coordinates: Group X file; CC BY 4.0."


def test_geometry_descriptors_reproduce_database_values():
    from structures_3d import geometry_descriptors, load_all_structure_blocks, molblock_atoms

    db = pd.read_csv("data/chalcogen_database_v1_master.csv.gz", low_memory=False)
    blocks = load_all_structure_blocks()
    cols = ["Planarity_Proxy_Z_Range", "Radius_of_Gyration", "Max_Interatomic_Distance",
            "Mean_Interatomic_Distance", "Std_Interatomic_Distance"]
    with_3d = db[db.Record_ID.map(lambda rid: structure_block_for(rid, blocks) is not None)]
    assert len(with_3d) == 3088 + 59
    sample = pd.concat([with_3d[with_3d.Record_ID.str.startswith("EROL_")].iloc[::300],
                        with_3d[with_3d.Record_ID.str.startswith("HAKAN_")]])
    for _, row in sample.iterrows():
        got = geometry_descriptors(molblock_atoms(structure_block_for(row.Record_ID, blocks)))
        for col in cols:
            assert abs(got[col] - row[col]) < 2e-4, (row.Record_ID, col)


def test_viewer_html_embeds_structure():
    import pytest
    from structures_3d import structure_viewer_html

    pytest.importorskip("py3Dmol")
    html = structure_viewer_html(SAMPLE)
    assert html and "3Dmol" in html and "DEV_0001" in html
