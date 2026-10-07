"""Access to the 3D structures distributed with ChalMolDB.

Development collection (data/structures_3d_v1.sdf.gz): the coordinates are the DFT-optimized structures (B3LYP/6-311+G(d), Jaguar) published by
Haciefendioglu and Yildirim, J. Chem. Inf. Model. 2025, 65, 5360-5369
(https://doi.org/10.1021/acs.jcim.5c00345), Supporting Information ci5c00345_si_002.zip,
under a CC BY 4.0 licence. They are redistributed unmodified; ChalMolDB adds only the
record identifier, InChIKey, source and licence fields to each SDF entry.

External collection (data/structures_3d_external.sdf.gz): B3LYP/LANL2DZ geometries supplied by the
H. Kayı group (Ankara University) for the records whose geometry file passes the identity, formula and
bond-length checks of scripts/build_external_structures.py.
"""

import gzip
from pathlib import Path

from public_labels import public_record_id

STRUCTURES_3D_PATH = Path("data/structures_3d_v1.sdf.gz")
EXTERNAL_STRUCTURES_3D_PATH = Path("data/structures_3d_external.sdf.gz")
STRUCTURE_BUNDLES = (STRUCTURES_3D_PATH, EXTERNAL_STRUCTURES_3D_PATH)
RECORD_ID_FIELD = "ChalMolDB_Record_ID"
COORDINATE_ATTRIBUTION = (
    "3D coordinates: Haciefendioglu & Yildirim, J. Chem. Inf. Model. 2025, 65, 5360–5369 "
    "(doi:10.1021/acs.jcim.5c00345), Supporting Information; B3LYP/6-311+G(d); "
    "CC BY 4.0; redistributed unmodified."
)


def parse_sdf_blocks(text):
    """Split multi-record SDF text into {public record ID: SDF block}."""
    blocks = {}
    for raw in text.split("$$$$"):
        block = raw.strip("\n")
        if not block.strip():
            continue
        lines = block.splitlines()
        record_id = None
        for i, line in enumerate(lines):
            if line.strip() == f"> <{RECORD_ID_FIELD}>" and i + 1 < len(lines):
                record_id = lines[i + 1].strip()
                break
        if record_id:
            blocks[record_id] = block + "\n$$$$\n"
    return blocks


def load_structure_blocks(path=STRUCTURES_3D_PATH):
    """Load all 3D structure blocks; returns an empty mapping if the bundle is absent."""
    path = Path(path)
    if not path.exists():
        return {}
    with gzip.open(path, "rt", encoding="utf-8") as handle:
        return parse_sdf_blocks(handle.read())


def load_all_structure_blocks(paths=STRUCTURE_BUNDLES):
    """Load and merge every available bundle (development and external)."""
    blocks = {}
    for path in paths:
        blocks.update(load_structure_blocks(path))
    return blocks


def sdf_field(block, name):
    """Return the value of an SDF data field in ``block``, or None."""
    lines = block.splitlines()
    for i, line in enumerate(lines):
        if line.strip() == f"> <{name}>" and i + 1 < len(lines):
            return lines[i + 1].strip()
    return None


def coordinate_attribution(block):
    """Attribution caption for one SDF block, built from its source and licence fields."""
    record_id = sdf_field(block, RECORD_ID_FIELD) or ""
    if record_id.startswith("DEV_"):
        return COORDINATE_ATTRIBUTION
    source = sdf_field(block, "Coordinate_Source")
    license_text = sdf_field(block, "Coordinate_License")
    if not source:
        return None
    return f"3D coordinates: {source}" + (f"; {license_text}." if license_text else ".")


def structure_block_for(record_id, blocks):
    """Return the SDF block for an internal or public record ID, or None."""
    if record_id is None:
        return None
    return blocks.get(public_record_id(str(record_id)))
