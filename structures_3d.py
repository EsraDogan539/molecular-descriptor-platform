"""Access to the 3D source structures of the ChalMolDB development collection.

The coordinates are the DFT-optimized structures (B3LYP/6-311+G(d), Jaguar) published by
Haciefendioglu and Yildirim, J. Chem. Inf. Model. 2025, 65, 5360-5369
(https://doi.org/10.1021/acs.jcim.5c00345), Supporting Information ci5c00345_si_002.zip,
under a CC BY 4.0 licence. They are redistributed unmodified; ChalMolDB adds only the
record identifier, InChIKey, source and licence fields to each SDF entry.
"""

import gzip
from pathlib import Path

from public_labels import public_record_id

STRUCTURES_3D_PATH = Path("data/structures_3d_v1.sdf.gz")
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


def structure_block_for(record_id, blocks):
    """Return the SDF block for an internal or public record ID, or None."""
    if record_id is None:
        return None
    return blocks.get(public_record_id(str(record_id)))
