"""Database v1.2: the polymer collection is added and the v1.1 records are unchanged."""

import gzip

import pandas as pd
import pytest

DB = "data/chalcogen_database_v1_master.csv.gz"
V11 = "data/releases/chalcogen_database_v1_1.csv.gz"
POLY_SDF = "data/structures_3d_polymer.sdf.gz"
SRC = "data/sources/liu_2025_nanoscale"


@pytest.fixture(scope="module")
def db():
    return pd.read_csv(DB, dtype=str, keep_default_na=False, na_values=[""])


@pytest.fixture(scope="module")
def pol(db):
    return db[db.Record_ID.str.startswith("POL_")].copy()


def test_v1_1_records_unchanged(db):
    v11 = pd.read_csv(V11, dtype=str, keep_default_na=False, na_values=[""])
    head = db.iloc[: len(v11)][v11.columns].reset_index(drop=True)
    pd.testing.assert_frame_equal(head, v11)


def test_polymer_ids_and_collection(pol):
    assert pol.Record_ID.tolist() == [f"POL_{i:04d}" for i in range(1, 1324)]
    assert set(pol.Split_Role) == {"Experimental Benchmark"}
    assert set(pol.Dataset_Owner) == {"Polymer collection"}
    assert pol.Source_URL.eq("https://doi.org/10.1039/d4nr03702b").all()


def test_polymer_values_match_source(pol):
    exp = pd.read_csv(f"{SRC}/conjugated-polymer-1096dp.csv")
    olig = pd.read_csv(f"{SRC}/oligomer-DFT-TDDFT-calculations-1096dp.csv")
    first = pol.iloc[:1096]
    assert (first.Source_Record_ID == ("1096:" + exp.idx.astype(str)).values).all()
    for col, src in (("Experimental_Eg_eV", exp["bandgap(eV)"]), ("Experimental_HOMO_eV", exp["HOMO(eV)"]),
                     ("Experimental_LUMO_eV", exp["LUMO(eV)"]), ("Eg_eV", olig["HOMO-LUMO-B3LYP(eV)"]),
                     ("HOMO_eV", olig["HOMO-B3LYP(eV)"]), ("TDDFT_S1_eV", olig["S1-state-B3LYP(eV)"])):
        assert (pd.to_numeric(first[col]).values - src.values).__abs__().max() < 1e-9, col


def test_every_polymer_record_has_experimental_gap_and_dft_gap(pol):
    assert pd.to_numeric(pol.Experimental_Eg_eV).notna().all()
    assert pd.to_numeric(pol.Eg_eV).notna().all()
    assert pol.InChIKey.notna().all() and pol.Source_Monomer_SMILES.notna().all()
    assert pol.Experimental_Eg_Source.notna().all()
    assert pd.to_numeric(pol.Experimental_HOMO_eV).notna().sum() == 1096


def test_polymer_geometries_linked_to_matching_records(pol):
    with gzip.open(POLY_SDF, "rt") as fh:
        blocks = [b for b in fh.read().split("$$$$\n") if b.strip()]
    keys = {}
    for b in blocks:
        lines = b.splitlines()
        f = {lines[i][3:-1]: lines[i + 1] for i in range(len(lines) - 1) if lines[i].startswith("> <")}
        keys[f["ChalMolDB_Record_ID"]] = f["ChalMolDB_InChIKey"]
    with_3d = pol[pol.Structure_Availability.str.contains("3D coordinates")]
    assert len(keys) == len(with_3d) == 1064
    assert all(keys[r] == k for r, k in zip(with_3d.Record_ID, with_3d.InChIKey))


def test_polymer_collection_has_no_tellurium_and_marks_controls(pol):
    assert pd.to_numeric(pol.Te_Count).fillna(0).eq(0).all()
    no_chalcogen = pd.to_numeric(pol.S_Count).eq(0) & pd.to_numeric(pol.Se_Count).eq(0)
    assert (pol.loc[no_chalcogen, "Scope_Flag"] != "Core_SSeTe").all()
