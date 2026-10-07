import pandas as pd
import pytest

DB = "data/chalcogen_database_v1_master.csv.gz"


@pytest.fixture(scope="module")
def db():
    return pd.read_csv(DB, low_memory=False)


def test_record_counts_unchanged(db):
    assert len(db) == 3360
    assert db.Record_ID.str.startswith("EROL_").sum() == 3088
    assert db.Record_ID.str.startswith("HAKAN_").sum() == 272
    assert db.Record_ID.is_unique


def test_development_donor_acceptor_ids_complete(db):
    dev = db[db.Record_ID.str.startswith("EROL_")]
    assert dev[["Donor_ID", "Acceptor_ID"]].notna().all().all()
    assert dev.Donor_ID.nunique() == 60
    assert dev.Acceptor_ID.nunique() == 52
    assert set(dev.loc[dev.Record_ID.isin(["EROL_1093", "EROL_1249"]), "Acceptor_ID"]) == {"A39"}
    assert dev.Source_URL.eq("https://doi.org/10.1021/acs.jcim.5c00345").all()


def test_no_person_named_collection_labels(db):
    for col in ("Dataset_Owner", "Dataset_Name"):
        values = " ".join(db[col].astype(str)).lower()
        assert "erol" not in values and "hakan" not in values


def test_duplicate_flag_is_shared_inchikey(db):
    shared = db.InChIKey.notna() & db.InChIKey.duplicated(keep=False)
    assert (db.Duplicate_Flag.astype(bool) == shared).all()


def test_external_structures_consistent_with_system_code(db):
    ext = db[db.Record_ID.str.startswith("HAKAN_") & db.InChIKey.notna()]
    assert len(ext) == 112
    for _, r in ext.iterrows():
        n = int(r.Oligomer_n)
        donor, acceptor = r.Donor_Chalcogen, r.Acceptor_Chalcogen
        counts = {e: 0 for e in ("O", "S", "Se", "Te")}
        counts[acceptor] += n
        counts[donor] += 2 * n
        if r.Family.startswith("ethylenedioxy"):
            counts["O"] += 4 * n
        for e, c in counts.items():
            assert int(r[f"{e}_Count"]) == c, (r.Record_ID, e)


def test_every_external_record_has_repeat_unit(db):
    ext = db[db.Record_ID.str.startswith("HAKAN_")]
    assert ext.Repeat_Unit_SMILES.notna().all()
    assert ext.Repeat_Unit_SMILES.str.count(r"\*").eq(2).all()


def test_unverified_b_system_records_are_flagged(db):
    b = db[db.Family.astype(str).str.contains("B system")]
    assert len(b) == 112
    no_structure = b[b.InChIKey.isna() & b.Oligomer_n.ne("polymer")]
    assert len(no_structure) == 80
    assert no_structure.Curation_Note.str.contains("unverified|could not be verified").all()
