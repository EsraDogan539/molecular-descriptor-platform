import pandas as pd
import pytest

DB = "data/chalcogen_database_v1_master.csv.gz"


@pytest.fixture(scope="module")
def db():
    return pd.read_csv(DB, low_memory=False)


def test_record_counts_unchanged(db):
    assert len(db) == 3248
    assert db.Record_ID.str.startswith("EROL_").sum() == 3088
    assert db.Record_ID.str.startswith("HAKAN_").sum() == 160
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
    assert len(ext) == 96
    for _, r in ext.iterrows():
        n = int(r.Oligomer_n)
        donor, acceptor = r.Donor_Chalcogen, r.Acceptor_Chalcogen
        counts = {e: 0 for e in ("O", "S", "Se", "Te")}
        counts[acceptor] += n
        counts[donor] += 2 * n
        for e, c in counts.items():
            assert int(r[f"{e}_Count"]) == c, (r.Record_ID, e)


def test_every_external_record_has_repeat_unit(db):
    ext = db[db.Record_ID.str.startswith("HAKAN_")]
    assert ext.Repeat_Unit_SMILES.notna().all()
    assert ext.Repeat_Unit_SMILES.str.count(r"\*").eq(2).all()


def test_b_system_removed_without_renumbering(db):
    assert not db.Family.astype(str).str.contains("B_system|B system").any()
    v1 = pd.read_csv("data/releases/chalcogen_database_v1.csv.gz", low_memory=False)
    kept = v1[~v1.Family.astype(str).str.startswith("B_system")]
    assert list(kept.Record_ID) == list(db.Record_ID)
    merged = kept.merge(db, on="Record_ID", suffixes=("_v1", "_v11"))
    for col in ("Eg_eV", "HOMO_eV", "LUMO_eV", "Experimental_Eg_eV"):
        a, b = merged[f"{col}_v1"], merged[f"{col}_v11"]
        if col == "Experimental_Eg_eV":
            assert (a.fillna("").astype(str) == b.fillna("").astype(str)).all(), col
        else:
            assert ((a - b).abs().lt(1e-9) | (a.isna() & b.isna())).all(), col


def test_source_annotations(db):
    ext = db[db.Record_ID.str.startswith("HAKAN_")]
    assert not db.Notes.astype(str).str.contains("Hakan|Erol").any()
    assert set(ext.Solvent_or_Conditions) == {"gas phase", "PCM (acetonitrile)"}
    assert (ext.Solvent_or_Conditions == "PCM (acetonitrile)").sum() == 16
    exp = db[db.Experimental_Eg_eV.notna()]
    assert len(exp) == 21
    assert exp.Experimental_Eg_Min_eV.notna().all() and exp.Experimental_Eg_Source.notna().all()
    assert (exp.Experimental_Eg_Min_eV <= exp.Experimental_Eg_Max_eV).all()
    no_value = db[db.Eg_eV.isna()]
    assert len(no_value) == 7
    assert no_value.Curation_Note.str.contains("6-31G\\(d\\) basis set is not defined for Te").all()


def test_annotation_is_idempotent(db):
    import sys
    sys.path.insert(0, "scripts")
    from annotate_external_sources import annotate

    again = annotate(db)
    assert again.to_csv(index=False) == db.to_csv(index=False)
