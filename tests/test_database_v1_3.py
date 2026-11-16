"""Database v1.3: polymers published 2015-2026 and literature model compounds are added; v1.2 values are unchanged."""

import pandas as pd
import pytest

DB = "data/chalcogen_database_v1_master.csv.gz"
V12 = "data/releases/chalcogen_database_v1_2.csv.gz"
SRC = "data/sources/dogan_kayi_2026_external"
UNVERIFIED = {"POL_1116", "POL_1117", "POL_1118", "POL_1119", "POL_1122", "POL_1123", "POL_1124", "POL_1125",
              "POL_1135", "POL_1141", "POL_1142", "POL_1143", "POL_1145", "POL_1147", "POL_1148"}


@pytest.fixture(scope="module")
def db():
    return pd.read_csv(DB, dtype=str, keep_default_na=False, na_values=[""])


@pytest.fixture(scope="module")
def v12():
    return pd.read_csv(V12, dtype=str, keep_default_na=False, na_values=[""])


def test_v1_2_records_unchanged_except_notes(db, v12):
    head = db.iloc[: len(v12)].reset_index(drop=True)
    assert head.Record_ID.tolist() == v12.Record_ID.tolist()
    diff = head.fillna("") != v12.fillna("")
    changed = {c: set(head.Record_ID[diff[c]]) for c in diff.columns if diff[c].any()}
    assert set(changed) <= {"Curation_Note", "Duplicate_Flag"}
    assert changed["Curation_Note"] == UNVERIFIED
    assert len(changed.get("Duplicate_Flag", set())) <= 5
    assert head.loc[head.Record_ID.isin(UNVERIFIED), "Curation_Note"].str.contains("could not be verified").all()


def test_added_records(db):
    assert len(db) == 5276 and db.Record_ID.is_unique
    new = db.iloc[5029:]
    pol = new[new.Record_ID.str.startswith("POL_")]
    lit = new[new.Record_ID.str.startswith("LIT_")]
    assert pol.Record_ID.tolist() == [f"POL_{i:04d}" for i in range(1324, 1521)]
    assert lit.Record_ID.tolist() == [f"LIT_{i:04d}" for i in range(219, 269)]
    assert pd.to_numeric(pol.Experimental_Eg_eV).notna().all() and pol.Eg_eV.isna().all()
    assert pol.Source_URL.str.startswith("https://doi.org/").all()
    assert pol.Solvent_or_Conditions.str.contains("film").sum() == 185
    assert pd.to_numeric(lit.Eg_eV).notna().all() and lit.Method.str.startswith("B3LYP").all()
    assert lit.Source_URL.str.startswith("https://doi.org/").all()
    assert new.InChIKey.notna().all()


def test_added_values_match_source(db):
    ext = pd.read_csv(f"{SRC}/external_polymers_2015_2026.csv")
    ext = ext[ext.backbone != "in ChalMolDB"]
    pol = db[db.Source_File.eq("external_polymers_2015_2026.csv")]
    assert pol.Source_Record_ID.tolist() == ext.id.tolist()
    assert (abs(pd.to_numeric(pol.Experimental_Eg_eV).values - ext.measured_optical_gap_eV.values) < 1e-9).all()
    lit = pd.read_csv(f"{SRC}/literature_model_compounds.csv")
    lit = lit[lit.used_in_model]
    rec = db[db.Source_File.eq("literature_model_compounds.csv")]
    assert rec.Source_Record_ID.tolist() == lit.id.tolist()
    assert (abs(pd.to_numeric(rec.Eg_eV).values - lit.gap_eV.values) < 1e-9).all()
