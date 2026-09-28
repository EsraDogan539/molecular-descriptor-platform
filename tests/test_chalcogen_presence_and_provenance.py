"""Regression tests for chalcogen presence filtering and provenance fields.

- chalcogen presence filtering covers metadata-level (external) records
- NaN-safe "Molecule / system" fallback
- DOI resolution from dedicated DOI/URL fields
"""

import math

import pandas as pd

from database_filters import apply_advanced_filters, chalcogen_presence_mask
from provenance_export import citation_ready_record, first_meaningful, record_doi


def _toy_df():
    nan = float("nan")
    return pd.DataFrame(
        {
            "Record_ID": ["DEV_1", "DEV_2", "EXT_1", "EXT_2", "EXT_3"],
            "S_Count": [2, 0, nan, nan, nan],
            "Se_Count": [1, 0, nan, nan, nan],
            "Te_Count": [0, 0, nan, nan, nan],
            "Donor_Chalcogen": [nan, nan, "Te", "S", "O"],
            "Acceptor_Chalcogen": [nan, nan, "Se", "S", "O"],
            "Chalcogen_Type": ["Mixed", nan, "Mixed", "S", nan],
        }
    )


def test_presence_uses_counts_and_metadata():
    df = _toy_df()
    assert chalcogen_presence_mask(df, "S").tolist() == [True, False, False, True, False]
    assert chalcogen_presence_mask(df, "Se").tolist() == [True, False, True, False, False]
    assert chalcogen_presence_mask(df, "Te").tolist() == [False, False, True, False, False]


def test_s_token_does_not_match_se():
    df = pd.DataFrame({"Donor_Chalcogen": ["Se"], "Acceptor_Chalcogen": ["Se"]})
    assert not chalcogen_presence_mask(df, "S").iloc[0]


def test_required_chalcogens_filter_combines_with_and():
    df = _toy_df()
    result = apply_advanced_filters(df, required_chalcogens=["Se", "Te"])
    assert result["Record_ID"].tolist() == ["EXT_1"]


def test_required_chalcogens_empty_is_noop():
    df = _toy_df()
    assert len(apply_advanced_filters(df, required_chalcogens=[])) == len(df)


def test_first_meaningful_skips_nan():
    assert first_meaningful(float("nan"), "SeTeSe") == "SeTeSe"
    assert first_meaningful(None, "  ", "x") == "x"
    assert first_meaningful(float("nan"), None) is None


def test_record_doi_falls_back_to_source_url():
    row = pd.Series(
        {
            "Reference": "Ozkilinc & Kayi, Journal of Molecular Modeling 2019, Table 4",
            "Source_URL": "https://doi.org/10.1007/s00894-019-4043-2",
        }
    )
    assert record_doi(row, row["Reference"]) == "10.1007/s00894-019-4043-2"


def test_citation_record_uses_system_code_and_source_url_doi():
    row = pd.Series(
        {
            "Record_ID": "HAKAN_0084",
            "Split_Role": "External Validation",
            "Molecule_Name": float("nan"),
            "System_Code": "SeTeSe",
            "Reference": "Kayi, Sen & Ozkilinc, Journal of Molecular Modeling 2024, Table 4",
            "Source_URL": "https://doi.org/10.1007/s00894-024-05985-2",
        }
    )
    table = citation_ready_record(row, "EXT_0084").set_index("Field")["Value"]
    assert table["Molecule / system"] == "SeTeSe"
    assert table["DOI"] == "10.1007/s00894-024-05985-2"
    assert table["DOI URL"] == "https://doi.org/10.1007/s00894-024-05985-2"


def test_release_database_chalcogen_presence_counts():
    """Filter counts must match the published statistics (2,907 / 644 / 119)."""
    df = pd.read_csv("data/chalcogen_database_v1_master.csv.gz", low_memory=False)
    assert int(chalcogen_presence_mask(df, "S").sum()) == 2907
    assert int(chalcogen_presence_mask(df, "Se").sum()) == 644
    assert int(chalcogen_presence_mask(df, "Te").sum()) == 119
    assert len(apply_advanced_filters(df, required_chalcogens=["Te"])) == 119
