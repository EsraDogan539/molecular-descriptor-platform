"""The headline numbers hard-coded in the app must match the database table."""
import re
from pathlib import Path

import pandas as pd

DB = pd.read_csv("data/chalcogen_database_v1_master.csv.gz", low_memory=False)
APP = (Path(__file__).resolve().parents[1] / "app_v08.py").read_text(encoding="utf-8")


def _shown(label):
    values = re.findall(r'metric-number">([\d,]+)</span><span class="metric-label">' + label, APP)
    assert values, label
    return {int(v.replace(",", "")) for v in values}


def test_headline_metrics_match_database():
    assert _shown("Records") == {len(DB)}
    assert _shown("Unique structures") == {DB.InChIKey.nunique()}
    assert _shown("S/Se/Te records") == {int((DB.Scope_Flag == "Core_SSeTe").sum())}
    assert _shown("Core S/Se/Te") == {int((DB.Scope_Flag == "Core_SSeTe").sum())}
    assert _shown("External records") == {int(DB.Record_ID.str.startswith("HAKAN_").sum())}
    assert _shown("Eg values") == {int(DB.Eg_eV.notna().sum())}
