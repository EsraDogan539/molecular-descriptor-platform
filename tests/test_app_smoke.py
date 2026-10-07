"""Headless smoke test of the Database page, with and without the molecule editor."""
from pathlib import Path

from streamlit.testing.v1 import AppTest

APP = Path(__file__).resolve().parents[1] / "app_v08.py"


def _database_page():
    at = AppTest.from_file(str(APP), default_timeout=120)
    at.query_params["page"] = "database"
    return at


def test_database_page_renders_without_exceptions():
    at = _database_page().run()
    assert not at.exception, at.exception


def test_draw_structure_toggle_renders_without_exceptions():
    at = _database_page().run()
    toggles = [t for t in at.toggle if t.key == "database_draw_structure"]
    assert toggles, "Draw structure toggle not rendered"
    toggles[0].set_value(True).run()
    assert not at.exception, at.exception
