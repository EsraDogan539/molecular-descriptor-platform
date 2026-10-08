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


def test_record_deep_link_selects_record_and_3d_viewer_renders():
    at = _database_page()
    at.query_params["record"] = "EXT_0148"
    at.run()
    assert not at.exception, at.exception
    assert at.session_state["database_text_query"] == "EXT_0148"
    viewers = [t for t in at.toggle if t.key == "structure_3d_view_EXT_0148"]
    assert viewers, "3D viewer toggle not rendered for a record with coordinates"
    viewers[0].set_value(True).run()
    assert not at.exception, at.exception


def test_documentation_page_renders_data_access():
    at = AppTest.from_file(str(APP), default_timeout=120)
    at.query_params["page"] = "documentation"
    at.run()
    assert not at.exception, at.exception
    assert any("Data access" in m.value for m in at.markdown)


def test_embedded_mode_links_keep_embed_flag():
    at = AppTest.from_file(str(APP), default_timeout=120)
    at.query_params["page"] = "about"
    at.query_params["embedded"] = "1"
    at.run()
    assert not at.exception, at.exception
    html = " ".join(m.value for m in at.markdown)
    assert "page=database&embedded=1&embed=true" in html
    assert "page=about&embedded=1&embed=true#cite" in html


def _result_count(at):
    import re
    for m in at.markdown:
        hit = re.search(r"([\d,]+) matching record", m.value)
        if hit:
            return int(hit.group(1).replace(",", ""))
    return None


def test_use_case_presets_load_results():
    from use_cases import PRESETS

    for name in PRESETS:
        at = _database_page()
        at.query_params["preset"] = name
        at.run()
        assert not at.exception, (name, at.exception)
        count = _result_count(at)
        assert count is not None and 0 < count < 3488, (name, count)


def test_home_shows_use_cases():
    at = AppTest.from_file(str(APP), default_timeout=120)
    at.query_params["page"] = "home"
    at.run()
    assert not at.exception, at.exception
    html = " ".join(m.value for m in at.markdown)
    assert "What you can do with ChalMolDB" in html and "preset=te_low_gap" in html
