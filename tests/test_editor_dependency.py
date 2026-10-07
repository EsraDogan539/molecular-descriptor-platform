def test_molecule_editor_component_is_installed():
    """The deployed app offers a drawing editor when streamlit-ketcher is importable."""
    from streamlit_ketcher import st_ketcher

    assert callable(st_ketcher)
