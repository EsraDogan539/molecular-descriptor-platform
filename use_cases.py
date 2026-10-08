"""Ready-made example queries for the "What can you do with ChalMolDB?" section.

Each preset is opened with a link such as ?page=database&preset=te_low_gap and fills the Database
workspace (search, structure search, advanced filters) before the widgets are drawn.
"""

PRESETS = {
    "te_low_gap": {
        "title": "Find low-band-gap tellurium systems",
        "text": "Records that contain Te and have a band gap below 1.5 eV.",
        "required_chalcogens": ["Te"],
        "eg_max": 1.5,
    },
    "chalcogen_swap": {
        "title": "See how swapping S, Se and Te changes the gap",
        "text": "One family of donor–acceptor–donor oligomers in which only the chalcogens change; "
                "pick two records in the comparison panel.",
        "text_query": "Molecular Modeling 2024",
    },
    "similar_structures": {
        "title": "Find structures similar to a molecule",
        "text": "Similarity search around 4,7-di(thiophen-2-yl)-2,1,3-benzothiadiazole; draw or paste your own "
                "molecule instead.",
        "structure_query": "c1csc(c1)-c1ccc(-c2cccs2)c2nsnc12",
        "structure_mode": "Similarity",
        "similarity": 0.40,
    },
}


def apply_preset(name, session_state, eg_bounds):
    """Write the widget state for preset `name`; returns True if the preset exists."""
    preset = PRESETS.get(name)
    if preset is None:
        return False
    session_state["database_text_query"] = preset.get("text_query", "")
    session_state["database_structure_query"] = preset.get("structure_query", "")
    if "structure_mode" in preset:
        session_state["database_structure_mode"] = preset["structure_mode"]
        session_state["database_structure_threshold"] = float(preset.get("similarity", 0.40))
    session_state["advanced_required_chalcogens"] = list(preset.get("required_chalcogens", []))
    lower, upper = eg_bounds
    if "eg_max" in preset and lower is not None:
        session_state["advanced_enable_Eg_eV"] = True
        session_state["advanced_min_Eg_eV"] = float(lower)
        session_state["advanced_max_Eg_eV"] = min(float(preset["eg_max"]), float(upper))
    else:
        session_state["advanced_enable_Eg_eV"] = False
    return True
