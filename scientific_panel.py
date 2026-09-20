import streamlit as st


SCIENTIFIC_COLUMNS = [
    "Molecule_ID",
    "Canonical SMILES",
    "InChIKey",
    "Chalcogen Type",
    "Target Chalcogen Count",
    "Target Chalcogen Fraction",
    "Aromatic Chalcogen Count",
    "NonAromatic Chalcogen Count",
    "Mixed Chalcogen Flag",
    "Ring Incorporated Chalcogen Count",
    "Chalcogen-C Bond Count",
    "Chalcogen-Heteroatom Bond Count",
    "Aromatic Neighbor Count",
    "Conjugated Bond Count",
    "Aromatic Bond Fraction",
    "Conjugated Atom Fraction",
    "Heteroaromatic Ring Count",
    "Duplicate Flag",
]


def display_scientific_core_panel(valid_df):
    st.subheader("Scientific Core")
    st.caption(
        "Chalcogen-focused identity, structural context, conjugation, "
        "and database-quality fields for publication-oriented analysis."
    )

    if valid_df.empty:
        st.info("No valid molecules are available for the descriptor summary.")
        return

    total = len(valid_df)
    s_count = int((valid_df["Sulfur Count"] > 0).sum()) if "Sulfur Count" in valid_df else 0
    se_count = int((valid_df["Selenium Count"] > 0).sum()) if "Selenium Count" in valid_df else 0
    te_count = int((valid_df["Tellurium Count"] > 0).sum()) if "Tellurium Count" in valid_df else 0
    duplicate_count = int(valid_df.get("Duplicate Flag", False).sum()) if "Duplicate Flag" in valid_df else 0

    c1, c2, c3, c4, c5 = st.columns(5)
    c1.metric("Molecules", total)
    c2.metric("Contains S", s_count)
    c3.metric("Contains Se", se_count)
    c4.metric("Contains Te", te_count)
    c5.metric("Repeated", duplicate_count)

    existing = [c for c in SCIENTIFIC_COLUMNS if c in valid_df.columns]
    display_df = valid_df[existing].rename(columns={
        "Molecule_ID": "Molecule ID",
        "NonAromatic Chalcogen Count": "Non-aromatic Chalcogen Count",
        "Mixed Chalcogen Flag": "Mixed Chalcogen",
        "Ring Incorporated Chalcogen Count": "Ring-incorporated Chalcogen Count",
        "Chalcogen-C Bond Count": "Chalcogen–C Bond Count",
        "Chalcogen-Heteroatom Bond Count": "Chalcogen–Heteroatom Bond Count",
        "Duplicate Flag": "Repeated Structure",
    })
    st.dataframe(display_df, use_container_width=True, hide_index=True)

    if "Chalcogen Type" in valid_df.columns:
        chalcogen_summary = (
            valid_df["Chalcogen Type"]
            .value_counts(dropna=False)
            .rename_axis("Chalcogen Type")
            .reset_index(name="Count")
        )
        st.markdown("#### Chalcogen distribution")
        st.dataframe(chalcogen_summary, use_container_width=True, hide_index=True)

    st.markdown("#### Quality notes")
    st.caption(
        "Repeated Structure is based on standardized structure identity (InChIKey when available, "
        "with canonical SMILES fallback). Oxygen is retained as a general elemental descriptor but "
        "is not included in the S/Se/Te-focused target chalcogen count."
    )
