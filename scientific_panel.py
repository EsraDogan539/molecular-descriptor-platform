import streamlit as st


SCIENTIFIC_COLUMNS = [
    "Molecule_ID",
    "Canonical SMILES",
    "InChIKey",
    "Chalcogen Type",
    "Heavy Chalcogen Count",
    "Chalcogen Fraction",
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
        st.info("Scientific Core için gösterilecek geçerli molekül bulunamadı.")
        return

    total = len(valid_df)
    s_count = int((valid_df["Chalcogen Type"] == "S").sum()) if "Chalcogen Type" in valid_df else 0
    se_count = int((valid_df["Chalcogen Type"] == "Se").sum()) if "Chalcogen Type" in valid_df else 0
    te_count = int((valid_df["Chalcogen Type"] == "Te").sum()) if "Chalcogen Type" in valid_df else 0
    duplicate_count = int(valid_df.get("Duplicate Flag", False).sum()) if "Duplicate Flag" in valid_df else 0

    c1, c2, c3, c4, c5 = st.columns(5)
    c1.metric("Molecules", total)
    c2.metric("S", s_count)
    c3.metric("Se", se_count)
    c4.metric("Te", te_count)
    c5.metric("Duplicates", duplicate_count)

    existing = [c for c in SCIENTIFIC_COLUMNS if c in valid_df.columns]
    st.dataframe(valid_df[existing], use_container_width=True, hide_index=True)

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
        "Duplicate Flag is based on canonical SMILES. Oxygen is retained as a general "
        "elemental descriptor but is not counted in the S/Se/Te-focused chalcogen total."
    )
