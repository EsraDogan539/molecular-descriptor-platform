import os
import pandas as pd
import streamlit as st


FULL_DATABASE_PATHS = [
    "data/chalcogen_database_v1_master.csv.gz",
    "data/chalcogen_database_v1_master.csv",
]
PREVIEW_DATABASE_PATH = "data/chalcogen_database_preview_v1.csv"


@st.cache_data(show_spinner=False)
def load_curated_database():
    for path in FULL_DATABASE_PATHS:
        if os.path.exists(path):
            compression = "gzip" if path.endswith(".gz") else "infer"
            return pd.read_csv(path, compression=compression), False, path

    if os.path.exists(PREVIEW_DATABASE_PATH):
        return pd.read_csv(PREVIEW_DATABASE_PATH), True, PREVIEW_DATABASE_PATH

    return pd.DataFrame(), True, None


def _safe_unique(df, column):
    if column not in df.columns:
        return []
    return sorted([str(v) for v in df[column].dropna().unique()])


def display_database_browser():
    st.subheader("Curated Chalcogen Database")
    st.caption(
        "Publication-oriented database browser. Development/training and independent "
        "external-validation records are kept as separate roles."
    )

    df, is_preview, source_path = load_curated_database()

    if df.empty:
        st.warning("Curated database file is not available in this build yet.")
        return

    if is_preview:
        st.info(
            "Bu geliştirme dalında veri tabanının küçük bir önizlemesi gösteriliyor. "
            "Tam v1 master dosyası bağlandığında aynı tarayıcı 3360 kaydın tamamını kullanacak."
        )

    total_expected = 3360
    core_expected = 3145
    dev_expected = 3088
    validation_expected = 272

    m1, m2, m3, m4 = st.columns(4)
    m1.metric("Database v1", f"{total_expected:,}")
    m2.metric("Core S/Se/Te", f"{core_expected:,}")
    m3.metric("Development", f"{dev_expected:,}")
    m4.metric("External validation", f"{validation_expected:,}")

    with st.expander("Database scope and split policy", expanded=False):
        st.markdown(
            "- **Development / Training:** Erol dataset\n"
            "- **External Validation Collection:** Hakan Kayı dataset\n"
            "- **Core scope:** records containing S, Se or Te\n"
            "- O-only and non-S/Se/Te records are retained as controls, not merged into the core benchmark.\n"
            "- External-validation records are not used for model fitting or feature selection."
        )

    st.markdown("#### Search and filter")
    f1, f2, f3, f4 = st.columns(4)

    split_options = ["All"] + _safe_unique(df, "Split_Role")
    scope_options = ["All"] + _safe_unique(df, "Scope_Flag")
    chalcogen_options = ["All"] + _safe_unique(df, "Chalcogen_Type")
    owner_options = ["All"] + _safe_unique(df, "Dataset_Owner")

    with f1:
        split_value = st.selectbox("Dataset role", split_options)
    with f2:
        scope_value = st.selectbox("Scope", scope_options)
    with f3:
        chalcogen_value = st.selectbox("Chalcogen type", chalcogen_options)
    with f4:
        owner_value = st.selectbox("Dataset source", owner_options)

    search_text = st.text_input(
        "Search Record ID, molecule/system code, donor/acceptor or InChIKey",
        value="",
    ).strip()

    filtered = df.copy()

    if split_value != "All" and "Split_Role" in filtered.columns:
        filtered = filtered[filtered["Split_Role"].astype(str) == split_value]
    if scope_value != "All" and "Scope_Flag" in filtered.columns:
        filtered = filtered[filtered["Scope_Flag"].astype(str) == scope_value]
    if chalcogen_value != "All" and "Chalcogen_Type" in filtered.columns:
        filtered = filtered[filtered["Chalcogen_Type"].astype(str) == chalcogen_value]
    if owner_value != "All" and "Dataset_Owner" in filtered.columns:
        filtered = filtered[filtered["Dataset_Owner"].astype(str) == owner_value]

    if search_text:
        searchable_columns = [
            c for c in [
                "Record_ID", "Molecule_Name", "Donor_ID", "Acceptor_ID",
                "System_Code", "InChIKey", "Canonical_SMILES"
            ] if c in filtered.columns
        ]
        if searchable_columns:
            mask = pd.Series(False, index=filtered.index)
            for column in searchable_columns:
                mask = mask | filtered[column].astype(str).str.contains(
                    search_text, case=False, na=False, regex=False
                )
            filtered = filtered[mask]

    st.caption(f"Gösterilen kayıt: {len(filtered)}")

    preferred_columns = [
        "Record_ID", "Dataset_Owner", "Split_Role", "Scope_Flag",
        "Molecule_Name", "Donor_ID", "Acceptor_ID", "System_Code",
        "Unit_Type", "Oligomer_n", "Chalcogen_Type",
        "HOMO_eV", "LUMO_eV", "Eg_eV", "Experimental_Eg_eV",
        "S_Count", "Se_Count", "Te_Count", "InChIKey",
        "Method", "Basis_Set", "Curation_Status",
    ]
    existing_columns = [c for c in preferred_columns if c in filtered.columns]

    st.dataframe(
        filtered[existing_columns],
        use_container_width=True,
        hide_index=True,
        height=520,
    )

    if "Split_Role" in filtered.columns:
        st.markdown("#### Current view summary")
        summary = (
            filtered.groupby([c for c in ["Split_Role", "Scope_Flag"] if c in filtered.columns], dropna=False)
            .size()
            .reset_index(name="Count")
        )
        st.dataframe(summary, use_container_width=True, hide_index=True)

    if not is_preview:
        st.download_button(
            "Filtered database CSV indir",
            filtered.to_csv(index=False).encode("utf-8"),
            "chalcogen_database_filtered.csv",
            "text/csv",
            use_container_width=True,
        )

    if source_path:
        st.caption(f"Database source in app build: {source_path}")
