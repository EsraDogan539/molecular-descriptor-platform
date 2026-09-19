import os
import pandas as pd
import streamlit as st

from rdkit import Chem
from rdkit.Chem import Draw


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


def _has_value(series):
    return series.notna() & series.astype(str).str.strip().ne("")


def _record_label(row):
    record_id = str(row.get("Record_ID", "Record"))
    system = row.get("System_Code")
    molecule = row.get("Molecule_Name")
    suffix = system if pd.notna(system) and str(system).strip() else molecule
    return f"{record_id} — {suffix}" if pd.notna(suffix) and str(suffix).strip() else record_id


def _render_record_detail(row):
    st.markdown("#### Record detail")

    left, right = st.columns([1, 2])

    with left:
        smiles = row.get("Canonical_SMILES")
        if pd.notna(smiles) and str(smiles).strip():
            mol = Chem.MolFromSmiles(str(smiles))
            if mol is not None:
                st.image(
                    Draw.MolToImage(mol, size=(440, 320)),
                    use_container_width=True,
                )
            else:
                st.info("Canonical SMILES is present but could not be rendered.")
        else:
            st.info("Exact standardized structure is not available for this record.")

        st.caption(f"Structure availability: {row.get('Structure_Availability', '—')}")
        st.caption(f"Curation status: {row.get('Curation_Status', '—')}")

    with right:
        id_cols = st.columns(3)
        id_cols[0].metric("Record ID", str(row.get("Record_ID", "—")))
        id_cols[1].metric("Dataset role", str(row.get("Split_Role", "—")))
        id_cols[2].metric("Scope", str(row.get("Scope_Flag", "—")))

        property_cols = st.columns(4)
        property_cols[0].metric("HOMO (eV)", row.get("HOMO_eV", "—"))
        property_cols[1].metric("LUMO (eV)", row.get("LUMO_eV", "—"))
        property_cols[2].metric("Eg (eV)", row.get("Eg_eV", "—"))
        property_cols[3].metric("Exp. Eg (eV)", row.get("Experimental_Eg_eV", "—"))

        details = {
            "Dataset source": row.get("Dataset_Owner"),
            "Molecule / system": row.get("Molecule_Name") or row.get("System_Code"),
            "Donor / Acceptor": " / ".join(
                [
                    str(v)
                    for v in [row.get("Donor_ID"), row.get("Acceptor_ID")]
                    if pd.notna(v) and str(v).strip()
                ]
            ) or "—",
            "Unit type": row.get("Unit_Type"),
            "Oligomer n": row.get("Oligomer_n"),
            "Chalcogen type": row.get("Chalcogen_Type"),
            "InChIKey": row.get("InChIKey"),
            "Method": row.get("Method"),
            "Basis set": row.get("Basis_Set"),
            "Conditions": row.get("Solvent_or_Conditions"),
            "Reference": row.get("Reference"),
        }
        detail_df = pd.DataFrame(
            [{"Field": key, "Value": value if pd.notna(value) else "—"} for key, value in details.items()]
        )
        st.dataframe(detail_df, use_container_width=True, hide_index=True)


def display_database_browser():
    st.subheader("Curated Chalcogen Database")
    st.caption(
        "Publication-oriented database browser with explicit dataset roles, "
        "structure availability, provenance and property completeness."
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
    m4.metric("External collection", f"{validation_expected:,}")

    completeness_cols = st.columns(3)
    structure_count = int(_has_value(df["Canonical_SMILES"]).sum()) if "Canonical_SMILES" in df.columns else 0
    eg_count = int(df["Eg_eV"].notna().sum()) if "Eg_eV" in df.columns else 0
    provenance_count = (
        int(_has_value(df["Reference"]).sum()) if "Reference" in df.columns else 0
    )
    completeness_cols[0].metric("Structures available", f"{structure_count:,}")
    completeness_cols[1].metric("Eg available", f"{eg_count:,}")
    completeness_cols[2].metric("Reference metadata", f"{provenance_count:,}")

    with st.expander("Database scope and curation policy", expanded=False):
        st.markdown(
            "- **Development / Training:** Erol dataset\n"
            "- **External validation collection:** Hakan Kayı dataset\n"
            "- **Core scope:** records containing S, Se or Te\n"
            "- O-only and non-S/Se/Te records are retained as controls.\n"
            "- Missing structures or properties remain explicit; they are not inferred in the curated release.\n"
            "- Repeated standardized structures are flagged rather than silently removed."
        )

    st.markdown("#### Search and filter")
    f1, f2, f3 = st.columns(3)

    split_options = ["All"] + _safe_unique(df, "Split_Role")
    scope_options = ["All"] + _safe_unique(df, "Scope_Flag")
    chalcogen_options = ["All"] + _safe_unique(df, "Chalcogen_Type")

    with f1:
        split_value = st.selectbox("Dataset role", split_options)
    with f2:
        scope_value = st.selectbox("Chemical scope", scope_options)
    with f3:
        chalcogen_value = st.selectbox("Chalcogen type", chalcogen_options)

    search_text = st.text_input(
        "Search by Record ID, molecule/system, donor/acceptor or InChIKey",
        value="",
    ).strip()

    filtered = df.copy()

    if split_value != "All" and "Split_Role" in filtered.columns:
        filtered = filtered[filtered["Split_Role"].astype(str) == split_value]
    if scope_value != "All" and "Scope_Flag" in filtered.columns:
        filtered = filtered[filtered["Scope_Flag"].astype(str) == scope_value]
    if chalcogen_value != "All" and "Chalcogen_Type" in filtered.columns:
        filtered = filtered[filtered["Chalcogen_Type"].astype(str) == chalcogen_value]

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

    st.caption(f"Showing {len(filtered):,} record(s)")

    preferred_columns = [
        "Record_ID", "Split_Role", "Scope_Flag", "Molecule_Name",
        "Donor_ID", "Acceptor_ID", "System_Code", "Unit_Type",
        "Chalcogen_Type", "HOMO_eV", "LUMO_eV", "Eg_eV",
        "Experimental_Eg_eV", "Structure_Availability",
        "Duplicate_Flag", "Curation_Status",
    ]
    existing_columns = [c for c in preferred_columns if c in filtered.columns]

    st.dataframe(
        filtered[existing_columns],
        use_container_width=True,
        hide_index=True,
        height=460,
    )

    if not filtered.empty:
        record_labels = [_record_label(row) for _, row in filtered.iterrows()]
        selected_label = st.selectbox(
            "Inspect a record",
            record_labels,
            index=0,
        )
        selected_position = record_labels.index(selected_label)
        selected_row = filtered.iloc[selected_position]
        _render_record_detail(selected_row)

    if "Split_Role" in filtered.columns:
        with st.expander("Current view summary", expanded=False):
            group_cols = [
                c for c in ["Split_Role", "Scope_Flag"]
                if c in filtered.columns
            ]
            summary = (
                filtered.groupby(group_cols, dropna=False)
                .size()
                .reset_index(name="Count")
            )
            st.dataframe(summary, use_container_width=True, hide_index=True)

    if not is_preview:
        st.download_button(
            "Download filtered database CSV",
            filtered.to_csv(index=False).encode("utf-8"),
            "chalcogen_database_filtered.csv",
            "text/csv",
            use_container_width=True,
        )

    if source_path:
        st.caption(f"Database source in app build: {source_path}")
