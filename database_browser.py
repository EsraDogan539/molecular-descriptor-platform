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


def _display_value(value, decimals=None):
    if value is None or pd.isna(value) or str(value).strip().lower() in {"", "nan", "none"}:
        return "—"
    if decimals is not None:
        try:
            return f"{float(value):.{decimals}f}"
        except (TypeError, ValueError):
            pass
    return str(value)


def _display_role(value):
    text = _display_value(value)
    if text == "Development/Training":
        return "Development"
    if text == "External Validation":
        return "External"
    return text


def _display_scope(value):
    mapping = {
        "Core_SSeTe": "Core S/Se/Te",
        "Control_NonSSeTe": "Non-core control",
        "O_Control": "O-only control",
    }
    return mapping.get(_display_value(value), _display_value(value))


def _display_collection(row):
    role = _display_role(row.get("Split_Role"))
    if role == "Development":
        return "Development collection"
    if role == "External":
        return "External validation collection"
    return role


def _display_record_id(value):
    text = _display_value(value)
    if text.startswith("EROL_"):
        return text.replace("EROL_", "DEV_", 1)
    if text.startswith("HAKAN_"):
        return text.replace("HAKAN_", "EXT_", 1)
    return text


def _record_label(row):
    record_id = _display_record_id(row.get("Record_ID", "Record"))
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
        identity_items = [
            ("Database ID", _display_record_id(row.get("Record_ID"))),
            ("Collection role", _display_role(row.get("Split_Role"))),
            ("Scope", _display_scope(row.get("Scope_Flag"))),
        ]
        for col, (label, value) in zip(id_cols, identity_items):
            with col:
                with st.container(border=True):
                    st.caption(label)
                    st.markdown(f"**{value}**")

        property_cols = st.columns(4)
        property_cols[0].metric("HOMO (eV)", _display_value(row.get("HOMO_eV"), 3))
        property_cols[1].metric("LUMO (eV)", _display_value(row.get("LUMO_eV"), 3))
        property_cols[2].metric("Eg (eV)", _display_value(row.get("Eg_eV"), 3))
        property_cols[3].metric("Exp. Eg (eV)", _display_value(row.get("Experimental_Eg_eV"), 3))

        details = {
            "Collection": _display_collection(row),
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
            [{"Field": key, "Value": _display_value(value)} for key, value in details.items()]
        )
        st.dataframe(
            detail_df,
            use_container_width=True,
            hide_index=True,
            height=360,
        )


def display_database_browser():
    st.subheader("Curated Chalcogen Database")
    st.caption(
        "Browse the versioned publication dataset by collection role, chemical scope, "
        "chalcogen class and molecular identity."
    )

    df, is_preview, source_path = load_curated_database()

    if df.empty:
        st.warning("Curated database file is not available in this build yet.")
        return

    if is_preview:
        st.info(
            "Preview data are loaded in this build. Publication metrics below refer to database v1."
        )

    total_records = len(df)
    core_count = int((df["Scope_Flag"].astype(str) == "Core_SSeTe").sum()) if "Scope_Flag" in df.columns else 0
    development_count = int((df["Split_Role"].astype(str) == "Development/Training").sum()) if "Split_Role" in df.columns else 0
    external_count = total_records - development_count if "Split_Role" in df.columns else 0

    st.caption("Database version: v1")
    m1, m2, m3, m4 = st.columns(4)
    m1.metric("Records", f"{total_records:,}")
    m2.metric("Core S/Se/Te", f"{core_count:,}")
    m3.metric("Development", f"{development_count:,}")
    m4.metric("External collection", f"{external_count:,}")

    structure_count = int(_has_value(df["Canonical_SMILES"]).sum()) if "Canonical_SMILES" in df.columns else 0
    eg_count = int(df["Eg_eV"].notna().sum()) if "Eg_eV" in df.columns else 0
    unique_structure_count = int(df["InChIKey"].dropna().astype(str).replace("", pd.NA).dropna().nunique()) if "InChIKey" in df.columns else 0

    completeness_cols = st.columns(3)
    completeness_cols[0].metric("Structures available", f"{structure_count:,}")
    completeness_cols[1].metric("Unique structures (development)", f"{unique_structure_count:,}")
    completeness_cols[2].metric("Eg available", f"{eg_count:,}")

    with st.expander("Database scope and curation policy", expanded=False):
        st.markdown(
            "- **Development collection:** structure-complete records used for database development and downstream modelling examples\n"
            "- **External validation collection:** separately curated literature/system-level records retained as an external collection\n"
            "- **Core scope:** records containing S, Se or Te\n"
            "- O-only and non-S/Se/Te records are retained as controls.\n"
            "- Missing structures and properties remain explicit.\n"
            "- Repeated standardized structures are retained with duplicate flags to preserve provenance."
        )

    st.markdown("#### Explore records")
    f1, f2, f3 = st.columns(3)

    split_options = ["All"] + _safe_unique(df, "Split_Role")
    scope_options = ["All"] + _safe_unique(df, "Scope_Flag")
    chalcogen_options = ["All"] + _safe_unique(df, "Chalcogen_Type")

    with f1:
        split_value = st.selectbox(
            "Collection role",
            split_options,
            format_func=lambda value: "All" if value == "All" else _display_role(value),
        )
    with f2:
        scope_value = st.selectbox(
            "Chemical scope",
            scope_options,
            format_func=lambda value: "All" if value == "All" else _display_scope(value),
        )
    with f3:
        chalcogen_value = st.selectbox("Chalcogen type", chalcogen_options)

    search_text = st.text_input(
        "Search by Database ID, molecule/system, donor/acceptor or InChIKey",
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
            if "Record_ID" in filtered.columns:
                public_ids = filtered["Record_ID"].map(_display_record_id)
                mask = mask | public_ids.str.contains(
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
    display_table = filtered[existing_columns].copy()
    display_table = display_table.rename(columns={
        "Record_ID": "Database ID",
        "Split_Role": "Collection role",
        "Scope_Flag": "Scope",
        "Molecule_Name": "Molecule / system",
        "Donor_ID": "Donor",
        "Acceptor_ID": "Acceptor",
        "System_Code": "System",
        "Unit_Type": "Unit",
        "Chalcogen_Type": "Chalcogen",
        "HOMO_eV": "HOMO (eV)",
        "LUMO_eV": "LUMO (eV)",
        "Eg_eV": "Eg (eV)",
        "Experimental_Eg_eV": "Exp. Eg (eV)",
        "Structure_Availability": "Structure",
        "Duplicate_Flag": "Repeated",
        "Curation_Status": "Curation status",
    })
    if "Database ID" in display_table.columns:
        display_table["Database ID"] = display_table["Database ID"].map(_display_record_id)
    if "Collection role" in display_table.columns:
        display_table["Collection role"] = display_table["Collection role"].map(_display_role)
    if "Scope" in display_table.columns:
        display_table["Scope"] = display_table["Scope"].map(_display_scope)

    for column in display_table.columns:
        display_table[column] = display_table[column].map(
            lambda value: _display_value(value)
        )

    st.dataframe(
        display_table,
        use_container_width=True,
        hide_index=True,
        height=420,
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
        with st.expander("Filtered view summary", expanded=False):
            group_cols = [
                c for c in ["Split_Role", "Scope_Flag"]
                if c in filtered.columns
            ]
            summary = (
                filtered.groupby(group_cols, dropna=False)
                .size()
                .reset_index(name="Count")
            )
            if "Split_Role" in summary.columns:
                summary["Split_Role"] = summary["Split_Role"].map(_display_role)
                summary = summary.rename(columns={"Split_Role": "Collection role"})
            if "Scope_Flag" in summary.columns:
                summary["Scope_Flag"] = summary["Scope_Flag"].map(_display_scope)
                summary = summary.rename(columns={"Scope_Flag": "Scope"})
            st.dataframe(summary, use_container_width=True, hide_index=True)

    if not is_preview:
        st.download_button(
            "Download filtered database CSV",
            filtered.to_csv(index=False).encode("utf-8"),
            "chalcogen_database_filtered.csv",
            "text/csv",
            use_container_width=True,
        )

    if is_preview and source_path:
        st.caption(f"Preview source: {source_path}")
