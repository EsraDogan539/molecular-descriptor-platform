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
        return "External collection"
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
    record_id = _display_record_id(row.get("Record_ID"))
    collection = _display_collection(row)
    chalcogen = _display_value(row.get("Chalcogen_Type"))

    st.markdown(f"### {record_id}")
    st.caption(f"{collection} · {chalcogen}")

    left, right = st.columns([1, 2.25], gap="large")

    with left:
        smiles = row.get("Canonical_SMILES")
        if pd.notna(smiles) and str(smiles).strip():
            mol = Chem.MolFromSmiles(str(smiles))
            if mol is not None:
                st.image(Draw.MolToImage(mol, size=(420, 320)), use_container_width=True)
            else:
                st.info("The standardized structure could not be rendered.")
        else:
            st.info("Exact standardized structure is not available for this record.")

        if pd.notna(smiles) and str(smiles).strip():
            st.caption("Canonical SMILES")
            st.code(str(smiles), language=None)

    with right:
        st.markdown("#### Identity")
        identity = {
            "Database ID": record_id,
            "Molecule / system": row.get("Molecule_Name") or row.get("System_Code"),
            "InChIKey": row.get("InChIKey"),
            "Structure availability": row.get("Structure_Availability"),
        }
        st.dataframe(
            pd.DataFrame(
                [{"Field": k, "Value": _display_value(v)} for k, v in identity.items()]
            ),
            use_container_width=True,
            hide_index=True,
        )

        st.markdown("#### Electronic properties")
        e1, e2, e3, e4 = st.columns(4)
        e1.metric("HOMO (eV)", _display_value(row.get("HOMO_eV"), 3))
        e2.metric("LUMO (eV)", _display_value(row.get("LUMO_eV"), 3))
        e3.metric("Eg (eV)", _display_value(row.get("Eg_eV"), 3))
        e4.metric("Exp. Eg (eV)", _display_value(row.get("Experimental_Eg_eV"), 3))

        chalcogen_fields = [
            ("O count", "O_Count"),
            ("S count", "S_Count"),
            ("Se count", "Se_Count"),
            ("Te count", "Te_Count"),
            ("Chalcogen type", "Chalcogen_Type"),
        ]
        available_chalcogen = [
            (label, column) for label, column in chalcogen_fields if column in row.index
        ]
        if available_chalcogen:
            st.markdown("#### Chalcogen environment")
            st.dataframe(
                pd.DataFrame([
                    {"Field": label, "Value": _display_value(row.get(column))}
                    for label, column in available_chalcogen
                ]),
                use_container_width=True,
                hide_index=True,
            )

        st.markdown("#### Provenance")
        provenance = {
            "Collection": collection,
            "Donor / Acceptor": " / ".join(
                [
                    str(v)
                    for v in [row.get("Donor_ID"), row.get("Acceptor_ID")]
                    if pd.notna(v) and str(v).strip()
                ]
            ) or "—",
            "Unit type": row.get("Unit_Type"),
            "Oligomer n": row.get("Oligomer_n"),
            "Method": row.get("Method"),
            "Basis set": row.get("Basis_Set"),
            "Conditions": row.get("Solvent_or_Conditions"),
            "Reference": row.get("Reference"),
            "Curation status": row.get("Curation_Status"),
        }
        st.dataframe(
            pd.DataFrame(
                [{"Field": k, "Value": _display_value(v)} for k, v in provenance.items()]
            ),
            use_container_width=True,
            hide_index=True,
        )


def _public_table(df):
    preferred_columns = [
        "Record_ID", "Molecule_Name", "System_Code", "Chalcogen_Type",
        "Eg_eV", "HOMO_eV", "LUMO_eV", "Split_Role",
    ]
    existing = [c for c in preferred_columns if c in df.columns]
    table = df[existing].copy()
    table = table.rename(columns={
        "Record_ID": "Database ID",
        "Molecule_Name": "Molecule / system",
        "System_Code": "System",
        "Chalcogen_Type": "Chalcogen",
        "Eg_eV": "Eg (eV)",
        "HOMO_eV": "HOMO (eV)",
        "LUMO_eV": "LUMO (eV)",
        "Split_Role": "Collection",
    })
    if "Database ID" in table.columns:
        table["Database ID"] = table["Database ID"].map(_display_record_id)
    if "Collection" in table.columns:
        table["Collection"] = table["Collection"].map(_display_role)

    if "Molecule / system" in table.columns and "System" in table.columns:
        table["Molecule / system"] = table["Molecule / system"].where(
            table["Molecule / system"].notna(), table["System"]
        )
        table = table.drop(columns=["System"])
    elif "System" in table.columns and "Molecule / system" not in table.columns:
        table = table.rename(columns={"System": "Molecule / system"})

    for col in table.columns:
        if col in {"Eg (eV)", "HOMO (eV)", "LUMO (eV)"}:
            table[col] = table[col].map(lambda v: _display_value(v, 3))
        else:
            table[col] = table[col].map(_display_value)
    return table


def _render_results(filtered, is_preview):
    st.caption(f"Showing {len(filtered):,} record(s)")
    st.dataframe(
        _public_table(filtered),
        use_container_width=True,
        hide_index=True,
        height=430,
    )

    if filtered.empty:
        return

    st.markdown("#### Record detail")
    labels = [_record_label(row) for _, row in filtered.iterrows()]
    selected = st.selectbox("Select a record", labels, label_visibility="collapsed")
    row = filtered.iloc[labels.index(selected)]
    _render_record_detail(row)

    if not is_preview:
        st.divider()
        st.download_button(
            "Download filtered records",
            filtered.to_csv(index=False).encode("utf-8"),
            "chalcogen_database_filtered.csv",
            "text/csv",
        )


def _render_statistics(df):
    st.markdown("### Database statistics")
    st.caption("Descriptive counts from database v1.")

    if "Split_Role" in df.columns:
        role = df["Split_Role"].map(_display_role).value_counts().rename("Records")
        st.markdown("#### Collection composition")
        st.bar_chart(role)

    if all(c in df.columns for c in ["S_Count", "Se_Count", "Te_Count"]):
        counts = pd.Series({
            "S": int((pd.to_numeric(df["S_Count"], errors="coerce").fillna(0) > 0).sum()),
            "Se": int((pd.to_numeric(df["Se_Count"], errors="coerce").fillna(0) > 0).sum()),
            "Te": int((pd.to_numeric(df["Te_Count"], errors="coerce").fillna(0) > 0).sum()),
        }, name="Records")
        st.markdown("#### Chalcogen coverage")
        st.bar_chart(counts)

    if "Eg_eV" in df.columns:
        eg = pd.to_numeric(df["Eg_eV"], errors="coerce").dropna()
        if not eg.empty:
            bins = pd.cut(eg, bins=16)
            hist = bins.value_counts(sort=False)
            hist.index = [f"{x.left:.2f}–{x.right:.2f}" for x in hist.index]
            st.markdown("#### Eg distribution")
            st.bar_chart(hist.rename("Records"))


def display_database_browser():
    st.header("Database")
    st.caption(
        "Browse curated records, search by molecular identity, and inspect "
        "structure, electronic properties and provenance."
    )

    df, is_preview, source_path = load_curated_database()
    if df.empty:
        st.warning("The curated database file is not available in this build.")
        return

    total = len(df)
    core = int((df["Scope_Flag"].astype(str) == "Core_SSeTe").sum()) if "Scope_Flag" in df.columns else 0
    unique = int(
        df["InChIKey"].dropna().astype(str).replace("", pd.NA).dropna().nunique()
    ) if "InChIKey" in df.columns else 0
    eg_count = int(df["Eg_eV"].notna().sum()) if "Eg_eV" in df.columns else 0

    st.caption(
        f"{total:,} records · {core:,} core S/Se/Te · "
        f"{unique:,} unique standardized structures · {eg_count:,} Eg values"
    )

    if is_preview:
        st.info("A preview dataset is loaded in this build.")

    browse_tab, search_tab, stats_tab = st.tabs(["Browse", "Search", "Statistics"])

    with browse_tab:
        f1, f2, f3, f4 = st.columns(4)

        split_options = ["All"] + _safe_unique(df, "Split_Role")
        scope_options = ["All"] + _safe_unique(df, "Scope_Flag")
        chalcogen_options = ["All"] + _safe_unique(df, "Chalcogen_Type")

        with f1:
            split_value = st.selectbox(
                "Collection",
                split_options,
                format_func=lambda v: "All" if v == "All" else _display_role(v),
                key="browse_collection",
            )
        with f2:
            scope_value = st.selectbox(
                "Scope",
                scope_options,
                format_func=lambda v: "All" if v == "All" else _display_scope(v),
                key="browse_scope",
            )
        with f3:
            chalcogen_value = st.selectbox(
                "Chalcogen",
                chalcogen_options,
                key="browse_chalcogen",
            )
        with f4:
            eg_value = st.selectbox(
                "Eg available",
                ["All", "Available", "Missing"],
                key="browse_eg",
            )

        filtered = df.copy()
        if split_value != "All" and "Split_Role" in filtered.columns:
            filtered = filtered[filtered["Split_Role"].astype(str) == split_value]
        if scope_value != "All" and "Scope_Flag" in filtered.columns:
            filtered = filtered[filtered["Scope_Flag"].astype(str) == scope_value]
        if chalcogen_value != "All" and "Chalcogen_Type" in filtered.columns:
            filtered = filtered[filtered["Chalcogen_Type"].astype(str) == chalcogen_value]
        if eg_value == "Available" and "Eg_eV" in filtered.columns:
            filtered = filtered[filtered["Eg_eV"].notna()]
        elif eg_value == "Missing" and "Eg_eV" in filtered.columns:
            filtered = filtered[filtered["Eg_eV"].isna()]

        _render_results(filtered, is_preview)

    with search_tab:
        query = st.text_input(
            "Search database",
            placeholder="Database ID, molecule/system, donor, acceptor, InChIKey or SMILES",
        ).strip()

        if not query:
            st.caption("Enter a molecular identifier or keyword.")
        else:
            searchable = [
                c for c in [
                    "Record_ID", "Molecule_Name", "Donor_ID", "Acceptor_ID",
                    "System_Code", "InChIKey", "Canonical_SMILES"
                ] if c in df.columns
            ]
            mask = pd.Series(False, index=df.index)
            for column in searchable:
                mask = mask | df[column].astype(str).str.contains(
                    query, case=False, na=False, regex=False
                )
            if "Record_ID" in df.columns:
                public_ids = df["Record_ID"].map(_display_record_id)
                mask = mask | public_ids.str.contains(
                    query, case=False, na=False, regex=False
                )
            _render_results(df[mask], is_preview)

    with stats_tab:
        _render_statistics(df)

    with st.expander("Curation notes", expanded=False):
        st.markdown(
            "- Development records contain exact standardized molecular structures.\n"
            "- External records are retained at the level supported by the source.\n"
            "- Missing structures and properties are left explicit.\n"
            "- Repeated standardized structures are retained to preserve record-level provenance."
        )

    if is_preview and source_path:
        st.caption(f"Preview source: {source_path}")
