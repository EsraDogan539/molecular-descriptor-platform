import os
import pandas as pd
import streamlit as st
import matplotlib.pyplot as plt

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


def _field_line(label, value):
    return f"**{label}:** {_display_value(value)}"


def _display_count(value):
    if value is None or pd.isna(value):
        return "—"
    try:
        return str(int(float(value)))
    except (TypeError, ValueError):
        return _display_value(value)


def _display_method(value):
    text = _display_value(value)
    if text == "—":
        return text
    lowered = text.lower()
    if "source sdf quantum-chemical properties" in lowered:
        return "Quantum-chemical calculation (exact level not stated)"
    return text


def _render_record_detail(row):
    record_id = _display_record_id(row.get("Record_ID"))
    collection = _display_collection(row)
    chalcogen = _display_value(row.get("Chalcogen_Type"))

    st.markdown(f"### {record_id}")
    st.caption(f"{collection} · {chalcogen}")

    left, right = st.columns([1.05, 2.25], gap="large")

    with left:
        smiles = row.get("Canonical_SMILES")
        if pd.notna(smiles) and str(smiles).strip():
            mol = Chem.MolFromSmiles(str(smiles))
            if mol is not None:
                st.image(Draw.MolToImage(mol, size=(420, 300)), use_container_width=True)
            else:
                st.info("The standardized structure could not be rendered.")

            st.caption("Canonical SMILES")
            st.code(str(smiles), language=None)
        else:
            st.info("Exact standardized structure is not available for this record.")

    with right:
        st.markdown("#### Identity")
        i1, i2, i3 = st.columns(3)
        i1.markdown(_field_line("Database ID", record_id))
        i2.markdown(_field_line(
            "Molecule / system",
            row.get("Molecule_Name") or row.get("System_Code"),
        ))
        i3.markdown(_field_line("InChIKey", row.get("InChIKey")))

        st.divider()

        st.markdown("#### Electronic properties")
        e1, e2, e3, e4 = st.columns(4)
        e1.metric("HOMO (eV)", _display_value(row.get("HOMO_eV"), 3))
        e2.metric("LUMO (eV)", _display_value(row.get("LUMO_eV"), 3))
        e3.metric("Eg (eV)", _display_value(row.get("Eg_eV"), 3))
        e4.metric("Exp. Eg (eV)", _display_value(row.get("Experimental_Eg_eV"), 3))

        st.divider()

        c1, c2 = st.columns(2)
        with c1:
            st.markdown("#### Chalcogen environment")
            st.markdown(_field_line("S count", _display_count(row.get("S_Count"))))
            st.markdown(_field_line("Se count", _display_count(row.get("Se_Count"))))
            st.markdown(_field_line("Te count", _display_count(row.get("Te_Count"))))
            st.markdown(_field_line("Chalcogen type", row.get("Chalcogen_Type")))

        with c2:
            st.markdown("#### Provenance")
            donor_acceptor = " / ".join(
                [
                    str(v)
                    for v in [row.get("Donor_ID"), row.get("Acceptor_ID")]
                    if pd.notna(v) and str(v).strip()
                ]
            ) or "—"
            st.markdown(_field_line("Donor / Acceptor", donor_acceptor))
            st.markdown(_field_line("Method", _display_method(row.get("Method"))))
            st.markdown(_field_line("Basis set", row.get("Basis_Set")))

        with st.expander("Additional record metadata", expanded=False):
            metadata = {
                "Collection": collection,
                "Scope": _display_scope(row.get("Scope_Flag")),
                "Unit type": row.get("Unit_Type"),
                "Oligomer n": row.get("Oligomer_n"),
                "Conditions": row.get("Solvent_or_Conditions"),
                "Structure availability": row.get("Structure_Availability"),
                "Curation status": row.get("Curation_Status"),
                "Repeated structure": row.get("Duplicate_Flag"),
                "Source method description": row.get("Method"),
                "Reference": row.get("Reference"),
            }
            st.dataframe(
                pd.DataFrame(
                    [{"Field": key, "Value": _display_value(value)} for key, value in metadata.items()]
                ),
                use_container_width=True,
                hide_index=True,
            )


def _public_table(df):
    preferred_columns = [
        "Record_ID", "Molecule_Name", "System_Code", "Chalcogen_Type",
        "Eg_eV", "HOMO_eV", "LUMO_eV",
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
    })

    if "Database ID" in table.columns:
        table["Database ID"] = table["Database ID"].map(_display_record_id)

    if "Molecule / system" in table.columns and "System" in table.columns:
        table["Molecule / system"] = table["Molecule / system"].where(
            table["Molecule / system"].notna(), table["System"]
        )
        table = table.drop(columns=["System"])
    elif "System" in table.columns and "Molecule / system" not in table.columns:
        table = table.rename(columns={"System": "Molecule / system"})

    for column in table.columns:
        if column in {"Eg (eV)", "HOMO (eV)", "LUMO (eV)"}:
            table[column] = table[column].map(lambda value: _display_value(value, 3))
        else:
            table[column] = table[column].map(_display_value)

    return table


def _apply_search(df, query):
    if not query:
        return df

    searchable = [
        c for c in [
            "Record_ID", "Molecule_Name", "Donor_ID", "Acceptor_ID",
            "System_Code", "InChIKey", "Canonical_SMILES"
        ]
        if c in df.columns
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

    return df[mask]


def _render_results(filtered, is_preview):
    st.markdown(
        f'<div class="result-count">{len(filtered):,} matching records</div>',
        unsafe_allow_html=True,
    )

    st.dataframe(
        _public_table(filtered),
        use_container_width=True,
        hide_index=True,
        height=420,
    )

    if filtered.empty:
        st.info("No records match the current search and filters.")
        return

    st.divider()
    st.markdown('<div class="section-rule-title">Inspect a record</div>', unsafe_allow_html=True)
    labels = [_record_label(row) for _, row in filtered.iterrows()]
    selected = st.selectbox(
        "Record",
        labels,
        label_visibility="collapsed",
    )
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


def _style_axes(ax, ylabel="Records"):
    ax.set_facecolor("white")
    ax.spines["top"].set_visible(False)
    ax.spines["right"].set_visible(False)
    ax.spines["left"].set_color("#C8D5E0")
    ax.spines["bottom"].set_color("#C8D5E0")
    ax.tick_params(axis="both", colors="#163A5B", labelsize=9)
    ax.set_ylabel(ylabel, color="#163A5B", fontsize=9)
    ax.grid(axis="y", linestyle=":", linewidth=0.6, alpha=0.5)
    ax.set_axisbelow(True)


def _render_statistics(df):
    st.markdown('<div class="section-rule-title">Distribution overview</div>', unsafe_allow_html=True)
    st.caption("Descriptive overview of database v1.")

    stat1, stat2 = st.columns(2, gap="large")

    with stat1:
        if "Split_Role" in df.columns:
            role = df["Split_Role"].map(_display_role).value_counts()
            fig, ax = plt.subplots(figsize=(5.2, 3.0))
            x = list(range(len(role)))
            role_colors = ["#1F4E79", "#5DA9E9"]
            bars = ax.bar(
                x,
                role.values,
                color=[role_colors[i % len(role_colors)] for i in range(len(role))],
                edgecolor="#163A5B",
                linewidth=0.8,
            )
            for i, bar in enumerate(bars):
                ax.text(
                    bar.get_x() + bar.get_width() / 2,
                    bar.get_height() + max(role.values) * 0.025,
                    f"{int(bar.get_height()):,}",
                    ha="center",
                    va="bottom",
                    fontsize=9,
                    color="#111827",
                )
            ax.set_xticks(x)
            ax.set_xticklabels(role.index)
            ax.set_title("Collection composition", loc="left", fontsize=11, fontweight="bold", pad=12)
            ax.margins(y=0.12)
            _style_axes(ax)
            fig.tight_layout(pad=1.4)
            st.pyplot(fig, use_container_width=True)
            plt.close(fig)

    with stat2:
        if all(c in df.columns for c in ["S_Count", "Se_Count", "Te_Count"]):
            def _element_presence(element, count_column):
                present = pd.to_numeric(
                    df[count_column], errors="coerce"
                ).fillna(0).gt(0)

                for metadata_column in ["Donor_Chalcogen", "Acceptor_Chalcogen", "Chalcogen_Type"]:
                    if metadata_column in df.columns:
                        metadata = df[metadata_column].fillna("").astype(str)
                        if element == "S":
                            match = metadata.str.contains(
                                r"(^|[^A-Za-z])S([^A-Za-z]|$)",
                                case=False,
                                regex=True,
                            )
                        elif element == "Se":
                            match = metadata.str.contains(
                                r"(^|[^A-Za-z])Se([^A-Za-z]|$)",
                                case=False,
                                regex=True,
                            )
                        else:
                            match = metadata.str.contains(
                                r"(^|[^A-Za-z])Te([^A-Za-z]|$)",
                                case=False,
                                regex=True,
                            )
                        present = present | match

                return int(present.sum())

            counts = pd.Series({
                "S": _element_presence("S", "S_Count"),
                "Se": _element_presence("Se", "Se_Count"),
                "Te": _element_presence("Te", "Te_Count"),
            })
            fig, ax = plt.subplots(figsize=(5.2, 3.0))
            x = list(range(len(counts)))
            element_colors = ["#D6A400", "#0F766E", "#6B4C8A"]
            bars = ax.bar(
                x,
                counts.values,
                color=element_colors,
                edgecolor="#163A5B",
                linewidth=0.8,
            )
            for i, bar in enumerate(bars):
                ax.text(
                    bar.get_x() + bar.get_width() / 2,
                    bar.get_height() + max(counts.values) * 0.025,
                    f"{int(bar.get_height()):,}",
                    ha="center",
                    va="bottom",
                    fontsize=9,
                    color="#111827",
                )
            ax.set_xticks(x)
            ax.set_xticklabels(counts.index)
            ax.set_title("Chalcogen coverage", loc="left", fontsize=11, fontweight="bold", pad=12)
            ax.margins(y=0.12)
            _style_axes(ax)
            fig.tight_layout(pad=1.4)
            st.pyplot(fig, use_container_width=True)
            plt.close(fig)

    if "Eg_eV" in df.columns:
        eg = pd.to_numeric(df["Eg_eV"], errors="coerce").dropna()
        if not eg.empty:
            fig, ax = plt.subplots(figsize=(10.6, 3.6))
            ax.hist(
                eg,
                bins=16,
                color="#5DA9E9",
                edgecolor="#163A5B",
                linewidth=0.8,
            )
            ax.set_xlabel("Eg (eV)", color="#163A5B", fontsize=9)
            ax.set_title("Eg distribution", loc="left", fontsize=11, fontweight="bold", pad=12)
            _style_axes(ax)
            fig.tight_layout(pad=1.4)
            st.pyplot(fig, use_container_width=True)
            plt.close(fig)

def display_database_browser():
    st.markdown(
        """
        <style>
        div[data-testid="stTextInput"] input {
            background: #FFFFFF !important;
            border: 1px solid #D9DEE7 !important;
            box-shadow: none !important;
        }
        div[data-testid="stTextInput"] input:focus {
            border-color: #1F4E79 !important;
            box-shadow: 0 0 0 1px #1F4E79 !important;
        }
        .result-count {
            color: #163A5B;
            font-size: .88rem;
            font-weight: 600;
            margin: .15rem 0 .45rem 0;
        }
        </style>
        """,
        unsafe_allow_html=True,
    )
    st.markdown(
        """
        <div class="page-intro">
          <div class="page-eyebrow">Curated records</div>
          <div class="page-title">Database</div>
          <div class="page-description">Browse curated records and inspect molecular identity, electronic properties, chalcogen context and record-level provenance.</div>
        </div>
        """,
        unsafe_allow_html=True,
    )

    df, is_preview, source_path = load_curated_database()
    if df.empty:
        st.warning("The curated database file is not available in this build.")
        return

    total = len(df)
    core = int(
        (df["Scope_Flag"].astype(str) == "Core_SSeTe").sum()
    ) if "Scope_Flag" in df.columns else 0
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

    query = st.text_input(
        "Search",
        placeholder="Database ID, molecule/system, donor, acceptor, InChIKey or SMILES",
    ).strip()

    f1, f2, f3, f4 = st.columns(4)

    split_options = ["All"] + _safe_unique(df, "Split_Role")
    scope_options = ["All"] + _safe_unique(df, "Scope_Flag")
    chalcogen_options = ["All"] + _safe_unique(df, "Chalcogen_Type")

    with f1:
        split_value = st.selectbox(
            "Collection",
            split_options,
            format_func=lambda value: "All" if value == "All" else _display_role(value),
            key="browse_collection",
        )
    with f2:
        scope_value = st.selectbox(
            "Scope",
            scope_options,
            format_func=lambda value: "All" if value == "All" else _display_scope(value),
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
            "Eg",
            ["All", "Available", "Missing"],
            key="browse_eg",
        )

    filtered = _apply_search(df.copy(), query)

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

    with st.expander("Curation notes", expanded=False):
        st.markdown(
            "- Development records contain exact standardized molecular structures.\n"
            "- External records are retained at the level supported by the source.\n"
            "- Missing structures and properties are left explicit.\n"
            "- Repeated standardized structures are retained to preserve record-level provenance."
        )

    if is_preview and source_path:
        st.caption(f"Preview source: {source_path}")


def display_database_statistics():
    """Render the publication-database statistics as a standalone page."""
    st.markdown(
        """
        <div class="page-intro">
          <div class="page-eyebrow">Database overview</div>
          <div class="page-title">Statistics</div>
          <div class="page-description">Database composition, chalcogen coverage and electronic-property availability across the current ChalMolDB release.</div>
        </div>
        """,
        unsafe_allow_html=True,
    )

    df, is_preview, source_path = load_curated_database()
    if df.empty:
        st.warning("The curated database file is not available in this build.")
        return

    total = len(df)
    core = int(
        (df["Scope_Flag"].astype(str) == "Core_SSeTe").sum()
    ) if "Scope_Flag" in df.columns else 0
    unique = int(
        df["InChIKey"].dropna().astype(str).replace("", pd.NA).dropna().nunique()
    ) if "InChIKey" in df.columns else 0
    eg_count = int(df["Eg_eV"].notna().sum()) if "Eg_eV" in df.columns else 0

    m1, m2, m3, m4 = st.columns(4)
    m1.metric("Records", f"{total:,}")
    m2.metric("Core S/Se/Te", f"{core:,}")
    m3.metric("Unique structures", f"{unique:,}")
    m4.metric("Eg values", f"{eg_count:,}")

    st.divider()
    _render_statistics(df)

    if is_preview:
        st.info("A preview dataset is loaded in this build.")
    if is_preview and source_path:
        st.caption(f"Preview source: {source_path}")
