import os
import html
import pandas as pd
import streamlit as st
import matplotlib.pyplot as plt

from rdkit import Chem
from rdkit.Chem import Draw

from database_quality import (
    available_reference_column,
    coverage_table,
    database_quality_summary,
    reference_mask,
)
from database_comparison import comparison_export, comparison_table
from structure_search import (
    MORGAN_N_BITS,
    MORGAN_RADIUS,
    SIMILARITY_METRIC,
    structure_search,
)


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
    smiles = row.get("Canonical_SMILES")

    st.markdown(
        f"""
        <div class="record-heading">
          <div class="record-id">{html.escape(str(record_id))}</div>
          <div class="record-meta">{html.escape(str(collection))} · {html.escape(str(chalcogen))}</div>
        </div>
        """,
        unsafe_allow_html=True,
    )

    with st.container(border=True):
        left, right = st.columns([0.88, 2.35], gap="large")

        with left:
            if pd.notna(smiles) and str(smiles).strip():
                mol = Chem.MolFromSmiles(str(smiles))
                if mol is not None:
                    st.image(Draw.MolToImage(mol, size=(360, 245)), use_container_width=True)
                else:
                    st.info("The standardized structure could not be rendered.")

                st.markdown('<div class="detail-label">Canonical SMILES</div>', unsafe_allow_html=True)
                st.code(str(smiles), language=None)
            else:
                st.info("Exact standardized structure is not available for this record.")

        with right:
            molecule_system = row.get("Molecule_Name") or row.get("System_Code")
            identity_html = f"""
            <div class="detail-section">
              <div class="detail-section-title">Identity</div>
              <div class="detail-field-grid detail-field-grid-3">
                <div class="detail-field"><span>Database ID</span><strong>{html.escape(str(_display_value(record_id)))}</strong></div>
                <div class="detail-field"><span>Molecule / system</span><strong>{html.escape(str(_display_value(molecule_system)))}</strong></div>
                <div class="detail-field"><span>InChIKey</span><strong>{html.escape(str(_display_value(row.get("InChIKey"))))}</strong></div>
              </div>
            </div>
            """
            st.markdown(identity_html, unsafe_allow_html=True)

            st.markdown('<div class="detail-section-title detail-metric-title">Electronic properties</div>', unsafe_allow_html=True)
            e1, e2, e3, e4 = st.columns(4)
            e1.metric("HOMO (eV)", _display_value(row.get("HOMO_eV"), 3))
            e2.metric("LUMO (eV)", _display_value(row.get("LUMO_eV"), 3))
            e3.metric("Eg (eV)", _display_value(row.get("Eg_eV"), 3))
            e4.metric("Exp. Eg (eV)", _display_value(row.get("Experimental_Eg_eV"), 3))

            donor_acceptor = " / ".join(
                [
                    str(v)
                    for v in [row.get("Donor_ID"), row.get("Acceptor_ID")]
                    if pd.notna(v) and str(v).strip()
                ]
            ) or "—"

            detail_html = f"""
            <div class="detail-bottom-grid">
              <div class="detail-subcard">
                <div class="detail-section-title">Chalcogen environment</div>
                <div class="compact-fields">
                  <div><span>S count</span><strong>{html.escape(_display_count(row.get("S_Count")))}</strong></div>
                  <div><span>Se count</span><strong>{html.escape(_display_count(row.get("Se_Count")))}</strong></div>
                  <div><span>Te count</span><strong>{html.escape(_display_count(row.get("Te_Count")))}</strong></div>
                  <div><span>Type</span><strong>{html.escape(str(_display_value(row.get("Chalcogen_Type"))))}</strong></div>
                </div>
              </div>
              <div class="detail-subcard">
                <div class="detail-section-title">Provenance</div>
                <div class="compact-fields provenance-fields">
                  <div><span>Dataset / owner</span><strong>{html.escape(str(_display_value(row.get("Dataset_Owner"))))}</strong></div>
                  <div><span>Donor / Acceptor</span><strong>{html.escape(donor_acceptor)}</strong></div>
                  <div><span>Method</span><strong>{html.escape(str(_display_method(row.get("Method"))))}</strong></div>
                  <div><span>Basis set</span><strong>{html.escape(str(_display_value(row.get("Basis_Set"))))}</strong></div>
                  <div><span>Reference / DOI</span><strong>{html.escape(str(_display_value(row.get(available_reference_column(row.to_frame().T)) if available_reference_column(row.to_frame().T) else None)))}</strong></div>
                </div>
              </div>
            </div>
            """
            st.markdown(detail_html, unsafe_allow_html=True)

            with st.expander("Additional record metadata", expanded=False):
                row_df = row.to_frame().T
                reference_column = available_reference_column(row_df)
                metadata = {
                    "Collection": collection,
                    "Dataset / owner": row.get("Dataset_Owner"),
                    "Scope": _display_scope(row.get("Scope_Flag")),
                    "Unit type": row.get("Unit_Type"),
                    "Oligomer n": row.get("Oligomer_n"),
                    "Conditions": row.get("Solvent_or_Conditions"),
                    "Structure availability": row.get("Structure_Availability"),
                    "Curation status": row.get("Curation_Status"),
                    "Repeated structure": row.get("Duplicate_Flag"),
                    "Source method description": row.get("Method"),
                    "Reference / DOI": row.get(reference_column) if reference_column else None,
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
        "Structure Similarity", "Eg_eV", "HOMO_eV", "LUMO_eV",
    ]
    existing = [c for c in preferred_columns if c in df.columns]
    table = df[existing].copy()

    table = table.rename(columns={
        "Record_ID": "Database ID",
        "Molecule_Name": "Molecule / system",
        "System_Code": "System",
        "Chalcogen_Type": "Chalcogen",
        "Structure Similarity": "Similarity",
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
        if column in {"Eg (eV)", "HOMO (eV)", "LUMO (eV)", "Similarity"}:
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
            "System_Code", "InChIKey", "Canonical_SMILES",
            "Dataset_Owner", "Method", "Basis_Set", "Curation_Status",
            "DOI_or_Reference", "DOI", "Reference", "Source_Reference", "Citation",
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


def _render_comparison_structure(row, title):
    st.markdown(f"**{title}**")
    record_id = _display_record_id(row.get("Record_ID"))
    molecule_system = row.get("Molecule_Name")
    if pd.isna(molecule_system) or not str(molecule_system).strip():
        molecule_system = row.get("System_Code")

    st.caption(
        f"{record_id} · {_display_value(molecule_system)} · "
        f"{_display_value(row.get('Chalcogen_Type'))}"
    )

    smiles = row.get("Canonical_SMILES")
    if pd.notna(smiles) and str(smiles).strip():
        mol = Chem.MolFromSmiles(str(smiles))
        if mol is not None:
            st.image(
                Draw.MolToImage(mol, size=(360, 235)),
                use_container_width=True,
            )
        st.code(str(smiles), language=None)
    else:
        st.info("Exact standardized structure is not available for this record.")

    p1, p2, p3 = st.columns(3)
    p1.metric("HOMO (eV)", _display_value(row.get("HOMO_eV"), 3))
    p2.metric("LUMO (eV)", _display_value(row.get("LUMO_eV"), 3))
    p3.metric("Eg (eV)", _display_value(row.get("Eg_eV"), 3))

    st.caption(
        f"Experimental Eg: {_display_value(row.get('Experimental_Eg_eV'))} · "
        f"Method: {_display_method(row.get('Method'))} · "
        f"Basis: {_display_value(row.get('Basis_Set'))}"
    )


def _render_record_comparison(filtered):
    if len(filtered) < 2:
        return

    record_lookup = {
        _record_label(row): row
        for _, row in filtered.iterrows()
    }
    labels = list(record_lookup.keys())

    with st.expander("Compare records", expanded=False):
        st.caption(
            "Compare two curated records descriptively. Property differences do not "
            "establish causal structure-property relationships."
        )
        c1, c2 = st.columns(2)
        with c1:
            reference_label = st.selectbox(
                "Reference record",
                labels,
                key="database_compare_reference",
            )
        with c2:
            candidate_options = [
                label for label in labels if label != reference_label
            ]
            candidate_label = st.selectbox(
                "Candidate record",
                candidate_options,
                key="database_compare_candidate",
            )

        reference_row = record_lookup[reference_label]
        candidate_row = record_lookup[candidate_label]

        left, right = st.columns(2, gap="large")
        with left:
            _render_comparison_structure(
                reference_row,
                f"Reference · {_display_record_id(reference_row.get('Record_ID'))}",
            )
        with right:
            _render_comparison_structure(
                candidate_row,
                f"Candidate · {_display_record_id(candidate_row.get('Record_ID'))}",
            )

        st.markdown(
            '<div class="section-rule-title">Property and provenance comparison</div>',
            unsafe_allow_html=True,
        )
        comparison_df = comparison_table(reference_row, candidate_row)
        if not comparison_df.empty:
            st.dataframe(
                comparison_df,
                use_container_width=True,
                hide_index=True,
            )

        reference_column = available_reference_column(filtered)
        export_df = comparison_export(
            reference_row=reference_row,
            candidate_row=candidate_row,
            reference_public_id=_display_record_id(reference_row.get("Record_ID")),
            candidate_public_id=_display_record_id(candidate_row.get("Record_ID")),
            reference_column=reference_column,
        )
        st.download_button(
            "Download comparison CSV",
            export_df.to_csv(index=False).encode("utf-8"),
            (
                f"{_display_record_id(reference_row.get('Record_ID'))}_vs_"
                f"{_display_record_id(candidate_row.get('Record_ID'))}_comparison.csv"
            ),
            "text/csv",
            key="download_database_comparison",
        )


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

    _render_record_comparison(filtered)

    st.divider()
    st.markdown('<div class="section-rule-title">Inspect a record</div>', unsafe_allow_html=True)
    record_lookup = {
        _record_label(row): row
        for _, row in filtered.iterrows()
    }
    labels = list(record_lookup.keys())
    selected = st.selectbox(
        "Record",
        labels,
        label_visibility="collapsed",
    )
    row = record_lookup[selected]
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
    ax.set_facecolor("#FFFFFF")
    ax.spines["top"].set_visible(False)
    ax.spines["right"].set_visible(False)
    ax.spines["left"].set_color("#CBD6DF")
    ax.spines["bottom"].set_color("#CBD6DF")
    ax.spines["left"].set_linewidth(0.8)
    ax.spines["bottom"].set_linewidth(0.8)
    ax.tick_params(axis="both", colors="#35516A", labelsize=8.5)
    ax.set_ylabel(ylabel, color="#50677B", fontsize=8.5, labelpad=7)
    ax.grid(axis="y", color="#DCE5EC", linestyle=":", linewidth=0.55, alpha=0.8)
    ax.set_axisbelow(True)


def _render_statistics(df):
    st.markdown('<div class="section-rule-title">Distribution overview</div>', unsafe_allow_html=True)
    st.caption("Descriptive overview of database v1.")

    stat1, stat2 = st.columns(2, gap="medium")

    with stat1:
        if "Split_Role" in df.columns:
            with st.container(border=True):
                role = df["Split_Role"].map(_display_role).value_counts()
                fig, ax = plt.subplots(figsize=(5.2, 2.55))
                fig.patch.set_facecolor("#FFFFFF")
                x = list(range(len(role)))
                role_colors = ["#1F4E79", "#5DA9E9"]
                bars = ax.bar(
                    x,
                    role.values,
                    color=[role_colors[i % len(role_colors)] for i in range(len(role))],
                    edgecolor="#315C7D",
                    linewidth=0.65,
                    width=0.68,
                )
                for bar in bars:
                    ax.text(
                        bar.get_x() + bar.get_width() / 2,
                        bar.get_height() + max(role.values) * 0.022,
                        f"{int(bar.get_height()):,}",
                        ha="center",
                        va="bottom",
                        fontsize=8.5,
                        color="#334155",
                    )
                ax.set_xticks(x)
                ax.set_xticklabels(role.index)
                ax.set_title("Collection composition", loc="left", fontsize=10, fontweight="semibold", color="#163A5B", pad=8)
                ax.margins(y=0.12)
                _style_axes(ax)
                fig.tight_layout(pad=0.9)
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

            with st.container(border=True):
                counts = pd.Series({
                    "S": _element_presence("S", "S_Count"),
                    "Se": _element_presence("Se", "Se_Count"),
                    "Te": _element_presence("Te", "Te_Count"),
                })
                fig, ax = plt.subplots(figsize=(5.2, 2.55))
                fig.patch.set_facecolor("#FFFFFF")
                x = list(range(len(counts)))
                element_colors = ["#D6A400", "#0F766E", "#6B4C8A"]
                bars = ax.bar(
                    x,
                    counts.values,
                    color=element_colors,
                    edgecolor="#315C7D",
                    linewidth=0.65,
                    width=0.68,
                )
                for bar in bars:
                    ax.text(
                        bar.get_x() + bar.get_width() / 2,
                        bar.get_height() + max(counts.values) * 0.022,
                        f"{int(bar.get_height()):,}",
                        ha="center",
                        va="bottom",
                        fontsize=8.5,
                        color="#334155",
                    )
                ax.set_xticks(x)
                ax.set_xticklabels(counts.index)
                ax.set_title("Chalcogen coverage", loc="left", fontsize=10, fontweight="semibold", color="#163A5B", pad=8)
                ax.margins(y=0.12)
                _style_axes(ax)
                fig.tight_layout(pad=0.9)
                st.pyplot(fig, use_container_width=True)
                plt.close(fig)

    if "Eg_eV" in df.columns:
        eg = pd.to_numeric(df["Eg_eV"], errors="coerce").dropna()
        if not eg.empty:
            with st.container(border=True):
                fig, ax = plt.subplots(figsize=(10.6, 2.95))
                fig.patch.set_facecolor("#FFFFFF")
                ax.hist(
                    eg,
                    bins=16,
                    color="#5DA9E9",
                    edgecolor="#315C7D",
                    linewidth=0.65,
                )
                ax.set_xlabel("Eg (eV)", color="#50677B", fontsize=8.5)
                ax.set_title("Eg distribution", loc="left", fontsize=10, fontweight="semibold", color="#163A5B", pad=8)
                _style_axes(ax)
                fig.tight_layout(pad=0.9)
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
        .record-heading {
            margin: .2rem 0 .65rem 0;
        }
        .record-id {
            color: #172033;
            font-size: 1.14rem;
            font-weight: 750;
            letter-spacing: -.01em;
        }
        .record-meta {
            color: #6B7280;
            font-size: .82rem;
            margin-top: .18rem;
        }
        .detail-label {
            color: #6B7280;
            font-size: .76rem;
            font-weight: 600;
            margin: .15rem 0 .35rem 0;
        }
        .detail-section {
            padding-bottom: .85rem;
            border-bottom: 1px solid #E5EAF0;
            margin-bottom: .85rem;
        }
        .detail-section-title {
            color: #163A5B;
            font-size: .92rem;
            font-weight: 750;
            letter-spacing: -.01em;
            margin-bottom: .58rem;
        }
        .detail-metric-title {
            margin-top: .05rem;
            margin-bottom: .45rem;
        }
        .detail-field-grid {
            display: grid;
            gap: .75rem;
        }
        .detail-field-grid-3 {
            grid-template-columns: repeat(3, minmax(0, 1fr));
        }
        .detail-field {
            min-width: 0;
        }
        .detail-field span,
        .compact-fields span {
            display: block;
            color: #788493;
            font-size: .72rem;
            margin-bottom: .13rem;
        }
        .detail-field strong,
        .compact-fields strong {
            display: block;
            color: #172033;
            font-size: .82rem;
            font-weight: 650;
            line-height: 1.38;
            overflow-wrap: anywhere;
        }
        .detail-bottom-grid {
            display: grid;
            grid-template-columns: .9fr 1.1fr;
            gap: .75rem;
            margin-top: .85rem;
        }
        .detail-subcard {
            border: 1px solid #E2E9EF;
            border-radius: 8px;
            background: #FBFCFD;
            padding: .8rem .9rem;
        }
        .compact-fields {
            display: grid;
            grid-template-columns: repeat(2, minmax(0, 1fr));
            gap: .58rem .8rem;
        }
        .provenance-fields {
            grid-template-columns: 1fr;
        }
        div[data-testid="stCode"] pre {
            white-space: pre !important;
            overflow-x: auto !important;
            font-size: .78rem !important;
        }
        @media (max-width: 760px) {
            .detail-field-grid-3,
            .detail-bottom-grid,
            .compact-fields {
                grid-template-columns: 1fr;
            }
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

    structure_query = ""
    structure_mode = "Exact"
    minimum_similarity = 0.40
    structure_top_n = 50

    with st.expander("Structure search", expanded=False):
        st.caption(
            "Search standardized curated structures by exact identity, graph substructure "
            "or Morgan/Tanimoto fingerprint similarity."
        )
        structure_query = st.text_input(
            "Query SMILES",
            placeholder="e.g. c1ccsc1",
            key="database_structure_query",
        ).strip()

        sc1, sc2 = st.columns([1, 1])
        with sc1:
            structure_mode = st.selectbox(
                "Search mode",
                ["Exact", "Substructure", "Similarity"],
                key="database_structure_mode",
            )
        with sc2:
            if structure_mode == "Similarity":
                structure_top_n = st.number_input(
                    "Maximum results",
                    min_value=1,
                    max_value=250,
                    value=50,
                    step=1,
                    key="database_structure_top_n",
                )

        if structure_mode == "Similarity":
            minimum_similarity = st.slider(
                "Minimum Tanimoto similarity",
                min_value=0.0,
                max_value=1.0,
                value=0.40,
                step=0.05,
                key="database_structure_threshold",
            )
            st.caption(
                f"Morgan fingerprint · radius {MORGAN_RADIUS} · "
                f"{MORGAN_N_BITS} bits · {SIMILARITY_METRIC}. "
                "Structural similarity does not imply equivalent electronic properties."
            )

    query = st.text_input(
        "Search",
        placeholder="ID, molecule/system, structure, method, source or reference",
    ).strip()

    f1, f2, f3, f4, f5 = st.columns(5)

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
    with f5:
        reference_column = available_reference_column(df)
        if reference_column:
            reference_value = st.selectbox(
                "Reference / DOI",
                ["All", "Available", "Missing"],
                key="browse_reference",
                help=f"Filtering uses the {reference_column} field in this database build.",
            )
        else:
            st.selectbox(
                "Reference / DOI",
                ["Not included in this release"],
                disabled=True,
                key="browse_reference_unavailable",
                help="No dedicated reference/DOI field is present in this database release.",
            )
            reference_value = "All"

    filtered = _apply_search(df, query)

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

    if reference_value != "All":
        ref_mask = reference_mask(filtered)
        filtered = filtered[ref_mask] if reference_value == "Available" else filtered[~ref_mask]

    if structure_query:
        try:
            filtered, standardized_query = structure_search(
                filtered,
                query_smiles=structure_query,
                mode=structure_mode,
                minimum_similarity=minimum_similarity,
                top_n=int(structure_top_n),
            )
            st.caption(
                f"Structure query · {structure_mode} · standardized as "
                f"{standardized_query['canonical_smiles']} · {len(filtered):,} records"
            )
        except ValueError as error:
            st.error(str(error))
            filtered = filtered.iloc[0:0].copy()

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

    st.markdown('<div class="section-rule-title">Data coverage</div>', unsafe_allow_html=True)
    coverage = coverage_table(df)
    if not coverage.empty:
        st.dataframe(
            coverage,
            use_container_width=True,
            hide_index=True,
            column_config={
                "Field": st.column_config.TextColumn(width="medium"),
                "Available": st.column_config.NumberColumn(format="%d"),
                "Missing": st.column_config.NumberColumn(format="%d"),
                "Coverage (%)": st.column_config.ProgressColumn(
                    min_value=0,
                    max_value=100,
                    format="%.1f%%",
                ),
            },
        )

    quality = database_quality_summary(df)
    st.caption(
        f"Quality checks · Missing record IDs: {quality['Missing Record IDs']:,} · "
        f"Duplicate record IDs: {quality['Duplicate Record IDs']:,} · "
        f"Records in repeated standardized structures: {quality['Repeated Standardized Structures']:,}"
    )

    st.divider()
    _render_statistics(df)

    if is_preview:
        st.info("A preview dataset is loaded in this build.")
    if is_preview and source_path:
        st.caption(f"Preview source: {source_path}")
