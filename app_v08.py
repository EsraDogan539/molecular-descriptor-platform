import os
import pandas as pd
import streamlit as st

from rdkit import Chem
from rdkit.Chem import Draw

from descriptor_engine import run_molecular_descriptor_platform
from database_metadata import add_database_export
from database_browser import display_database_browser
from similarity_search import display_similarity_search_panel
from scientific_panel import display_scientific_core_panel


st.set_page_config(
    page_title="Chalcogen Molecular Database",
    page_icon=None,
    layout="wide",
    initial_sidebar_state="collapsed",
)

st.markdown(
    """
    <style>
    :root {
        --ink: #111827;
        --muted: #667085;
        --line: #E5E7EB;
        --soft: #F8FAFC;
        --accent: #315F86;
    }
    .block-container {
        max-width: 1240px;
        padding-top: 2.6rem;
        padding-bottom: 2.5rem;
    }
    h1, h2, h3, h4 {
        color: var(--ink);
        letter-spacing: -0.01em;
    }
    h1 { font-size: 2.15rem !important; font-weight: 650 !important; }
    h2 { font-size: 1.45rem !important; font-weight: 650 !important; }
    h3 { font-size: 1.10rem !important; font-weight: 650 !important; }
    p, label, .stCaption { color: var(--ink); }
    .academic-header {
        display: flex;
        align-items: center;
        justify-content: space-between;
        border-bottom: 1px solid var(--line);
        padding: .2rem 0 .75rem 0;
        margin-bottom: .7rem;
        gap: 1.5rem;
    }
    .academic-brand {
        color: var(--ink) !important;
        text-decoration: none !important;
        font-size: 1.02rem;
        font-weight: 700;
        white-space: nowrap;
    }
    .academic-nav {
        display: flex;
        gap: 1.35rem;
        align-items: center;
        font-size: .90rem;
    }
    .academic-nav a {
        color: #374151 !important;
        text-decoration: none !important;
        padding: .2rem 0 .45rem 0;
        border-bottom: 2px solid transparent;
    }
    .academic-nav a:hover {
        color: var(--ink) !important;
        border-bottom-color: #AAB4C0;
    }
    .academic-nav a.active {
        color: var(--ink) !important;
        font-weight: 650;
        border-bottom-color: var(--accent);
    }
    .home-intro {
        max-width: 720px;
        margin: 1.95rem auto 1.15rem auto;
        text-align: center;
    }
    .home-intro h1 {
        margin-bottom: .55rem;
    }
    .home-intro p {
        color: var(--muted);
        font-size: 1.02rem;
        line-height: 1.6;
    }
    .metric-row {
        display: flex;
        justify-content: center;
        align-items: stretch;
        gap: 0;
        margin: 1.15rem auto 1.5rem auto;
        max-width: 820px;
    }
    .metric-item {
        min-width: 155px;
        padding: .05rem 1.15rem;
        text-align: center;
        border-right: 1px solid var(--line);
    }
    .metric-item:last-child { border-right: 0; }
    .metric-number {
        display: block;
        color: var(--ink);
        font-size: 1.13rem;
        font-weight: 700;
        line-height: 1.15;
    }
    .metric-label {
        display: block;
        color: var(--muted);
        font-size: .80rem;
        margin-top: .28rem;
        white-space: nowrap;
    }
    .quiet-note {
        color: var(--muted);
        font-size: .82rem;
        text-align: center;
        margin-top: 1.25rem;
    }
    .site-footer {
        border-top: 1px solid var(--line);
        margin-top: 2.2rem;
        padding-top: .9rem;
        color: var(--muted);
        font-size: .80rem;
        display: flex;
        justify-content: space-between;
        gap: 1rem;
        flex-wrap: wrap;
    }
    .site-footer a {
        color: #4B5563 !important;
        text-decoration: none !important;
        margin-left: .85rem;
    }
    .site-footer a:hover { text-decoration: underline !important; }
    @media (max-width: 760px) {
        .academic-header { align-items: flex-start; flex-direction: column; }
        .academic-nav { gap: .9rem; flex-wrap: wrap; }
        .metric-row { flex-wrap: wrap; }
        .metric-item { min-width: 50%; margin-bottom: .8rem; }
    }
    div[data-testid="stMetric"] {
        background: transparent;
        border: 0;
        padding: .25rem 0;
    }
    div[data-testid="stMetric"] [data-testid="stMetricValue"] {
        color: var(--ink);
        font-size: 1.25rem;
    }
    div[data-testid="stMetric"] [data-testid="stMetricLabel"] {
        color: var(--muted);
        font-size: .82rem;
    }
    div[data-testid="stDataFrame"] {
        border: 1px solid var(--line);
        border-radius: 4px;
        overflow: hidden;
    }
    div[data-testid="stFileUploader"] {
        border: 1px solid var(--line);
        border-radius: 4px;
        padding: .6rem;
    }
    div[data-testid="stExpander"] {
        border: 1px solid var(--line);
        border-radius: 4px;
    }
    .stButton > button, .stDownloadButton > button {
        border-radius: 4px !important;
        font-weight: 600 !important;
        box-shadow: none !important;
    }
    .stButton > button[kind="primary"] {
        background: var(--accent) !important;
        border-color: var(--accent) !important;
        color: white !important;
    }
    a[data-testid="stBaseLinkButton-primary"],
    a[data-testid="stLinkButton"] {
        background: var(--accent) !important;
        border-color: var(--accent) !important;
        color: white !important;
    }
    a[data-testid="stBaseLinkButton-primary"]:hover,
    a[data-testid="stLinkButton"]:hover {
        background: #274F70 !important;
        border-color: #274F70 !important;
        color: white !important;
    }
    button[data-baseweb="tab"] {
        font-size: .92rem;
    }
    hr {
        border-color: var(--line);
    }
    footer { visibility: hidden; }
    </style>
    """,
    unsafe_allow_html=True,
)


DESCRIPTOR_GROUPS = {
    "Basic": [
        "Molecule_ID", "Original SMILES", "Canonical SMILES", "InChIKey",
        "Molecular Formula", "Molecular Weight", "Exact Molecular Weight",
        "LogP", "TPSA", "H-Bond Donors", "H-Bond Acceptors",
    ],
    "Structural": [
        "Molecule_ID", "Canonical SMILES", "Rotatable Bonds", "Ring Count",
        "Aromatic Ring Count", "Aliphatic Ring Count", "Heavy Atom Count",
        "Heteroatom Count", "Fraction Csp3", "Formal Charge",
    ],
    "Topological": [
        "Molecule_ID", "Canonical SMILES", "Bertz CT", "Balaban J",
    ],
    "Element Counts": [
        "Molecule_ID", "Canonical SMILES", "Carbon Count", "Nitrogen Count",
        "Oxygen Count", "Sulfur Count", "Selenium Count", "Tellurium Count",
        "Phosphorus Count", "Halogen Count",
    ],
    "Chalcogen Core": [
        "Molecule_ID", "Canonical SMILES", "Chalcogen Type",
        "Target Chalcogen Count", "Target Chalcogen Fraction",
        "Aromatic Chalcogen Count", "NonAromatic Chalcogen Count",
        "Mixed Chalcogen Flag", "Ring Incorporated Chalcogen Count",
        "Chalcogen-C Bond Count", "Chalcogen-Heteroatom Bond Count",
        "Aromatic Neighbor Count", "Conjugated Bond Count",
        "Aromatic Bond Fraction", "Conjugated Atom Fraction",
        "Heteroaromatic Ring Count", "Duplicate Flag",
    ],
}


def filter_descriptor_columns(dataframe, selected_groups):
    selected = []
    for group in selected_groups:
        selected.extend(DESCRIPTOR_GROUPS.get(group, []))
    selected = list(dict.fromkeys(selected))
    existing = [c for c in selected if c in dataframe.columns]
    return dataframe[existing]


def display_molecule_gallery(valid_df, max_molecules=12):
    if valid_df.empty:
        st.info("No valid molecules are available to display.")
        return

    count = min(len(valid_df), max_molecules)
    for start in range(0, count, 4):
        cols = st.columns(4)
        for j in range(4):
            idx = start + j
            if idx >= count:
                break
            row = valid_df.iloc[idx]
            with cols[j]:
                st.markdown(f"**{row['Molecule_ID']}**")
                mol = Chem.MolFromSmiles(row["Canonical SMILES"])
                if mol is not None:
                    st.image(Draw.MolToImage(mol, size=(300, 220)), use_container_width=True)
                st.caption(
                    f"{row.get('Molecular Formula', '—')} · "
                    f"{row.get('Chalcogen Type', '—')}"
                )
                st.code(row["Canonical SMILES"], language=None)


sample_df = pd.DataFrame({
    "Molecule_ID": ["S_001", "SE_001", "TE_001", "MIX_001"],
    "SMILES": ["c1ccsc1", "c1cc[se]c1", "c1cc[te]c1", "C[Te]c1ccsc1"],
    "HOMO_eV": [-5.20, -5.10, -4.95, -5.00],
    "LUMO_eV": [-2.80, -2.85, -2.90, -2.95],
    "Eg_eV": [2.40, 2.25, 2.05, 2.05],
    "Property_Source": ["DFT", "DFT", "DFT", "DFT"],
    "DOI_or_Reference": ["example", "example", "example", "example"],
})
sample_csv = sample_df.to_csv(index=False).encode("utf-8")


page = str(st.query_params.get("page", "home")).lower()
if page not in {"home", "database", "analyze", "documentation", "about"}:
    page = "home"

nav_items = [
    ("Database", "database"),
    ("Analyze", "analyze"),
    ("Documentation", "documentation"),
    ("About", "about"),
]
nav_html = "".join(
    f'<a class="{"active" if page == key else ""}" href="?page={key}">{label}</a>'
    for label, key in nav_items
)

st.markdown(
    f"""
    <div class="academic-header">
      <a class="academic-brand" href="?page=home">Chalcogen Molecular Database</a>
      <nav class="academic-nav">{nav_html}</nav>
    </div>
    """,
    unsafe_allow_html=True,
)

if page == "home":
    st.markdown(
        """
        <div class="home-intro">
          <h1>Chalcogen Molecular Database</h1>
          <p>
            Curated molecular records with standardized structures
            and S/Se/Te-specific annotations.
          </p>
        </div>
        <div class="metric-row">
          <div class="metric-item"><span class="metric-number">3,360</span><span class="metric-label">Records</span></div>
          <div class="metric-item"><span class="metric-number">3,145</span><span class="metric-label">Core S/Se/Te</span></div>
          <div class="metric-item"><span class="metric-number">2,983</span><span class="metric-label">Unique structures</span></div>
          <div class="metric-item"><span class="metric-number">3,353</span><span class="metric-label">Eg values</span></div>
        </div>
        """,
        unsafe_allow_html=True,
    )

    left, c1, c2, right = st.columns([1.7, 1, 1, 1.7])
    with c1:
        st.link_button("Browse database", "?page=database", type="primary", use_container_width=True)
    with c2:
        st.link_button("Analyze your dataset", "?page=analyze", use_container_width=True)

    st.markdown(
        """
        <div class="quiet-note">
          Database v1 &nbsp;&middot;&nbsp; Structure standardization with RDKit
          &nbsp;&middot;&nbsp; Record-level provenance retained
        </div>
        <div class="site-footer">
          <span>Chalcogen Molecular Database · v1</span>
          <span>
            <a href="?page=documentation">Documentation</a>
            <a href="https://github.com/EsraDogan539/molecular-descriptor-platform" target="_blank">GitHub</a>
            <a href="?page=about">About</a>
          </span>
        </div>
        """,
        unsafe_allow_html=True,
    )
    st.stop()


if page == "database":
    display_database_browser()
    st.divider()
    st.caption(
        "Database v1 · Scientific Core v0.8 · "
        "Curated records are read-only in the public browser."
    )
    st.stop()


if page == "documentation":
    st.header("Documentation")
    st.caption("Database scope, structure handling and user-analysis workflow.")

    st.markdown("### Database scope")
    st.write(
        "The public database contains curated development and external records. "
        "Exact standardized structures are displayed only where supported by the source."
    )

    st.markdown("### Molecular identity")
    st.write(
        "Structure-complete records are standardized with RDKit and represented by "
        "canonical SMILES, InChI and InChIKey. Repeated standardized structures are "
        "retained to preserve record-level provenance."
    )

    st.markdown("### Analyze your dataset")
    st.write(
        "Upload a CSV containing Molecule_ID and SMILES. The platform validates the "
        "structures, calculates general and S/Se/Te-aware descriptors, and keeps user "
        "data separate from the curated publication database."
    )

    st.markdown("### Citation and data release")
    st.write(
        "The recommended citation and permanent dataset DOI will be added with the "
        "archived publication release."
    )
    st.stop()


if page == "about":
    st.header("About")
    st.write(
        "Chalcogen Molecular Database is a research resource for curated molecular "
        "records, electronic-property data and interpretable S/Se/Te structural annotations."
    )
    st.write(
        "The platform is designed for transparent database exploration and reproducible "
        "downstream cheminformatics workflows."
    )
    st.caption("Database v1 · Scientific Core v0.8")
    st.stop()


st.header("Analyze Your Dataset")
st.caption(
    "Upload molecular data to validate structures and calculate general and "
    "chalcogen-aware descriptors. User uploads remain separate from the curated database."
)

upload_left, upload_right = st.columns([2.1, 1])
with upload_left:
    uploaded_file = st.file_uploader(
        "Choose a CSV file",
        type=["csv"],
        help="Required columns: Molecule_ID and SMILES.",
    )
with upload_right:
    st.markdown("**Input format**")
    st.caption("Required: Molecule_ID, SMILES")
    st.caption("Optional property and provenance fields are retained.")
    st.download_button(
        "Download example CSV",
        sample_csv,
        "sample_chalcogen_database.csv",
        "text/csv",
        use_container_width=True,
    )

if uploaded_file is None:
    st.info("Choose a CSV file to begin. Analysis options appear after the file is loaded.")
    st.stop()

try:
    input_df = pd.read_csv(uploaded_file)
    st.markdown("### Input preview")
    st.dataframe(input_df.head(12), use_container_width=True, hide_index=True)

    with st.expander("Analysis settings", expanded=False):
        project_name = st.text_input("Project name", value="chalcogen_project")
        selected_groups = st.multiselect(
            "Descriptor groups",
            options=list(DESCRIPTOR_GROUPS.keys()),
            default=["Basic", "Structural", "Chalcogen Core"],
        )
        max_molecule_cards = st.slider("Maximum structures to display", 4, 24, 12, 4)
    show_molecule_cards = True

    if not selected_groups:
        st.warning("Select at least one descriptor group.")

    if st.button(
        "Run analysis",
        type="primary",
        disabled=not selected_groups,
    ):
        with st.spinner("Validating structures and calculating descriptors..."):
            results = run_molecular_descriptor_platform(
                input_df=input_df,
                project_name=project_name,
            )
            results = add_database_export(
                input_df=input_df,
                results=results,
                project_name=project_name,
            )

        valid_df = results["valid_df"]
        invalid_df = results["invalid_df"]
        database_df = results["database_df"]
        summary_df = results["summary_df"]
        filtered_valid_df = filter_descriptor_columns(valid_df, selected_groups)

        total_records = int(summary_df.loc[0, "Total Records"])
        valid_count = int(summary_df.loc[0, "Valid Molecules"])
        invalid_count = int(summary_df.loc[0, "Invalid Molecules"])
        duplicate_count = int(summary_df.loc[0, "Duplicate Molecules"])

        st.caption(
            f"{total_records:,} submitted · {valid_count:,} valid · "
            f"{invalid_count:,} invalid · {duplicate_count:,} repeated structures"
        )

        tabs = st.tabs([
            "Records",
            "Descriptors",
            "Structures",
            "Similarity",
            "Export",
        ])

        with tabs[0]:
            st.markdown("### Records")
            st.caption(
                "Standardized identifiers, retained source fields and quality flags."
            )
            st.dataframe(database_df, use_container_width=True, hide_index=True)
            if not invalid_df.empty:
                st.markdown("#### Invalid structures")
                st.dataframe(invalid_df, use_container_width=True, hide_index=True)

        with tabs[1]:
            st.markdown("### Descriptors")
            chalcogen_columns = [
                c for c in DESCRIPTOR_GROUPS["Chalcogen Core"]
                if c in valid_df.columns
            ]
            general_columns = [
                c for c in filtered_valid_df.columns
                if c not in chalcogen_columns
            ]
            descriptor_tabs = st.tabs(["Chalcogen-aware", "General"])
            with descriptor_tabs[0]:
                st.dataframe(
                    valid_df[chalcogen_columns],
                    use_container_width=True,
                    hide_index=True,
                )
            with descriptor_tabs[1]:
                if general_columns:
                    st.dataframe(
                        filtered_valid_df[general_columns],
                        use_container_width=True,
                        hide_index=True,
                    )
                else:
                    st.info("No general descriptor groups are selected.")
            with st.expander("Descriptor summary", expanded=False):
                display_scientific_core_panel(valid_df)

        with tabs[2]:
            st.markdown("### Structures")
            display_molecule_gallery(valid_df, max_molecule_cards)

        with tabs[3]:
            st.markdown("### Similarity Search")
            st.caption(
                "Exploratory fingerprint-based molecular similarity. "
                "Similarity does not imply equivalent electronic behavior."
            )
            display_similarity_search_panel(valid_df)

        with tabs[4]:
            st.markdown("### Export")
            selected_csv = filtered_valid_df.to_csv(index=False).encode("utf-8")
            database_csv = database_df.to_csv(index=False).encode("utf-8")

            r1c1, r1c2 = st.columns([3, 1])
            with r1c1:
                st.markdown("**Processed records**")
                st.caption("Standardized records with retained provenance and quality fields.")
            with r1c2:
                st.download_button(
                    "Download CSV",
                    database_csv,
                    f"{project_name}_processed_database.csv",
                    "text/csv",
                    use_container_width=True,
                    key="download_processed",
                )

            st.divider()
            r2c1, r2c2 = st.columns([3, 1])
            with r2c1:
                st.markdown("**Descriptor table**")
                st.caption("Descriptor groups selected for this analysis.")
            with r2c2:
                st.download_button(
                    "Download CSV",
                    selected_csv,
                    f"{project_name}_selected_descriptors.csv",
                    "text/csv",
                    use_container_width=True,
                    key="download_descriptors",
                )

            st.divider()
            with open(results["zip_file"], "rb") as file:
                zip_data = file.read()
            r3c1, r3c2 = st.columns([3, 1])
            with r3c1:
                st.markdown("**Complete package**")
                st.caption("Processed data, descriptors, fingerprints and run summary.")
            with r3c2:
                st.download_button(
                    "Download ZIP",
                    zip_data,
                    os.path.basename(results["zip_file"]),
                    "application/zip",
                    use_container_width=True,
                    key="download_package",
                )

except Exception as error:
    st.error(f"Dataset could not be processed: {error}")

st.divider()
st.caption(
    "Chalcogen Molecular Database · Database v1 · Scientific Core v0.8"
)
