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
    page_title="Molecular Descriptor Platform v0.8",
    page_icon="🧬",
    layout="wide",
    initial_sidebar_state="expanded",
)

st.markdown(
    """
    <style>
    .block-container {max-width: 1380px; padding-top: 2rem; padding-bottom: 3rem;}
    .hero {padding: 2rem; border: 1px solid #E5E9F2; border-radius: 24px;
           background: linear-gradient(135deg,#FFFFFF,#F5F7FF); margin-bottom: 1.5rem;}
    .hero h1 {margin:0 0 .4rem 0; font-size:2.6rem;}
    .hero p {color:#667085; margin:0; font-size:1.03rem;}
    div[data-testid="stMetric"] {border:1px solid #E7EAF2; border-radius:18px; padding:1rem; background:white;}
    div[data-testid="stDataFrame"] {border:1px solid #E7EAF2; border-radius:16px; overflow:hidden;}
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
    for start in range(0, count, 3):
        cols = st.columns(3)
        for j in range(3):
            idx = start + j
            if idx >= count:
                break
            row = valid_df.iloc[idx]
            with cols[j]:
                with st.container(border=True):
                    st.markdown(f"### {row['Molecule_ID']}")
                    mol = Chem.MolFromSmiles(row["Canonical SMILES"])
                    if mol is not None:
                        st.image(
                            Draw.MolToImage(mol, size=(360, 280)),
                            use_container_width=True,
                        )
                    st.write(f"**Molecular formula:** {row['Molecular Formula']}")
                    st.write(f"**Chalcogen class:** {row.get('Chalcogen Type', '—')}")
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


with st.sidebar:
    st.title("🧬 Molecular Descriptor Platform")
    platform_mode = st.radio(
        "Mode",
        ["Our Curated Database", "Analyze Your Dataset"],
        index=0,
    )

    st.divider()

    if platform_mode == "Analyze Your Dataset":
        project_name = st.text_input("Project name", value="chalcogen_project")
        selected_groups = st.multiselect(
            "Descriptor groups",
            options=list(DESCRIPTOR_GROUPS.keys()),
            default=["Basic", "Structural", "Chalcogen Core"],
        )
        max_molecule_cards = st.slider(
            "Molecule cards",
            3,
            30,
            12,
            3,
        )
        show_molecule_cards = st.checkbox(
            "Show molecule structures",
            value=True,
        )

        st.divider()
        st.download_button(
            "Download example CSV",
            sample_csv,
            "sample_chalcogen_database.csv",
            "text/csv",
            use_container_width=True,
        )
        st.caption(
            "Uploaded data never enter the curated publication database automatically. "
            "Required columns: Molecule_ID and SMILES."
        )
    else:
        project_name = "chalcogen_project"
        selected_groups = ["Basic", "Structural", "Chalcogen Core"]
        max_molecule_cards = 12
        show_molecule_cards = True
        st.caption(
            "Database v1 is read-only in this browser. Use Analyze Your Dataset "
            "for user-supplied structures and descriptor calculations."
        )


st.markdown(
    """
    <div class="hero">
      <h1>🧬 Molecular Descriptor Platform</h1>
      <p>A curated chalcogen-focused molecular database with interpretable structural annotations and user-side descriptor analysis.</p>
    </div>
    """,
    unsafe_allow_html=True,
)


if platform_mode == "Our Curated Database":
    display_database_browser()
    st.caption("Molecular Descriptor Platform — Database v1 / Scientific Core v0.8")
    st.stop()


st.subheader("Analyze Your Dataset")
st.caption(
    "Upload a CSV with Molecule_ID and SMILES. The platform validates structures, "
    "calculates general and chalcogen-aware descriptors, and returns a quality-controlled "
    "analysis package. Uploaded records never modify the curated publication database."
)

uploaded_file = st.file_uploader("Upload molecular CSV", type=["csv"])

if uploaded_file is None:
    c1, c2, c3 = st.columns(3)
    c1.info("1. Upload and validate molecular structures")
    c2.info("2. Calculate selected descriptor layers")
    c3.info("3. Review quality flags and export results")

if uploaded_file is not None:
    try:
        input_df = pd.read_csv(uploaded_file)
        st.subheader("Uploaded data")
        st.dataframe(input_df.head(20), use_container_width=True, hide_index=True)

        if not selected_groups:
            st.warning("Select at least one descriptor group.")

        if st.button(
            "Run descriptor analysis",
            type="primary",
            disabled=not selected_groups,
            use_container_width=True,
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
            success_rate = float(summary_df.loc[0, "Success Rate (%)"])

            m1, m2, m3, m4, m5 = st.columns(5)
            m1.metric("Records", total_records)
            m2.metric("Valid", valid_count)
            m3.metric("Invalid", invalid_count)
            m4.metric("Repeated structures", duplicate_count)
            m5.metric("Validation rate", f"{success_rate:.2f}%")

            tabs = st.tabs([
                "🗂️ Quality & Records",
                "🧪 Descriptor Layer",
                "🧬 Structures",
                "🔎 Similarity Search",
                "⬇️ Export",
            ])

            with tabs[0]:
                st.subheader("Quality-controlled records")
                st.caption(
                    "Canonical identifiers, user-supplied property fields, provenance metadata "
                    "and quality flags are retained together."
                )
                st.dataframe(database_df, use_container_width=True, hide_index=True)
                if not invalid_df.empty:
                    st.markdown("#### Invalid or missing structures")
                    st.dataframe(invalid_df, use_container_width=True, hide_index=True)

            with tabs[1]:
                st.subheader("Selected descriptor layer")
                st.caption(
                    "General molecular descriptors and chalcogen-aware annotations are calculated "
                    "from the standardized structure."
                )

                chalcogen_columns = [
                    c for c in DESCRIPTOR_GROUPS["Chalcogen Core"]
                    if c in valid_df.columns
                ]
                general_columns = [
                    c for c in filtered_valid_df.columns
                    if c not in chalcogen_columns
                ]

                descriptor_tabs = st.tabs([
                    "Chalcogen-aware descriptors",
                    "General descriptors",
                ])

                with descriptor_tabs[0]:
                    st.caption(
                        "S/Se/Te-focused structural annotations used in the publication descriptor layer."
                    )
                    st.dataframe(
                        valid_df[chalcogen_columns],
                        use_container_width=True,
                        hide_index=True,
                    )

                with descriptor_tabs[1]:
                    st.caption(
                        "General molecular descriptors selected from the sidebar."
                    )
                    if general_columns:
                        st.dataframe(
                            filtered_valid_df[general_columns],
                            use_container_width=True,
                            hide_index=True,
                        )
                    else:
                        st.info("No general descriptor groups are selected.")

                with st.expander("Scientific Core summary", expanded=False):
                    display_scientific_core_panel(valid_df)

            with tabs[2]:
                if show_molecule_cards:
                    display_molecule_gallery(valid_df, max_molecule_cards)
                else:
                    st.info("Molecule structure cards are disabled in the sidebar.")

            with tabs[3]:
                st.caption(
                    "Fingerprint-based similarity is provided as an exploratory tool and does not "
                    "change the curated database."
                )
                display_similarity_search_panel(valid_df)

            with tabs[4]:
                selected_csv = filtered_valid_df.to_csv(index=False).encode("utf-8")
                database_csv = database_df.to_csv(index=False).encode("utf-8")

                st.download_button(
                    "Download processed records CSV",
                    database_csv,
                    f"{project_name}_processed_database.csv",
                    "text/csv",
                    use_container_width=True,
                )
                st.download_button(
                    "Download selected descriptors CSV",
                    selected_csv,
                    f"{project_name}_selected_descriptors.csv",
                    "text/csv",
                    use_container_width=True,
                )
                with open(results["zip_file"], "rb") as file:
                    zip_data = file.read()
                st.download_button(
                    "Download complete analysis package (ZIP)",
                    zip_data,
                    os.path.basename(results["zip_file"]),
                    "application/zip",
                    use_container_width=True,
                )

    except Exception as error:
        st.error(f"Dataset could not be processed: {error}")


st.caption("Molecular Descriptor Platform — Database v1 / Scientific Core v0.8")
