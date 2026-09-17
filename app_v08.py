import os
import pandas as pd
import streamlit as st

from rdkit import Chem
from rdkit.Chem import Draw

from descriptor_engine import run_molecular_descriptor_platform
from analysis_panel import display_analysis_panel
from similarity_panel import display_similarity_panel
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
        "Heavy Chalcogen Count", "Chalcogen Fraction",
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
        st.info("Gösterilecek geçerli molekül bulunamadı.")
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
                        st.image(Draw.MolToImage(mol, size=(360, 280)), use_container_width=True)
                    st.write(f"**Formül:** {row['Molecular Formula']}")
                    if "Chalcogen Type" in row:
                        st.write(f"**Chalcogen:** {row['Chalcogen Type']}")
                    st.code(row["Canonical SMILES"], language=None)


sample_df = pd.DataFrame({
    "Molecule_ID": ["S_001", "SE_001", "TE_001", "MIX_001"],
    "SMILES": ["c1ccsc1", "c1cc[se]c1", "c1cc[te]c1", "C[Te]c1ccsc1"],
})
sample_csv = sample_df.to_csv(index=False).encode("utf-8")


with st.sidebar:
    st.title("🧬 Platform Ayarları")
    project_name = st.text_input("Proje adı", value="chalcogen_project")
    selected_groups = st.multiselect(
        "Deskriptör grupları",
        options=list(DESCRIPTOR_GROUPS.keys()),
        default=["Basic", "Structural", "Chalcogen Core"],
    )
    max_molecule_cards = st.slider("Gösterilecek molekül sayısı", 3, 30, 12, 3)
    show_molecule_cards = st.checkbox("Molekül yapılarını göster", value=True)
    st.download_button(
        "Örnek CSV indir", sample_csv, "sample_chalcogen_molecules.csv", "text/csv",
        use_container_width=True,
    )
    st.caption("Girdi dosyası Molecule_ID ve SMILES sütunlarını içermelidir.")


st.markdown(
    """
    <div class="hero">
      <h1>🧬 Molecular Descriptor Platform</h1>
      <p>Chalcogen-focused molecular database preparation, interpretable descriptors, fingerprints and similarity analysis.</p>
    </div>
    """,
    unsafe_allow_html=True,
)

uploaded_file = st.file_uploader("Molekül CSV dosyanızı yükleyin", type=["csv"])

if uploaded_file is None:
    c1, c2, c3 = st.columns(3)
    c1.info("1. CSV yükle ve yapıları doğrula")
    c2.info("2. Scientific Core descriptor'larını hesapla")
    c3.info("3. Analiz, karşılaştırma ve export")

if uploaded_file is not None:
    try:
        input_df = pd.read_csv(uploaded_file)
        st.subheader("Yüklenen Veri")
        st.dataframe(input_df.head(20), use_container_width=True, hide_index=True)

        if not selected_groups:
            st.warning("En az bir deskriptör grubu seçmelisiniz.")

        if st.button(
            "🧪 Scientific Core Analizini Çalıştır",
            type="primary",
            disabled=not selected_groups,
            use_container_width=True,
        ):
            with st.spinner("Moleküller işleniyor..."):
                results = run_molecular_descriptor_platform(
                    input_df=input_df,
                    project_name=project_name,
                )

            valid_df = results["valid_df"]
            invalid_df = results["invalid_df"]
            summary_df = results["summary_df"]
            filtered_valid_df = filter_descriptor_columns(valid_df, selected_groups)

            total_records = int(summary_df.loc[0, "Total Records"])
            valid_count = int(summary_df.loc[0, "Valid Molecules"])
            invalid_count = int(summary_df.loc[0, "Invalid Molecules"])
            duplicate_count = int(summary_df.loc[0, "Duplicate Molecules"])
            success_rate = float(summary_df.loc[0, "Success Rate (%)"])

            m1, m2, m3, m4, m5 = st.columns(5)
            m1.metric("Toplam Kayıt", total_records)
            m2.metric("Geçerli", valid_count)
            m3.metric("Geçersiz", invalid_count)
            m4.metric("Duplicate", duplicate_count)
            m5.metric("Başarı", f"{success_rate:.2f}%")

            tabs = st.tabs([
                "📋 Genel Bakış",
                "🧪 Scientific Core",
                "🧬 Molekül Yapıları",
                "📊 Veri Analizi",
                "⚖️ Karşılaştırma",
                "🔎 Benzer Molekül Arama",
                "⬇️ İndirmeler",
            ])

            with tabs[0]:
                st.subheader("Seçili Moleküler Deskriptörler")
                st.dataframe(filtered_valid_df, use_container_width=True, hide_index=True)
                if not invalid_df.empty:
                    st.subheader("Geçersiz SMILES Kayıtları")
                    st.dataframe(invalid_df, use_container_width=True, hide_index=True)

            with tabs[1]:
                display_scientific_core_panel(valid_df)

            with tabs[2]:
                if show_molecule_cards:
                    display_molecule_gallery(valid_df, max_molecule_cards)
                else:
                    st.info("Molekül yapıları sol menüden kapatıldı.")

            with tabs[3]:
                display_analysis_panel(valid_df)

            with tabs[4]:
                display_similarity_panel(valid_df)

            with tabs[5]:
                display_similarity_search_panel(valid_df)

            with tabs[6]:
                selected_csv = filtered_valid_df.to_csv(index=False).encode("utf-8")
                st.download_button(
                    "Seçili Deskriptörleri CSV Olarak İndir",
                    selected_csv,
                    f"{project_name}_selected_descriptors.csv",
                    "text/csv",
                    use_container_width=True,
                )
                with open(results["zip_file"], "rb") as file:
                    zip_data = file.read()
                st.download_button(
                    "Tüm Sonuçları ZIP Olarak İndir",
                    zip_data,
                    os.path.basename(results["zip_file"]),
                    "application/zip",
                    use_container_width=True,
                )

    except Exception as error:
        st.error(f"Dosya işlenemedi: {error}")

st.caption("Molecular Descriptor Platform — v0.8 Scientific Core development branch")
