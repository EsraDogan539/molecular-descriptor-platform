
import os
import pandas as pd
import streamlit as st

from rdkit import Chem
from rdkit.Chem import Draw

from descriptor_engine import run_molecular_descriptor_platform
from analysis_panel import display_analysis_panel
from similarity_panel import display_similarity_panel
from similarity_search import display_similarity_search_panel


st.set_page_config(
    page_title="Molecular Descriptor Platform",
    page_icon="🧬",
    layout="wide",
    initial_sidebar_state="expanded"
)


st.markdown(
    """
    


<style>
@import url('https://fonts.googleapis.com/css2?family=Inter:wght@400;500;600;700;800&family=Manrope:wght@600;700;800&display=swap');

:root {
    --bg-main: #F6F8FC;
    --bg-soft: #EEF2FF;
    --surface: rgba(255, 255, 255, 0.88);
    --surface-solid: #FFFFFF;

    --primary: #5B5FEF;
    --primary-dark: #4046C9;
    --primary-soft: #E9ECFF;

    --cyan: #35C6E8;
    --teal: #2DD4BF;

    --text-main: #1F2740;
    --text-secondary: #5F6B85;
    --text-muted: #8C96AD;

    --border: #DCE3F0;
    --border-soft: #E8ECF5;

    --success-bg: #DDF5E8;
    --success-text: #1E8A52;

    --warning-bg: #FFF2D6;
    --warning-text: #A66B00;

    --error-bg: #FCE2E2;
    --error-text: #B83232;

    --shadow-soft: 0 10px 30px rgba(37, 49, 90, 0.07);
    --shadow-card: 0 14px 36px rgba(37, 49, 90, 0.09);
}


/* --------------------------------------------------
   GLOBAL
-------------------------------------------------- */

html,
body,
[class*="css"] {
    font-family: "Inter", sans-serif;
}

.stApp {
    background:
        radial-gradient(
            circle at 8% 0%,
            rgba(91, 95, 239, 0.13),
            transparent 28%
        ),
        radial-gradient(
            circle at 92% 5%,
            rgba(53, 198, 232, 0.10),
            transparent 24%
        ),
        linear-gradient(
            180deg,
            #FAFBFF 0%,
            var(--bg-main) 48%,
            #F8FAFF 100%
        );

    color: var(--text-main);
}

.block-container {
    max-width: 1380px;
    padding-top: 2.1rem;
    padding-bottom: 3rem;
}

h1,
h2,
h3,
h4 {
    font-family: "Manrope", sans-serif;
    color: var(--text-main);
    letter-spacing: -0.025em;
}

p,
label,
span {
    color: inherit;
}


/* --------------------------------------------------
   HERO
-------------------------------------------------- */

.main-title {
    font-family: "Manrope", sans-serif;
    font-size: clamp(2.4rem, 5vw, 3.5rem);
    line-height: 1.05;
    font-weight: 800;
    letter-spacing: -0.055em;
    margin-bottom: 0.55rem;

    background: linear-gradient(
        90deg,
        #252B48 0%,
        #4C51D8 55%,
        #258DAD 100%
    );

    -webkit-background-clip: text;
    -webkit-text-fill-color: transparent;
    background-clip: text;
}

.subtitle {
    max-width: 850px;
    font-size: 1.08rem;
    line-height: 1.75;
    color: var(--text-secondary);
    margin-bottom: 1.35rem;
}

.hero-card {
    position: relative;
    overflow: hidden;

    padding: 2.4rem;
    margin-bottom: 1.8rem;

    border-radius: 30px;
    border: 1px solid rgba(210, 219, 247, 0.95);

    background:
        linear-gradient(
            135deg,
            rgba(255, 255, 255, 0.94),
            rgba(240, 244, 255, 0.84)
        );

    box-shadow:
        0 22px 52px rgba(49, 61, 120, 0.10),
        inset 0 1px 0 rgba(255, 255, 255, 0.85);

    backdrop-filter: blur(18px);
    -webkit-backdrop-filter: blur(18px);
}

.hero-card::before {
    content: "";
    position: absolute;
    width: 280px;
    height: 280px;
    right: -80px;
    top: -105px;
    border-radius: 50%;

    background: radial-gradient(
        circle,
        rgba(91, 95, 239, 0.20),
        rgba(53, 198, 232, 0.05) 60%,
        transparent 72%
    );
}

.hero-card::after {
    content: "";
    position: absolute;
    width: 180px;
    height: 180px;
    right: 130px;
    bottom: -110px;
    border-radius: 50%;

    background: radial-gradient(
        circle,
        rgba(45, 212, 191, 0.16),
        transparent 70%
    );
}


/* --------------------------------------------------
   SIDEBAR
-------------------------------------------------- */

section[data-testid="stSidebar"] {
    background:
        linear-gradient(
            180deg,
            rgba(255, 255, 255, 0.98),
            rgba(245, 247, 253, 0.98)
        );

    border-right: 1px solid var(--border-soft);
}

section[data-testid="stSidebar"] > div {
    padding-top: 1.25rem;
}

section[data-testid="stSidebar"] h1,
section[data-testid="stSidebar"] h2,
section[data-testid="stSidebar"] h3 {
    color: var(--text-main);
}

section[data-testid="stSidebar"] hr {
    border-color: var(--border-soft);
}


/* --------------------------------------------------
   FORM FIELDS
-------------------------------------------------- */

div[data-baseweb="input"] > div,
div[data-baseweb="select"] > div {
    border-radius: 14px !important;
    border-color: var(--border) !important;
    background: rgba(255, 255, 255, 0.90) !important;
}

div[data-baseweb="input"]:focus-within > div,
div[data-baseweb="select"]:focus-within > div {
    border-color: var(--primary) !important;

    box-shadow:
        0 0 0 3px rgba(91, 95, 239, 0.12) !important;
}


/* Multiselect etiketleri */

span[data-baseweb="tag"] {
    color: #FFFFFF !important;

    background: linear-gradient(
        135deg,
        var(--primary),
        #7181FF
    ) !important;

    border-radius: 999px !important;
    border: none !important;
}


/* Slider */

div[data-testid="stSlider"] div[role="slider"] {
    background: var(--primary) !important;
    border-color: var(--primary) !important;

    box-shadow:
        0 0 0 5px rgba(91, 95, 239, 0.12);
}

div[data-testid="stSlider"] div[data-baseweb="slider"] > div > div {
    background: linear-gradient(
        90deg,
        var(--primary),
        var(--cyan)
    ) !important;
}


/* Checkbox */

div[data-testid="stCheckbox"] label span:first-child {
    border-radius: 6px !important;
}

div[data-testid="stCheckbox"]
input:checked + div {
    background-color: var(--primary) !important;
    border-color: var(--primary) !important;
}


/* --------------------------------------------------
   FILE UPLOADER
-------------------------------------------------- */

div[data-testid="stFileUploader"] > section {
    padding: 1.25rem !important;

    border-radius: 20px !important;
    border: 1.5px dashed #C9D5F3 !important;

    background:
        linear-gradient(
            135deg,
            rgba(255, 255, 255, 0.92),
            rgba(243, 247, 255, 0.86)
        );

    transition:
        border-color 0.18s ease,
        box-shadow 0.18s ease,
        transform 0.18s ease;
}

div[data-testid="stFileUploader"] > section:hover {
    border-color: var(--primary) !important;

    box-shadow:
        0 12px 28px rgba(91, 95, 239, 0.11);

    transform: translateY(-1px);
}


/* --------------------------------------------------
   PRIMARY BUTTON
-------------------------------------------------- */

.stButton > button {
    width: 100%;
    min-height: 52px;

    border: none;
    border-radius: 16px;

    color: #FFFFFF;
    font-weight: 700;
    font-size: 0.98rem;

    background: linear-gradient(
        90deg,
        var(--primary) 0%,
        #6A7CFF 56%,
        var(--cyan) 100%
    );

    box-shadow:
        0 12px 28px rgba(91, 95, 239, 0.28);

    transition:
        transform 0.18s ease,
        box-shadow 0.18s ease,
        filter 0.18s ease;
}

.stButton > button:hover {
    transform: translateY(-2px);

    box-shadow:
        0 18px 34px rgba(91, 95, 239, 0.34);

    filter: saturate(1.08);
}

.stButton > button:active {
    transform: translateY(0);
}


/* --------------------------------------------------
   DOWNLOAD BUTTON
-------------------------------------------------- */

.stDownloadButton > button {
    min-height: 48px;

    border-radius: 15px;
    border: 1px solid var(--border);

    background: rgba(255, 255, 255, 0.94);

    color: #35405D;
    font-weight: 650;

    box-shadow:
        0 7px 20px rgba(32, 44, 84, 0.05);

    transition:
        border-color 0.18s ease,
        transform 0.18s ease,
        box-shadow 0.18s ease;
}

.stDownloadButton > button:hover {
    transform: translateY(-1px);

    border-color: rgba(91, 95, 239, 0.48);

    box-shadow:
        0 11px 24px rgba(91, 95, 239, 0.10);
}


/* --------------------------------------------------
   INFO CARDS
-------------------------------------------------- */

.info-card {
    min-height: 180px;
    padding: 1.45rem;

    border-radius: 22px;
    border: 1px solid var(--border-soft);

    background: rgba(255, 255, 255, 0.90);

    box-shadow:
        0 12px 30px rgba(37, 49, 90, 0.065);

    transition:
        transform 0.18s ease,
        box-shadow 0.18s ease,
        border-color 0.18s ease;
}

.info-card:hover {
    transform: translateY(-3px);

    border-color: rgba(91, 95, 239, 0.28);

    box-shadow:
        0 18px 38px rgba(37, 49, 90, 0.10);
}


/* --------------------------------------------------
   METRIC CARDS
-------------------------------------------------- */

div[data-testid="stMetric"] {
    position: relative;
    overflow: hidden;

    min-height: 130px;
    padding: 1.2rem;

    border-radius: 22px;
    border: 1px solid var(--border-soft);

    background: rgba(255, 255, 255, 0.92);

    box-shadow:
        0 13px 30px rgba(37, 49, 90, 0.065);

    transition:
        transform 0.18s ease,
        box-shadow 0.18s ease;
}

div[data-testid="stMetric"]::after {
    content: "";

    position: absolute;
    width: 105px;
    height: 105px;

    right: -42px;
    top: -50px;

    border-radius: 50%;

    background: radial-gradient(
        circle,
        rgba(91, 95, 239, 0.15),
        transparent 72%
    );
}

div[data-testid="stMetric"]:hover {
    transform: translateY(-2px);

    box-shadow:
        0 18px 38px rgba(37, 49, 90, 0.09);
}

div[data-testid="stMetricLabel"] {
    color: var(--text-secondary);
    font-weight: 600;
}

div[data-testid="stMetricValue"] {
    color: var(--text-main);

    font-family: "Manrope", sans-serif;
    font-size: 2.25rem;
    font-weight: 800;
}


/* --------------------------------------------------
   TABS
-------------------------------------------------- */

.stTabs [data-baseweb="tab-list"] {
    display: flex;
    flex-wrap: wrap;

    gap: 9px;

    padding: 0.35rem 0;
    margin-bottom: 1rem;

    border-bottom: none;
}

.stTabs [data-baseweb="tab"] {
    min-height: 43px;

    padding: 0 17px;

    border-radius: 999px;
    border: 1px solid var(--border-soft);

    color: var(--text-secondary);
    font-weight: 600;

    background: rgba(255, 255, 255, 0.78);

    transition:
        background 0.18s ease,
        color 0.18s ease,
        transform 0.18s ease,
        box-shadow 0.18s ease;
}

.stTabs [data-baseweb="tab"]:hover {
    color: var(--primary);

    transform: translateY(-1px);

    border-color: rgba(91, 95, 239, 0.28);

    background: rgba(239, 242, 255, 0.88);
}

.stTabs [aria-selected="true"] {
    color: #FFFFFF !important;

    border-color: transparent !important;

    background: linear-gradient(
        90deg,
        var(--primary),
        #7181FF
    ) !important;

    box-shadow:
        0 8px 20px rgba(91, 95, 239, 0.24);
}


/* --------------------------------------------------
   DATAFRAME / TABLE
-------------------------------------------------- */

div[data-testid="stDataFrame"] {
    overflow: hidden;

    border-radius: 19px;
    border: 1px solid var(--border-soft);

    background: var(--surface-solid);

    box-shadow:
        0 11px 28px rgba(37, 49, 90, 0.055);
}


/* --------------------------------------------------
   MOLECULE CARDS / BORDERED CONTAINERS
-------------------------------------------------- */

div[data-testid="stVerticalBlockBorderWrapper"] {
    overflow: hidden;

    border-radius: 22px !important;
    border-color: var(--border-soft) !important;

    background: rgba(255, 255, 255, 0.90);

    box-shadow:
        0 12px 30px rgba(37, 49, 90, 0.06);

    transition:
        transform 0.18s ease,
        box-shadow 0.18s ease,
        border-color 0.18s ease;
}

div[data-testid="stVerticalBlockBorderWrapper"]:hover {
    transform: translateY(-3px);

    border-color: rgba(91, 95, 239, 0.24) !important;

    box-shadow:
        0 18px 40px rgba(37, 49, 90, 0.10);
}


/* --------------------------------------------------
   CODE / SMILES BOX
-------------------------------------------------- */

div[data-testid="stCodeBlock"] {
    border-radius: 13px;

    border: 1px solid var(--border-soft);

    background: #F7F9FE;
}


/* --------------------------------------------------
   ALERTS
-------------------------------------------------- */

div[data-testid="stAlert"] {
    border-radius: 18px;
    border-width: 1px;

    box-shadow:
        0 8px 20px rgba(37, 49, 90, 0.045);
}


/* --------------------------------------------------
   DIVIDER
-------------------------------------------------- */

hr {
    margin-top: 1.1rem;
    margin-bottom: 1.1rem;

    border: none;
    border-top: 1px solid var(--border-soft);
}


/* --------------------------------------------------
   FOOTER
-------------------------------------------------- */

.footer {
    padding-top: 2.8rem;
    padding-bottom: 1.2rem;

    text-align: center;

    color: var(--text-muted);
    font-size: 0.94rem;
    line-height: 1.7;
}


/* --------------------------------------------------
   RESPONSIVE
-------------------------------------------------- */

@media (max-width: 900px) {
    .block-container {
        padding-top: 1.25rem;
    }

    .hero-card {
        padding: 1.55rem;
        border-radius: 22px;
    }

    .main-title {
        font-size: 2.2rem;
    }

    .stTabs [data-baseweb="tab"] {
        padding: 0 12px;
        font-size: 0.88rem;
    }
}


/* --------------------------------------------------
   CLEANUP
-------------------------------------------------- */

#MainMenu {
    visibility: hidden;
}

footer {
    visibility: hidden;
}

header[data-testid="stHeader"] {
    background: transparent;
}
</style>



    """,
    unsafe_allow_html=True
)


DESCRIPTOR_GROUPS = {
    "Basic": [
        "Molecule_ID",
        "Original SMILES",
        "Canonical SMILES",
        "Molecular Formula",
        "Molecular Weight",
        "Exact Molecular Weight",
        "LogP",
        "TPSA",
        "H-Bond Donors",
        "H-Bond Acceptors"
    ],

    "Structural": [
        "Molecule_ID",
        "Original SMILES",
        "Canonical SMILES",
        "Rotatable Bonds",
        "Ring Count",
        "Aromatic Ring Count",
        "Aliphatic Ring Count",
        "Heavy Atom Count",
        "Heteroatom Count",
        "Fraction Csp3",
        "Formal Charge"
    ],

    "Topological": [
        "Molecule_ID",
        "Original SMILES",
        "Canonical SMILES",
        "Bertz CT",
        "Balaban J"
    ],

    "Element Counts": [
        "Molecule_ID",
        "Original SMILES",
        "Canonical SMILES",
        "Carbon Count",
        "Nitrogen Count",
        "Oxygen Count",
        "Sulfur Count",
        "Selenium Count",
        "Tellurium Count",
        "Phosphorus Count",
        "Halogen Count"
    ],

    "Chalcogen": [
        "Molecule_ID",
        "Original SMILES",
        "Canonical SMILES",
        "Oxygen Count",
        "Sulfur Count",
        "Selenium Count",
        "Tellurium Count",
        "Total Chalcogen Count",
        "Heavy Chalcogen Count",
        "Contains S",
        "Contains Se",
        "Contains Te"
    ]
}


def filter_descriptor_columns(dataframe, selected_groups):
    selected_columns = []

    for group in selected_groups:
        selected_columns.extend(
            DESCRIPTOR_GROUPS.get(group, [])
        )

    selected_columns = list(
        dict.fromkeys(selected_columns)
    )

    existing_columns = [
        column
        for column in selected_columns
        if column in dataframe.columns
    ]

    return dataframe[existing_columns]


def display_molecule_gallery(valid_df, max_molecules=12):
    if valid_df.empty:
        st.info("Gösterilecek geçerli molekül bulunamadı.")
        return

    molecule_count = min(
        len(valid_df),
        max_molecules
    )

    st.caption(
        f"İlk {molecule_count} geçerli molekül gösteriliyor."
    )

    cards_per_row = 3

    for start_index in range(
        0,
        molecule_count,
        cards_per_row
    ):
        columns = st.columns(cards_per_row)

        for column_index in range(cards_per_row):
            molecule_index = start_index + column_index

            if molecule_index >= molecule_count:
                break

            row = valid_df.iloc[molecule_index]

            with columns[column_index]:
                with st.container(border=True):
                    st.markdown(
                        f"### {row['Molecule_ID']}"
                    )

                    mol = Chem.MolFromSmiles(
                        row["Canonical SMILES"]
                    )

                    if mol is not None:
                        image = Draw.MolToImage(
                            mol,
                            size=(360, 280)
                        )

                        st.image(
                            image,
                            use_container_width=True
                        )

                    st.write(
                        f"**Formül:** "
                        f"{row['Molecular Formula']}"
                    )

                    st.code(
                        row["Canonical SMILES"],
                        language=None
                    )


sample_df = pd.DataFrame({
    "Molecule_ID": [
        "MOL_001",
        "MOL_002",
        "MOL_003",
        "MOL_004"
    ],
    "SMILES": [
        "c1ccccc1",
        "CCO",
        "CC(=O)O",
        "c1ccsc1"
    ]
})

sample_csv = sample_df.to_csv(
    index=False
).encode("utf-8")


with st.sidebar:
    st.title("🧬 Platform Ayarları")

    project_name = st.text_input(
        "Proje adı",
        value="molecular_project"
    )

    st.divider()

    selected_groups = st.multiselect(
        "Deskriptör grupları",
        options=list(DESCRIPTOR_GROUPS.keys()),
        default=[
            "Basic",
            "Structural",
            "Chalcogen"
        ]
    )

    st.divider()

    max_molecule_cards = st.slider(
        "Gösterilecek molekül sayısı",
        min_value=3,
        max_value=30,
        value=12,
        step=3
    )

    show_molecule_cards = st.checkbox(
        "Molekül yapılarını göster",
        value=True
    )

    st.divider()

    st.download_button(
        label="Örnek CSV indir",
        data=sample_csv,
        file_name="sample_molecules.csv",
        mime="text/csv",
        use_container_width=True
    )

    st.caption(
        "Girdi dosyası Molecule_ID ve SMILES "
        "sütunlarını içermelidir."
    )


st.markdown(
    """
    <div class="hero-card">
        <div class="main-title">
            🧬 Molecular Descriptor Platform
        </div>
        <div class="subtitle">
            Moleküler deskriptör, kalkojen analizi,
            fingerprint ve benzerlik tarama platformu
        </div>
        <b>CSV yükle → Hesapla → Analiz et → Karşılaştır → İndir</b>
    </div>
    """,
    unsafe_allow_html=True
)


uploaded_file = st.file_uploader(
    "Molekül CSV dosyanızı yükleyin",
    type=["csv"]
)


if uploaded_file is None:
    col1, col2, col3 = st.columns(3)

    with col1:
        st.markdown(
            """
            <div class="info-card">
                <h3>1. Veri Yükleme</h3>
                <p>Molecule_ID ve SMILES sütunlarını içeren CSV yükleyin.</p>
            </div>
            """,
            unsafe_allow_html=True
        )

    with col2:
        st.markdown(
            """
            <div class="info-card">
                <h3>2. Moleküler Analiz</h3>
                <p>Deskriptör, kalkojen ve fingerprint özelliklerini hesaplayın.</p>
            </div>
            """,
            unsafe_allow_html=True
        )

    with col3:
        st.markdown(
            """
            <div class="info-card">
                <h3>3. Aday Tarama</h3>
                <p>Molekülleri karşılaştırın ve en benzer adayları bulun.</p>
            </div>
            """,
            unsafe_allow_html=True
        )


if uploaded_file is not None:

    try:
        input_df = pd.read_csv(uploaded_file)

        st.subheader("Yüklenen Veri")

        st.dataframe(
            input_df.head(20),
            use_container_width=True,
            hide_index=True
        )

        if not selected_groups:
            st.warning(
                "En az bir deskriptör grubu seçmelisiniz."
            )

        calculate_button = st.button(
            "🧪 Deskriptörleri Hesapla",
            type="primary",
            disabled=not selected_groups,
            use_container_width=True
        )

        if calculate_button:

            with st.spinner(
                "Moleküller işleniyor..."
            ):
                results = (
                    run_molecular_descriptor_platform(
                        input_df=input_df,
                        project_name=project_name
                    )
                )

            st.success(
                "Hesaplama başarıyla tamamlandı."
            )

            valid_df = results["valid_df"]
            invalid_df = results["invalid_df"]
            summary_df = results["summary_df"]

            filtered_valid_df = filter_descriptor_columns(
                valid_df,
                selected_groups
            )

            total_records = int(
                summary_df.loc[0, "Total Records"]
            )

            valid_count = int(
                summary_df.loc[0, "Valid Molecules"]
            )

            invalid_count = int(
                summary_df.loc[0, "Invalid Molecules"]
            )

            success_rate = float(
                summary_df.loc[
                    0,
                    "Success Rate (%)"
                ]
            )

            metric_col1, metric_col2, metric_col3, metric_col4 = (
                st.columns(4)
            )

            metric_col1.metric(
                "Toplam Kayıt",
                total_records
            )

            metric_col2.metric(
                "Geçerli Molekül",
                valid_count
            )

            metric_col3.metric(
                "Geçersiz Molekül",
                invalid_count
            )

            metric_col4.metric(
                "Başarı Oranı",
                f"{success_rate:.2f}%"
            )

            tabs = st.tabs([
                "📋 Genel Bakış",
                "🧬 Molekül Yapıları",
                "📊 Veri Analizi",
                "⚖️ Karşılaştırma",
                "🔎 Benzer Molekül Arama",
                "⬇️ İndirmeler"
            ])

            with tabs[0]:
                st.subheader(
                    "Seçili Moleküler Deskriptörler"
                )

                st.caption(
                    "Gösterilen gruplar: "
                    + ", ".join(selected_groups)
                )

                st.dataframe(
                    filtered_valid_df,
                    use_container_width=True,
                    hide_index=True
                )

                if not invalid_df.empty:
                    st.subheader(
                        "Geçersiz SMILES Kayıtları"
                    )

                    st.dataframe(
                        invalid_df,
                        use_container_width=True,
                        hide_index=True
                    )

            with tabs[1]:
                if show_molecule_cards:
                    display_molecule_gallery(
                        valid_df,
                        max_molecule_cards
                    )
                else:
                    st.info(
                        "Molekül yapıları sol menüden kapatıldı."
                    )

            with tabs[2]:
                display_analysis_panel(valid_df)

            with tabs[3]:
                display_similarity_panel(valid_df)

            with tabs[4]:
                display_similarity_search_panel(valid_df)

            with tabs[5]:
                st.subheader("Sonuç Dosyaları")

                selected_csv = (
                    filtered_valid_df
                    .to_csv(index=False)
                    .encode("utf-8")
                )

                st.download_button(
                    label="Seçili Deskriptörleri CSV Olarak İndir",
                    data=selected_csv,
                    file_name=(
                        f"{project_name}_"
                        "selected_descriptors.csv"
                    ),
                    mime="text/csv",
                    use_container_width=True
                )

                zip_file = results["zip_file"]

                with open(zip_file, "rb") as file:
                    zip_data = file.read()

                st.download_button(
                    label="Tüm Sonuçları ZIP Olarak İndir",
                    data=zip_data,
                    file_name=os.path.basename(zip_file),
                    mime="application/zip",
                    use_container_width=True
                )

    except Exception as error:
        st.error(
            f"Dosya işlenemedi: {error}"
        )


st.markdown(
    """
    <div class="footer">
        Molecular Descriptor Platform — MVP v0.7<br>
        Developed by Esra Nur Doğan |
        Computational Chemical Engineering & AI
    </div>
    """,
    unsafe_allow_html=True
)
