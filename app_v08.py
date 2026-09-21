import os
import platform
import pandas as pd
import numpy as np
import streamlit as st

from rdkit import Chem
from rdkit import rdBase
from rdkit.Chem import Draw

from descriptor_engine import run_molecular_descriptor_platform
from database_metadata import add_database_export
from database_browser import display_database_browser, display_database_statistics
from similarity_search import display_similarity_search_panel
from scientific_panel import display_scientific_core_panel

CHALMOLDB_ICON_SVG = """
<svg viewBox="0 0 120 120" aria-hidden="true" focusable="false">
  <path d="M60 18 L92 36 L92 74 L60 102 L28 74 L28 36 Z"
        fill="none" stroke="#1F4E79" stroke-width="5"
        stroke-linecap="round" stroke-linejoin="round"/>
  <g font-family="Arial,Helvetica,sans-serif" font-weight="700"
     text-anchor="middle" dominant-baseline="central">
    <circle cx="60" cy="18" r="16" fill="#D6A400"/>
    <text x="60" y="18" font-size="20" fill="#FFFFFF">S</text>
    <circle cx="28" cy="74" r="17" fill="#0F766E"/>
    <text x="28" y="74" font-size="18" fill="#FFFFFF">Se</text>
    <circle cx="92" cy="74" r="17" fill="#6B4C8A"/>
    <text x="92" y="74" font-size="18" fill="#FFFFFF">Te</text>
    <circle cx="28" cy="36" r="8" fill="#5DA9E9"/>
    <circle cx="92" cy="36" r="8" fill="#5DA9E9"/>
    <circle cx="60" cy="102" r="8" fill="#5DA9E9"/>
  </g>
</svg>
"""

CHALMOLDB_LOGO_URL = "https://raw.githubusercontent.com/EsraDogan539/molecular-descriptor-platform/v0.8-scientific-core/assets/chalmoldb_logo.svg"
CHALMOLDB_ICON_URL = "https://raw.githubusercontent.com/EsraDogan539/molecular-descriptor-platform/v0.8-scientific-core/assets/chalmoldb_icon.svg"


st.set_page_config(
    page_title="ChalMolDB | Chalcogen Molecular Database",
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
        --accent: #1F4E79;
        --accent-dark: #163A5B;
        --accent-light: #5DA9E9;
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
        display: inline-flex;
        align-items: center;
        gap: .62rem;
        white-space: nowrap;
    }
    .brand-logo-img {
        width: 244px;
        height: auto;
        display: block;
    }
    .hero-logo-img {
        width: min(560px, 92%);
        height: auto;
        display: block;
        margin-bottom: 1.35rem;
    }
    .brand-copy { display: none; }
    .brand-copy strong {
        color: var(--accent-dark);
        font-size: 1.03rem;
        font-weight: 750;
        letter-spacing: -0.015em;
    }
    .brand-copy small {
        color: var(--muted);
        font-size: .67rem;
        font-weight: 500;
        margin-top: .19rem;
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
        max-width: 760px;
        margin: 2.0rem auto 1.15rem auto;
        text-align: center;
    }
    .hero-brand {
        display: inline-flex;
        align-items: center;
        justify-content: center;
        gap: 1rem;
        margin-bottom: .72rem;
    }
    .hero-title-wrap {
        text-align: left;
    }
    .hero-title-wrap h1 {
        color: var(--accent-dark);
        margin: 0 !important;
        font-size: 2.45rem !important;
        line-height: 1 !important;
    }
    .hero-subbrand {
        color: var(--muted);
        font-size: .92rem;
        margin-top: .36rem;
        letter-spacing: .01em;
    }
    .home-intro p {
        color: var(--muted);
        font-size: 1.02rem;
        line-height: 1.6;
    }
    .hero-shell {
        display: grid;
        grid-template-columns: 1.22fr .78fr;
        gap: 0;
        min-height: 420px;
        margin: 1.35rem 0 1.2rem 0;
        border: 1px solid #DCE6EF;
        border-radius: 12px;
        overflow: hidden;
        background: linear-gradient(135deg, #F7FBFE 0%, #FFFFFF 68%);
    }
    .hero-copy {
        padding: 3.0rem 3.15rem 2.75rem 3.15rem;
        display: flex;
        flex-direction: column;
        justify-content: center;
    }
    .hero-kicker {
        color: var(--accent);
        font-size: .78rem;
        font-weight: 700;
        letter-spacing: .11em;
        text-transform: uppercase;
        margin-bottom: .8rem;
    }
    .hero-heading {
        color: var(--accent-dark);
        font-size: 3.15rem;
        font-weight: 760;
        letter-spacing: -0.035em;
        line-height: 1;
        margin: 0;
    }
    .hero-name {
        color: var(--muted);
        font-size: 1.03rem;
        margin-top: .55rem;
    }
    .hero-copy p {
        color: #44546A;
        font-size: 1.10rem;
        line-height: 1.65;
        max-width: 650px;
        margin: 1.35rem 0 0 0;
    }
    .hero-visual {
        position: relative;
        min-height: 420px;
        background:
          radial-gradient(circle at 76% 19%, rgba(93,169,233,.28), transparent 26%),
          radial-gradient(circle at 27% 82%, rgba(31,78,121,.13), transparent 24%),
          linear-gradient(145deg, #EAF4FA 0%, #F9FCFE 100%);
        overflow: hidden;
    }
    .hero-visual::before {
        content: "";
        position: absolute;
        inset: 0;
        background:
          radial-gradient(circle at 78% 23%, rgba(255,255,255,.96) 0 20px, rgba(93,169,233,.20) 21px 31px, transparent 32px),
          radial-gradient(circle at 67% 54%, rgba(255,255,255,.94) 0 15px, rgba(31,78,121,.20) 16px 25px, transparent 26px),
          radial-gradient(circle at 43% 67%, rgba(255,255,255,.95) 0 19px, rgba(93,169,233,.18) 20px 30px, transparent 31px);
        opacity: .95;
    }
    .molecule-stage {
        position: absolute;
        inset: 0;
        display: flex;
        align-items: center;
        justify-content: center;
    }
    .molecule-network {
        position: relative;
        width: 340px;
        height: 290px;
        transform: rotate(-4deg);
    }
    .atom {
        position: absolute;
        width: 42px;
        height: 42px;
        border-radius: 50%;
        border: 2px solid rgba(22,58,91,.22);
        background: rgba(255,255,255,.92);
        box-shadow: 0 8px 22px rgba(22,58,91,.10);
        z-index: 2;
    }
    .atom.core {
        width: 54px;
        height: 54px;
        background: #FFFFFF;
        border-color: rgba(31,78,121,.35);
    }
    .atom.s { background: #F3CF58; border-color: #D6A400; }
    .atom.se { background: #3BA79C; border-color: #0F766E; }
    .atom.te { background: #80639A; border-color: #6B4C8A; }
    .atom.a1 { left: 34px; top: 92px; }
    .atom.a2 { left: 106px; top: 36px; }
    .atom.a3 { left: 184px; top: 70px; }
    .atom.a4 { left: 210px; top: 150px; }
    .atom.a5 { left: 126px; top: 174px; }
    .atom.a6 { left: 58px; top: 154px; }
    .atom.a7 { left: 124px; top: 104px; }
    .mol-bond {
        position: absolute;
        height: 4px;
        background: rgba(31,78,121,.34);
        transform-origin: left center;
        border-radius: 8px;
        z-index: 1;
    }
    .mol-bond.b1 { width: 92px; left: 62px; top: 106px; transform: rotate(-38deg); }
    .mol-bond.b2 { width: 88px; left: 132px; top: 58px; transform: rotate(24deg); }
    .mol-bond.b3 { width: 88px; left: 201px; top: 95px; transform: rotate(72deg); }
    .mol-bond.b4 { width: 92px; left: 150px; top: 190px; transform: rotate(-17deg); }
    .mol-bond.b5 { width: 82px; left: 79px; top: 176px; transform: rotate(15deg); }
    .mol-bond.b6 { width: 74px; left: 72px; top: 132px; transform: rotate(55deg); }
    .mol-bond.b7 { width: 72px; left: 148px; top: 128px; transform: rotate(-25deg); }
    .element-tag {
        position: absolute;
        z-index: 3;
        color: white;
        font-size: .70rem;
        font-weight: 800;
        line-height: 42px;
        text-align: center;
        width: 42px;
        height: 42px;
        pointer-events: none;
    }
    .element-tag.t1 { left: 106px; top: 36px; }
    .element-tag.t2 { left: 58px; top: 154px; }
    .element-tag.t3 { left: 210px; top: 150px; }
    .feature-grid {
        display: grid;
        grid-template-columns: repeat(4, 1fr);
        gap: 1rem;
        margin: 1.15rem 0 1.35rem 0;
    }
    .feature-card {
        border: 1px solid #DFE7EF;
        border-radius: 10px;
        padding: 1.15rem 1.15rem 1.05rem 1.15rem;
        background: #FFFFFF;
        min-height: 135px;
    }
    .feature-icon {
        width: 38px;
        height: 38px;
        border-radius: 50%;
        background: #F4F9FC;
        color: var(--accent);
        display: flex;
        align-items: center;
        justify-content: center;
        font-weight: 800;
        font-size: .74rem;
        margin-bottom: .85rem;
        border: 1px solid #CFE0EC;
        box-shadow: inset 0 0 0 3px #FFFFFF;
    }
    .feature-card strong {
        display: block;
        color: var(--accent-dark);
        font-size: .94rem;
        margin-bottom: .35rem;
    }
    .feature-card span {
        display: block;
        color: var(--muted);
        font-size: .81rem;
        line-height: 1.5;
    }
    .section-label {
        color: var(--accent-dark);
        font-size: .82rem;
        font-weight: 700;
        letter-spacing: .08em;
        text-transform: uppercase;
        margin: .35rem 0 .65rem 0;
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
        margin-top: 2.2rem;
        padding: 1.45rem 1.55rem;
        background: var(--accent-dark);
        color: #D8E4EE;
        font-size: .80rem;
        display: flex;
        justify-content: space-between;
        gap: 1rem;
        flex-wrap: wrap;
        border-radius: 8px;
    }
    .site-footer strong { color: #FFFFFF; }
    .site-footer a {
        color: #E7EFF6 !important;
        text-decoration: none !important;
        margin-left: .85rem;
    }
    .site-footer a:hover { text-decoration: underline !important; }
    @media (max-width: 760px) {
        .academic-header { align-items: flex-start; flex-direction: column; }
        .academic-nav { gap: .9rem; flex-wrap: wrap; }
        .metric-row { flex-wrap: wrap; }
        .metric-item { min-width: 50%; margin-bottom: .8rem; }
        .hero-shell { grid-template-columns: 1fr; }
        .hero-visual { min-height: 220px; }
        .hero-copy { padding: 2rem 1.5rem; }
        .hero-heading { font-size: 2.45rem; }
        .feature-grid { grid-template-columns: 1fr 1fr; }
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
        padding: .35rem;
        background: #FFFFFF;
    }
    div[data-testid="stFileUploader"] section {
        background: #FFFFFF !important;
        border: 1px dashed #D9DEE7 !important;
        border-radius: 4px !important;
    }
    .analyze-note {
        color: var(--muted);
        font-size: .88rem;
        margin-top: .8rem;
        padding-top: .75rem;
        border-top: 1px solid var(--line);
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
        background: var(--accent-dark) !important;
        border-color: var(--accent-dark) !important;
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


DISPLAY_LABELS = {
    "Molecule_ID": "Molecule ID",
    "Original SMILES": "Original SMILES",
    "Canonical SMILES": "Canonical SMILES",
    "InChIKey": "InChIKey",
    "Molecular Formula": "Molecular Formula",
    "Molecular Weight": "Molecular Weight",
    "Exact Molecular Weight": "Exact Molecular Weight",
    "HOMO_eV": "HOMO (eV)",
    "LUMO_eV": "LUMO (eV)",
    "Eg_eV": "Eg (eV)",
    "Property_Source": "Property Source",
    "DOI_or_Reference": "DOI / Reference",
    "Chalcogen Type": "Chalcogen Type",
    "Target Chalcogen Count": "Target Chalcogen Count",
    "Target Chalcogen Fraction": "Target Chalcogen Fraction",
    "Aromatic Chalcogen Count": "Aromatic Chalcogen Count",
    "NonAromatic Chalcogen Count": "Non-aromatic Chalcogen Count",
    "Mixed Chalcogen Flag": "Mixed Chalcogen",
    "Ring Incorporated Chalcogen Count": "Ring-incorporated Chalcogen Count",
    "Chalcogen-C Bond Count": "Chalcogen–C Bond Count",
    "Chalcogen-Heteroatom Bond Count": "Chalcogen–Heteroatom Bond Count",
    "Aromatic Neighbor Count": "Aromatic Neighbor Count",
    "Conjugated Bond Count": "Conjugated Bond Count",
    "Aromatic Bond Fraction": "Aromatic Bond Fraction",
    "Conjugated Atom Fraction": "Conjugated Atom Fraction",
    "Heteroaromatic Ring Count": "Heteroaromatic Ring Count",
    "Duplicate Flag": "Repeated Structure",
}


def display_table(dataframe):
    table = dataframe.copy()
    table = table.rename(columns={c: DISPLAY_LABELS.get(c, c.replace("_", " ")) for c in table.columns})
    return table


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
if page not in {"home", "database", "statistics", "analyze", "documentation", "about"}:
    page = "home"

nav_items = [
    ("Database", "database"),
    ("Statistics", "statistics"),
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
      <a class="academic-brand" href="?page=home">
        <img class="brand-logo-img" src="https://raw.githubusercontent.com/EsraDogan539/molecular-descriptor-platform/v0.8-scientific-core/assets/chalmoldb_logo.svg" alt="ChalMolDB — Chalcogen Molecular Database">
      </a>
      <nav class="academic-nav">{nav_html}</nav>
    </div>
    """,
    unsafe_allow_html=True,
)

if page == "home":
    st.markdown(
        """
        <div class="hero-shell">
          <div class="hero-copy">
            <div class="hero-kicker">Curated molecular data · S / Se / Te</div>
            <img class="hero-logo-img" src="https://raw.githubusercontent.com/EsraDogan539/molecular-descriptor-platform/v0.8-scientific-core/assets/chalmoldb_logo.svg" alt="ChalMolDB — Chalcogen Molecular Database">
            <p>
              A curated molecular database for chalcogen-focused electronic property studies,
              combining standardized identity, provenance and interpretable structural annotations.
            </p>
          </div>
          <div class="hero-visual" aria-hidden="true">
            <div class="molecule-stage">
              <div class="molecule-network">
                <div class="mol-bond b1"></div>
                <div class="mol-bond b2"></div>
                <div class="mol-bond b3"></div>
                <div class="mol-bond b4"></div>
                <div class="mol-bond b5"></div>
                <div class="mol-bond b6"></div>
                <div class="mol-bond b7"></div>
                <div class="atom a1"></div>
                <div class="atom s a2"></div>
                <div class="atom core a3"></div>
                <div class="atom te a4"></div>
                <div class="atom core a5"></div>
                <div class="atom se a6"></div>
                <div class="atom core a7"></div>
                <div class="element-tag t1">S</div>
                <div class="element-tag t2">Se</div>
                <div class="element-tag t3">Te</div>
              </div>
            </div>
          </div>
        </div>
        """,
        unsafe_allow_html=True,
    )

    left, c1, c2, c3, right = st.columns([.55, 1, 1, 1, .55])
    with c1:
        st.link_button("Explore database", "?page=database", type="primary", use_container_width=True)
    with c2:
        st.link_button("View statistics", "?page=statistics", use_container_width=True)
    with c3:
        st.link_button("Analyze molecules", "?page=analyze", use_container_width=True)

    st.markdown(
        """
        <div class="section-label">What ChalMolDB provides</div>
        <div class="feature-grid">
          <div class="feature-card">
            <div class="feature-icon">▤</div>
            <strong>Curated Molecular Records</strong>
            <span>Standardized molecular identity, electronic properties and record-level provenance.</span>
          </div>
          <div class="feature-card">
            <div class="feature-icon">⬡</div>
            <strong>Chalcogen-Aware Descriptors</strong>
            <span>Interpretable annotations describing sulfur, selenium and tellurium environments.</span>
          </div>
          <div class="feature-card">
            <div class="feature-icon">◎</div>
            <strong>Similarity Search</strong>
            <span>Structure-based exploratory comparison using molecular fingerprints.</span>
          </div>
          <div class="feature-card">
            <div class="feature-icon">⇩</div>
            <strong>Export & Reuse</strong>
            <span>Downloadable processed records, descriptor tables and reproducible analysis outputs.</span>
          </div>
        </div>
        """,
        unsafe_allow_html=True,
    )

    st.markdown(
        """
        <div class="section-label">Database at a glance</div>
        <div class="metric-row">
          <div class="metric-item"><span class="metric-number">3,360</span><span class="metric-label">Records</span></div>
          <div class="metric-item"><span class="metric-number">2,983</span><span class="metric-label">Unique structures</span></div>
          <div class="metric-item"><span class="metric-number">3,145</span><span class="metric-label">S/Se/Te records</span></div>
          <div class="metric-item"><span class="metric-number">272</span><span class="metric-label">External records</span></div>
        </div>
        <div class="quiet-note">
          Database v1 &nbsp;&middot;&nbsp; Structure standardization with RDKit
          &nbsp;&middot;&nbsp; Record-level provenance retained
        </div>
        <div class="site-footer">
          <span><strong>ChalMolDB</strong> · Chalcogen Molecular Database · Database v1</span>
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
        "ChalMolDB · Database v1 · Scientific Core v0.8 · "
        "Curated records are read-only in the public browser."
    )
    st.stop()


if page == "statistics":
    display_database_statistics()
    st.divider()
    st.caption("ChalMolDB · Database v1 · Descriptive statistics")
    st.stop()


if page == "documentation":
    st.header("Documentation")
    st.caption("Database scope, molecular identity, curation rules and user-analysis workflow.")

    st.markdown(
        """
        <div class="metric-row" style="margin:1.1rem 0 1.8rem 0; justify-content:flex-start;">
          <div class="metric-item" style="padding-left:0;"><span class="metric-number">3,360</span><span class="metric-label">Records</span></div>
          <div class="metric-item"><span class="metric-number">3,145</span><span class="metric-label">Core S/Se/Te</span></div>
          <div class="metric-item"><span class="metric-number">2,983</span><span class="metric-label">Unique structures</span></div>
          <div class="metric-item"><span class="metric-number">3,353</span><span class="metric-label">Eg values</span></div>
        </div>
        """,
        unsafe_allow_html=True,
    )

    st.markdown("### Database scope")
    st.write(
        "Database v1 contains curated development and external records. Exact standardized "
        "structures are displayed only where the source supports structure-level identity."
    )

    st.markdown("### Molecular identity")
    st.write(
        "Structure-complete records are standardized with RDKit and represented by canonical "
        "SMILES, InChI and InChIKey. Repeated standardized structures are retained to preserve "
        "record-level provenance."
    )

    st.markdown("### Curation principles")
    st.markdown(
        "- Missing structures or properties are left explicit; unavailable values are not imputed.\n"
        "- S/Se/Te-focused annotations are reported separately from general molecular descriptors.\n"
        "- Record identity and standardized structure identity are treated as distinct concepts."
    )

    st.markdown("### Analyze your dataset")
    st.write(
        "Upload a CSV containing Molecule_ID and SMILES. The platform validates structures, "
        "calculates general and S/Se/Te-aware descriptors, and keeps user-supplied data separate "
        "from the curated publication database."
    )

    st.markdown("### Software environment")
    with st.expander("Runtime versions", expanded=False):
        st.code(
            "\n".join([
                f"Python {platform.python_version()}",
                f"Streamlit {st.__version__}",
                f"RDKit {rdBase.rdkitVersion}",
                f"pandas {pd.__version__}",
                f"NumPy {np.__version__}",
            ]),
            language=None,
        )

    st.markdown("### Citation and data release")
    st.write(
        "The recommended citation and permanent dataset DOI will be added with the archived "
        "publication release."
    )
    st.stop()


if page == "about":
    st.header("About")
    st.caption("ChalMolDB — a curated research resource for chalcogen-focused molecular data.")

    st.write(
        "ChalMolDB (Chalcogen Molecular Database) brings together standardized molecular identity, "
        "electronic-property data and interpretable S/Se/Te structural annotations in a "
        "single searchable resource."
    )

    st.write(
        "The platform supports transparent record inspection, structure-aware data curation "
        "and reproducible downstream cheminformatics workflows. Predictive modelling is treated "
        "as a secondary use case rather than the primary purpose of the database."
    )

    st.markdown("### Current release")
    st.markdown(
        "**Database v1** · **Scientific Core v0.8**  \\n"
        "3,360 records · 2,983 unique standardized development structures"
    )

    st.markdown("### Resources")
    st.markdown(
        "[Documentation](?page=documentation) · "
        "[GitHub repository](https://github.com/EsraDogan539/molecular-descriptor-platform)"
    )
    st.stop()


st.header("Analyze Your Dataset")
st.caption(
    "Upload molecular structures to validate records and calculate general and "
    "S/Se/Te-aware descriptors."
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
    st.caption("Required columns: Molecule_ID, SMILES")
    st.caption("Additional property or provenance columns are preserved.")
    st.download_button(
        "Download example CSV",
        sample_csv,
        "sample_chalcogen_database.csv",
        "text/csv",
        use_container_width=True,
    )

if uploaded_file is None:
    st.markdown(
        '<div class="analyze-note">Analysis options appear after a CSV file is loaded. '
        'User uploads are processed separately from the curated database.</div>',
        unsafe_allow_html=True,
    )
    st.stop()

try:
    input_df = pd.read_csv(uploaded_file)
    st.markdown("### Input preview")
    st.dataframe(display_table(input_df.head(12)), use_container_width=True, hide_index=True)

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
            st.dataframe(display_table(database_df), use_container_width=True, hide_index=True)
            if not invalid_df.empty:
                st.markdown("#### Invalid structures")
                st.dataframe(display_table(invalid_df), use_container_width=True, hide_index=True)

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
                    display_table(valid_df[chalcogen_columns]),
                    use_container_width=True,
                    hide_index=True,
                )
            with descriptor_tabs[1]:
                if general_columns:
                    st.dataframe(
                        display_table(filtered_valid_df[general_columns]),
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
                st.caption("Standardized records with retained source fields and quality flags.")
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
                st.caption("Selected descriptor groups for the current analysis.")
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
                st.caption("Processed records, descriptors, fingerprints and run summary.")
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
    "ChalMolDB · Database v1 · Scientific Core v0.8"
)
