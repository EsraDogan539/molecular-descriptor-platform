import os
import platform
import pandas as pd
import numpy as np
import streamlit as st

from rdkit import Chem
from rdkit import rdBase
from rdkit.Chem import Draw

from descriptor_engine import run_molecular_descriptor_platform, validate_input_dataframe
from database_metadata import add_database_export
from database_browser import display_database_browser, display_database_statistics
from similarity_search import display_similarity_search_panel
from scientific_panel import display_scientific_core_panel
from descriptor_dictionary import (
    DESCRIPTOR_DICTIONARY_VERSION,
    descriptor_dictionary_dataframe,
)
from release_metadata import (
    ASSET_BASE_URL,
    DATABASE_VERSION,
    RELEASE_LABEL,
    SCIENTIFIC_CORE_DISPLAY_VERSION,
)

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

CHALMOLDB_LOGO_URL = f"{ASSET_BASE_URL}/chalmoldb_logo.svg"
CHALMOLDB_ICON_URL = f"{ASSET_BASE_URL}/chalmoldb_icon.svg"


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
        gap: .72rem;
        white-space: nowrap;
    }
    .brand-mark-img {
        width: 66px;
        height: 58px;
        object-fit: contain;
        display: block;
        flex: 0 0 auto;
    }
    .brand-copy {
        display: flex;
        flex-direction: column;
        justify-content: center;
        line-height: 1;
    }
    .brand-copy strong {
        color: var(--accent-dark);
        font-size: 1.58rem;
        font-weight: 800;
        letter-spacing: -0.035em;
    }
    .brand-copy small {
        color: var(--muted);
        font-size: .69rem;
        font-weight: 500;
        margin-top: .30rem;
        letter-spacing: .01em;
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
        display: flex;
        align-items: center;
        gap: .88rem;
        margin: .1rem 0 1.25rem 0;
    }
    .hero-mark-img {
        width: 116px;
        height: 98px;
        object-fit: contain;
        display: block;
        flex: 0 0 auto;
    }
    .hero-title-wrap {
        text-align: left;
    }
    .hero-title-wrap h1 {
        color: var(--accent-dark);
        margin: 0 !important;
        font-size: 3.35rem !important;
        font-weight: 800 !important;
        letter-spacing: -0.045em;
        line-height: .94 !important;
    }
    .hero-subbrand {
        color: #657286;
        font-size: 1.11rem;
        font-weight: 500;
        margin-top: .30rem;
        letter-spacing: .005em;
    }
    .home-intro p {
        color: var(--muted);
        font-size: 1.02rem;
        line-height: 1.6;
    }
    .hero-shell {
        display: grid;
        grid-template-columns: 1.16fr .84fr;
        gap: 0;
        min-height: 388px;
        margin: 1.35rem 0 1.2rem 0;
        border: 1px solid #DCE6EF;
        border-radius: 12px;
        overflow: hidden;
        background: linear-gradient(135deg, #F7FBFE 0%, #FFFFFF 68%);
    }
    .hero-copy {
        padding: 2.95rem 3.15rem 2.75rem 3.15rem;
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
        margin-bottom: 1.05rem;
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
        font-size: 1.08rem;
        line-height: 1.67;
        max-width: 610px;
        margin: 1.15rem 0 0 0;
    }
    .hero-visual {
        position: relative;
        min-height: 388px;
        background:
          radial-gradient(circle at 72% 22%, rgba(255,255,255,.78), transparent 29%),
          radial-gradient(circle at 34% 78%, rgba(111,160,194,.12), transparent 28%),
          linear-gradient(145deg, #E6F1F8 0%, #F6FAFD 100%);
        overflow: hidden;
        display: flex;
        align-items: center;
        justify-content: center;
    }
    .hero-visual::after {
        content: "";
        position: absolute;
        inset: 0;
        background: linear-gradient(90deg, rgba(255,255,255,.18), transparent 24%);
        pointer-events: none;
    }
    .hero-molecule-svg {
        width: 119%;
        max-width: none;
        height: auto;
        display: block;
        transform: translate(7%, 1%);
        filter: drop-shadow(0 20px 26px rgba(22,58,91,.13));
    }

    .page-intro {
        margin: 1.45rem 0 1.65rem 0;
        padding: 1.55rem 1.7rem 1.5rem 1.7rem;
        border: 1px solid #DFE8EF;
        border-radius: 10px;
        background: linear-gradient(135deg, #F8FBFD 0%, #FFFFFF 72%);
    }
    .page-eyebrow {
        color: var(--accent);
        font-size: .74rem;
        font-weight: 750;
        letter-spacing: .11em;
        text-transform: uppercase;
        margin-bottom: .48rem;
    }
    .page-title {
        color: var(--accent-dark);
        font-size: 2rem;
        line-height: 1.08;
        font-weight: 780;
        letter-spacing: -0.025em;
        margin: 0;
    }
    .page-description {
        color: var(--muted);
        font-size: .96rem;
        line-height: 1.58;
        max-width: 790px;
        margin-top: .6rem;
    }
    .info-grid {
        display: grid;
        grid-template-columns: repeat(2, minmax(0, 1fr));
        gap: 1rem;
        margin: 1.1rem 0 1.4rem 0;
    }
    .info-card {
        border: 1px solid #DFE7EF;
        border-radius: 9px;
        padding: 1.15rem 1.2rem;
        background: #FFFFFF;
    }
    .info-card strong {
        color: var(--accent-dark);
        display: block;
        font-size: .93rem;
        margin-bottom: .35rem;
    }
    .info-card span {
        color: var(--muted);
        font-size: .84rem;
        line-height: 1.5;
        display: block;
    }
    .release-card {
        border: 1px solid #DFE7EF;
        border-radius: 9px;
        background: #FBFCFD;
        padding: 1rem 1.1rem;
        max-width: 620px;
        margin-top: .55rem;
    }
    .release-card strong {
        display: block;
        color: var(--accent-dark);
        font-size: .96rem;
        margin-bottom: .28rem;
    }
    .release-card span {
        display: block;
        color: var(--muted);
        font-size: .84rem;
        line-height: 1.45;
    }
    .resource-links {
        margin-top: .45rem;
        font-size: .88rem;
    }
    .resource-links a {
        color: var(--accent) !important;
        text-decoration: none !important;
        margin-right: 1rem;
        border-bottom: 1px solid rgba(31,78,121,.25);
    }
    .resource-links a:hover {
        border-bottom-color: var(--accent);
    }
    .section-rule-title {
        color: var(--accent-dark);
        font-size: .78rem;
        font-weight: 750;
        letter-spacing: .09em;
        text-transform: uppercase;
        margin: 1.55rem 0 .65rem 0;
    }

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
        width: 42px;
        height: 42px;
        display: flex;
        align-items: center;
        justify-content: center;
        margin-bottom: .9rem;
        color: var(--accent);
    }
    .feature-icon svg {
        width: 32px;
        height: 32px;
        fill: none;
        stroke: currentColor;
        stroke-width: 1.8;
        stroke-linecap: round;
        stroke-linejoin: round;
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
    .about-intro {
        color: var(--ink);
        font-size: .96rem;
        line-height: 1.65;
        max-width: 980px;
        margin: .2rem 0 1rem 0;
    }
    .about-grid {
        display: grid;
        grid-template-columns: repeat(2, minmax(0, 1fr));
        gap: 1rem;
        margin: .5rem 0 1.35rem 0;
    }
    .about-card {
        border: 1px solid #DFE7EF;
        border-radius: 10px;
        padding: 1rem 1.1rem;
        background: #FFFFFF;
        min-height: 150px;
    }
    .about-card strong {
        display: block;
        color: var(--accent-dark);
        font-size: .96rem;
        margin-bottom: .38rem;
    }
    .about-card p {
        margin: 0;
        color: var(--muted);
        font-size: .82rem;
        line-height: 1.55;
    }
    .citation-box {
        border-left: 3px solid var(--accent);
        background: #F7FAFC;
        padding: .9rem 1rem;
        margin: .25rem 0 1.4rem 0;
        border-radius: 0 8px 8px 0;
        color: var(--muted);
        font-size: .82rem;
        line-height: 1.55;
    }
    .citation-box strong {
        color: var(--accent-dark);
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
        margin-top: 2.4rem;
        padding: 2rem 2.15rem;
        background: #0D263B;
        color: #BFD0DE;
        font-size: .81rem;
        display: flex;
        justify-content: space-between;
        align-items: center;
        gap: 1.4rem;
        flex-wrap: wrap;
        border-radius: 10px;
        box-shadow: 0 14px 32px rgba(13,38,59,.10);
        border-top: 1px solid rgba(255,255,255,.08);
    }
    .site-footer strong {
        color: #FFFFFF;
        font-size: .94rem;
        letter-spacing: .01em;
    }
    .site-footer a {
        color: #DDE8F0 !important;
        text-decoration: none !important;
        margin-left: 1rem;
    }
    .site-footer a:hover { text-decoration: underline !important; }
    @media (max-width: 760px) {
        .academic-header { align-items: flex-start; flex-direction: column; }
        .academic-nav { gap: .9rem; flex-wrap: wrap; }
        .metric-row { flex-wrap: wrap; }
        .metric-item { min-width: 50%; margin-bottom: .8rem; }
        .hero-shell { grid-template-columns: 1fr; }
        .hero-visual { min-height: 210px; }
        .hero-copy { padding: 2rem 1.5rem; }
        .hero-brand { gap: .85rem; }
        .hero-mark-img { width: 88px; height: 76px; }
        .hero-title-wrap h1 { font-size: 2.45rem !important; }
        .brand-mark-img { width: 58px; height: 50px; }
        .brand-copy strong { font-size: 1.38rem; }
        .hero-molecule-svg { width: 108%; transform: translateX(4%); }
        .feature-grid { grid-template-columns: 1fr 1fr; }
        .info-grid { grid-template-columns: 1fr; }
        .page-title { font-size: 1.7rem; }
    }
    div[data-testid="stMetric"] {
        background: #FFFFFF;
        border: 1px solid #E1E8EF;
        border-radius: 8px;
        padding: .78rem .9rem;
        min-height: 84px;
    }
    div[data-testid="stMetric"] [data-testid="stMetricValue"] {
        color: var(--accent-dark);
        font-size: 1.28rem;
        font-weight: 700;
    }
    div[data-testid="stMetric"] [data-testid="stMetricLabel"] {
        color: var(--muted);
        font-size: .82rem;
    }
    div[data-testid="stDataFrame"] {
        border: 1px solid #DFE7EF;
        border-radius: 8px;
        overflow: hidden;
        box-shadow: 0 4px 14px rgba(22,58,91,.035);
    }
    .analyze-toolbar {
        display: grid;
        grid-template-columns: 1.65fr .85fr;
        gap: 1rem;
        align-items: stretch;
        margin-bottom: 1.2rem;
    }
    .upload-side-card {
        border: 1px solid #DFE7EF;
        border-radius: 9px;
        background: #FBFCFD;
        padding: 1rem 1.05rem;
        min-height: 118px;
        display: flex;
        flex-direction: column;
        justify-content: center;
    }
    .upload-side-card strong {
        display: block;
        color: var(--accent-dark);
        font-size: .88rem;
        margin-bottom: .45rem;
    }
    .upload-side-card span {
        display: block;
        color: var(--muted);
        font-size: .79rem;
        line-height: 1.48;
        margin-bottom: .22rem;
    }
    .analysis-summary {
        display: grid;
        grid-template-columns: repeat(4, minmax(0,1fr));
        gap: .75rem;
        margin: .95rem 0 1.1rem 0;
    }
    .analysis-summary-item {
        border: 1px solid #E0E7EE;
        border-radius: 8px;
        padding: .72rem .82rem;
        background: #FFFFFF;
    }
    .analysis-summary-item span {
        display: block;
        color: #7A8794;
        font-size: .70rem;
        text-transform: uppercase;
        letter-spacing: .06em;
    }
    .analysis-summary-item strong {
        display: block;
        color: var(--accent-dark);
        font-size: 1.05rem;
        margin-top: .14rem;
    }
    .workspace-title {
        color: var(--accent-dark);
        font-size: 1.06rem;
        font-weight: 750;
        margin: 1rem 0 .2rem 0;
    }
    .workspace-caption {
        color: var(--muted);
        font-size: .82rem;
        margin-bottom: .75rem;
    }
    .export-row {
        border: 1px solid #E0E7EE;
        border-radius: 8px;
        padding: .85rem .95rem;
        margin-bottom: .7rem;
        background: #FFFFFF;
    }
    .export-row strong {
        color: var(--accent-dark);
    }
    @media (max-width: 760px) {
        .analyze-toolbar,
        .analysis-summary {
            grid-template-columns: 1fr;
        }
    }

    div[data-testid="stFileUploader"] {
        border: 1px solid var(--line);
        border-radius: 8px;
        padding: .42rem;
        background: #FFFFFF;
    }
    div[data-testid="stFileUploader"] section {
        background: #FFFFFF !important;
        border: 1px dashed #D9DEE7 !important;
        border-radius: 7px !important;
    }
    .analyze-note {
        color: var(--muted);
        font-size: .88rem;
        margin-top: .8rem;
        padding-top: .75rem;
        border-top: 1px solid var(--line);
    }
    div[data-testid="stExpander"] {
        border: 1px solid #DFE7EF;
        border-radius: 8px;
        background: #FFFFFF;
    }
    div[data-testid="stExpander"] details[open] > div {
        padding-top: .45rem !important;
        padding-bottom: .65rem !important;
    }
    div[data-testid="stExpander"] label {
        margin-bottom: .18rem !important;
    }
    div[data-baseweb="select"] > div {
        border-radius: 7px !important;
        border-color: #D9E2EA !important;
    }
    div[data-testid="stTextInput"] input,
    div[data-testid="stNumberInput"] input {
        border-radius: 7px !important;
    }
    div[data-testid="stTabs"] [data-baseweb="tab-list"] {
        gap: 1.1rem;
        border-bottom: 1px solid #E3E9EF;
    }
    div[data-testid="stTabs"] button[data-baseweb="tab"] {
        padding-left: .1rem;
        padding-right: .1rem;
    }
    .stButton > button, .stDownloadButton > button {
        border-radius: 7px !important;
        font-weight: 600 !important;
        box-shadow: none !important;
    }
    .analyze-run-anchor + div[data-testid="stButton"] {
        width: 210px;
    }
    .analyze-run-anchor + div[data-testid="stButton"] > button {
        width: 210px !important;
        min-height: 44px;
        font-weight: 700 !important;
        letter-spacing: .01em;
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
    f'<a class="{"active" if page == key else ""}" href="?page={key}" target="_self">{label}</a>'
    for label, key in nav_items
)

st.markdown(
    f"""
    <div class="academic-header">
      <a class="academic-brand" href="?page=home" target="_self">
        <img class="brand-mark-img" src="{CHALMOLDB_ICON_URL}" alt="">
        <span class="brand-copy">
          <strong>ChalMolDB</strong>
          <small>Chalcogen Molecular Database</small>
        </span>
      </a>
      <nav class="academic-nav">{nav_html}</nav>
    </div>
    """,
    unsafe_allow_html=True,
)

if page == "home":
    st.markdown(
        f"""
        <div class="hero-shell">
          <div class="hero-copy">
            <div class="hero-kicker">Curated molecular data · S / Se / Te</div>
            <div class="hero-brand">
              <img class="hero-mark-img" src="{CHALMOLDB_ICON_URL}" alt="">
              <div class="hero-title-wrap">
                <h1>ChalMolDB</h1>
                <div class="hero-subbrand">Chalcogen Molecular Database</div>
              </div>
            </div>
            <p>
              A curated molecular database for chalcogen-focused electronic property studies,
              combining standardized identity, provenance and interpretable structural annotations.
            </p>
          </div>
          <div class="hero-visual" aria-hidden="true">
            <svg class="hero-molecule-svg" viewBox="0 0 760 540" xmlns="http://www.w3.org/2000/svg">
              <defs>
                <radialGradient id="heroBg" cx="68%" cy="24%" r="74%"><stop stop-color="#ffffff"/><stop offset=".52" stop-color="#eef6fb"/><stop offset="1" stop-color="#dfeaf2"/></radialGradient>
                <radialGradient id="heroC" cx="32%" cy="25%" r="72%"><stop stop-color="#ffffff"/><stop offset=".48" stop-color="#dbe1e7"/><stop offset="1" stop-color="#8b96a1"/></radialGradient>
                <radialGradient id="heroS" cx="30%" cy="24%" r="72%"><stop stop-color="#fff5b6"/><stop offset=".42" stop-color="#e7bd39"/><stop offset="1" stop-color="#ad7e00"/></radialGradient>
                <radialGradient id="heroSe" cx="30%" cy="24%" r="72%"><stop stop-color="#c6fff6"/><stop offset=".42" stop-color="#39aa9e"/><stop offset="1" stop-color="#0d6b63"/></radialGradient>
                <radialGradient id="heroTe" cx="30%" cy="24%" r="72%"><stop stop-color="#eee6f7"/><stop offset=".42" stop-color="#846b9f"/><stop offset="1" stop-color="#55406d"/></radialGradient>
                <linearGradient id="heroBond" x1="0" x2="1"><stop stop-color="#99a4af"/><stop offset=".5" stop-color="#cbd2d9"/><stop offset="1" stop-color="#7f8a95"/></linearGradient>
                <filter id="heroShadow" x="-40%" y="-40%" width="180%" height="180%"><feDropShadow dx="0" dy="9" stdDeviation="10" flood-color="#173b5b" flood-opacity=".20"/></filter>
              </defs>
              <g opacity=".17" stroke="#7D95A7" stroke-width="2" fill="none">
                <path d="M560 46 650 98 650 202 560 254 470 202 470 98Z"/>
                <path d="M642 360 686 386 686 437 642 463 598 437 598 386Z"/>
              </g>
              <g filter="url(#heroShadow)" stroke-linecap="round">
                <g stroke="url(#heroBond)" stroke-width="14">
                  <path d="M182 257 282 176"/><path d="M282 176 398 207"/><path d="M398 207 459 310"/>
                  <path d="M459 310 351 367"/><path d="M351 367 233 348"/><path d="M233 348 182 257"/>
                  <path d="M282 176 310 274"/><path d="M398 207 310 274"/><path d="M310 274 351 367"/>
                  <path d="M459 310 568 267"/><path d="M568 267 632 341"/>
                </g>
                <g stroke="#f9fbfc" stroke-width="3">
                  <circle cx="182" cy="257" r="31" fill="url(#heroC)"/><circle cx="282" cy="176" r="39" fill="url(#heroS)"/>
                  <circle cx="398" cy="207" r="31" fill="url(#heroC)"/><circle cx="459" cy="310" r="43" fill="url(#heroTe)"/>
                  <circle cx="351" cy="367" r="31" fill="url(#heroC)"/><circle cx="233" cy="348" r="43" fill="url(#heroSe)"/>
                  <circle cx="310" cy="274" r="33" fill="url(#heroC)"/><circle cx="568" cy="267" r="27" fill="url(#heroC)"/><circle cx="632" cy="341" r="24" fill="url(#heroC)"/>
                </g>
                <g font-family="Arial,Helvetica,sans-serif" font-weight="800" text-anchor="middle" dominant-baseline="central" fill="#fff">
                  <text x="282" y="176" font-size="23">S</text><text x="233" y="348" font-size="22">Se</text><text x="459" y="310" font-size="22">Te</text>
                </g>
              </g>
            </svg>
          </div>
        </div>
        """,
        unsafe_allow_html=True,
    )

    def _navigate_home(target_page: str) -> None:
        st.query_params["page"] = target_page
        st.rerun()

    left, c1, c2, c3, right = st.columns([.55, 1, 1, 1, .55])
    with c1:
        if st.button("Explore database", type="primary", use_container_width=True, key="home_explore_database"):
            _navigate_home("database")
    with c2:
        if st.button("View statistics", use_container_width=True, key="home_view_statistics"):
            _navigate_home("statistics")
    with c3:
        if st.button("Analyze molecules", use_container_width=True, key="home_analyze_molecules"):
            _navigate_home("analyze")

    st.markdown(
        """
        <div class="section-label">About ChalMolDB</div>
        <div class="about-intro">
          <strong>Find structures, compare electronic properties, trace provenance, and reproduce your queries in one chalcogen-focused research environment.</strong><br>
          ChalMolDB is a curated research platform for sulfur-, selenium- and tellurium-containing molecular systems.
        </div>
        <div class="about-grid">
          <div class="about-card">
            <strong>What is ChalMolDB?</strong>
            <p>A searchable, structure-aware database for systematic exploration of chalcogen-focused molecular data.</p>
          </div>
          <div class="about-card">
            <strong>Who is it for?</strong>
            <p>Researchers in computational chemistry, cheminformatics, organic electronics and materials informatics working with S/Se/Te-containing molecules.</p>
          </div>
          <div class="about-card">
            <strong>Where does the data come from?</strong>
            <p>Records are curated from documented computational and literature-derived sources, with original provenance retained at record level.</p>
          </div>
          <div class="about-card">
            <strong>Why is it useful?</strong>
            <p>Search molecular motifs, compare structures, filter HOMO/LUMO/Eg values, inspect provenance, and replay saved queries.</p>
          </div>
        </div>
        <div class="citation-box">
          <strong>How to cite ChalMolDB</strong><br>
          Current archived release: Database v1 · Scientific Core v0.11.0<br>
          Version DOI: 10.5281/zenodo.22903550 · Concept DOI: 10.5281/zenodo.22903549
        </div>
        """,
        unsafe_allow_html=True,
    )

    st.markdown(
        """
        <div class="section-label">What ChalMolDB provides</div>
        <div class="feature-grid">
          <div class="feature-card">
            <div class="feature-icon"><svg viewBox="0 0 32 32" aria-hidden="true"><ellipse cx="16" cy="7" rx="9" ry="4"/><path d="M7 7v9c0 2.2 4 4 9 4s9-1.8 9-4V7"/><path d="M7 16v8c0 2.2 4 4 9 4s9-1.8 9-4v-8"/></svg></div>
            <strong>Curated Molecular Records</strong>
            <span>Standardized molecular identity, electronic properties and record-level provenance.</span>
          </div>
          <div class="feature-card">
            <div class="feature-icon"><svg viewBox="0 0 32 32" aria-hidden="true"><path d="M10 7 16 3l6 4v7l-6 4-6-4Z"/><path d="m22 14 6 4v7l-6 4-6-4v-7"/><circle cx="10" cy="7" r="1.6"/><circle cx="22" cy="14" r="1.6"/></svg></div>
            <strong>Chalcogen-Aware Descriptors</strong>
            <span>Interpretable annotations describing sulfur, selenium and tellurium environments.</span>
          </div>
          <div class="feature-card">
            <div class="feature-icon"><svg viewBox="0 0 32 32" aria-hidden="true"><circle cx="13" cy="13" r="7"/><path d="m18 18 8 8"/><circle cx="13" cy="13" r="2.5"/></svg></div>
            <strong>Similarity Search</strong>
            <span>Structure-based exploratory comparison using molecular fingerprints.</span>
          </div>
          <div class="feature-card">
            <div class="feature-icon"><svg viewBox="0 0 32 32" aria-hidden="true"><path d="M16 4v15"/><path d="m10 14 6 6 6-6"/><path d="M7 24v4h18v-4"/></svg></div>
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
          Database v{DATABASE_VERSION} &nbsp;&middot;&nbsp; Structure standardization with RDKit
          &nbsp;&middot;&nbsp; Record-level provenance retained
        </div>
        <div class="site-footer">
          <span><strong>ChalMolDB</strong> · Chalcogen Molecular Database · Database v{DATABASE_VERSION}</span>
          <span>
            <a href="?page=documentation" target="_self">Documentation</a>
            <a href="https://github.com/EsraDogan539/molecular-descriptor-platform" target="_blank">GitHub</a>
            <a href="?page=about" target="_self">About</a>
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
        f"ChalMolDB · {RELEASE_LABEL} · "
        "Curated records are read-only in the public browser."
    )
    st.stop()


if page == "statistics":
    display_database_statistics()
    st.divider()
    st.caption(
        f"ChalMolDB · Database v{DATABASE_VERSION} · Descriptive statistics"
    )
    st.stop()


if page == "documentation":
    st.markdown(
        """
        <div class="page-intro">
          <div class="page-eyebrow">Reference guide</div>
          <div class="page-title">Documentation</div>
          <div class="page-description">Database scope, structure-first search, molecular comparison, provenance, reproducibility and the user-analysis workflow for ChalMolDB.</div>
        </div>
        """,
        unsafe_allow_html=True,
    )

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

    st.markdown('<div class="section-rule-title">Database scope</div>', unsafe_allow_html=True)
    st.write(
        f"Database v{DATABASE_VERSION} contains curated development and external records. Exact standardized "
        "structures are displayed only where the source supports structure-level identity."
    )

    st.markdown('<div class="section-rule-title">Molecular identity</div>', unsafe_allow_html=True)
    st.write(
        "Structure-complete records are standardized with RDKit and represented by canonical "
        "SMILES, InChI and InChIKey. Repeated standardized structures are retained to preserve "
        "record-level provenance."
    )

    st.markdown('<div class="section-rule-title">Structure-first database search</div>', unsafe_allow_html=True)
    st.write(
        "The Database workspace supports Exact standardized-identity search, RDKit graph "
        "substructure matching and Morgan/Tanimoto similarity search. Similarity uses a Morgan "
        "fingerprint with radius 2 and 2048 bits. Structural similarity is descriptive and does "
        "not establish equivalent electronic properties."
    )

    st.markdown('<div class="section-rule-title">Advanced research filters</div>', unsafe_allow_html=True)
    st.write(
        "Optional filters can constrain HOMO, LUMO and Eg ranges; minimum S, Se and Te counts; "
        "and method, basis-set and curation metadata. Numeric filters are applied only when "
        "explicitly enabled, and missing values do not silently satisfy an active numeric filter."
    )

    st.markdown('<div class="section-rule-title">Curated record comparison</div>', unsafe_allow_html=True)
    st.write(
        "Filtered database records can be compared side by side. The comparison retains molecular "
        "identity, HOMO, LUMO, Eg, experimental Eg where available, chalcogen context, calculation "
        "method and basis-set metadata. Numeric deltas are reported as candidate minus reference; "
        "multi-valued experimental measurements are preserved rather than silently collapsed."
    )

    st.markdown('<div class="section-rule-title">Curation principles</div>', unsafe_allow_html=True)
    st.markdown(
        "- Missing structures or properties are left explicit; unavailable values are not imputed.\n"
        "- S/Se/Te-focused annotations are reported separately from general molecular descriptors.\n"
        "- Record identity and standardized structure identity are treated as distinct concepts."
    )

    st.markdown('<div class="section-rule-title">Descriptor dictionary</div>', unsafe_allow_html=True)
    st.write(
        "Descriptor definitions are versioned with the scientific core so that the terminology "
        "shown in the interface and the exported analysis context remains explicit."
    )
    descriptor_reference_df = descriptor_dictionary_dataframe()
    descriptor_group_filter = st.multiselect(
        "Filter descriptor dictionary",
        options=["Basic", "Structural", "Topological", "Element Counts", "Chalcogen Core", "Quality"],
        default=[],
        key="documentation_descriptor_filter",
        placeholder="All descriptor groups",
    )
    if descriptor_group_filter:
        descriptor_reference_df = descriptor_dictionary_dataframe(descriptor_group_filter)
    st.dataframe(
        descriptor_reference_df,
        use_container_width=True,
        hide_index=True,
        column_config={
            "Descriptor": st.column_config.TextColumn(width="medium"),
            "Group": st.column_config.TextColumn(width="small"),
            "Unit": st.column_config.TextColumn(width="small"),
            "Definition": st.column_config.TextColumn(width="large"),
            "Interpretation": st.column_config.TextColumn(width="large"),
        },
    )
    st.caption(f"Descriptor Dictionary v{DESCRIPTOR_DICTIONARY_VERSION}")

    st.markdown('<div class="section-rule-title">Provenance and data coverage</div>', unsafe_allow_html=True)
    st.write(
        "Record-level provenance is retained internally from the source material, while the public "
        "interface uses neutral collection labels. Calculation method, basis set, curation status "
        "and source/reference text are shown where available. Missing provenance is kept explicit "
        "rather than inferred."
    )
    st.write(
        "The Statistics page reports field-level availability for electronic properties, "
        "standardized structures, method metadata and references. Coverage describes what is "
        "present in the release and should not be interpreted as a quality score."
    )

    st.write(
        "Valid DOI values are normalized to https://doi.org links in record detail views. "
        "Citation-ready record CSV exports preserve record identity, standardized structure, "
        "electronic properties, method/basis metadata, source reference and release versions. "
        "Non-DOI references remain explicit source text."
    )

    st.markdown('<div class="section-rule-title">Release quality checks</div>', unsafe_allow_html=True)
    st.write(
        "The release workflow automatically audits missing or duplicate record IDs, "
        "populated Canonical SMILES that cannot be parsed, non-numeric populated electronic-property "
        "values and SMILES/InChIKey consistency before release."
    )

    st.markdown('<div class="section-rule-title">Reproducibility manifests</div>', unsafe_allow_html=True)
    st.write(
        "Every complete user-analysis package includes a machine-readable run_manifest.json file. "
        "It records the normalized project name, UTC creation time, input row/column structure, "
        "a SHA-256 checksum of the submitted table, analysis counts, Descriptor Dictionary version, "
        "Morgan/MACCS fingerprint settings, software versions and the files included in the package."
    )
    st.write(
        "Database searches can also export a query manifest JSON. It records the raw and standardized "
        "structure query, search mode, active filters, Morgan/Tanimoto settings where applicable, "
        "result count, release versions, a deterministic query-configuration SHA-256 checksum and "
        "an ordered result-record-set SHA-256 checksum."
    )
    st.write(
        "Saved query manifests can be replayed in the Database workspace. ChalMolDB validates the "
        "manifest schema, restores compatible search and filter settings, constrains restored numeric "
        "values to the current release bounds and re-runs the query against the currently loaded "
        "database. A release-difference notice is shown when the saved manifest and current release "
        "do not match."
    )

    st.markdown('<div class="section-rule-title">Analyze your dataset</div>', unsafe_allow_html=True)
    st.write(
        "Upload a CSV containing Molecule_ID and SMILES. The platform validates structures, "
        "calculates general and S/Se/Te-aware descriptors, and keeps user-supplied data separate "
        "from the curated publication database."
    )

    st.markdown('<div class="section-rule-title">Software environment</div>', unsafe_allow_html=True)
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

    st.markdown('<div class="section-rule-title">Citation and data release</div>', unsafe_allow_html=True)
    st.write(
        "The recommended citation and permanent dataset DOI will be added with the archived "
        "publication release."
    )
    st.stop()


if page == "about":
    st.markdown(
        """
        <div class="page-intro">
          <div class="page-eyebrow">Research resource</div>
          <div class="page-title">About ChalMolDB</div>
          <div class="page-description">A curated research resource for chalcogen-focused molecular data, designed around transparent identity, provenance and reusable scientific descriptors.</div>
        </div>
        """,
        unsafe_allow_html=True,
    )

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

    st.markdown('<div class="section-rule-title">Current release</div>', unsafe_allow_html=True)
    st.markdown(
        f"""
        <div class="release-card">
          <strong>{RELEASE_LABEL}</strong>
          <span>3,360 records · 2,983 unique standardized development structures</span>
        </div>
        """,
        unsafe_allow_html=True,
    )

    st.markdown('<div class="section-rule-title">Resources</div>', unsafe_allow_html=True)
    st.markdown(
        """
        <div class="resource-links">
          <a href="?page=documentation" target="_self">Documentation</a>
          <a href="https://github.com/EsraDogan539/molecular-descriptor-platform" target="_blank">GitHub repository</a>
        </div>
        """,
        unsafe_allow_html=True,
    )
    st.stop()


st.markdown(
    """
    <div class="page-intro">
      <div class="page-eyebrow">User analysis workspace</div>
      <div class="page-title">Analyze Your Dataset</div>
      <div class="page-description">Upload molecular structures to validate records and calculate general and S/Se/Te-aware descriptors while keeping user data separate from the curated publication database.</div>
    </div>
    """,
    unsafe_allow_html=True,
)

st.markdown('<div class="section-rule-title">Upload dataset</div>', unsafe_allow_html=True)
upload_left, upload_right = st.columns([1.75, 1], gap="large")
with upload_left:
    st.markdown('<div class="upload-align-anchor"></div>', unsafe_allow_html=True)
    uploaded_file = st.file_uploader(
        "Choose a CSV file",
        type=["csv"],
        help="Required columns: Molecule_ID and SMILES.",
    )
with upload_right:
    st.markdown(
        """
        <div class="upload-side-card">
          <strong>Input format</strong>
          <span>Required: Molecule_ID, SMILES</span>
          <span>Optional property and provenance columns are preserved.</span>
          <span>User uploads remain separate from the curated ChalMolDB release.</span>
        </div>
        """,
        unsafe_allow_html=True,
    )
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
    validate_input_dataframe(input_df)

    if len(input_df) > 5000:
        st.warning(
            f"This file contains {len(input_df):,} rows. Descriptor and fingerprint "
            "generation may take longer for large datasets."
        )

    st.markdown('<div class="section-rule-title">Input preview</div>', unsafe_allow_html=True)
    st.dataframe(display_table(input_df.head(12)), use_container_width=True, hide_index=True)

    with st.expander("Analysis settings", expanded=True):
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

    st.markdown('<div class="analyze-run-anchor"></div>', unsafe_allow_html=True)
    if st.button(
        "Run analysis",
        type="primary",
        disabled=not selected_groups,
    ):
        with st.status("Running molecular analysis...", expanded=True) as status:
            st.write("Validating structures and standardizing molecular identity...")
            results = run_molecular_descriptor_platform(
                input_df=input_df,
                project_name=project_name,
            )
            st.write("Preparing curated metadata and downloadable outputs...")
            results = add_database_export(
                input_df=input_df,
                results=results,
                project_name=project_name,
            )
            status.update(label="Analysis complete", state="complete", expanded=False)

        valid_df = results["valid_df"]
        invalid_df = results["invalid_df"]
        database_df = results["database_df"]
        summary_df = results["summary_df"]
        filtered_valid_df = filter_descriptor_columns(valid_df, selected_groups)

        total_records = int(summary_df.loc[0, "Total Records"])
        valid_count = int(summary_df.loc[0, "Valid Molecules"])
        invalid_count = int(summary_df.loc[0, "Invalid Molecules"])
        duplicate_count = int(summary_df.loc[0, "Duplicate Molecules"])

        st.markdown(
            f"""
            <div class="analysis-summary">
              <div class="analysis-summary-item"><span>Submitted</span><strong>{total_records:,}</strong></div>
              <div class="analysis-summary-item"><span>Valid</span><strong>{valid_count:,}</strong></div>
              <div class="analysis-summary-item"><span>Invalid</span><strong>{invalid_count:,}</strong></div>
              <div class="analysis-summary-item"><span>Repeated</span><strong>{duplicate_count:,}</strong></div>
            </div>
            """,
            unsafe_allow_html=True,
        )

        tabs = st.tabs([
            "Records",
            "Descriptors",
            "Structures",
            "Similarity",
            "Export",
        ])

        with tabs[0]:
            st.markdown('<div class="workspace-title">Records</div>', unsafe_allow_html=True)
            st.caption(
                "Standardized identifiers, retained source fields and quality flags."
            )
            st.dataframe(display_table(database_df), use_container_width=True, hide_index=True)
            if not invalid_df.empty:
                st.markdown("#### Invalid structures")
                st.dataframe(display_table(invalid_df), use_container_width=True, hide_index=True)

        with tabs[1]:
            st.markdown('<div class="workspace-title">Descriptors</div>', unsafe_allow_html=True)
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
            with st.expander("Descriptor dictionary", expanded=False):
                dictionary_groups = [
                    group for group in selected_groups
                    if group in {"Basic", "Structural", "Topological", "Element Counts", "Chalcogen Core"}
                ]
                dictionary_df = descriptor_dictionary_dataframe(dictionary_groups)
                st.caption(
                    f"Scientific definitions and interpretation notes · "
                    f"Dictionary v{DESCRIPTOR_DICTIONARY_VERSION}"
                )
                st.dataframe(
                    dictionary_df,
                    use_container_width=True,
                    hide_index=True,
                    column_config={
                        "Descriptor": st.column_config.TextColumn(width="medium"),
                        "Group": st.column_config.TextColumn(width="small"),
                        "Unit": st.column_config.TextColumn(width="small"),
                        "Definition": st.column_config.TextColumn(width="large"),
                        "Interpretation": st.column_config.TextColumn(width="large"),
                    },
                )
                st.download_button(
                    "Download descriptor dictionary",
                    dictionary_df.to_csv(index=False).encode("utf-8"),
                    f"chalmoldb_descriptor_dictionary_v{DESCRIPTOR_DICTIONARY_VERSION}.csv",
                    "text/csv",
                    key="download_descriptor_dictionary",
                )

            with st.expander("Descriptor summary", expanded=False):
                display_scientific_core_panel(valid_df)

        with tabs[2]:
            st.markdown('<div class="workspace-title">Structures</div>', unsafe_allow_html=True)
            display_molecule_gallery(valid_df, max_molecule_cards)

        with tabs[3]:
            st.markdown('<div class="workspace-title">Similarity Search</div>', unsafe_allow_html=True)
            st.caption(
                "Exploratory fingerprint-based molecular similarity. "
                "Similarity does not imply equivalent electronic behavior."
            )
            display_similarity_search_panel(valid_df)

        with tabs[4]:
            st.markdown('<div class="workspace-title">Export</div>', unsafe_allow_html=True)
            selected_csv = filtered_valid_df.to_csv(index=False).encode("utf-8")
            database_csv = database_df.to_csv(index=False).encode("utf-8")

            with st.container(border=True):
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

            with st.container(border=True):
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

            with open(results["zip_file"], "rb") as file:
                zip_data = file.read()
            with st.container(border=True):
                r3c1, r3c2 = st.columns([3, 1])
                with r3c1:
                    st.markdown("**Complete package**")
                    st.caption(
                        "Processed records, descriptors, fingerprints, curated metadata, "
                        "descriptor dictionary and reproducibility manifest."
                    )
                with r3c2:
                    st.download_button(
                        "Download ZIP",
                        zip_data,
                        os.path.basename(results["zip_file"]),
                        "application/zip",
                        use_container_width=True,
                        key="download_package",
                    )

except pd.errors.EmptyDataError:
    st.error("The uploaded CSV is empty or could not be read as a table.")
except pd.errors.ParserError:
    st.error("The uploaded file could not be parsed as CSV. Check the delimiter and quoting.")
except UnicodeDecodeError:
    st.error("The uploaded CSV could not be decoded. Please save it as UTF-8 and try again.")
except ValueError as error:
    st.error(str(error))
except Exception as error:
    st.error(
        "The dataset could not be processed because of an unexpected error. "
        f"Details: {error}"
    )

st.divider()
st.caption(
    f"ChalMolDB · {RELEASE_LABEL}"
)
