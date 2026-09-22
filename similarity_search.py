import pandas as pd
import streamlit as st

from rdkit import Chem, DataStructs
from rdkit.Chem import Draw, rdFingerprintGenerator


MORGAN_RADIUS = 2
MORGAN_N_BITS = 2048
SIMILARITY_METRIC = "Tanimoto"
FINGERPRINT_METHOD = "Morgan"

morgan_generator = rdFingerprintGenerator.GetMorganGenerator(
    radius=MORGAN_RADIUS,
    fpSize=MORGAN_N_BITS,
)


def calculate_morgan_fingerprint(smiles):
    mol = Chem.MolFromSmiles(str(smiles))

    if mol is None:
        return None

    return morgan_generator.GetFingerprint(mol)


def _safe_float(value):
    try:
        if pd.isna(value):
            return None
        return float(value)
    except (TypeError, ValueError):
        return None


def _absolute_difference(target, reference):
    target_value = _safe_float(target)
    reference_value = _safe_float(reference)
    if target_value is None or reference_value is None:
        return None
    return round(abs(target_value - reference_value), 4)


def find_similar_molecules(
    valid_df,
    reference_id,
    top_n=10,
    minimum_similarity=0.0,
):
    reference_rows = valid_df[
        valid_df["Molecule_ID"].astype(str) == str(reference_id)
    ]

    if reference_rows.empty:
        raise ValueError("Reference molecule was not found.")

    reference_row = reference_rows.iloc[0]
    reference_smiles = reference_row["Canonical SMILES"]
    reference_fp = calculate_morgan_fingerprint(reference_smiles)

    if reference_fp is None:
        raise ValueError(
            "A fingerprint could not be generated for the reference molecule."
        )

    comparison_columns = [
        "Molecule_ID",
        "Canonical SMILES",
        "Molecular Formula",
        "Molecular Weight",
        "LogP",
        "TPSA",
        "Chalcogen Type",
        "Target Chalcogen Count",
        "Target Chalcogen Fraction",
        "Conjugated Atom Fraction",
        "Aromatic Bond Fraction",
    ]
    available_columns = [
        column for column in comparison_columns
        if column in valid_df.columns
    ]

    target_rows = []
    target_fingerprints = []

    for _, row in valid_df[available_columns].iterrows():
        molecule_id = str(row["Molecule_ID"])

        if molecule_id == str(reference_id):
            continue

        target_fp = calculate_morgan_fingerprint(row["Canonical SMILES"])
        if target_fp is None:
            continue

        target_rows.append(row)
        target_fingerprints.append(target_fp)

    if not target_fingerprints:
        return pd.DataFrame()

    similarities = DataStructs.BulkTanimotoSimilarity(
        reference_fp,
        target_fingerprints,
    )

    results = []
    for row, similarity in zip(target_rows, similarities):
        if similarity < minimum_similarity:
            continue

        result = {
            "Molecule_ID": str(row["Molecule_ID"]),
            "Canonical SMILES": row["Canonical SMILES"],
            "Molecular Formula": row.get("Molecular Formula"),
            "Morgan Similarity": round(float(similarity), 4),
            "Molecular Weight": row.get("Molecular Weight"),
            "LogP": row.get("LogP"),
            "TPSA": row.get("TPSA"),
            "MW Difference": _absolute_difference(
                row.get("Molecular Weight"),
                reference_row.get("Molecular Weight"),
            ),
            "LogP Difference": _absolute_difference(
                row.get("LogP"),
                reference_row.get("LogP"),
            ),
            "TPSA Difference": _absolute_difference(
                row.get("TPSA"),
                reference_row.get("TPSA"),
            ),
        }

        for column in [
            "Chalcogen Type",
            "Target Chalcogen Count",
            "Target Chalcogen Fraction",
            "Conjugated Atom Fraction",
            "Aromatic Bond Fraction",
        ]:
            if column in row.index:
                result[column] = row.get(column)

        results.append(result)

    results_df = pd.DataFrame(results)

    if results_df.empty:
        return results_df

    results_df = (
        results_df
        .sort_values(by="Morgan Similarity", ascending=False)
        .head(top_n)
        .reset_index(drop=True)
    )
    results_df.insert(0, "Rank", range(1, len(results_df) + 1))
    return results_df


def prepare_similarity_export(
    similar_df,
    reference_id,
    minimum_similarity,
):
    export_df = similar_df.copy()
    export_df.insert(1, "Reference Molecule_ID", str(reference_id))
    export_df["Fingerprint Method"] = FINGERPRINT_METHOD
    export_df["Morgan Radius"] = MORGAN_RADIUS
    export_df["Fingerprint Bits"] = MORGAN_N_BITS
    export_df["Similarity Metric"] = SIMILARITY_METRIC
    export_df["Minimum Similarity"] = minimum_similarity
    return export_df


def _format_value(value, decimals=4):
    if value is None or pd.isna(value):
        return "—"
    if isinstance(value, (int, float)):
        return f"{float(value):.{decimals}f}"
    return str(value)


def _comparison_table(reference_row, target_row):
    fields = [
        ("Molecular Weight", "g/mol"),
        ("LogP", ""),
        ("TPSA", "Å²"),
        ("Target Chalcogen Count", "count"),
        ("Target Chalcogen Fraction", "fraction"),
        ("Conjugated Atom Fraction", "fraction"),
        ("Aromatic Bond Fraction", "fraction"),
    ]

    rows = []
    for field, unit in fields:
        if field not in reference_row.index or field not in target_row.index:
            continue

        reference_value = reference_row.get(field)
        target_value = target_row.get(field)
        delta = _absolute_difference(target_value, reference_value)

        rows.append({
            "Descriptor": field,
            "Unit": unit or "—",
            "Reference": reference_value,
            "Candidate": target_value,
            "Absolute difference": delta,
        })

    return pd.DataFrame(rows)


def _render_molecule_summary(row, title):
    st.markdown(f"**{title}**")
    mol = Chem.MolFromSmiles(str(row["Canonical SMILES"]))
    if mol is not None:
        st.image(
            Draw.MolToImage(mol, size=(320, 225)),
            use_container_width=True,
        )

    st.caption(
        f"{row.get('Molecular Formula', '—')} · "
        f"{row.get('Chalcogen Type', '—')}"
    )
    st.markdown(
        f"**MW:** {_format_value(row.get('Molecular Weight'))} &nbsp;&nbsp; "
        f"**LogP:** {_format_value(row.get('LogP'))} &nbsp;&nbsp; "
        f"**TPSA:** {_format_value(row.get('TPSA'))}"
    )
    st.code(str(row["Canonical SMILES"]), language=None)


def display_similarity_search_panel(valid_df):
    """Rank and compare molecules by Morgan/Tanimoto structural similarity."""

    if len(valid_df) < 2:
        st.info(
            "At least two valid molecules are required for similarity search."
        )
        return

    molecule_options = valid_df["Molecule_ID"].astype(str).tolist()
    control_col_1, control_col_2, control_col_3 = st.columns(
        [1.35, 0.75, 1.1]
    )

    with control_col_1:
        reference_id = st.selectbox(
            "Reference molecule",
            options=molecule_options,
            key="similarity_search_reference",
        )

    maximum_result_count = min(20, len(valid_df) - 1)

    with control_col_2:
        top_n = st.number_input(
            "Number of results",
            min_value=1,
            max_value=maximum_result_count,
            value=min(5, maximum_result_count),
            step=1,
        )

    with control_col_3:
        minimum_similarity = st.slider(
            "Minimum Tanimoto similarity",
            min_value=0.0,
            max_value=1.0,
            value=0.0,
            step=0.05,
        )

    reference_row = valid_df[
        valid_df["Molecule_ID"].astype(str) == reference_id
    ].iloc[0]

    st.caption(
        "Structural similarity · Morgan fingerprint, radius 2, 2048 bits · "
        "Tanimoto coefficient. A high fingerprint similarity does not establish "
        "equivalent electronic properties."
    )

    similar_df = find_similar_molecules(
        valid_df=valid_df,
        reference_id=reference_id,
        top_n=int(top_n),
        minimum_similarity=minimum_similarity,
    )

    st.markdown(
        '<div class="section-rule-title">Similarity results</div>',
        unsafe_allow_html=True,
    )

    if similar_df.empty:
        st.warning(
            "No molecules meet the selected similarity threshold."
        )
        return

    compact_columns = [
        "Rank",
        "Molecule_ID",
        "Molecular Formula",
        "Chalcogen Type",
        "Morgan Similarity",
        "MW Difference",
        "LogP Difference",
        "TPSA Difference",
    ]
    compact_columns = [
        column for column in compact_columns
        if column in similar_df.columns
    ]

    display_df = similar_df[compact_columns].rename(columns={
        "Molecule_ID": "Molecule ID",
        "Molecular Formula": "Molecular Formula",
        "Chalcogen Type": "Chalcogen",
        "Morgan Similarity": "Similarity",
        "MW Difference": "ΔMW",
        "LogP Difference": "ΔLogP",
        "TPSA Difference": "ΔTPSA",
    })

    st.dataframe(
        display_df,
        hide_index=True,
        use_container_width=True,
        column_config={
            "Similarity": st.column_config.ProgressColumn(
                min_value=0.0,
                max_value=1.0,
                format="%.3f",
            ),
        },
    )

    st.markdown(
        '<div class="section-rule-title">Pairwise comparison</div>',
        unsafe_allow_html=True,
    )
    candidate_id = st.selectbox(
        "Candidate molecule",
        options=similar_df["Molecule_ID"].astype(str).tolist(),
        key="similarity_comparison_candidate",
        help="Select one ranked result for a direct comparison with the reference molecule.",
    )
    candidate_row = valid_df[
        valid_df["Molecule_ID"].astype(str) == candidate_id
    ].iloc[0]
    similarity_value = float(
        similar_df.loc[
            similar_df["Molecule_ID"].astype(str) == candidate_id,
            "Morgan Similarity",
        ].iloc[0]
    )

    ref_col, candidate_col = st.columns(2, gap="large")
    with ref_col:
        _render_molecule_summary(
            reference_row,
            f"Reference · {reference_id}",
        )
    with candidate_col:
        _render_molecule_summary(
            candidate_row,
            f"Candidate · {candidate_id}",
        )

    st.metric("Morgan / Tanimoto similarity", f"{similarity_value:.3f}")
    st.caption(
        "Descriptor differences below are descriptive context only. "
        "They are not components of the Tanimoto similarity score."
    )

    comparison_df = _comparison_table(reference_row, candidate_row)
    if not comparison_df.empty:
        st.dataframe(
            comparison_df,
            hide_index=True,
            use_container_width=True,
        )

    export_df = prepare_similarity_export(
        similar_df,
        reference_id=reference_id,
        minimum_similarity=minimum_similarity,
    )
    similarity_csv = export_df.to_csv(index=False).encode("utf-8")

    st.download_button(
        label="Download similarity results CSV",
        data=similarity_csv,
        file_name=f"{reference_id}_similarity_search.csv",
        mime="text/csv",
    )
