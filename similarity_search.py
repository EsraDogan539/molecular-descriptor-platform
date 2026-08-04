
import pandas as pd
import streamlit as st

from rdkit import Chem, DataStructs
from rdkit.Chem import Draw, rdFingerprintGenerator


MORGAN_RADIUS = 2
MORGAN_N_BITS = 2048

morgan_generator = rdFingerprintGenerator.GetMorganGenerator(
    radius=MORGAN_RADIUS,
    fpSize=MORGAN_N_BITS
)


def calculate_morgan_fingerprint(smiles):
    mol = Chem.MolFromSmiles(smiles)

    if mol is None:
        return None

    return morgan_generator.GetFingerprint(mol)


def find_similar_molecules(
    valid_df,
    reference_id,
    top_n=10,
    minimum_similarity=0.0
):
    reference_rows = valid_df[
        valid_df["Molecule_ID"].astype(str) == str(reference_id)
    ]

    if reference_rows.empty:
        raise ValueError("Referans molekül bulunamadı.")

    reference_row = reference_rows.iloc[0]
    reference_smiles = reference_row["Canonical SMILES"]

    reference_fp = calculate_morgan_fingerprint(
        reference_smiles
    )

    if reference_fp is None:
        raise ValueError(
            "Referans molekülün fingerprint'i oluşturulamadı."
        )

    results = []

    for _, row in valid_df.iterrows():

        molecule_id = str(row["Molecule_ID"])

        if molecule_id == str(reference_id):
            continue

        target_fp = calculate_morgan_fingerprint(
            row["Canonical SMILES"]
        )

        if target_fp is None:
            continue

        similarity = DataStructs.TanimotoSimilarity(
            reference_fp,
            target_fp
        )

        if similarity < minimum_similarity:
            continue

        results.append({
            "Molecule_ID": molecule_id,
            "Canonical SMILES": row["Canonical SMILES"],
            "Molecular Formula": row["Molecular Formula"],
            "Morgan Similarity": round(similarity, 4),
            "Molecular Weight": row["Molecular Weight"],
            "LogP": row["LogP"],
            "TPSA": row["TPSA"],
            "MW Difference": round(
                abs(
                    float(row["Molecular Weight"])
                    - float(reference_row["Molecular Weight"])
                ),
                4
            ),
            "LogP Difference": round(
                abs(
                    float(row["LogP"])
                    - float(reference_row["LogP"])
                ),
                4
            ),
            "TPSA Difference": round(
                abs(
                    float(row["TPSA"])
                    - float(reference_row["TPSA"])
                ),
                4
            )
        })

    results_df = pd.DataFrame(results)

    if results_df.empty:
        return results_df

    results_df = (
        results_df
        .sort_values(
            by="Morgan Similarity",
            ascending=False
        )
        .head(top_n)
        .reset_index(drop=True)
    )

    results_df.insert(
        0,
        "Rank",
        range(1, len(results_df) + 1)
    )

    return results_df


def display_similarity_search_panel(valid_df):
    """
    Kullanıcının bir referans molekül seçerek veri setindeki
    en benzer molekülleri sıralamasını sağlar.
    """

    st.subheader("Benzer Molekül Arama")

    if len(valid_df) < 2:
        st.info(
            "Benzer molekül araması için en az iki geçerli "
            "molekül gereklidir."
        )
        return

    molecule_options = (
        valid_df["Molecule_ID"]
        .astype(str)
        .tolist()
    )

    control_col_1, control_col_2, control_col_3 = st.columns(3)

    with control_col_1:
        reference_id = st.selectbox(
            "Referans molekül",
            options=molecule_options,
            key="similarity_search_reference"
        )

    maximum_result_count = min(
        20,
        len(valid_df) - 1
    )

    with control_col_2:
        top_n = st.number_input(
            "Gösterilecek sonuç sayısı",
            min_value=1,
            max_value=maximum_result_count,
            value=min(5, maximum_result_count),
            step=1
        )

    with control_col_3:
        minimum_similarity = st.slider(
            "Minimum benzerlik",
            min_value=0.0,
            max_value=1.0,
            value=0.0,
            step=0.05
        )

    reference_row = valid_df[
        valid_df["Molecule_ID"].astype(str) == reference_id
    ].iloc[0]

    st.markdown("### Referans Molekül")

    reference_col_1, reference_col_2 = st.columns([1, 2])

    with reference_col_1:
        reference_mol = Chem.MolFromSmiles(
            reference_row["Canonical SMILES"]
        )

        if reference_mol is not None:
            reference_image = Draw.MolToImage(
                reference_mol,
                size=(450, 320)
            )

            st.image(
                reference_image,
                use_container_width=True
            )

    with reference_col_2:
        st.write(
            f"**Molecule ID:** {reference_row['Molecule_ID']}"
        )

        st.write(
            f"**Moleküler formül:** "
            f"{reference_row['Molecular Formula']}"
        )

        st.write(
            f"**Moleküler ağırlık:** "
            f"{reference_row['Molecular Weight']}"
        )

        st.write(
            f"**LogP:** {reference_row['LogP']}"
        )

        st.write(
            f"**TPSA:** {reference_row['TPSA']}"
        )

        st.code(
            reference_row["Canonical SMILES"],
            language=None
        )

    similar_df = find_similar_molecules(
        valid_df=valid_df,
        reference_id=reference_id,
        top_n=int(top_n),
        minimum_similarity=minimum_similarity
    )

    st.markdown("### Benzerlik Sonuçları")

    if similar_df.empty:
        st.warning(
            "Seçilen benzerlik eşiğine uygun molekül bulunamadı."
        )
        return

    st.dataframe(
        similar_df,
        hide_index=True,
        use_container_width=True
    )

    st.markdown("### En Benzer Molekül Kartları")

    cards_per_row = 3

    for start_index in range(
        0,
        len(similar_df),
        cards_per_row
    ):
        card_columns = st.columns(cards_per_row)

        for column_index in range(cards_per_row):
            result_index = start_index + column_index

            if result_index >= len(similar_df):
                break

            result_row = similar_df.iloc[result_index]

            with card_columns[column_index]:
                with st.container(border=True):

                    st.markdown(
                        f"### #{result_row['Rank']} "
                        f"{result_row['Molecule_ID']}"
                    )

                    mol = Chem.MolFromSmiles(
                        result_row["Canonical SMILES"]
                    )

                    if mol is not None:
                        molecule_image = Draw.MolToImage(
                            mol,
                            size=(400, 300)
                        )

                        st.image(
                            molecule_image,
                            use_container_width=True
                        )

                    st.metric(
                        "Morgan Benzerliği",
                        f"{result_row['Morgan Similarity']:.3f}"
                    )

                    st.write(
                        f"**Formül:** "
                        f"{result_row['Molecular Formula']}"
                    )

                    st.write(
                        f"**MW farkı:** "
                        f"{result_row['MW Difference']}"
                    )

                    st.write(
                        f"**LogP farkı:** "
                        f"{result_row['LogP Difference']}"
                    )

                    st.write(
                        f"**TPSA farkı:** "
                        f"{result_row['TPSA Difference']}"
                    )

                    st.code(
                        result_row["Canonical SMILES"],
                        language=None
                    )

    similarity_csv = similar_df.to_csv(
        index=False
    ).encode("utf-8")

    st.download_button(
        label="Benzerlik Sonuçlarını CSV Olarak İndir",
        data=similarity_csv,
        file_name=(
            f"{reference_id}_similarity_search.csv"
        ),
        mime="text/csv"
    )
