
import pandas as pd
import streamlit as st

from rdkit import Chem, DataStructs
from rdkit.Chem import Draw, MACCSkeys, rdFingerprintGenerator


MORGAN_RADIUS = 2
MORGAN_N_BITS = 2048

morgan_generator = rdFingerprintGenerator.GetMorganGenerator(
    radius=MORGAN_RADIUS,
    fpSize=MORGAN_N_BITS
)


def calculate_similarity(smiles_1, smiles_2):
    """
    İki molekül için Morgan ve MACCS Tanimoto
    benzerlik skorlarını hesaplar.
    """

    mol_1 = Chem.MolFromSmiles(smiles_1)
    mol_2 = Chem.MolFromSmiles(smiles_2)

    if mol_1 is None or mol_2 is None:
        raise ValueError(
            "Moleküllerden biri geçerli bir SMILES içermiyor."
        )

    morgan_fp_1 = morgan_generator.GetFingerprint(mol_1)
    morgan_fp_2 = morgan_generator.GetFingerprint(mol_2)

    maccs_fp_1 = MACCSkeys.GenMACCSKeys(mol_1)
    maccs_fp_2 = MACCSkeys.GenMACCSKeys(mol_2)

    morgan_similarity = DataStructs.TanimotoSimilarity(
        morgan_fp_1,
        morgan_fp_2
    )

    maccs_similarity = DataStructs.TanimotoSimilarity(
        maccs_fp_1,
        maccs_fp_2
    )

    return {
        "Morgan Tanimoto": morgan_similarity,
        "MACCS Tanimoto": maccs_similarity,
        "Average Similarity": (
            morgan_similarity + maccs_similarity
        ) / 2
    }


def classify_similarity(score):
    """
    Benzerlik skorunu yorumlayan kısa sınıflandırma.
    """

    if score >= 0.85:
        return "Çok yüksek benzerlik"

    if score >= 0.65:
        return "Yüksek benzerlik"

    if score >= 0.40:
        return "Orta düzey benzerlik"

    if score >= 0.20:
        return "Düşük benzerlik"

    return "Çok düşük benzerlik"


def create_descriptor_comparison(row_1, row_2):
    """
    Seçilen iki molekülün temel deskriptörlerini
    karşılaştıran tabloyu oluşturur.
    """

    descriptor_definitions = [
        ("Molecular Weight", "Moleküler Ağırlık"),
        ("Exact Molecular Weight", "Kesin Moleküler Ağırlık"),
        ("LogP", "LogP"),
        ("TPSA", "TPSA"),
        ("H-Bond Donors", "H-Bond Verici"),
        ("H-Bond Acceptors", "H-Bond Alıcı"),
        ("Rotatable Bonds", "Dönebilir Bağ"),
        ("Ring Count", "Halka Sayısı"),
        ("Aromatic Ring Count", "Aromatik Halka"),
        ("Heavy Atom Count", "Ağır Atom"),
        ("Heteroatom Count", "Heteroatom"),
        ("Total Chalcogen Count", "Toplam Kalkojen"),
        ("Heavy Chalcogen Count", "Ağır Kalkojen")
    ]

    comparison_rows = []

    for column_name, display_name in descriptor_definitions:

        if (
            column_name not in row_1.index
            or column_name not in row_2.index
        ):
            continue

        value_1 = row_1[column_name]
        value_2 = row_2[column_name]

        try:
            difference = round(
                float(value_1) - float(value_2),
                4
            )
        except (TypeError, ValueError):
            difference = None

        comparison_rows.append({
            "Deskriptör": display_name,
            str(row_1["Molecule_ID"]): value_1,
            str(row_2["Molecule_ID"]): value_2,
            "Fark (1 - 2)": difference
        })

    return pd.DataFrame(comparison_rows)


def display_molecule_structure(row):
    """
    Tek molekülün yapı kartını gösterir.
    """

    smiles = row["Canonical SMILES"]
    mol = Chem.MolFromSmiles(smiles)

    st.markdown(
        f"### {row['Molecule_ID']}"
    )

    if mol is not None:
        image = Draw.MolToImage(
            mol,
            size=(500, 350)
        )

        st.image(
            image,
            use_container_width=True
        )

    st.write(
        f"**Moleküler formül:** "
        f"{row['Molecular Formula']}"
    )

    st.code(
        smiles,
        language=None
    )


def display_similarity_panel(valid_df):
    """
    Molekül seçimi ve karşılaştırma arayüzünü gösterir.
    """

    st.subheader(
        "Molekül Benzerlik ve Karşılaştırma Paneli"
    )

    if len(valid_df) < 2:
        st.info(
            "Benzerlik analizi için en az iki geçerli "
            "molekül gereklidir."
        )
        return

    molecule_options = valid_df[
        "Molecule_ID"
    ].astype(str).tolist()

    selector_col_1, selector_col_2 = st.columns(2)

    with selector_col_1:
        selected_id_1 = st.selectbox(
            "Birinci molekül",
            options=molecule_options,
            index=0,
            key="similarity_molecule_1"
        )

    second_default_index = (
        1 if len(molecule_options) > 1 else 0
    )

    with selector_col_2:
        selected_id_2 = st.selectbox(
            "İkinci molekül",
            options=molecule_options,
            index=second_default_index,
            key="similarity_molecule_2"
        )

    if selected_id_1 == selected_id_2:
        st.warning(
            "Karşılaştırmak için iki farklı molekül seçin."
        )
        return

    row_1 = valid_df[
        valid_df["Molecule_ID"].astype(str)
        == selected_id_1
    ].iloc[0]

    row_2 = valid_df[
        valid_df["Molecule_ID"].astype(str)
        == selected_id_2
    ].iloc[0]

    similarity_results = calculate_similarity(
        row_1["Canonical SMILES"],
        row_2["Canonical SMILES"]
    )

    average_similarity = similarity_results[
        "Average Similarity"
    ]

    metric_col_1, metric_col_2, metric_col_3 = (
        st.columns(3)
    )

    metric_col_1.metric(
        "Morgan Tanimoto",
        f"{similarity_results['Morgan Tanimoto']:.3f}"
    )

    metric_col_2.metric(
        "MACCS Tanimoto",
        f"{similarity_results['MACCS Tanimoto']:.3f}"
    )

    metric_col_3.metric(
        "Ortalama Benzerlik",
        f"{average_similarity:.3f}"
    )

    st.info(
        f"Yorum: **{classify_similarity(average_similarity)}**"
    )

    structure_col_1, structure_col_2 = st.columns(2)

    with structure_col_1:
        with st.container(border=True):
            display_molecule_structure(row_1)

    with structure_col_2:
        with st.container(border=True):
            display_molecule_structure(row_2)

    st.markdown("### Deskriptör Karşılaştırması")

    comparison_df = create_descriptor_comparison(
        row_1,
        row_2
    )

    st.dataframe(
        comparison_df,
        hide_index=True,
        use_container_width=True
    )

    comparison_csv = comparison_df.to_csv(
        index=False
    ).encode("utf-8")

    st.download_button(
        label="Karşılaştırma Tablosunu CSV Olarak İndir",
        data=comparison_csv,
        file_name=(
            f"{selected_id_1}_vs_"
            f"{selected_id_2}_comparison.csv"
        ),
        mime="text/csv"
    )
