
import numpy as np
import pandas as pd
import streamlit as st


def create_histogram_dataframe(series, bins=10):
    """
    Sayısal bir pandas serisinden histogram tablosu oluşturur.
    """

    clean_series = pd.to_numeric(
        series,
        errors="coerce"
    ).dropna()

    if clean_series.empty:
        return pd.DataFrame(
            columns=["Aralık", "Molekül Sayısı"]
        )

    unique_count = clean_series.nunique()

    if unique_count == 1:
        value = clean_series.iloc[0]

        return pd.DataFrame({
            "Aralık": [str(round(value, 3))],
            "Molekül Sayısı": [len(clean_series)]
        })

    selected_bins = min(
        bins,
        max(3, unique_count)
    )

    counts, bin_edges = np.histogram(
        clean_series,
        bins=selected_bins
    )

    interval_labels = [
        f"{bin_edges[index]:.2f}–{bin_edges[index + 1]:.2f}"
        for index in range(len(bin_edges) - 1)
    ]

    return pd.DataFrame({
        "Aralık": interval_labels,
        "Molekül Sayısı": counts
    })


def display_numeric_summary(valid_df):
    """
    Temel sayısal deskriptörlerin özet istatistiklerini gösterir.
    """

    st.markdown("### Temel İstatistikler")

    descriptors = [
        ("Molecular Weight", "Moleküler Ağırlık"),
        ("LogP", "LogP"),
        ("TPSA", "TPSA")
    ]

    available_descriptors = [
        item
        for item in descriptors
        if item[0] in valid_df.columns
    ]

    if not available_descriptors:
        st.info(
            "Özetlenecek sayısal deskriptör bulunamadı."
        )
        return

    summary_columns = st.columns(
        len(available_descriptors)
    )

    for column, (descriptor, label) in zip(
        summary_columns,
        available_descriptors
    ):
        numeric_series = pd.to_numeric(
            valid_df[descriptor],
            errors="coerce"
        ).dropna()

        if numeric_series.empty:
            continue

        with column:
            with st.container(border=True):
                st.markdown(f"#### {label}")

                st.metric(
                    "Ortalama",
                    f"{numeric_series.mean():.3f}"
                )

                stat_col1, stat_col2 = st.columns(2)

                stat_col1.metric(
                    "Minimum",
                    f"{numeric_series.min():.3f}"
                )

                stat_col2.metric(
                    "Maksimum",
                    f"{numeric_series.max():.3f}"
                )


def display_descriptor_distributions(valid_df):
    """
    Moleküler ağırlık, LogP ve TPSA histogramlarını gösterir.
    """

    st.markdown("### Deskriptör Dağılımları")

    chart_definitions = [
        (
            "Molecular Weight",
            "Moleküler Ağırlık Dağılımı"
        ),
        (
            "LogP",
            "LogP Dağılımı"
        ),
        (
            "TPSA",
            "TPSA Dağılımı"
        )
    ]

    for descriptor, title in chart_definitions:

        if descriptor not in valid_df.columns:
            continue

        histogram_df = create_histogram_dataframe(
            valid_df[descriptor],
            bins=10
        )

        if histogram_df.empty:
            continue

        st.markdown(f"#### {title}")

        chart_df = histogram_df.set_index("Aralık")

        st.bar_chart(
            chart_df,
            y="Molekül Sayısı",
            use_container_width=True
        )


def display_chalcogen_analysis(valid_df):
    """
    S, Se ve Te içeren molekül sayılarını gösterir.
    """

    st.markdown("### Kalkojen Dağılımı")

    chalcogen_definitions = [
        ("Contains S", "S içeren"),
        ("Contains Se", "Se içeren"),
        ("Contains Te", "Te içeren")
    ]

    chalcogen_counts = {}

    for column_name, label in chalcogen_definitions:

        if column_name in valid_df.columns:
            values = pd.to_numeric(
                valid_df[column_name],
                errors="coerce"
            ).fillna(0)

            chalcogen_counts[label] = int(
                (values > 0).sum()
            )

    if not chalcogen_counts:
        st.info(
            "Kalkojen analizi için gerekli sütunlar bulunamadı."
        )
        return

    metric_columns = st.columns(
        len(chalcogen_counts)
    )

    for metric_column, (label, count) in zip(
        metric_columns,
        chalcogen_counts.items()
    ):
        metric_column.metric(
            label,
            count
        )

    chalcogen_chart_df = pd.DataFrame({
        "Kalkojen": list(chalcogen_counts.keys()),
        "Molekül Sayısı": list(
            chalcogen_counts.values()
        )
    }).set_index("Kalkojen")

    st.bar_chart(
        chalcogen_chart_df,
        y="Molekül Sayısı",
        use_container_width=True
    )


def display_formula_frequency(valid_df, top_n=10):
    """
    En sık görülen moleküler formülleri gösterir.
    """

    if "Molecular Formula" not in valid_df.columns:
        return

    formula_counts = (
        valid_df["Molecular Formula"]
        .dropna()
        .astype(str)
        .value_counts()
        .head(top_n)
        .rename_axis("Moleküler Formül")
        .reset_index(name="Molekül Sayısı")
    )

    if formula_counts.empty:
        return

    st.markdown("### En Sık Görülen Moleküler Formüller")

    col1, col2 = st.columns([1, 2])

    with col1:
        st.dataframe(
            formula_counts,
            hide_index=True,
            use_container_width=True
        )

    with col2:
        st.bar_chart(
            formula_counts.set_index(
                "Moleküler Formül"
            ),
            y="Molekül Sayısı",
            use_container_width=True
        )


def display_analysis_panel(valid_df):
    """
    Moleküler veri seti için genel analiz panelini gösterir.
    """

    st.subheader("Veri Seti Analiz Paneli")

    if valid_df.empty:
        st.info(
            "Analiz edilecek geçerli molekül bulunamadı."
        )
        return

    st.caption(
        f"Analiz edilen geçerli molekül sayısı: {len(valid_df)}"
    )

    display_numeric_summary(valid_df)

    st.markdown("---")

    display_descriptor_distributions(valid_df)

    st.markdown("---")

    display_chalcogen_analysis(valid_df)

    st.markdown("---")

    display_formula_frequency(valid_df)
