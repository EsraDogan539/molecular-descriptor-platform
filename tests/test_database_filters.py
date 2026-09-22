import pandas as pd

from database_filters import apply_advanced_filters, numeric_bounds


def sample_df():
    return pd.DataFrame(
        {
            "Record_ID": ["A", "B", "C", "D"],
            "HOMO_eV": [-5.8, -5.2, -4.9, None],
            "LUMO_eV": [-3.2, -2.7, -2.1, None],
            "Eg_eV": [2.6, 2.5, 2.8, None],
            "S_Count": [2, 0, 1, None],
            "Se_Count": [0, 2, 1, None],
            "Te_Count": [0, 0, 1, None],
            "Dataset_Owner": ["Development", "Development", "External", "External"],
            "Method": ["DFT", "DFT", "Experiment", None],
            "Basis_Set": ["B3LYP", "B3LYP", "NA", None],
            "Curation_Status": ["Ready", "Ready", "Reviewed", "Reviewed"],
        }
    )


def test_numeric_bounds_ignores_missing_values():
    low, high = numeric_bounds(sample_df(), "HOMO_eV")
    assert low == -5.8
    assert high == -4.9


def test_no_advanced_filters_preserves_all_rows():
    df = sample_df()
    result = apply_advanced_filters(df)
    assert result["Record_ID"].tolist() == ["A", "B", "C", "D"]


def test_numeric_range_excludes_missing_and_outside_values():
    result = apply_advanced_filters(
        sample_df(),
        numeric_ranges={"HOMO_eV": (-5.3, -4.8)},
    )
    assert result["Record_ID"].tolist() == ["B", "C"]


def test_minimum_chalcogen_count_is_explicit():
    result = apply_advanced_filters(
        sample_df(),
        minimum_counts={"Se_Count": 1},
    )
    assert result["Record_ID"].tolist() == ["B", "C"]


def test_categorical_filters_use_exact_matching():
    result = apply_advanced_filters(
        sample_df(),
        categorical_filters={
            "Dataset_Owner": "External",
            "Curation_Status": "Reviewed",
        },
    )
    assert result["Record_ID"].tolist() == ["C", "D"]


def test_filters_compose_across_categories():
    result = apply_advanced_filters(
        sample_df(),
        numeric_ranges={"Eg_eV": (2.7, 3.0)},
        minimum_counts={"Te_Count": 1},
        categorical_filters={"Dataset_Owner": "External"},
    )
    assert result["Record_ID"].tolist() == ["C"]
