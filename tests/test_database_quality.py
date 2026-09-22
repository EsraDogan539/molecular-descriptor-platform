import unittest

import pandas as pd

from database_quality import (
    available_reference_column,
    classify_numeric_measurement,
    coverage_table,
    database_quality_summary,
    provenance_label,
    reference_mask,
)


class DatabaseQualityTests(unittest.TestCase):
    def setUp(self):
        self.df = pd.DataFrame({
            "Record_ID": ["A", "B", "C", "C"],
            "Canonical_SMILES": ["c1ccsc1", "", None, "c1cc[se]c1"],
            "InChIKey": ["KEY1", "KEY2", "", "KEY1"],
            "HOMO_eV": [-5.1, None, -4.9, -5.0],
            "LUMO_eV": [-2.8, -2.7, None, -2.9],
            "Eg_eV": [2.3, 2.1, None, 2.1],
            "Experimental_Eg_eV": [None, None, None, 2.0],
            "Dataset_Owner": ["Set A", "Set A", "Set B", "Set B"],
            "Method": ["DFT", "", None, "DFT"],
            "Reference": ["10.1000/a", "", None, "Paper B"],
        })

    def test_reference_detection_is_schema_aware(self):
        self.assertEqual(available_reference_column(self.df), "Reference")
        self.assertEqual(reference_mask(self.df).tolist(), [True, False, False, True])

    def test_coverage_table_reports_available_and_missing(self):
        coverage = coverage_table(self.df).set_index("Field")
        self.assertEqual(int(coverage.loc["Eg", "Available"]), 3)
        self.assertEqual(int(coverage.loc["Eg", "Missing"]), 1)
        self.assertEqual(int(coverage.loc["Standardized structure", "Available"]), 2)
        self.assertEqual(int(coverage.loc["Reference / source", "Available"]), 2)

    def test_quality_summary_flags_record_id_duplicates(self):
        summary = database_quality_summary(self.df)
        self.assertEqual(summary["Missing Record IDs"], 0)
        self.assertEqual(summary["Duplicate Record IDs"], 2)
        self.assertEqual(summary["Repeated Standardized Structures"], 2)

    def test_provenance_label_prefers_reference_then_method(self):
        self.assertEqual(provenance_label(self.df.iloc[0]), "Reference available")
        self.assertEqual(provenance_label(self.df.iloc[1]), "Dataset provenance")
        self.assertEqual(provenance_label(self.df.iloc[3]), "Reference available")

    def test_reference_helpers_handle_missing_reference_column(self):
        no_ref = self.df.drop(columns=["Reference"])
        self.assertIsNone(available_reference_column(no_ref))
        self.assertFalse(reference_mask(no_ref).any())

    def test_multi_valued_experimental_measurements_are_preserved(self):
        self.assertEqual(classify_numeric_measurement("1.80; 1.50"), "multi-valued")
        self.assertEqual(classify_numeric_measurement("1.46; 0.83"), "multi-valued")
        self.assertEqual(classify_numeric_measurement("2.10"), "scalar")
        self.assertEqual(classify_numeric_measurement("reported"), "non-numeric")

        df = pd.DataFrame({
            "Experimental_Eg_eV": ["1.80; 1.50", "2.10", None, "1.46; 0.83"],
        })
        coverage = coverage_table(df).set_index("Field")
        self.assertEqual(int(coverage.loc["Experimental Eg", "Available"]), 3)
        self.assertEqual(
            int(coverage.loc["Experimental Eg (multi-valued records)", "Available"]),
            2,
        )


if __name__ == "__main__":
    unittest.main()
