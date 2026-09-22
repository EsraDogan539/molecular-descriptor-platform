import unittest

import pandas as pd

from database_quality import (
    available_reference_column,
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
        self.assertEqual(int(coverage.loc["Reference / DOI", "Available"]), 2)

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


if __name__ == "__main__":
    unittest.main()
