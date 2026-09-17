import unittest
import pandas as pd

from descriptor_engine import (
    calculate_single_molecule_descriptors,
    process_molecular_dataset,
)


class ScientificCoreTests(unittest.TestCase):
    def test_sulfur_classification_and_aromatic_environment(self):
        result, error = calculate_single_molecule_descriptors(
            "S_001", "c1ccsc1"
        )
        self.assertIsNone(error)
        self.assertEqual(result["Chalcogen Type"], "S")
        self.assertEqual(result["Heavy Chalcogen Count"], 1)
        self.assertEqual(result["Aromatic Chalcogen Count"], 1)
        self.assertEqual(result["Ring Incorporated Chalcogen Count"], 1)
        self.assertEqual(result["Heteroaromatic Ring Count"], 1)
        self.assertTrue(result["InChI"])
        self.assertTrue(result["InChIKey"])

    def test_selenium_classification(self):
        result, error = calculate_single_molecule_descriptors(
            "SE_001", "c1cc[se]c1"
        )
        self.assertIsNone(error)
        self.assertEqual(result["Chalcogen Type"], "Se")
        self.assertEqual(result["Selenium Count"], 1)
        self.assertEqual(result["Total Chalcogen Count"], 1)

    def test_tellurium_classification(self):
        result, error = calculate_single_molecule_descriptors(
            "TE_001", "c1cc[te]c1"
        )
        self.assertIsNone(error)
        self.assertEqual(result["Chalcogen Type"], "Te")
        self.assertEqual(result["Tellurium Count"], 1)
        self.assertEqual(result["Total Chalcogen Count"], 1)

    def test_oxygen_is_not_counted_as_target_chalcogen(self):
        result, error = calculate_single_molecule_descriptors(
            "O_001", "c1ccoc1"
        )
        self.assertIsNone(error)
        self.assertEqual(result["Chalcogen Type"], "None")
        self.assertEqual(result["Total Chalcogen Count"], 0)
        self.assertEqual(result["Oxygen Count"], 1)

    def test_mixed_chalcogen_detection(self):
        result, error = calculate_single_molecule_descriptors(
            "MIX_001", "C[Te]c1ccsc1"
        )
        self.assertIsNone(error)
        self.assertEqual(result["Chalcogen Type"], "Mixed")
        self.assertEqual(result["Mixed Chalcogen Flag"], 1)
        self.assertEqual(result["Total Chalcogen Count"], 2)

    def test_invalid_smiles(self):
        result, error = calculate_single_molecule_descriptors(
            "BAD_001", "not_a_smiles"
        )
        self.assertIsNone(result)
        self.assertEqual(error["Status"], "Invalid SMILES")

    def test_duplicate_detection_uses_canonical_smiles(self):
        df = pd.DataFrame({
            "Molecule_ID": ["A", "B", "C"],
            "SMILES": ["c1ccsc1", "c1sccc1", "c1cc[se]c1"],
        })
        valid_df, invalid_df, summary = process_molecular_dataset(df)
        self.assertTrue(invalid_df.empty)
        self.assertEqual(summary["Duplicate Molecules"], 2)
        duplicate_rows = valid_df[valid_df["Duplicate Flag"]]
        self.assertEqual(len(duplicate_rows), 2)


if __name__ == "__main__":
    unittest.main()
