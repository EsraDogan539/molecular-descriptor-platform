import unittest
import pandas as pd

from descriptor_engine import process_molecular_dataset
from database_metadata import build_curated_database


class DatabaseMetadataTests(unittest.TestCase):
    def test_metadata_aliases_are_normalized_and_preserved(self):
        input_df = pd.DataFrame({
            "Molecule_ID": ["S_001"],
            "SMILES": ["c1ccsc1"],
            "HOMO": [-5.2],
            "LUMO": [-2.8],
            "Eg": [2.4],
            "Source": ["DFT"],
            "DOI": ["10.0000/example"],
            "Method": ["B3LYP"],
        })

        valid_df, _, _ = process_molecular_dataset(input_df)
        database_df = build_curated_database(input_df, valid_df)

        self.assertAlmostEqual(database_df.loc[0, "HOMO_eV"], -5.2)
        self.assertAlmostEqual(database_df.loc[0, "LUMO_eV"], -2.8)
        self.assertAlmostEqual(database_df.loc[0, "Eg_eV"], 2.4)
        self.assertEqual(database_df.loc[0, "Property_Source"], "DFT")
        self.assertEqual(database_df.loc[0, "DOI_or_Reference"], "10.0000/example")
        self.assertEqual(database_df.loc[0, "Calculation_Method"], "B3LYP")
        self.assertEqual(database_df.loc[0, "Data_Quality_Flag"], "Review")

    def test_missing_optional_metadata_does_not_fail(self):
        input_df = pd.DataFrame({
            "Molecule_ID": ["SE_001"],
            "SMILES": ["c1cc[se]c1"],
        })

        valid_df, _, _ = process_molecular_dataset(input_df)
        database_df = build_curated_database(input_df, valid_df)

        self.assertEqual(database_df.loc[0, "Molecule_ID"], "SE_001")
        self.assertEqual(database_df.loc[0, "Data_Quality_Flag"], "Review")


if __name__ == "__main__":
    unittest.main()
