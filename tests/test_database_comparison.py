import unittest

import pandas as pd

from database_comparison import comparison_export, comparison_table
from release_metadata import DATABASE_VERSION, SCIENTIFIC_CORE_VERSION


class DatabaseComparisonTests(unittest.TestCase):
    def setUp(self):
        self.reference = pd.Series({
            "Record_ID": "EROL_0001",
            "Canonical_SMILES": "c1ccsc1",
            "InChIKey": "REFKEY",
            "Dataset_Owner": "Dataset A",
            "HOMO_eV": -5.50,
            "LUMO_eV": -2.50,
            "Eg_eV": 3.00,
            "Experimental_Eg_eV": "1.80; 1.50",
            "S_Count": 1,
            "Se_Count": 0,
            "Te_Count": 0,
            "Chalcogen_Type": "S",
            "Method": "DFT",
            "Basis_Set": "6-31G*",
            "DOI_or_Reference": "10.1000/ref",
        })
        self.candidate = pd.Series({
            "Record_ID": "EROL_0002",
            "Canonical_SMILES": "c1cc[se]c1",
            "InChIKey": "CANDKEY",
            "Dataset_Owner": "Dataset B",
            "HOMO_eV": -5.20,
            "LUMO_eV": -2.30,
            "Eg_eV": 2.90,
            "Experimental_Eg_eV": "1.46; 0.83",
            "S_Count": 0,
            "Se_Count": 1,
            "Te_Count": 0,
            "Chalcogen_Type": "Se",
            "Method": "DFT",
            "Basis_Set": "def2-SVP",
            "DOI_or_Reference": "10.1000/cand",
        })

    def test_numeric_property_deltas_are_signed_candidate_minus_reference(self):
        table = comparison_table(self.reference, self.candidate)
        homo = table.loc[table["Field"] == "HOMO (eV)"].iloc[0]
        eg = table.loc[table["Field"] == "Eg (eV)"].iloc[0]

        self.assertAlmostEqual(homo["Candidate - Reference"], 0.3, places=4)
        self.assertAlmostEqual(eg["Candidate - Reference"], -0.1, places=4)

    def test_multi_valued_experimental_measurements_are_preserved(self):
        table = comparison_table(self.reference, self.candidate)
        experimental = table.loc[
            table["Field"] == "Experimental Eg (eV)"
        ].iloc[0]

        self.assertEqual(experimental["Reference"], "1.80; 1.50")
        self.assertEqual(experimental["Candidate"], "1.46; 0.83")
        self.assertTrue(pd.isna(experimental["Candidate - Reference"]))

    def test_export_contains_identity_provenance_and_release_versions(self):
        export = comparison_export(
            reference_row=self.reference,
            candidate_row=self.candidate,
            reference_public_id="DEV_0001",
            candidate_public_id="DEV_0002",
            reference_column="DOI_or_Reference",
        )

        values = export.set_index("Field")
        self.assertEqual(
            values.loc["Database record ID", "Reference"],
            "DEV_0001",
        )
        self.assertEqual(
            values.loc["Reference / DOI", "Candidate"],
            "10.1000/cand",
        )
        self.assertEqual(
            str(values.loc["Database version", "Reference"]),
            DATABASE_VERSION,
        )
        self.assertEqual(
            values.loc["Scientific core version", "Candidate"],
            SCIENTIFIC_CORE_VERSION,
        )


if __name__ == "__main__":
    unittest.main()
