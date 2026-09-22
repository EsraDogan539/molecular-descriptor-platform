import unittest

import pandas as pd

from provenance_export import (
    citation_ready_record,
    doi_url,
    extract_doi,
)
from release_metadata import DATABASE_VERSION, SCIENTIFIC_CORE_VERSION


class ProvenanceExportTests(unittest.TestCase):
    def test_extracts_plain_and_prefixed_doi(self):
        self.assertEqual(
            extract_doi("doi:10.1234/ABC.567"),
            "10.1234/ABC.567",
        )
        self.assertEqual(
            extract_doi("https://doi.org/10.5555/test-value"),
            "10.5555/test-value",
        )

    def test_returns_none_for_non_doi_reference(self):
        self.assertIsNone(
            extract_doi(
                "Band Gap and Reorganization Energy Prediction of Conducting Polymers"
            )
        )

    def test_builds_https_doi_url(self):
        self.assertEqual(
            doi_url("10.1000/example"),
            "https://doi.org/10.1000/example",
        )

    def test_citation_ready_export_preserves_record_and_versions(self):
        row = pd.Series({
            "Record_ID": "EROL_0001",
            "Dataset_Owner": "Example dataset",
            "Split_Role": "Development/Training",
            "Scope_Flag": "Core_SSeTe",
            "Molecule_Name": "Molecule A",
            "System_Code": "A1",
            "Canonical_SMILES": "c1ccsc1",
            "InChIKey": "EXAMPLEKEY",
            "Chalcogen_Type": "S",
            "S_Count": 1,
            "Se_Count": 0,
            "Te_Count": 0,
            "HOMO_eV": -5.5,
            "LUMO_eV": -2.5,
            "Eg_eV": 3.0,
            "Experimental_Eg_eV": "1.80; 1.50",
            "Method": "DFT",
            "Basis_Set": "6-31G*",
            "DOI_or_Reference": "10.1000/example",
            "Curation_Status": "Ready",
        })

        export = citation_ready_record(
            row=row,
            public_record_id="DEV_0001",
        ).set_index("Field")["Value"]

        self.assertEqual(export["Database record ID"], "DEV_0001")
        self.assertEqual(str(export["Database version"]), DATABASE_VERSION)
        self.assertEqual(
            export["Scientific core version"],
            SCIENTIFIC_CORE_VERSION,
        )
        self.assertEqual(export["Canonical SMILES"], "c1ccsc1")
        self.assertEqual(export["Collection source"], "Development collection")
        self.assertNotIn("dataset", str(export["Collection source"]).lower())
        self.assertEqual(export["Experimental Eg (eV)"], "1.80; 1.50")
        self.assertEqual(export["DOI"], "10.1000/example")
        self.assertEqual(
            export["DOI URL"],
            "https://doi.org/10.1000/example",
        )


if __name__ == "__main__":
    unittest.main()
