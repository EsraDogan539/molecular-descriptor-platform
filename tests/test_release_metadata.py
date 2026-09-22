import unittest

import pandas as pd

from release_metadata import (
    DATABASE_VERSION,
    DESCRIPTOR_DICTIONARY_VERSION,
    RELEASE_LABEL,
    SCIENTIFIC_CORE_DISPLAY_VERSION,
    SCIENTIFIC_CORE_VERSION,
)
from run_manifest import build_run_manifest


class ReleaseMetadataTests(unittest.TestCase):
    def test_release_label_is_composed_from_central_versions(self):
        self.assertIn(f"Database v{DATABASE_VERSION}", RELEASE_LABEL)
        self.assertIn(
            f"Scientific Core {SCIENTIFIC_CORE_DISPLAY_VERSION}",
            RELEASE_LABEL,
        )

    def test_run_manifest_uses_central_release_versions(self):
        df = pd.DataFrame({
            "Molecule_ID": ["S_001"],
            "SMILES": ["c1ccsc1"],
        })
        manifest = build_run_manifest(
            input_df=df,
            project_name="version_test",
            summary={
                "Total Records": 1,
                "Valid Molecules": 1,
                "Invalid Molecules": 0,
                "Duplicate Molecules": 0,
                "Success Rate (%)": 100.0,
            },
            descriptor_dictionary_version=DESCRIPTOR_DICTIONARY_VERSION,
            output_files=[],
            fingerprint_settings={},
        )

        self.assertEqual(manifest["database_version"], DATABASE_VERSION)
        self.assertEqual(
            manifest["scientific_core_version"],
            SCIENTIFIC_CORE_VERSION,
        )
        self.assertEqual(
            manifest["descriptor_dictionary_version"],
            DESCRIPTOR_DICTIONARY_VERSION,
        )


if __name__ == "__main__":
    unittest.main()
