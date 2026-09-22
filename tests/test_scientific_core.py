import unittest
import zipfile
import pandas as pd

from descriptor_engine import (
    calculate_single_molecule_descriptors,
    process_molecular_dataset,
    process_molecular_dataset_with_fingerprints,
    create_fingerprint_dataset,
    sanitize_project_name,
    validate_input_dataframe,
    run_molecular_descriptor_platform,
)

from descriptor_dictionary import (
    DESCRIPTOR_DEFINITIONS,
    DESCRIPTOR_DICTIONARY_VERSION,
    descriptor_dictionary_dataframe,
)


class ScientificCoreTests(unittest.TestCase):
    def test_sulfur_classification_and_aromatic_environment(self):
        result, error = calculate_single_molecule_descriptors(
            "S_001", "c1ccsc1"
        )
        self.assertIsNone(error)
        self.assertEqual(result["Chalcogen Type"], "S")
        self.assertEqual(result["Target Chalcogen Count"], 1)
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
        self.assertEqual(result["Target Chalcogen Count"], 1)

    def test_tellurium_classification(self):
        result, error = calculate_single_molecule_descriptors(
            "TE_001", "c1cc[te]c1"
        )
        self.assertIsNone(error)
        self.assertEqual(result["Chalcogen Type"], "Te")
        self.assertEqual(result["Tellurium Count"], 1)
        self.assertEqual(result["Target Chalcogen Count"], 1)

    def test_oxygen_is_not_counted_as_target_chalcogen(self):
        result, error = calculate_single_molecule_descriptors(
            "O_001", "c1ccoc1"
        )
        self.assertIsNone(error)
        self.assertEqual(result["Chalcogen Type"], "None")
        self.assertEqual(result["Target Chalcogen Count"], 0)
        self.assertEqual(result["Oxygen Count"], 1)

    def test_mixed_chalcogen_detection(self):
        result, error = calculate_single_molecule_descriptors(
            "MIX_001", "C[Te]c1ccsc1"
        )
        self.assertIsNone(error)
        self.assertEqual(result["Chalcogen Type"], "Mixed")
        self.assertEqual(result["Mixed Chalcogen Flag"], 1)
        self.assertEqual(result["Target Chalcogen Count"], 2)

    def test_invalid_smiles(self):
        result, error = calculate_single_molecule_descriptors(
            "BAD_001", "not_a_smiles"
        )
        self.assertIsNone(result)
        self.assertEqual(error["Status"], "Invalid SMILES")

    def test_duplicate_detection_uses_standardized_structure_identity(self):
        df = pd.DataFrame({
            "Molecule_ID": ["A", "B", "C"],
            "SMILES": ["c1ccsc1", "c1sccc1", "c1cc[se]c1"],
        })
        valid_df, invalid_df, summary = process_molecular_dataset(df)
        self.assertTrue(invalid_df.empty)
        self.assertEqual(summary["Duplicate Molecules"], 2)
        duplicate_rows = valid_df[valid_df["Duplicate Flag"]]
        self.assertEqual(len(duplicate_rows), 2)

    def test_validation_rejects_duplicate_molecule_ids(self):
        df = pd.DataFrame({
            "Molecule_ID": ["A", "A"],
            "SMILES": ["c1ccsc1", "c1cc[se]c1"],
        })
        with self.assertRaisesRegex(ValueError, "must be unique"):
            validate_input_dataframe(df)

    def test_validation_rejects_missing_molecule_id(self):
        df = pd.DataFrame({
            "Molecule_ID": ["A", ""],
            "SMILES": ["c1ccsc1", "c1cc[se]c1"],
        })
        with self.assertRaisesRegex(ValueError, "missing or empty"):
            validate_input_dataframe(df)

    def test_project_name_is_sanitized_for_output_files(self):
        self.assertEqual(
            sanitize_project_name("../../my project:01"),
            "my_project_01",
        )
        self.assertEqual(sanitize_project_name(""), "chalcogen_project")

    def test_combined_pipeline_matches_separate_outputs(self):
        df = pd.DataFrame({
            "Molecule_ID": ["S1", "SE1", "BAD"],
            "SMILES": ["c1ccsc1", "c1cc[se]c1", "not_a_smiles"],
        })

        valid_old, invalid_old, summary_old = process_molecular_dataset(df)
        fp_old, _ = create_fingerprint_dataset(df)

        valid_new, invalid_new, fp_new, summary_new = (
            process_molecular_dataset_with_fingerprints(df)
        )

        pd.testing.assert_frame_equal(
            valid_old.reset_index(drop=True),
            valid_new.reset_index(drop=True),
            check_dtype=False,
        )
        pd.testing.assert_frame_equal(
            invalid_old.reset_index(drop=True),
            invalid_new.reset_index(drop=True),
            check_dtype=False,
        )
        pd.testing.assert_frame_equal(
            fp_old.reset_index(drop=True),
            fp_new.reset_index(drop=True),
            check_dtype=False,
        )
        self.assertEqual(summary_old, summary_new)

    def test_descriptor_dictionary_is_unique_and_versioned(self):
        dictionary_df = descriptor_dictionary_dataframe()
        self.assertFalse(dictionary_df.empty)
        self.assertEqual(DESCRIPTOR_DICTIONARY_VERSION, "0.9.0")
        self.assertEqual(
            len(dictionary_df["Descriptor"]),
            dictionary_df["Descriptor"].nunique(),
        )
        self.assertTrue(
            {"Descriptor", "Group", "Unit", "Definition", "Interpretation"}
            .issubset(dictionary_df.columns)
        )

    def test_descriptor_dictionary_covers_scientific_core_outputs(self):
        result, invalid = calculate_single_molecule_descriptors(
            "DICT_001",
            "c1ccsc1",
        )
        self.assertIsNone(invalid)
        documented = {item["Descriptor"] for item in DESCRIPTOR_DEFINITIONS}
        excluded = {
            "Molecule_ID", "Original SMILES", "Canonical SMILES", "InChI",
            "InChIKey", "Valid SMILES", "Molecular Formula", "Status",
            "Contains S", "Contains Se", "Contains Te",
        }
        scientific_outputs = set(result) - excluded
        self.assertTrue(scientific_outputs.issubset(documented))

    def test_analysis_export_records_descriptor_dictionary_version(self):
        df = pd.DataFrame({
            "Molecule_ID": ["S1"],
            "SMILES": ["c1ccsc1"],
        })
        results = run_molecular_descriptor_platform(
            df,
            project_name="dictionary_test",
        )

        self.assertEqual(
            results["descriptor_dictionary_version"],
            DESCRIPTOR_DICTIONARY_VERSION,
        )
        self.assertEqual(
            results["summary_df"].iloc[0]["Descriptor Dictionary Version"],
            DESCRIPTOR_DICTIONARY_VERSION,
        )

        with zipfile.ZipFile(results["zip_file"]) as archive:
            names = archive.namelist()
        self.assertTrue(
            any("descriptor_dictionary_v0.9.0.csv" in name for name in names)
        )


if __name__ == "__main__":
    unittest.main()
