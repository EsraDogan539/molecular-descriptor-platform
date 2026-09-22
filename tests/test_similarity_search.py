import unittest

import pandas as pd

from descriptor_engine import process_molecular_dataset
from similarity_search import (
    FINGERPRINT_METHOD,
    MORGAN_N_BITS,
    MORGAN_RADIUS,
    SIMILARITY_METRIC,
    find_similar_molecules,
    prepare_similarity_export,
)


class SimilaritySearchTests(unittest.TestCase):
    def setUp(self):
        raw = pd.DataFrame({
            "Molecule_ID": ["REF", "SAME", "SE", "TE"],
            "SMILES": [
                "c1ccsc1",
                "c1sccc1",
                "c1cc[se]c1",
                "c1cc[te]c1",
            ],
        })
        self.valid_df, _, _ = process_molecular_dataset(raw)

    def test_equivalent_structure_ranks_first(self):
        results = find_similar_molecules(
            self.valid_df,
            reference_id="REF",
            top_n=3,
            minimum_similarity=0.0,
        )
        self.assertFalse(results.empty)
        self.assertEqual(results.iloc[0]["Molecule_ID"], "SAME")
        self.assertAlmostEqual(
            float(results.iloc[0]["Morgan Similarity"]),
            1.0,
            places=4,
        )

    def test_minimum_similarity_filters_results(self):
        results = find_similar_molecules(
            self.valid_df,
            reference_id="REF",
            top_n=3,
            minimum_similarity=1.0,
        )
        self.assertEqual(results["Molecule_ID"].tolist(), ["SAME"])

    def test_similarity_export_contains_method_metadata(self):
        results = find_similar_molecules(
            self.valid_df,
            reference_id="REF",
            top_n=2,
            minimum_similarity=0.2,
        )
        exported = prepare_similarity_export(
            results,
            reference_id="REF",
            minimum_similarity=0.2,
        )
        self.assertTrue((exported["Reference Molecule_ID"] == "REF").all())
        self.assertTrue((exported["Fingerprint Method"] == FINGERPRINT_METHOD).all())
        self.assertTrue((exported["Morgan Radius"] == MORGAN_RADIUS).all())
        self.assertTrue((exported["Fingerprint Bits"] == MORGAN_N_BITS).all())
        self.assertTrue((exported["Similarity Metric"] == SIMILARITY_METRIC).all())
        self.assertTrue((exported["Minimum Similarity"] == 0.2).all())


if __name__ == "__main__":
    unittest.main()
