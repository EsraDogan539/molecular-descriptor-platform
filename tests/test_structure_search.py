import unittest

import pandas as pd

from structure_search import (
    MORGAN_N_BITS,
    MORGAN_RADIUS,
    SIMILARITY_METRIC,
    standardize_query_smiles,
    structure_search,
)


class StructureSearchTests(unittest.TestCase):
    def setUp(self):
        self.df = pd.DataFrame({
            "Record_ID": ["A", "B", "C", "D"],
            "Canonical_SMILES": [
                "c1ccsc1",
                "c1sccc1",
                "c1cc[se]c1",
                "CCO",
            ],
            "InChIKey": ["", "", "", ""],
            "Eg_eV": [2.4, 2.4, 2.1, 5.0],
        })

    def test_query_smiles_is_canonicalized(self):
        query = standardize_query_smiles("c1sccc1")
        self.assertEqual(query["canonical_smiles"], "c1ccsc1")
        self.assertTrue(query["inchikey"])

    def test_exact_search_matches_equivalent_smiles(self):
        result, query = structure_search(
            self.df,
            query_smiles="c1sccc1",
            mode="Exact",
        )
        self.assertEqual(set(result["Record_ID"]), {"A", "B"})
        self.assertTrue((result["Structure Search Mode"] == "Exact").all())
        self.assertEqual(query["canonical_smiles"], "c1ccsc1")

    def test_substructure_search_uses_molecular_graph(self):
        result, _ = structure_search(
            self.df,
            query_smiles="c1ccsc1",
            mode="Substructure",
        )
        self.assertEqual(set(result["Record_ID"]), {"A", "B"})

    def test_similarity_search_ranks_identical_structures_first(self):
        result, _ = structure_search(
            self.df,
            query_smiles="c1ccsc1",
            mode="Similarity",
            minimum_similarity=0.0,
            top_n=4,
        )
        self.assertEqual(result.iloc[0]["Record_ID"], "A")
        self.assertAlmostEqual(
            float(result.iloc[0]["Structure Similarity"]),
            1.0,
            places=4,
        )
        self.assertEqual(result.iloc[1]["Record_ID"], "B")
        self.assertAlmostEqual(
            float(result.iloc[1]["Structure Similarity"]),
            1.0,
            places=4,
        )

    def test_similarity_settings_are_scientifically_explicit(self):
        self.assertEqual(MORGAN_RADIUS, 2)
        self.assertEqual(MORGAN_N_BITS, 2048)
        self.assertEqual(SIMILARITY_METRIC, "Tanimoto")

    def test_invalid_query_smiles_is_rejected(self):
        with self.assertRaisesRegex(ValueError, "not a valid SMILES"):
            structure_search(
                self.df,
                query_smiles="not_a_smiles",
                mode="Exact",
            )

    def test_empty_query_smiles_is_rejected(self):
        with self.assertRaisesRegex(ValueError, "Enter a SMILES query"):
            structure_search(
                self.df,
                query_smiles="",
                mode="Exact",
            )


if __name__ == "__main__":
    unittest.main()
