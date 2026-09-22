import unittest

from query_manifest import (
    build_query_manifest,
    result_set_sha256,
    sha256_json,
)


class QueryManifestTests(unittest.TestCase):
    def _build(self, **overrides):
        params = {
            "text_query": "",
            "structure_query": "c1ccsc1",
            "standardized_structure_query": "c1ccsc1",
            "structure_mode": "Similarity",
            "minimum_similarity": 0.4,
            "maximum_results": 50,
            "fingerprint_method": "Morgan",
            "fingerprint_radius": 2,
            "fingerprint_bits": 2048,
            "similarity_metric": "Tanimoto",
            "filters": {
                "collection": "All",
                "scope": "All",
                "chalcogen": "All",
                "eg": "All",
                "reference_or_doi": "All",
            },
            "result_record_ids": ["EROL_0001", "EROL_0002"],
        }
        params.update(overrides)
        return build_query_manifest(**params)

    def test_query_checksum_is_deterministic(self):
        first = self._build()
        second = self._build()
        self.assertEqual(
            first["query_checksum_sha256"],
            second["query_checksum_sha256"],
        )

    def test_result_set_checksum_is_deterministic_and_order_sensitive(self):
        self.assertEqual(
            result_set_sha256(["A", "B"]),
            result_set_sha256(["A", "B"]),
        )
        self.assertNotEqual(
            result_set_sha256(["A", "B"]),
            result_set_sha256(["B", "A"]),
        )

    def test_similarity_settings_are_recorded(self):
        manifest = self._build()
        fingerprint = manifest["query"]["fingerprint"]

        self.assertEqual(fingerprint["method"], "Morgan")
        self.assertEqual(fingerprint["radius"], 2)
        self.assertEqual(fingerprint["bits"], 2048)
        self.assertEqual(fingerprint["metric"], "Tanimoto")
        self.assertEqual(fingerprint["minimum_similarity"], 0.4)
        self.assertEqual(fingerprint["maximum_results"], 50)
        self.assertEqual(manifest["result_count"], 2)

    def test_non_similarity_query_omits_fingerprint_block(self):
        manifest = self._build(
            structure_mode="Exact",
        )
        self.assertIsNone(manifest["query"]["fingerprint"])

    def test_filter_change_changes_query_checksum(self):
        first = self._build()
        second = self._build(
            filters={
                "collection": "Development/Training",
                "scope": "All",
                "chalcogen": "All",
                "eg": "All",
                "reference_or_doi": "All",
            }
        )
        self.assertNotEqual(
            first["query_checksum_sha256"],
            second["query_checksum_sha256"],
        )

    def test_sha256_json_is_key_order_independent(self):
        self.assertEqual(
            sha256_json({"a": 1, "b": 2}),
            sha256_json({"b": 2, "a": 1}),
        )


if __name__ == "__main__":
    unittest.main()
