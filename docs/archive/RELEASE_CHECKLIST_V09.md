# ChalMolDB v0.9 Release Checklist

Use this checklist before tagging v0.9.0.

## Scientific core

- [x] Full pytest suite passes.
- [x] Descriptor Dictionary v0.9.0 covers all exposed scientific outputs.
- [x] Combined descriptor/fingerprint pipeline matches validated outputs.
- [x] Similarity ranking and thresholds pass regression tests.
- [x] Similarity exports record fingerprint method, radius, bit count and metric.

## Database quality

- [x] Database release audit passes.
- [x] Missing Record IDs = 0.
- [x] Duplicate Record IDs = 0.
- [x] Populated Canonical SMILES parse successfully.
- [x] SMILES/InChIKey mismatches are reviewed.
- [x] Numeric-property audit contains no unexplained non-numeric values.
- [x] Multi-valued experimental measurements are preserved explicitly and not silently collapsed.

## Reproducibility

- [x] Complete ZIP contains descriptor dictionary.
- [x] Complete ZIP contains run_manifest.json.
- [x] Manifest includes input SHA-256 checksum.
- [x] Manifest records software versions.
- [x] Manifest records fingerprint settings.
- [x] Manifest output inventory matches the final ZIP package.
- [x] Curated database export is reflected in the refreshed manifest.

## Performance

- [x] CI performance baseline runs successfully.
- [x] Descriptor/fingerprint pipeline remains within the documented v0.9 performance range.
- [x] No regression reintroduces duplicate RDKit parsing in the main analysis path.

## Automated smoke tests

- [x] Valid S/Se/Te/Mixed dataset passes.
- [x] Missing SMILES-column validation passes.
- [x] Mixed valid/invalid SMILES processing passes.
- [x] Complete export package smoke check passes.
- [x] Similarity search returns expected result count.

## UI review

- [x] Home page desktop review.
- [x] Database search/filter review.
- [x] Statistics coverage table review.
- [x] Analyze upload/settings review.
- [x] Descriptor Dictionary review.
- [x] Similarity pairwise comparison review.
- [x] Export/package review.
- [x] Documentation review.
- [x] Mobile review at approximately 400 px.
- [x] No page-level horizontal overflow.

## Release preparation

- [x] Update release notes.
- [x] Update roadmap completion state.
- [x] Confirm CI green on the final release commit.
- [x] Confirm Streamlit deployment uses the intended release branch/commit.
- [x] Create v0.9.0 tag only after the final UI smoke review.
