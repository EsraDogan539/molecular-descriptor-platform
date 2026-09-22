# ChalMolDB v0.9 Release Checklist

Use this checklist before tagging v0.9.0.

## Scientific core

- [ ] Full pytest suite passes.
- [ ] Descriptor Dictionary v0.9.0 covers all exposed scientific outputs.
- [ ] Combined descriptor/fingerprint pipeline matches validated outputs.
- [ ] Similarity ranking and thresholds pass regression tests.
- [ ] Similarity exports record fingerprint method, radius, bit count and metric.

## Database quality

- [ ] Database release audit passes.
- [ ] Missing Record IDs = 0.
- [ ] Duplicate Record IDs = 0.
- [ ] Populated Canonical SMILES parse successfully.
- [ ] SMILES/InChIKey mismatches are reviewed.
- [ ] Numeric-property audit contains no unexplained non-numeric values.
- [ ] Multi-valued experimental measurements are preserved explicitly and not silently collapsed.

## Reproducibility

- [ ] Complete ZIP contains descriptor dictionary.
- [ ] Complete ZIP contains run_manifest.json.
- [ ] Manifest includes input SHA-256 checksum.
- [ ] Manifest records software versions.
- [ ] Manifest records fingerprint settings.
- [ ] Manifest output inventory matches the final ZIP package.
- [ ] Curated database export is reflected in the refreshed manifest.

## Performance

- [ ] CI performance baseline runs successfully.
- [ ] Descriptor/fingerprint pipeline remains within the documented v0.9 performance range.
- [ ] No regression reintroduces duplicate RDKit parsing in the main analysis path.

## Automated smoke tests

- [ ] Valid S/Se/Te/Mixed dataset passes.
- [ ] Missing SMILES-column validation passes.
- [ ] Mixed valid/invalid SMILES processing passes.
- [ ] Complete export package smoke check passes.
- [ ] Similarity search returns expected result count.

## UI review

- [ ] Home page desktop review.
- [ ] Database search/filter review.
- [ ] Statistics coverage table review.
- [ ] Analyze upload/settings review.
- [ ] Descriptor Dictionary review.
- [ ] Similarity pairwise comparison review.
- [ ] Export/package review.
- [ ] Documentation review.
- [ ] Mobile review at approximately 400 px.
- [ ] No page-level horizontal overflow.

## Release preparation

- [ ] Update release notes.
- [ ] Update roadmap completion state.
- [ ] Confirm CI green on the final release commit.
- [ ] Confirm Streamlit deployment uses the intended release branch/commit.
- [ ] Create v0.9.0 tag only after the final UI smoke review.
