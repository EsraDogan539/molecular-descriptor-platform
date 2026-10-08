# ChalMolDB v0.10.0 — Structure-First Research Workflow

ChalMolDB v0.10.0 extends the Scientific Core with a structure-first curated-database workflow, record-level comparison, provenance-aware citation export and reproducible database-query manifests.

## Highlights

- Exact standardized-structure search in the curated database
- RDKit graph substructure search
- Morgan/Tanimoto structural similarity search
- Explicit Morgan radius 2 / 2048-bit fingerprint settings
- Side-by-side curated record comparison
- HOMO, LUMO, Eg and experimental Eg comparison context
- Preservation of multi-valued experimental measurements
- Record-level provenance with DOI normalization
- Citation-ready record CSV export
- Reproducible database-query JSON manifest
- SHA-256 query-configuration checksum
- SHA-256 ordered result-record-set checksum
- Centralized release/version metadata
- v0.10-specific release smoke tests
- Structure-search performance benchmarking

## Scientific safeguards

- Structural similarity is not presented as evidence of equivalent electronic properties.
- Missing structures, properties and provenance remain explicit.
- Non-DOI references are preserved as source text rather than converted or inferred.
- Multi-valued experimental measurements are not silently collapsed.
- Record comparison is descriptive; candidate-minus-reference deltas are not causal claims.
- Query manifests record the search state and do not imply that future database releases will return the same records.

## Reproducibility

Database query manifests record:
- input and standardized SMILES,
- Exact / Substructure / Similarity mode,
- active database filters,
- Morgan/Tanimoto settings where applicable,
- result count,
- Database and Scientific Core versions,
- deterministic query checksum,
- ordered result-set checksum.

Record-level citation-ready exports include molecular identity, electronic properties, chalcogen context, method/basis metadata, source/reference information and release versions.

## Validation

Release validation includes:
- full pytest suite,
- curated-database quality audit,
- v0.9 regression smoke checks,
- v0.10 structure-search smoke checks,
- v0.10 comparison/provenance smoke checks,
- v0.10 query-manifest smoke checks,
- CI performance baseline,
- desktop and 400 px live UI review before tagging.

Release: Database v1 · Scientific Core v0.10
