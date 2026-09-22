# ChalMolDB v0.10 Release Hardening

## Release theme
**Structure-first research workflow**

v0.10 extends the v0.9 scientific core with curated-database structure search, record comparison, provenance-aware citation export and reproducible query manifests.

## New scientific workflow

### Structure-first database search
The Database workspace supports three explicit structure-search modes:
- Exact standardized molecular identity
- Graph substructure matching
- Morgan/Tanimoto similarity

Similarity settings are fixed and disclosed as Morgan radius 2, 2048 bits and Tanimoto similarity. A structural-similarity result is descriptive and is not presented as evidence of equivalent electronic behavior.

### Curated record comparison
Two filtered database records can be compared side by side. The comparison includes molecular structures, HOMO, LUMO, Eg, experimental Eg where available, chalcogen context, method and basis-set metadata. Numeric deltas are descriptive candidate-minus-reference values.

Multi-valued experimental measurements remain source-preserving strings and are not silently collapsed.

### Provenance and citation-ready export
Record detail pages preserve dataset/source metadata. Valid DOI values are normalized to https://doi.org links. Non-DOI references remain explicit text.

Citation-ready CSV export includes record identity, standardized structure, InChIKey, electronic properties, chalcogen counts, method/basis metadata, reference/DOI, curation status, Database version and Scientific Core version.

### Reproducible query manifest
Database search state can be exported as a machine-readable JSON manifest containing:
- raw and standardized structure query,
- exact/substructure/similarity mode,
- active text and database filters,
- Morgan/Tanimoto settings when applicable,
- result count,
- Database and Scientific Core versions,
- deterministic SHA-256 checksum of the query configuration,
- SHA-256 checksum of the ordered result-record set.

The manifest records the search state; it does not claim future database releases will return the same records.

## Validation
v0.10 validation includes:
- unit tests for structure search, record comparison, DOI/provenance export and query manifests,
- database release audit,
- legacy v0.9 release smoke checks,
- v0.10 release smoke checks,
- performance baseline with exact, substructure and structure-similarity timings,
- live Streamlit review.

## Scope
v0.10 does not add predictive ML for unseen molecules, generative design, automated molecule ranking or automated scientific conclusions.
