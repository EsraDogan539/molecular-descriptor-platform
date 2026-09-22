# ChalMolDB v0.11.0 — Evidence-Aware Molecular Exploration

ChalMolDB v0.11.0 extends the structure-first v0.10 workflow with publication-friendly database statistics, research-oriented filtering, replayable query manifests and neutral public provenance labels.

## Highlights

- Publication-friendly Statistics page with release snapshot and coverage overview
- Explicit distinction between reference-backed records and parsed DOI records
- HOMO, LUMO and Eg range filtering
- Minimum S, Se and Te count filtering
- Method, basis-set and curation-status filtering
- Advanced filter state included in reproducible query manifests
- Query manifest upload and replay
- Release-difference notices for replayed queries
- Neutral public collection labels and DEV/EXT record identifiers
- Person-linked internal dataset labels excluded from public-facing exports
- v0.11-specific automated smoke checks

## Scientific safeguards

- Missing values do not satisfy an actively enabled numeric filter.
- Query replay re-runs the saved configuration against the current database release; it does not claim identical results across releases.
- Restored numeric filter values are constrained to the current release bounds.
- Structural similarity remains descriptive and is not treated as electronic-property equivalence.
- Public-facing labels are neutralized without deleting internal provenance from the curated source data.
- Reference-backed records are not presented as DOI-backed unless a DOI can actually be parsed.

## Reproducibility

v0.11 query manifests retain the v1.0 schema and record:
- text and structure queries,
- exact/substructure/similarity mode,
- basic and advanced filters,
- Morgan/Tanimoto settings where applicable,
- result count,
- Database and Scientific Core versions,
- query-configuration SHA-256,
- ordered result-record-set SHA-256.

Saved manifests can be replayed through the Database workspace after schema validation.

## Validation

Release validation includes:
- full pytest suite,
- curated-database quality audit,
- v0.9 regression smoke checks,
- v0.10 release smoke checks,
- v0.11 advanced-filter/public-label/query-replay smoke checks,
- CI performance baseline,
- desktop and 400 px live UI review before tagging.

Release target: Database v1 · Scientific Core v0.11.0
