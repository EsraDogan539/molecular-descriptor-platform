# ChalMolDB v0.11 Release Hardening

## Release theme
**Evidence-aware molecular exploration**

v0.11 keeps the v0.10 structure-first scientific core stable while improving database interpretation, filtering, reproducibility and publication-facing presentation.

## Publication-friendly statistics

The Statistics workspace now follows a release-snapshot-first structure:
1. database size and identity summary,
2. collection/chalcogen/Eg distributions,
3. field-level data availability,
4. compact quality-control details.

Repeated-structure groups are computed from standardized InChIKey identity. Reference coverage and parsed DOI coverage are reported separately.

## Advanced research filters

The Database workspace supports optional filters for:
- HOMO, LUMO and Eg ranges,
- minimum S, Se and Te counts,
- calculation method,
- basis set,
- curation status.

Numeric filters only act when explicitly enabled. Missing values are excluded from an actively filtered numeric field rather than imputed.

## Neutral public provenance

Internal source provenance remains in the curated data, while public-facing labels are neutral:
- Development collection
- External literature collection
- DEV_* public record IDs
- EXT_* public record IDs

Public record views, downloadable database output, citation-ready exports and comparison exports use the neutral public layer.

## Query replay

Schema-v1.0 query manifests can be uploaded and replayed. The replay workflow:
- validates the manifest,
- restores compatible text/structure queries,
- restores basic and advanced filters,
- restores similarity settings,
- constrains numeric values to current database bounds,
- warns when the saved and current release identities differ,
- re-runs the query against the currently loaded release.

Replay is a reproducibility aid, not a guarantee of identical result sets after database changes.

## Validation

v0.11 validation includes:
- unit tests for advanced filters,
- public-label sanitization tests,
- query replay parser tests,
- legacy v0.9/v0.10 smoke coverage,
- v0.11 release smoke checks,
- database release audit,
- performance baseline,
- live desktop/mobile acceptance.

## Scope

v0.11 does not add predictive ML for unseen molecules, generative design, automated molecule ranking or automated scientific conclusions. Those remain outside this publication release.
