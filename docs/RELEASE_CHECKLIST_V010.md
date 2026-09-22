# ChalMolDB v0.10 Release Checklist

Use this checklist before tagging v0.10.0.

## Release identity
- [x] Central release metadata module is used by UI and manifests.
- [x] Database version and Scientific Core version are exported.
- [x] Stale hard-coded v0.9 UI labels are removed.

## Structure-first search
- [x] Exact standardized-structure search passes tests.
- [x] Substructure search passes tests.
- [x] Morgan/Tanimoto similarity search passes tests.
- [x] Invalid and empty SMILES validation passes.
- [x] Similarity UI states that structural similarity does not imply electronic-property equivalence.
- [x] Structure-search live desktop review completed.

## Curated record comparison
- [x] Two curated records can be compared side by side.
- [x] HOMO/LUMO/Eg and experimental Eg are displayed.
- [x] Multi-valued experimental measurements are preserved.
- [x] Comparison export includes release/provenance metadata.
- [x] Explicit comparison-submit flow prevents stale record rendering.
- [x] Live comparison mapping review completed.

## Provenance and citation workflow
- [x] DOI extraction/normalization tests pass.
- [x] DOI links use https://doi.org when a valid DOI is present.
- [x] Non-DOI references remain explicit text.
- [x] Citation-ready record CSV includes identity, properties, provenance and release versions.
- [x] Citation-ready record export live review completed.

## Reproducible database queries
- [x] Query manifest records structure input and standardized structure.
- [x] Query manifest records search mode and active filters.
- [x] Similarity query manifest records Morgan radius, bit count, Tanimoto metric, threshold and result limit.
- [x] Query configuration SHA-256 is deterministic.
- [x] Result-record-set SHA-256 is emitted.
- [x] Query manifest live export review completed.

## Performance
- [x] Existing v0.9 performance baseline remains in CI.
- [x] Exact structure search is benchmarked.
- [x] Substructure search is benchmarked.
- [x] Database structure-similarity search is benchmarked.
- [x] Record final v0.10 CI benchmark timings in documentation.

## Automated validation
- [x] Full pytest suite passes on v0.10 development branch.
- [x] Database release audit remains enabled.
- [x] v0.9 regression smoke checks remain enabled.
- [x] v0.10 structure-search smoke check is enabled.
- [x] v0.10 comparison/provenance smoke check is enabled.
- [x] v0.10 query-manifest smoke check is enabled.
- [ ] Final release-candidate CI is green.

## UI review
- [x] Database structure search desktop review.
- [x] Curated record comparison desktop review.
- [x] Citation-ready export desktop review.
- [x] Query manifest desktop review.
- [ ] Documentation desktop review.
- [ ] 400 px mobile review for Database structure search.
- [ ] 400 px mobile review for record comparison.
- [ ] 400 px mobile review for query reproducibility.
- [ ] No page-level horizontal overflow.

## Release preparation
- [ ] Update v0.10 roadmap completion state.
- [ ] Prepare v0.10.0 release notes.
- [ ] Confirm Streamlit deployment uses final v0.10 release commit.
- [ ] Tag v0.10.0 only after final UI/mobile review.
