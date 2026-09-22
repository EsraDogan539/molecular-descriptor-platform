# ChalMolDB v0.10 Roadmap

## Release theme

**Structure-first research workflow**

ChalMolDB v0.10 will turn the v0.9 scientific core into a more direct research workflow: start from a molecular structure, search the curated database, compare related structures and properties, inspect provenance, and export a reproducible result set.

The release should remain focused, publication-oriented and conservative about scientific interpretation.

## 1. Version and release identity

- Centralize visible version strings in one module.
- Define explicit constants for Database version, Scientific Core version and Descriptor Dictionary version.
- Remove stale hard-coded version labels from UI components.
- Add version metadata to analysis exports and the run manifest.
- Add regression tests for release/version consistency.

## 2. Structure-first database search

- Accept a query SMILES directly in the Database workspace.
- Standardize the query before searching.
- Support:
  - exact standardized-structure match,
  - substructure search,
  - Morgan/Tanimoto similarity search against curated structures.
- Keep text/provenance filters available alongside structure search.
- Clearly distinguish exact, substructure and fingerprint-similarity results.
- Never imply that structural similarity establishes equivalent electronic behavior.

## 3. Curated-database molecular comparison

- Allow a selected database record to be compared directly with another curated record.
- Show side-by-side structures and a compact descriptor/property delta table.
- Include available HOMO, LUMO, Eg and experimental Eg context.
- Preserve multi-valued experimental measurements without collapsing them silently.
- Add CSV export for comparison results with method metadata.

## 4. Provenance and citation workflow

- Improve record-level provenance presentation.
- Render DOI/reference as clickable links when valid.
- Add a citation-ready record export containing:
  - Record ID,
  - standardized identity,
  - electronic properties,
  - source/method/basis metadata,
  - DOI/reference,
  - database/scientific-core version.
- Keep unavailable provenance explicit rather than inferred.

## 5. Reproducible query exports

- Add machine-readable search configuration to exported result packages.
- Record:
  - query structure,
  - standardized query identity,
  - search mode,
  - fingerprint settings,
  - thresholds,
  - active filters,
  - result count,
  - software/database versions.
- Add deterministic query checksum where practical.
- Ensure exported CSV/JSON inventories match the final ZIP package.

## 6. Scientific UX

- Keep the existing restrained academic visual language.
- Add a clear structure-search entry point without cluttering the Database page.
- Use progressive disclosure for advanced search controls.
- Keep desktop and 400 px mobile layouts usable.
- Make scientific caveats visible but concise.

## 7. Performance and scalability

- Cache standardized curated structures/fingerprints where safe.
- Avoid reparsing the complete database for every query.
- Benchmark exact, substructure and similarity search separately.
- Document expected performance for the current Database v1 size.
- Add regression thresholds only where stable enough to be meaningful.

## 8. Validation and release process

- Add unit tests for exact, substructure and similarity search.
- Add tests for invalid/empty query SMILES.
- Add provenance/citation export tests.
- Add reproducible query-manifest tests.
- Extend the release smoke workflow.
- Complete desktop and 400 px mobile live reviews before tagging v0.10.0.

## Scope guardrails

The following are intentionally out of scope for v0.10 unless required by the core workflow:

- predictive machine-learning models,
- automatic property prediction for unseen molecules,
- generative molecule design,
- user accounts or persistent cloud storage,
- large-scale external database ingestion,
- automated scientific conclusions or ranking of molecules.

## Definition of done

v0.10 is ready when a user can:

1. enter a molecular structure,
2. find exact/substructure/similar curated records,
3. inspect structure, electronic-property and provenance context,
4. compare selected records,
5. export the search/comparison with enough metadata to reproduce the workflow,
6. use the complete flow on desktop and mobile without layout failure.
