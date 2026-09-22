# ChalMolDB v0.9 Development Roadmap

ChalMolDB v0.9 builds on the stable v0.8 Scientific Core release. The goal is to improve speed, scientific usability, provenance visibility, and analysis depth without changing the validated scientific definitions introduced in v0.8.

## Release principles

- Keep v0.8.0 as the stable reference release.
- Preserve publication-aligned descriptor definitions unless a versioned scientific change is explicitly documented.
- Keep curated database records separate from user-uploaded analysis data.
- Prefer measurable improvements over feature volume.
- Maintain responsive behavior and clear error states.
- Add tests for every new data-processing rule.

## 1. Performance and scalability

- Profile database loading, filtering, descriptor generation, structure rendering, and similarity search.
- Cache stable curated-database transformations where safe.
- Reduce unnecessary dataframe copies and repeated RDKit parsing.
- Avoid rendering molecular structures until they are needed.
- Add performance benchmarks for representative datasets.
- Define practical row-count guidance for interactive analysis.

### Exit criteria
- Database browsing remains responsive on the full curated release.
- Repeated navigation does not trigger avoidable expensive recomputation.
- Analysis performance is measured and documented.

## 2. Descriptor experience

- Add an in-app descriptor dictionary linked to the versioned scientific definitions.
- Improve descriptor grouping and discoverability.
- Add concise scientific descriptions, units, and interpretation notes.
- Make chalcogen-aware and general descriptor layers easier to compare.
- Add export metadata that records descriptor version and calculation context.

### Exit criteria
- Every exposed descriptor has a visible definition and unit where applicable.
- UI terminology stays synchronized with the descriptor dictionary.

## 3. Similarity and molecular comparison

- Improve reference-molecule comparison layout.
- Add side-by-side molecular identity and descriptor-difference views.
- Add configurable fingerprint settings only where scientifically justified.
- Preserve the distinction between structural similarity and electronic-property equivalence.
- Add downloadable similarity result tables with method metadata.

### Exit criteria
- Similarity results clearly state method and threshold.
- Users can inspect why two records are similar or different without leaving the workflow.

## 4. Database quality and provenance

- Strengthen provenance visibility at record level.
- Improve DOI/reference discovery and filtering.
- Add clearer source/method/quality summaries.
- Add coverage diagnostics for missing electronic-property fields.
- Add database-quality checks for duplicate IDs, malformed structures, and inconsistent metadata.
- Document release-level database coverage.

### Exit criteria
- Provenance is inspectable for every curated record where source metadata exists.
- Quality checks run automatically before a database release.

## 5. Export and reproducibility

- Add version metadata to CSV/ZIP exports.
- Include a machine-readable run manifest in analysis packages.
- Record descriptor groups, fingerprint settings, timestamps, and software version.
- Improve export naming while preserving safe file handling.

### Exit criteria
- An exported analysis package contains enough metadata to reproduce the calculation settings.

## 6. Reliability and release process

- Keep full CI green on the development branch.
- Expand tests for validation, metadata handling, similarity, and export behavior.
- Add a v0.9 smoke-test checklist before release.
- Re-run desktop and mobile review before tagging v0.9.0.
- Keep release notes synchronized with implemented changes.

## Suggested implementation order

1. Performance profiling and baseline measurements
2. Descriptor dictionary UI
3. Provenance and DOI/reference improvements
4. Similarity comparison refinement
5. Reproducibility manifest and export metadata
6. Full regression, smoke, and mobile testing

## Non-goals for v0.9

- No silent imputation of missing scientific values.
- No automatic ingestion of user uploads into the curated database.
- No major redesign of the validated v0.8 visual identity.
- No opaque predictive scoring layer inside ChalMolDB without a separately documented scientific workflow.

## Branch

Development branch: `v0.9-development`

Stable baseline: `v0.8.0`
