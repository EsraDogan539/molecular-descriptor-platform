# ChalMolDB v0.11 Release Checklist

Use this checklist before tagging v0.11.0.

## Release identity
- [x] v0.11 development identity is active on the development branch.
- [x] Set Scientific Core version to 0.11.0 only after final acceptance.
- [x] Confirm final release label is shown consistently in UI and manifests.

## Publication-friendly statistics
- [x] Release snapshot shows total records, unique standardized structures, core S/Se/Te records, Eg coverage, reference-backed records and repeated-structure groups.
- [x] Repeated structure groups are computed from standardized InChIKey identity.
- [x] Reference-backed and parsed-DOI counts are distinguished.
- [x] Distribution overview includes collection composition, chalcogen coverage and Eg distribution.
- [x] Data-availability table and quality-control summary are retained.
- [x] Desktop live review completed.

## Advanced database filtering
- [x] Optional HOMO/LUMO/Eg ranges are implemented.
- [x] Optional minimum S/Se/Te counts are implemented.
- [x] Method, basis-set and curation-status filters are implemented.
- [x] Active advanced filters are written to the query manifest.
- [x] Missing values do not silently satisfy active numeric filters.
- [x] Unit tests cover range, count, categorical and composed filtering.

## Public-facing provenance
- [x] Person-linked dataset names are removed from public record detail views.
- [x] Public record IDs map internal source prefixes to DEV/EXT identifiers.
- [x] Public CSV export uses neutral collection labels.
- [x] Citation-ready export uses neutral collection labels.
- [x] Comparison export uses neutral collection labels.
- [x] Automated tests reject person-linked labels in public dataframe output.

## Query replay
- [x] Query manifests can be uploaded and validated.
- [x] Text query, structure query and search mode can be restored.
- [x] Similarity threshold and result limit can be restored.
- [x] Core and advanced filters can be restored.
- [x] Restored numeric values are constrained to current release bounds.
- [x] Release differences produce an explicit notice.
- [x] Compatible v0.10 schema-1.0 manifests remain replayable.

## Documentation
- [x] Advanced filtering is documented.
- [x] Query replay is documented.
- [x] Public provenance wording is aligned with neutral collection labels.
- [ ] Add permanent DOI and recommended citation after archival release.

## Automated validation
- [x] v0.9 regression smoke checks remain enabled.
- [x] v0.10 release smoke checks remain enabled.
- [x] v0.11 release smoke checks are implemented.
- [x] Full pytest suite passes.
- [x] Syntax checks pass.
- [x] Database release audit passes.
- [x] Performance baseline passes.
- [x] Final release-candidate CI is green.

## UI acceptance
- [x] Statistics desktop review completed.
- [x] Advanced filters desktop review completed.
- [x] Query replay desktop review completed.
- [x] 400 px mobile review completed for Database filters and replay.
- [x] No page-level horizontal overflow.

## Release preparation
- [x] Freeze v0.11 scientific behavior.
- [x] Set release metadata to v0.11.0.
- [x] Prepare final v0.11.0 release notes.
- [ ] Confirm Streamlit deployment uses final release commit.
- [ ] Tag v0.11.0.
- [ ] Archive v0.11.0 and obtain DOI.
- [ ] Add DOI/citation metadata to platform and manuscript.
