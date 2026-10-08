# ChalMolDB v0.8 Scientific Core — Release Checklist

## Functional smoke test
- [ ] Home, Database, Statistics, Analyze, Documentation and About load without exceptions.
- [ ] Navigation preserves the active page state.
- [ ] Curated database search and all four filters return expected records.
- [ ] Record inspection renders molecular structure, identity, electronic properties and provenance.
- [ ] Statistics charts render with the full database and with the preview fallback.
- [ ] Example CSV uploads successfully.
- [ ] Analyze rejects missing required columns with a clear message.
- [ ] Analyze rejects empty or duplicate Molecule_ID values with a clear message.
- [ ] Invalid SMILES are reported without stopping valid records from processing.
- [ ] Records, Descriptors, Structures, Similarity and Export tabs render after analysis.
- [ ] Processed CSV, descriptor CSV and complete ZIP downloads work.

## Scientific integrity
- [ ] S, Se and Te counts match representative hand-checked structures.
- [ ] Target Chalcogen Count excludes oxygen.
- [ ] Duplicate Flag uses InChIKey with canonical SMILES fallback.
- [ ] Missing values remain explicit and are not silently imputed.
- [ ] Record-level provenance is retained in curated exports.
- [ ] Similarity language states that structural similarity does not imply equivalent electronic behavior.

## Reliability and privacy
- [ ] Each analysis writes to an isolated temporary output directory.
- [ ] User project names cannot create arbitrary file-system paths.
- [ ] Uploaded user records do not enter the curated database.
- [ ] Public upload limit is appropriate for Streamlit Cloud resources.
- [ ] No secrets, credentials or private source files are committed.

## UI and responsive review
- [ ] Desktop review at approximately 1440 px width.
- [ ] Tablet review around 768–1024 px width.
- [ ] Mobile review around 390–430 px width.
- [ ] Header/navigation remains usable on mobile.
- [ ] Dataframes remain horizontally scrollable where necessary.
- [ ] Long SMILES/InChIKey values do not break layout.
- [ ] Loading, warning, error and empty states are readable.

## CI and release
- [ ] Full test suite passes on the release branch.
- [ ] Python syntax check passes for app and scientific-core modules.
- [ ] Streamlit deploy completes successfully.
- [ ] Live app smoke test completed after deploy.
- [ ] Release commit/tag recorded.
- [ ] Database version, Scientific Core version and documentation are synchronized.
