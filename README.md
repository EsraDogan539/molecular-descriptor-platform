
# ChalMolDB — Chalcogen Molecular Database

A curated chalcogen-focused molecular database and Streamlit research interface for standardized molecular identity, electronic-property data, interpretable S/Se/Te descriptors, structure-first exploration, provenance-aware inspection, and reproducible query workflows.

**Current release:** Database v1.1 · Scientific Core v0.11.0 (tag `db-v1.1`): 3,248 records, 3,079 unique standardized structures  
**Concept DOI (all versions, resolves to the latest):** https://doi.org/10.5281/zenodo.22903549  
**Database v1 (v0.11.0) version DOI:** https://doi.org/10.5281/zenodo.22903550  
**Web platform:** https://chalmoldb.streamlit.app

Changes between releases: [`docs/DATABASE_CHANGELOG.md`](docs/DATABASE_CHANGELOG.md).
Programmatic access (CSV, 3D SDF, checksums, Python/R examples): [`docs/DATA_ACCESS.md`](docs/DATA_ACCESS.md).

## Citation

Doğan, E. N. (2026). *ChalMolDB Database v1.1 · Scientific Core v0.11.0 — Chalcogen Molecular Database* [Computer software]. Zenodo. https://doi.org/10.5281/zenodo.22903549

The concept DOI covers all versions; cite the version DOI of the release you used once it is listed on Zenodo.

## Features

- CSV upload with `Molecule_ID` and `SMILES` columns
- SMILES validation
- Molecular descriptor calculation
- Chalcogen-specific descriptors for S, Se and Te
- Morgan and MACCS fingerprints
- Molecular structure cards
- Dataset analysis panel
- Pairwise molecule comparison
- Similar molecule search with Morgan Tanimoto similarity
- CSV and ZIP output files

## Input Format

```csv
Molecule_ID,SMILES
MOL_001,c1ccccc1
MOL_002,CCO

## Scientific Core terminology

The publication-oriented scientific core uses a fixed S/Se/Te-focused descriptor vocabulary. The current release is Database v1.1 with Scientific Core v0.11.0.

- **Target Chalcogen Count** = `Sulfur Count + Selenium Count + Tellurium Count`
- **Target Chalcogen Fraction** = `Target Chalcogen Count / Heavy Atom Count`
- Oxygen is retained as a general elemental descriptor but is not included in the target S/Se/Te count.
- Repeated structures are identified by **InChIKey** when available, with canonical SMILES used as a fallback.
- Repeated structures are flagged rather than silently removed so that record-level provenance can be retained.

The exact implementation-aligned definitions are documented in
`docs/chalcogen_descriptor_dictionary_v1.csv`. These names are intended to remain synchronized with the platform UI and the manuscript descriptor table.

## Publication workflow

The platform separates two workflows:

1. **Our Curated Database** — browse and filter the versioned publication database.
2. **Analyze Your Dataset** — validate user-supplied structures and calculate the same general and chalcogen-aware descriptor layers.

Uploaded user records do not enter the curated publication database automatically.
