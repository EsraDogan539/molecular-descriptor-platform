
# ChalMolDB — Chalcogen Molecular Database

A curated chalcogen-focused molecular database and Streamlit research interface for standardized molecular identity, electronic-property data, interpretable S/Se/Te descriptors, structure-first exploration, provenance-aware inspection, and reproducible query workflows.

**Current release:** Database v1.1 · Scientific Core v0.11.0 (tag `db-v1.1`): 3,488 records, 3,173 unique standardized structures  
**Database v1.1 version DOI:** https://doi.org/10.5281/zenodo.23238884  
**Concept DOI (all versions, resolves to the latest):** https://doi.org/10.5281/zenodo.22903549  
**Database v1 (v0.11.0) version DOI:** https://doi.org/10.5281/zenodo.22903550  
**Web platform:** https://chalmoldb.streamlit.app

Changes between releases: [`docs/DATABASE_CHANGELOG.md`](docs/DATABASE_CHANGELOG.md).
Programmatic access (CSV, 3D SDF, checksums, Python/R examples): [`docs/DATA_ACCESS.md`](docs/DATA_ACCESS.md).
Embedding the platform in another website: [`docs/EMBEDDING.md`](docs/EMBEDDING.md).

## Citation

Doğan, E. N.; Kayı, H. (2026). *ChalMolDB Database v1.1 · Scientific Core v0.11.0 — Chalcogen Molecular Database* [Computer software]. Zenodo. https://doi.org/10.5281/zenodo.23238884

The concept DOI (10.5281/zenodo.22903549) covers all versions; cite the version DOI of the release you used.

## What ChalMolDB offers

**Curated database (Database v1.1)**
- 3,488 records: 3,088 donor–acceptor molecules with DFT-optimized structures and 400 oligomer and polymer
  records of two chalcogen-substituted donor–acceptor–donor families
- Standardized identity (RDKit canonical SMILES, InChI, InChIKey); repeated structures flagged, not removed
- Record-level provenance: source, level of theory, solvation, curation notes
- 3D structures for 3,147 records, with an interactive in-browser viewer and SDF download

**Web platform** (https://chalmoldb.streamlit.app)
- Text, metadata and S/Se/Te-content filters; HOMO/LUMO/Eg ranges
- Structure search by SMILES or by drawing (exact, substructure, Morgan/Tanimoto similarity)
- Side-by-side record comparison, statistics view, citation-ready record export
- Replayable query manifests with SHA-256 checksums
- Permanent record links: `https://chalmoldb.streamlit.app/?page=database&record=EXT_0162`

**Analyze your own molecules**
- Upload a CSV with `Molecule_ID` and `SMILES`; structures are validated and standardized
- General and chalcogen-aware descriptors, Morgan and MACCS fingerprints, similarity search
- Results exported as CSV/ZIP; uploads never enter the curated database

```csv
Molecule_ID,SMILES
MOL_001,c1ccsc1
MOL_002,c1cc[se]c1
```

## Data access and licences

Records, 3D structures and checksums can be read directly by scripts; see [`docs/DATA_ACCESS.md`](docs/DATA_ACCESS.md).
Database content: CC BY 4.0 ([`LICENSE-DATA.md`](LICENSE-DATA.md)). Software: MIT ([`LICENSE`](LICENSE)).

## Team

- **Esra Nur Doğan** (ORCID 0000-0001-5755-7596): database design, data curation, software development
- **Hakan Kayı** (ORCID 0000-0001-7300-0325): scientific supervision, source of the external DFT data;
  corresponding author (hkayi@ankara.edu.tr)

Department of Chemical Engineering, Faculty of Engineering, Ankara University, Ankara, Türkiye.

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

1. **Curated Database** — browse and filter the versioned publication database.
2. **Analyze Your Dataset** — validate user-supplied structures and calculate the same general and chalcogen-aware descriptor layers.

Uploaded user records do not enter the curated publication database automatically.
