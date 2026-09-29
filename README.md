
# ChalMolDB — Chalcogen Molecular Database

A curated chalcogen-focused molecular database and Streamlit research interface for standardized molecular identity, electronic-property data, interpretable S/Se/Te descriptors, structure-first exploration, provenance-aware inspection, and reproducible query workflows.

**Current archived release:** Database v1 · Scientific Core v0.11.0  
**Version DOI:** https://doi.org/10.5281/zenodo.22903550  
**Concept DOI (all versions):** https://doi.org/10.5281/zenodo.22903549  
**GitHub release:** https://github.com/EsraDogan539/molecular-descriptor-platform/releases/tag/v0.11.0

## Citation

Doğan, E. N. (2026). *ChalMolDB v0.11.0 — Chalcogen Molecular Database* (Version 0.11.0) [Computer software]. Zenodo. https://doi.org/10.5281/zenodo.22903550

For references to the evolving ChalMolDB software record across versions, use the concept DOI: https://doi.org/10.5281/zenodo.22903549

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

The publication-oriented scientific core uses a fixed S/Se/Te-focused descriptor vocabulary. The archived publication release is `v0.11.0`.

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
