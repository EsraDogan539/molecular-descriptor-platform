
# Molecular Descriptor Platform

A Streamlit-based molecular descriptor and similarity analysis platform developed for cheminformatics and materials informatics studies.

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

The publication-oriented `v0.8-scientific-core` branch uses a fixed S/Se/Te-focused descriptor vocabulary.

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
