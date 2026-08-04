
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
