# Programmatic access to ChalMolDB

ChalMolDB is distributed as static, versioned files. No API key or account is needed. The files can be read
straight from the GitHub repository or from the Zenodo archive.

| File | Content |
|---|---|
| `data/public/chalmoldb_records.csv` | All records of the current release, with public `DEV_`/`EXT_` identifiers (UTF-8 CSV) |
| `data/structures_3d_v1.sdf.gz` | 3D structures of the development collection (B3LYP/6-311+G(d); CC BY 4.0, ref. Haciefendioglu & Yildirim 2025) |
| `data/structures_3d_external.sdf.gz` | 3D structures of the external oligomers (B3LYP/LANL2DZ; CC BY 4.0, H. Kayı group) |
| `docs/chalcogen_descriptor_dictionary_v1.csv` | Definitions of the chalcogen-aware descriptors |
| `data/public/MANIFEST.json` | Release version, field list and SHA-256 checksum of every file above |

The field definitions are given in Table S1 of the Supporting Information of the ChalMolDB article and in
`docs/DATABASE_CHANGELOG.md`. Each SDF entry carries the field `ChalMolDB_Record_ID`, which links it to the
`Record_ID` column of the CSV.

## Base URLs

- Latest release: `https://raw.githubusercontent.com/EsraDogan539/molecular-descriptor-platform/main/`
- A fixed release (recommended for reproducible work): replace `main` with the release tag, or use the
  Zenodo archive (concept DOI for all versions: https://doi.org/10.5281/zenodo.22903549).

## Python

```python
import gzip, io, json, urllib.request
import pandas as pd

BASE = "https://raw.githubusercontent.com/EsraDogan539/molecular-descriptor-platform/main/"

records = pd.read_csv(BASE + "data/public/chalmoldb_records.csv", low_memory=False)
manifest = json.load(urllib.request.urlopen(BASE + "data/public/MANIFEST.json"))
print(manifest["database_version"], len(records))

# Tellurium-containing records with an exact structure
te = records[(records["Te_Count"] > 0) & records["InChIKey"].notna()]

# 3D structures (RDKit), keyed by public record ID
from rdkit import Chem
raw = gzip.decompress(urllib.request.urlopen(BASE + "data/structures_3d_external.sdf.gz").read())
mols = {m.GetProp("ChalMolDB_Record_ID"): m
        for m in Chem.ForwardSDMolSupplier(io.BytesIO(raw), removeHs=False) if m is not None}
```

## R

```r
records <- read.csv("https://raw.githubusercontent.com/EsraDogan539/molecular-descriptor-platform/main/data/public/chalmoldb_records.csv")
```

## Links to single records

Every record has a permanent link in the web platform:
`https://chalmoldb.streamlit.app/?page=database&record=DEV_0001` (or `EXT_0148`, ...).

## Checking a download

```python
import hashlib
digest = hashlib.sha256(open("chalmoldb_records.csv", "rb").read()).hexdigest()
assert digest == manifest["files"]["data/public/chalmoldb_records.csv"]["sha256"]
```
