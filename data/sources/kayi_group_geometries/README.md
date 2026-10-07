# Kayı group geometry files (external collection)

Gaussian input files with the geometries of the donor–acceptor–donor (D–A–D) systems of the
external collection, provided by the authors (H. Kayı group, Department of Chemical Engineering, Ankara
University) and released with ChalMolDB with their permission. Files are included unmodified; only the folder
names were changed (spaces replaced by underscores).

| Folder | Systems | Family in ChalMolDB |
|---|---|---|
| `A-DAD_Monomer`, `A-DAD_Hexamer` | 16 + 16 | benzochalcogenadiazole (Ozkilinc & Kayi 2019) — repeat unit of the polymer records |
| `B-DAD_Monomer`, `B-DAD_Hexamer` | 16 + 16 | B_system / benzochalcogenadiazole-like DAD (3,4-ethylenedioxy-chalcogenophene donors) |
| `C-DAD_Monomer`, `C-DAD_Hexamer` | 16 + 16 | chalcogendiazoloquinoxaline (Kayi, Sen & Ozkilinc 2024) |

Each system code (e.g. `SeTeSe`) gives the donor, acceptor and donor chalcogen. The files are TD-DFT
(B3LYP/LANL2DZ) inputs. According to the authors, geometries were optimized at B3LYP/LANL2DZ; for the quinoxaline
family this workflow (B3LYP/LANL2DZ geometries, LC-BLYP/LANL2DZ electronic properties) is stated in Kayi, Sen &
Ozkilinc, J. Mol. Model. 2024, 30, 179.

Structures are derived by `scripts/build_external_structures.py` (results in `data/build/external_structures.csv`).
A file's chemical identity is accepted when the intended family structure has the same full atom/bond graph as
the coordinates or as the file's Gaussian connectivity block. Coordinates are released only when, in addition,
the formula matches and every ring chalcogen bond length is plausible for that element.

## Findings of the automated check (to be confirmed with the authors)

| Set | Identity confirmed | Coordinates released | Issue |
|---|---|---|---|
| C monomers, C hexamers | 32/32 | 32/32 | none |
| A hexamers | 16/16 | 15/16 | `A-OSO6`: S–N 1.897 Å, identical to the Se–N length of `A-OSeO6` (relabelled Se geometry) |
| A monomers | 16/16 | 3/16 | ring bond lengths do not change with the element (e.g. X–N 1.767 Å for O, S, Se and Te): geometries were not re-optimized after the chalcogen was changed |
| B monomers | 16/16 | 1/16 (`OOO`) | all other files keep the `OOO` bond lengths (X–C 1.40 Å, X–N 1.43 Å for S, Se and Te) |
| B hexamers | 0/16 | 0/16 | every file has 2 H fewer than the intended structure (C108H60 instead of C108H62); one benzo ring is missing both H atoms and has a 1.265 Å C–C bond |
