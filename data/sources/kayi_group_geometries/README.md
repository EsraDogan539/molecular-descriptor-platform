# Kayı group geometry files (external collection)

Gaussian input files with the geometries of the donor–acceptor–donor (D–A–D) systems of the
external collection, provided by the authors (H. Kayı group, Department of Chemical Engineering, Ankara
University) and released with ChalMolDB with their permission. Files are included unmodified; only the folder
names were changed (spaces replaced by underscores).

| Folder | Systems | Family in ChalMolDB |
|---|---|---|
| `A-DAD_Monomer`, `A-DAD_Hexamer` | 16 + 16 | benzochalcogenadiazole (Ozkilinc & Kayi 2019) — repeat unit of the polymer records |
| `C-DAD_Monomer`, `C-DAD_Hexamer` | 16 + 16 | chalcogendiazoloquinoxaline (Kayi, Sen & Ozkilinc 2024) |

Each system code (e.g. `SeTeSe`) gives the donor, acceptor and donor chalcogen. The files are TD-DFT
(B3LYP/LANL2DZ) inputs. According to the authors, geometries were optimized at B3LYP/LANL2DZ; for the quinoxaline
family this workflow (B3LYP/LANL2DZ geometries, LC-BLYP/LANL2DZ electronic properties) is stated in Kayi, Sen &
Ozkilinc, J. Mol. Model. 2024, 30, 179.

Structures are derived by `scripts/build_external_structures.py` (results in `data/build/external_structures.csv`).
A file's chemical identity is accepted when the intended family structure has the same full atom/bond graph as
the coordinates or as the file's Gaussian connectivity block. Coordinates are released only when, in addition,
the formula matches and every ring chalcogen bond length is plausible for that element.

## Results of the automated check

| Set | Identity confirmed | Coordinates released | Note |
|---|---|---|---|
| C monomers, C hexamers | 32/32 | 32/32 | used for the chalcogendiazoloquinoxaline records |
| A hexamers | 16/16 | 15/16 | `A-OSO6`: S–N 1.897 Å, identical to the Se–N length of `A-OSeO6` (relabelled Se geometry) |
| A monomers | 16/16 | 3/16 | ring bond lengths do not change with the element (e.g. X–N 1.767 Å for O, S, Se and Te) |

The A family appears in ChalMolDB only as polymer records, so the A files serve as structure references
for the repeat unit; no A coordinates are attached to records.

The geometry files of a third family (B system, 3,4-ethylenedioxychalcogenophene donors) were also checked.
The hexamer files contain two H atoms fewer than the intended structures, and 15 of 16 monomer files were not
re-optimized after the chalcogen was changed. The 112 B-system records were therefore removed in Database
v1.1, and their geometry files are not distributed.
