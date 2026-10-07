# Kayı group geometry files (external collection)

Gaussian input files containing the optimized geometries of the donor–acceptor–donor (D–A–D) systems of the
external collection, provided by the authors (H. Kayı group, Department of Chemical Engineering, Ankara
University) and released with ChalMolDB with their permission. Files are included unmodified; only the folder
names were changed (spaces replaced by underscores).

| Folder | Systems | Family in ChalMolDB |
|---|---|---|
| `A-DAD_Monomer`, `A-DAD_Hexamer` | 16 + 16 | benzochalcogenadiazole (Ozkilinc & Kayi 2019) — repeat unit of the polymer records |
| `B-DAD_Monomer`, `B-DAD_Hexamer` | 16 + 16 | B_system / benzochalcogenadiazole-like DAD (3,4-ethylenedioxy-chalcogenophene donors) |
| `C-DAD_Monomer`, `C-DAD_Hexamer` | 16 + 16 | chalcogendiazoloquinoxaline (Kayi, Sen & Ozkilinc 2024) |

Each system code (e.g. `SeTeSe`) gives the donor, acceptor and donor chalcogen. The files are TD-DFT
(B3LYP/LANL2DZ) inputs; the coordinates are the B3LYP/LANL2DZ-optimized geometries. For the quinoxaline family
this workflow (B3LYP/LANL2DZ geometries, LC-BLYP/LANL2DZ electronic properties) is stated in Kayi, Sen &
Ozkilinc, J. Mol. Model. 2024, 30, 179.

Structures are derived by `scripts/build_external_structures.py`, which accepts a structure only if its full
atom/bond graph matches the graph perceived from these coordinates.
