# 3D structures of the development collection

`structures_3d_v1.sdf.gz` contains the DFT-optimized 3D structures (B3LYP/6-311+G(d), Jaguar) of the 3,088
development records of ChalMolDB Database v1, one SDF entry per record.

**Source:** T. Haciefendioglu, E. Yildirim, *Band Gap and Reorganization Energy Prediction of Conducting Polymers
by the Integration of Machine Learning and Density Functional Theory*, J. Chem. Inf. Model. 2025, 65, 5360–5369,
https://doi.org/10.1021/acs.jcim.5c00345 — Supporting Information file `ci5c00345_si_002.zip`
(obtained via Europe PMC, PMC12152970).

**Licence:** CC BY 4.0 (https://creativecommons.org/licenses/by/4.0/). Coordinates and original SDF fields are
unmodified. ChalMolDB adds four fields to each entry: `ChalMolDB_Record_ID` (public DEV_ identifier),
`ChalMolDB_InChIKey`, `Coordinate_Source` and `Coordinate_License`.

**Mapping check:** each source file `Structures_N.sdf` was matched to the record with `Source_Record_ID = N`;
for all 3,088 records the molecule title, atom and bond counts, S/Se/O counts and the HOMO–LUMO gap (Hartree)
agree with the curated database.

# 3D structures of the external collection (Database v1.1)

`structures_3d_external.sdf.gz` contains 33 B3LYP/LANL2DZ geometries supplied by the H. Kayı group (Ankara
University) and released with ChalMolDB under CC BY 4.0: the chalcogendiazoloquinoxaline monomers and hexamers
(32) and the B-system OOO monomer (1). Bond orders come from the family template; coordinates are those of the
geometry files (`data/sources/kayi_group_geometries`). A structure is included only if its atom/bond graph,
formula and ring chalcogen bond lengths pass the checks in `scripts/build_external_structures.py`
(full results: `data/build/external_structures.csv`). Each entry carries `ChalMolDB_Record_ID` (public EXT_
identifier), `ChalMolDB_InChIKey`, `Coordinate_Source` and `Coordinate_License`.
