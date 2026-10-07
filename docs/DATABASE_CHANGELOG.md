# ChalMolDB database changelog

## Database v1.1 (in preparation; not yet released)

Built from Database v1 by `scripts/build_database_v1_1.py`. Record count and record identifiers unchanged
(3,360 records: 3,088 development, 272 external). Property values are unchanged.

**Collection metadata**
- `Dataset_Owner`, `Dataset_Name` and `Paper_Use` use neutral collection labels.

**Development collection**
- `Donor_ID` / `Acceptor_ID` filled for all 3,088 records (previously 52), parsed from the source titles.
  DEV_1093 and DEV_1249 carry only a donor in the source title; their acceptor is assigned as A39
  (the acceptor fragment used by all other A39 records, and the only acceptor otherwise missing for donors D21/D24).
- DEV_0316 and DEV_0317 are a duplicate in the source Supporting Information (D6A4 twice, D6A6 absent);
  both are kept and annotated.
- Source DOI (10.1021/acs.jcim.5c00345), full reference and DFT level (B3LYP/6-311+G(d), Jaguar) added.

**External collection**
- New field `Repeat_Unit_SMILES` for all 272 records (repeat unit with two `*` attachment points).
- Exact structures (SMILES, InChI, InChIKey, element and atom counts) for 112 records:
  - chalcogendiazoloquinoxaline monomers and hexamers (32): confirmed against the authors' geometry files,
    3D coordinates included;
  - chalcogendiazoloquinoxaline dimers to pentamers (64): from the family template verified on the
    monomer and hexamer files;
  - B-system monomers (16): identity confirmed from the geometry files; 3D coordinates included for OOO only.
- B-system records: family renamed to *ethylenedioxychalcogenophene-benzochalcogenadiazole DAD (B system)*;
  reference changed to the Kayı group calculations. The B-system hexamer geometry files contain two H atoms
  fewer than the intended structures, and 15 of 16 monomer files were not re-optimized after the chalcogen
  was changed. B-system property values are therefore kept as supplied but marked unverified in
  `Curation_Note`, and no structure is assigned to B-system dimer to hexamer records.
- New field `Curation_Note` with record-level curation remarks.

**Duplicate flag**
- `Duplicate_Flag` now means: the record's InChIKey is shared with at least one other record
  (210 development records; no overlap between the collections). External records were previously flagged
  on equal system code and band gap; that flag is removed.

**3D structures**
- `data/structures_3d_external.sdf.gz`: 33 external structures (CC BY 4.0, H. Kayı group).
