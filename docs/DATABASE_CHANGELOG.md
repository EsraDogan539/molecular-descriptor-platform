# ChalMolDB database changelog

## Database v1.1 (in preparation; not yet released)

Built from the archived Database v1 table (`data/releases/chalcogen_database_v1.csv.gz`) by
`scripts/build_database_v1_1.py`. 3,488 records (3,088 development, 400 external); property values of the retained
records are unchanged, record identifiers are not renumbered, and new records continue the numbering (EXT_0273–EXT_0512).

**Added: benzochalcogenadiazole oligomers (240 external records, EXT_0273–EXT_0512)**
- Monomer to hexamer HOMO, LUMO and Eg values of the 16 systems from Table 4 of Ozkilinc & Kayi,
  J. Mol. Model. 2019, 25, 167, at B3LYP/6-31G(d) (systems without Te), B3LYP/LANL2DZ and B3LYP/LANL2DZ with PCM
  (acetonitrile); transcribed in `data/sources/ozkilinc_kayi_2019_table4_oligomers.csv` and checked record by record.
- Exact structures for all 240 records (monomer and hexamer identity confirmed against the authors' geometry files,
  dimer to pentamer from the verified repeat unit); 3D coordinates for 27 records whose gas-phase geometry file
  passes all checks. Experimental monomer gaps (27 records) with their literature references.
- Three internal inconsistencies of the source table were corrected by withholding HOMO and LUMO (Eg kept) for
  7 records: the OSeO monomer to pentamer at B3LYP/LANL2DZ (LUMO - HOMO exceeds Eg by about 0.2 eV), the SeTeSe
  hexamer with PCM (HOMO/LUMO copied from the pentamer) and the SSeS hexamer at B3LYP/LANL2DZ (HOMO/LUMO copied
  from the SSS hexamer). In each case only the published Eg values reproduce the published polymer gap by the
  article's 1/n extrapolation (OSeO: 1.37 eV from Eg, 1.54 eV from LUMO - HOMO). Reasons are given in `Curation_Note`.

**Removed: B-system family (112 external records, EXT_0001–EXT_0112)**
- The authors' geometry files for this family are not consistent with the reported systems: every hexamer file
  has two H atoms fewer than the intended structure (C108H60 instead of C108H62), and 15 of 16 monomer files keep
  the bond lengths of the OOO system for S, Se and Te (not re-optimized after the chalcogen was changed).
  The property values could not be traced to valid geometries, so the family is removed. The identifiers are
  retired, not reused.

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
- New field `Repeat_Unit_SMILES` for all 160 records (repeat unit with two `*` attachment points).
- Exact structures (SMILES, InChI, InChIKey, element and atom counts) for the 96 chalcogendiazoloquinoxaline
  oligomer records: monomers and hexamers (32) confirmed against the authors' geometry files with 3D coordinates
  included; dimers to pentamers (64) from the family template verified on those files.
- New field `Curation_Note` with record-level curation remarks.

**Duplicate flag**
- `Duplicate_Flag` now means: the record's InChIKey is shared with at least one other record
  (210 development records and 240 external records, the latter because the same oligomer is reported at several levels of theory; no overlap between the collections). External records were previously flagged
  on equal system code and band gap; that flag is removed.

**3D structures**
- `data/structures_3d_external.sdf.gz`: 59 external structures (CC BY 4.0, H. Kayı group).

**Source annotations**
- External notes cite the source articles; solvation recorded as gas phase or PCM (acetonitrile); numeric
  Experimental_Eg_Min/Max_eV and Experimental_Eg_Source; explanations for the seven B3LYP/6-31G(d) Te records
  without values (6-31G(d) is not defined for Te).

**Summary statistics (v1 → v1.1)**

| | v1 | v1.1 |
|---|---|---|
| Records | 3,360 | 3,488 |
| External records | 272 | 400 |
| S/Se/Te core scope / control or non-core | 3,145 / 215 | 3,262 / 226 |
| Records containing S / Se / Te | 2,907 / 644 / 119 | 2,969 / 706 / 154 |
| Records with exact structure | 3,088 | 3,424 |
| Records with HOMO and LUMO | 3,280 | 3,417 |
| Unique InChIKeys | 2,983 | 3,173 |
| Te records with exact structure | 0 | 126 |
| Records with 3D coordinates | 3,088 | 3,147 |
| Experimental Eg values | 26 | 48 |
