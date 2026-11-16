# ChalMolDB database changelog

## Database v1.3 (2026-11-16, tag `db-v1.3`)

Built from the Database v1.2 table (`data/releases/chalcogen_database_v1_2.csv.gz`) by
`scripts/build_database_v1_3.py`. 5,276 records (3,088 development, 400 external, 1,520 polymer, 268 literature).
All v1.2 values are unchanged (checked cell by cell in `tests/test_database_v1_2.py`). Two kinds of cells of
v1.2 records change: the curation note of 15 polymers (below) and the derived `Duplicate_Flag` of one record whose
structure now also occurs among the added records.

**Added: 197 polymers published 2015–2026 (polymer collection, POL_1324–POL_1520)**
- Linear donor–acceptor polymers from articles published 2015–2026, with the optical gap measured in the source
  article (185 film values, 12 solution values marked in `Solvent_or_Conditions`) and the source DOI in every record
  (`data/sources/dogan_kayi_2026_external/external_polymers_2015_2026.csv`).
- The record structure is the methyl-truncated oligomer built from the repeat unit with the convention of the
  existing polymer records; the repeat unit with full side chains is kept in `Source_Monomer_SMILES`. These records
  have no DFT values.
- Polymers already in the collection (same name with the same backbone, or the same source article; 23) are not
  added again.

**Added: 50 model compounds of literature polymers (literature collection, LIT_0219–LIT_0268)**
- B3LYP HOMO–LUMO gaps (HOMO and LUMO where given) of oligomer models reported in the source articles: gas phase,
  basis set and number of repeat units stated (`data/sources/dogan_kayi_2026_external/literature_model_compounds.csv`).
- The structure is built from the repeat unit with the ChalMolDB oligomer convention; the model of the source
  article may truncate side chains differently (noted in `Curation_Note`).

**Curation notes: 15 polymers of the 227-polymer set**
- POL_1116–1119, POL_1122–1125, POL_1135, POL_1141–1143, POL_1145, POL_1147, POL_1148: the source article given
  in the compilation could not be verified. The records and their values are kept as published in the compilation.

**Summary statistics (v1.2 → v1.3)**

| | v1.2 | v1.3 |
|---|---|---|
| Records | 5,029 | 5,276 |
| Unique InChIKeys | 4,229 | 4,443 |
| S/Se/Te core scope | 4,781 | 5,019 |
| Records containing S / Se / Te | 4,344 / 847 / 227 | 4,579 / 874 / 227 |
| Experimental optical gaps | 1,371 | 1,568 |

## Database v1.2 (2026-10-22, tag `db-v1.2`, DOI 10.5281/zenodo.23558298)

Built from the Database v1.1 table (`data/releases/chalcogen_database_v1_1.csv.gz`) by
`scripts/build_database_v1_2.py`. 5,029 records (3,088 development, 400 external, 1,323 polymer, 218 literature). All v1.1 records are
unchanged (checked cell by cell in `tests/test_database_v1_2.py`).

**Added: polymer collection (1,323 records, POL_0001–POL_1323)**
- Source: B. Liu, Y. Yan, M. Liu, *Nanoscale* 2025, 17, 7865 (doi:10.1039/d4nr03702b), data and geometries from the
  authors' GitHub repository (MIT License; copied unmodified to `data/sources/liu_2025_nanoscale/`).
- POL_0001–POL_1096: 1,096 conjugated polymers for organic solar cells with experimental HOMO, LUMO (cyclic
  voltammetry) and optical gap, compiled by S. Nagasawa, E. Al-Naamani, A. Saeki, *J. Phys. Chem. Lett.* 2018, 9, 2639.
  POL_1097–POL_1323: 227 polymers reported 2011–2024 with experimental optical gap and the DOI of the source article.
- For every polymer: B3LYP-D3/6-31G* (Gaussian 16) HOMO, LUMO and HOMO–LUMO gap of the oligomer built from the
  methyl-truncated repeat unit (`HOMO_eV`, `LUMO_eV`, `Eg_eV`); TD-B3LYP lowest singlet excitation (`TDDFT_S1_eV`)
  for the 1,096 set. The record structure (`Canonical_SMILES`, InChI, InChIKey) is the computed oligomer; the
  source monomer with full side chains is kept in `Source_Monomer_SMILES`.
- 3D coordinates for 1,064 oligomers (`data/structures_3d_polymer.sdf.gz`). The geometry file is attached only when its
  InChIKey equals that of the oligomer SMILES; 32 files are not attached (formula or connectivity differs from
  the SMILES; reason in `Curation_Note`).
- Chalcogen content: 1,296 records contain S and 73 contain Se; no Te (the source excluded its single Te polymer).
- 154 InChIKey groups (440 records) share an oligomer structure, mostly polymers that differ only in their side
  chains, which are truncated for the calculation; they are flagged by `Duplicate_Flag`, not removed.

**Added: literature collection (218 records, LIT_0001–LIT_0218)**
- 214 records of thia-, seleno- and telluro[n]helicenes from R. Kumar, R. Kar, D. K. Maity, *J. Comput. Chem.* 2025,
  46, e70175 (CC BY 4.0), transcribed from the Supporting Information (`data/sources/kumar_2025_helicenes/`):
  HOMO/LUMO for 13 functional settings (thia n = 1–5, seleno and telluro n = 1–4; 169 records), HSE06 (12),
  B2PLYP-D3(BJ) and ωB97XD gaps (20) and CCSD(T) fundamental gaps IP − EA (13), the reference of the source
  benchmark. The records are the unsubstituted helicenes (the CCSD(T) total energies correspond to the bare ring
  systems); the ring-fusion pattern is taken from the source geometries of the telluro[n]helicenes.
- 4 Se/Te fused bicyclic monomers (SeSe, SeTe, TeSe, TeTe; B3PW91) from Table 2 of Z. Büyükmumcu, F. Selampinar,
  *ACS Omega* 2026, 11, 15694 (CC BY 4.0). The polymer band gaps of that article are given only in figures and are
  not transcribed.

**New fields:** `Experimental_HOMO_eV`, `Experimental_LUMO_eV`, `TDDFT_S1_eV`, `Source_Monomer_SMILES`.

**Summary statistics (v1.1 → v1.2)**

| | v1.1 | v1.2 |
|---|---|---|
| Records | 3,488 | 5,029 |
| Unique InChIKeys | 3,173 | 4,229 |
| S/Se/Te core scope | 3,262 | 4,781 |
| Records containing S / Se / Te | 2,969 / 706 / 154 | 4,344 / 847 / 227 |
| Records with 3D coordinates | 3,147 | 4,211 |
| Experimental optical gaps | 48 | 1,371 |
| Experimental HOMO/LUMO | 0 | 1,096 |

## Database v1.1 (2026-10-08, tag `db-v1.1`, DOI 10.5281/zenodo.23238884)

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

**Removed: one unpublished system family (112 external records, EXT_0001–EXT_0112)**
- These records are not part of a published source and are not included. The identifiers are retired, not reused.

**Collection metadata**
- `Dataset_Owner`, `Dataset_Name` and `Paper_Use` use neutral collection labels.

**Development collection**
- `Donor_ID` / `Acceptor_ID` filled for all 3,088 records (previously 52), parsed from the source titles.
  DEV_1093 and DEV_1249 carry only a donor in the source title; their acceptor is assigned as A39
  (the acceptor fragment used by all other A39 records, and the only acceptor otherwise missing for donors D21/D24).
- DEV_0316 and DEV_0317 are a duplicate in the source Supporting Information (D6A4 twice, D6A6 absent);
  both are kept and annotated.
- Source DOI (10.1021/acs.jcim.5c00345), full reference and DFT level (B3LYP/6-311+G(d), Jaguar) added.
- One label format: the 52 source labels written as "D58,A1" are written as "D58A1" like the other records;
  the source label is kept in `Curation_Note`.

**External collection**
- New field `Repeat_Unit_SMILES` for all 400 records (repeat unit with two `*` attachment points).
- Exact structures (SMILES, InChI, InChIKey, element and atom counts) for the 96 chalcogendiazoloquinoxaline
  oligomer records: monomers and hexamers (32) confirmed against the authors' geometry files with 3D coordinates
  included; dimers to pentamers (64) from the family template verified on those files.
- New field `Curation_Note` with record-level curation remarks.

**Duplicate flag**
- `Duplicate_Flag` now means: the record's InChIKey is shared with at least one other record
  (210 development records and 240 external records, the latter because the same oligomer is reported at several levels of theory; no overlap between the collections). External records were previously flagged
  on equal system code and band gap; that flag is removed.

**3D structures**
- `data/structures_3d_external.sdf.gz`: 59 record-linked external structures from 50 validated geometry files
  (18 benzochalcogenadiazole, 32 chalcogendiazoloquinoxaline; CC BY 4.0, H. Kayı group).

**Source annotations**
- External notes cite the source articles; solvation recorded as gas phase or PCM (acetonitrile); numeric
  Experimental_Eg_Min/Max_eV and Experimental_Eg_Source; explanations for the seven B3LYP/6-31G(d) Te records
  without values (6-31G(d) is not defined for Te).

- No empty category fields: `Chalcogen_Type` is "None (control)" for the 226 control records, and
  `Solvent_or_Conditions` is "not stated in source" for the development collection.

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
