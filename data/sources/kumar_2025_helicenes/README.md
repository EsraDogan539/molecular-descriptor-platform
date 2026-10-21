# Source data: thia-, seleno- and telluro[n]helicenes (J. Comput. Chem. 2025)

> R. Kumar, R. Kar, D. K. Maity, "Accurate Prediction of HOMO–LUMO Gap Using DFT Functional and Application to
> Next-Generation Organic Telluro[n]Helicenes Materials", *J. Comput. Chem.* **2025**, *46* (18), e70175.
> doi:10.1002/jcc.70175 (CC BY 4.0; full text and Supporting Information via Europe PMC, PMC12223471).

| File | Source |
|---|---|
| `helicene_orbital_energies_13_functionals.csv` | SI Tables S3/S4 (thia[n], n = 1–5), S6/S7 (seleno[n], n = 1–4), S9/S10 (telluro[n], n = 1–4): HOMO and LUMO in eV for 13 functional settings. Set 1: geometry optimized with the same functional; Set 2: single point on the B3LYP geometry. Basis: 6-311++G(d,p), LANL2DZ for Te (as stated in the article). Gas phase. |
| `helicene_hse06_table_s12.csv` | SI Table S12: HSE06 HOMO, LUMO (a.u.) and gap (eV), n = 1–4 |
| `helicene_gaps_table_s13.csv` | SI Table S13: B2PLYP-D3(BJ) gaps n = 1–5 and ωB97XD gaps of telluro[n]helicenes (def2-TZVPP for Te) |
| `helicene_ccsdt_tables_s5_s8_s11.csv` | SI Tables S5, S8, S11: CCSD(T) energies, ionization potential, electron affinity and IP − EA |
| `telluro_n_helicenes_end_substituted_B3LYP-D_DCM.xyz` | SI section 1: coordinates of the end-substituted (Br, TMS) telluro[n]helicenes, n = 1–10 |

Values were parsed from the SI document and checked: every eV value equals the a.u. value × 27.2114 within 0.01 eV
(Tables S3–S10), and IP − EA equals the tabulated gap (Tables S5, S8, S11).

Structure assignment in ChalMolDB: the benchmark systems are the unsubstituted helicenes (the CCSD(T) total
energies correspond to the bare ring systems, e.g. −552.11 a.u. for [1]TH, thiophene); the ring-fusion pattern is
taken from the telluro[n]helicene geometries of SI section 1 with the Br and TMS end groups replaced by H, and the
S and Se analogues have the same pattern with Te replaced. The geometries themselves belong to the substituted
molecules and are therefore not attached as 3D coordinates.
