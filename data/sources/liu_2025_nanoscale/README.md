# Source data: conjugated polymers with experimental optical gaps (Liu, Yan & Liu, Nanoscale 2025)

Files in this folder are unmodified copies from the public GitHub repository that accompanies

> B. Liu, Y. Yan, M. Liu, "Harnessing DFT and machine learning for accurate optical gap prediction in conjugated
> polymers", *Nanoscale* **2025**, *17*, 7865. https://doi.org/10.1039/d4nr03702b

Repository: https://github.com/Liu-Group-UF/Machine-Learning-for-Accurate-Optical-Gap-Prediction-in-Conjugated-Polymers
(commit 678ccffb4abf1e9cf961b2601866cf6deb278433). The repository is released under the MIT License; the licence
text is kept in `LICENSE` in this folder.

| File | Content |
|---|---|
| `conjugated-polymer-1096dp.csv` | 1,096 conjugated polymers: experimental HOMO, LUMO (cyclic voltammetry) and optical gap, source reference number, monomer and oligomer SMILES |
| `monomer-DFT-calculations-1096dp.csv` | B3LYP HOMO, LUMO and HOMO–LUMO gap of the monomers |
| `oligomer-DFT-TDDFT-calculations-1096dp.csv` | Oligomer HOMO, LUMO and gap at B3LYP, ωB97XD, CAM-B3LYP and PBE; lowest singlet excitation (TD-B3LYP) |
| `oligomer-xyz_opt_1096dp.zip` | Optimized oligomer geometries (`donor-<idx>_opt.xyz`) |
| `new_227/new-polymer-*.csv` | 227 polymers reported 2011–2024: experimental optical gap, year, source DOI, monomer and oligomer SMILES |
| `new_227/oligomer-DFT-*.csv` | B3LYP oligomer HOMO, LUMO and gap for those polymers |

The experimental values of the 1,096 polymers were compiled by S. Nagasawa, E. Al-Naamani and A. Saeki,
*J. Phys. Chem. Lett.* **2018**, *9*, 2639 (https://doi.org/10.1021/acs.jpclett.8b00635); `Ref.No` refers to the
reference list of that compilation. According to the Nanoscale article, geometries and orbital energies were
computed with Gaussian 16 at B3LYP-D3/6-31G*, with long alkyl side chains truncated to methyl groups and the
monomer repeated to build the oligomer.

ChalMolDB imports these files with `scripts/build_database_v1_2.py` as the polymer collection (records POL_0001–).
