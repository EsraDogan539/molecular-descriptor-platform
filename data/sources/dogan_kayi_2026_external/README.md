# Polymers published 2015–2026 and model compounds of literature polymers

Curated by E. N. Doğan and H. Kayı (Ankara University) for Database v1.3. Released under CC BY 4.0.
Every row gives the source article (reference and DOI); the values are those published in that article.

`external_polymers_2015_2026.csv`: linear donor–acceptor polymers from articles published 2015–2026.
- `repeat_unit_smiles`: repeat unit with full side chains, `*` at the two attachment points.
- `model_oligomer_smiles`, `model_oligomer_n`: methyl-truncated oligomer built from the repeat unit with the
  convention of the polymer collection (the record structure in the database).
- `measured_optical_gap_eV`, `gap_state`: optical gap measured in the source article (film or solution).
- `backbone`: `in ChalMolDB` marks polymers that are already in the polymer collection (same name with the same
  backbone, or the same source article); these rows are not added again.
- `used_in_external_test`, `exclusion_reason`: use as the external test set of the accompanying study
  (film values of the polymers that are not in ChalMolDB).

`literature_model_compounds.csv`: DFT HOMO–LUMO gaps of model compounds of literature polymers.
- `model_note`, `n_units`, `functional`, `basis`, `location_in_source`: the model and level of theory as described
  in the source article.
- `used_in_model`: gas-phase B3LYP models with a stated basis set and number of repeat units (added to the
  database); the other rows (periodic models, solvent models, enforced coplanarity, unstated basis set or model
  size) are listed for completeness only.
