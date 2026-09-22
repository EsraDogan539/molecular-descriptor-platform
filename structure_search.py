"""Structure-first search helpers for the curated ChalMolDB database."""

import pandas as pd

from rdkit import Chem, DataStructs
from rdkit.Chem import rdFingerprintGenerator


MORGAN_RADIUS = 2
MORGAN_N_BITS = 2048
SIMILARITY_METRIC = "Tanimoto"

_morgan_generator = rdFingerprintGenerator.GetMorganGenerator(
    radius=MORGAN_RADIUS,
    fpSize=MORGAN_N_BITS,
)


def standardize_query_smiles(smiles):
    """Parse and canonicalize a query SMILES, returning identity metadata."""
    if smiles is None or not str(smiles).strip():
        raise ValueError("Enter a SMILES query before running structure search.")

    text = str(smiles).strip()
    mol = Chem.MolFromSmiles(text)
    if mol is None:
        raise ValueError("The structure query is not a valid SMILES string.")

    canonical = Chem.MolToSmiles(mol, canonical=True)

    try:
        inchi = Chem.MolToInchi(mol)
        inchikey = Chem.InchiToInchiKey(inchi) if inchi else ""
    except Exception:
        inchi = ""
        inchikey = ""

    return {
        "query_smiles": text,
        "canonical_smiles": canonical,
        "inchi": inchi,
        "inchikey": inchikey,
        "mol": mol,
    }


def _valid_structure_rows(df):
    if "Canonical_SMILES" not in df.columns:
        return []

    rows = []
    for index, smiles in df["Canonical_SMILES"].items():
        if pd.isna(smiles) or not str(smiles).strip():
            continue
        mol = Chem.MolFromSmiles(str(smiles).strip())
        if mol is not None:
            rows.append((index, mol))
    return rows


def exact_structure_search(df, query):
    """Return records with the same standardized molecular identity."""
    canonical = query["canonical_smiles"]
    inchikey = query.get("inchikey", "")

    if inchikey and "InChIKey" in df.columns:
        keys = df["InChIKey"].fillna("").astype(str).str.strip()
        matched = df[keys == inchikey].copy()
        if not matched.empty:
            matched["Structure Search Mode"] = "Exact"
            return matched

    if "Canonical_SMILES" not in df.columns:
        return df.iloc[0:0].copy()

    canonical_targets = []
    for value in df["Canonical_SMILES"]:
        if pd.isna(value) or not str(value).strip():
            canonical_targets.append("")
            continue
        mol = Chem.MolFromSmiles(str(value).strip())
        canonical_targets.append(
            Chem.MolToSmiles(mol, canonical=True) if mol is not None else ""
        )

    mask = pd.Series(canonical_targets, index=df.index).eq(canonical)
    matched = df[mask].copy()
    matched["Structure Search Mode"] = "Exact"
    return matched


def substructure_search(df, query):
    """Return records whose standardized structure contains the query graph."""
    query_mol = query["mol"]
    matched_indices = []

    for index, target_mol in _valid_structure_rows(df):
        if target_mol.HasSubstructMatch(query_mol):
            matched_indices.append(index)

    matched = df.loc[matched_indices].copy()
    matched["Structure Search Mode"] = "Substructure"
    return matched


def similarity_structure_search(
    df,
    query,
    minimum_similarity=0.4,
    top_n=50,
):
    """Rank curated structures by Morgan/Tanimoto fingerprint similarity."""
    query_fp = _morgan_generator.GetFingerprint(query["mol"])

    indices = []
    fingerprints = []
    for index, target_mol in _valid_structure_rows(df):
        indices.append(index)
        fingerprints.append(_morgan_generator.GetFingerprint(target_mol))

    if not fingerprints:
        result = df.iloc[0:0].copy()
        result["Structure Similarity"] = pd.Series(dtype=float)
        result["Structure Search Mode"] = pd.Series(dtype=str)
        return result

    similarities = DataStructs.BulkTanimotoSimilarity(
        query_fp,
        fingerprints,
    )

    ranked = [
        (index, float(similarity))
        for index, similarity in zip(indices, similarities)
        if similarity >= minimum_similarity
    ]
    ranked.sort(key=lambda item: item[1], reverse=True)
    ranked = ranked[: int(top_n)]

    if not ranked:
        result = df.iloc[0:0].copy()
        result["Structure Similarity"] = pd.Series(dtype=float)
        result["Structure Search Mode"] = pd.Series(dtype=str)
        return result

    result = df.loc[[index for index, _ in ranked]].copy()
    similarity_map = dict(ranked)
    result["Structure Similarity"] = [
        round(similarity_map[index], 4)
        for index in result.index
    ]
    result["Structure Search Mode"] = "Similarity"
    return result


def structure_search(
    df,
    query_smiles,
    mode,
    minimum_similarity=0.4,
    top_n=50,
):
    """Standardize the query and dispatch to the selected search mode."""
    query = standardize_query_smiles(query_smiles)

    normalized_mode = str(mode).strip().lower()
    if normalized_mode == "exact":
        result = exact_structure_search(df, query)
    elif normalized_mode == "substructure":
        result = substructure_search(df, query)
    elif normalized_mode == "similarity":
        result = similarity_structure_search(
            df,
            query,
            minimum_similarity=minimum_similarity,
            top_n=top_n,
        )
    else:
        raise ValueError(
            "Structure search mode must be Exact, Substructure or Similarity."
        )

    return result, query
