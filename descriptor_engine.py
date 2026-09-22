import os
import re
import tempfile
import zipfile
import numpy as np
import pandas as pd

from rdkit import Chem, DataStructs, RDLogger
from rdkit.Chem import (
    Descriptors,
    rdMolDescriptors,
    rdFingerprintGenerator,
    MACCSkeys
)

RDLogger.DisableLog("rdApp.error")


MORGAN_RADIUS = 2
MORGAN_N_BITS = 2048
CHALCOGEN_SYMBOLS = ("S", "Se", "Te")

morgan_generator = rdFingerprintGenerator.GetMorganGenerator(
    radius=MORGAN_RADIUS,
    fpSize=MORGAN_N_BITS
)


def validate_input_dataframe(input_df):
    required_columns = ["Molecule_ID", "SMILES"]

    if input_df is None or input_df.empty:
        raise ValueError("The uploaded CSV does not contain any data rows.")

    missing_columns = [
        column
        for column in required_columns
        if column not in input_df.columns
    ]

    if missing_columns:
        raise ValueError(
            "Missing required columns: " + ", ".join(missing_columns)
            + ". Required columns are Molecule_ID and SMILES."
        )

    molecule_ids = input_df["Molecule_ID"]
    missing_id_mask = molecule_ids.isna() | molecule_ids.astype(str).str.strip().eq("")
    if missing_id_mask.any():
        rows = (input_df.index[missing_id_mask] + 2).tolist()
        preview = ", ".join(map(str, rows[:8]))
        suffix = "..." if len(rows) > 8 else ""
        raise ValueError(
            f"Molecule_ID is missing or empty in CSV row(s): {preview}{suffix}"
        )

    normalized_ids = molecule_ids.astype(str).str.strip()
    duplicate_ids = normalized_ids[normalized_ids.duplicated(keep=False)].unique().tolist()
    if duplicate_ids:
        preview = ", ".join(map(str, duplicate_ids[:8]))
        suffix = "..." if len(duplicate_ids) > 8 else ""
        raise ValueError(
            "Molecule_ID values must be unique. Duplicate ID(s): "
            f"{preview}{suffix}"
        )


def sanitize_project_name(project_name):
    name = str(project_name or "").strip()
    if not name:
        return "chalcogen_project"
    name = re.sub(r"[^A-Za-z0-9._-]+", "_", name)
    name = name.strip("._-")
    return name[:80] or "chalcogen_project"


def count_selected_atoms(mol):
    atom_symbols = [atom.GetSymbol() for atom in mol.GetAtoms()]

    return {
        "Carbon Count": atom_symbols.count("C"),
        "Nitrogen Count": atom_symbols.count("N"),
        "Oxygen Count": atom_symbols.count("O"),
        "Sulfur Count": atom_symbols.count("S"),
        "Selenium Count": atom_symbols.count("Se"),
        "Tellurium Count": atom_symbols.count("Te"),
        "Phosphorus Count": atom_symbols.count("P"),
        "Halogen Count": sum(
            atom_symbols.count(symbol)
            for symbol in ["F", "Cl", "Br", "I"]
        )
    }


def classify_chalcogen_type(atom_counts):
    present = []
    if atom_counts["Sulfur Count"] > 0:
        present.append("S")
    if atom_counts["Selenium Count"] > 0:
        present.append("Se")
    if atom_counts["Tellurium Count"] > 0:
        present.append("Te")

    if not present:
        return "None"
    if len(present) == 1:
        return present[0]
    return "Mixed"


def calculate_chalcogen_aware_descriptors(mol):
    chalcogen_atoms = [
        atom for atom in mol.GetAtoms()
        if atom.GetSymbol() in CHALCOGEN_SYMBOLS
    ]

    heavy_atom_count = max(Descriptors.HeavyAtomCount(mol), 1)

    aromatic_chalcogen_count = sum(
        int(atom.GetIsAromatic()) for atom in chalcogen_atoms
    )
    non_aromatic_chalcogen_count = (
        len(chalcogen_atoms) - aromatic_chalcogen_count
    )

    chalcogen_c_bond_count = 0
    chalcogen_heteroatom_bond_count = 0
    aromatic_neighbor_count = 0
    ring_incorporated_count = 0

    seen_bonds = set()
    for atom in chalcogen_atoms:
        if atom.IsInRing():
            ring_incorporated_count += 1

        for bond in atom.GetBonds():
            bond_idx = bond.GetIdx()
            if bond_idx in seen_bonds:
                continue
            seen_bonds.add(bond_idx)

            neighbor = bond.GetOtherAtom(atom)
            symbol = neighbor.GetSymbol()

            if symbol == "C":
                chalcogen_c_bond_count += 1
            elif symbol != "H":
                chalcogen_heteroatom_bond_count += 1

            if neighbor.GetIsAromatic():
                aromatic_neighbor_count += 1

    conjugated_bonds = [
        bond for bond in mol.GetBonds()
        if bond.GetIsConjugated()
    ]
    aromatic_bonds = [
        bond for bond in mol.GetBonds()
        if bond.GetIsAromatic()
    ]
    total_bonds = max(mol.GetNumBonds(), 1)

    conjugated_atom_indices = set()
    for bond in conjugated_bonds:
        conjugated_atom_indices.add(bond.GetBeginAtomIdx())
        conjugated_atom_indices.add(bond.GetEndAtomIdx())

    heteroaromatic_rings = 0
    ring_info = mol.GetRingInfo()
    for ring in ring_info.AtomRings():
        ring_atoms = [mol.GetAtomWithIdx(idx) for idx in ring]
        if ring_atoms and all(atom.GetIsAromatic() for atom in ring_atoms):
            if any(atom.GetSymbol() not in ("C", "H") for atom in ring_atoms):
                heteroaromatic_rings += 1

    total_chalcogens = len(chalcogen_atoms)

    return {
        "Target Chalcogen Fraction": round(
            total_chalcogens / heavy_atom_count,
            4
        ),
        "Aromatic Chalcogen Count": aromatic_chalcogen_count,
        "NonAromatic Chalcogen Count": non_aromatic_chalcogen_count,
        "Mixed Chalcogen Flag": int(
            len({a.GetSymbol() for a in chalcogen_atoms}) > 1
        ),
        "Ring Incorporated Chalcogen Count": ring_incorporated_count,
        "Chalcogen-C Bond Count": chalcogen_c_bond_count,
        "Chalcogen-Heteroatom Bond Count": chalcogen_heteroatom_bond_count,
        "Aromatic Neighbor Count": aromatic_neighbor_count,
        "Conjugated Bond Count": len(conjugated_bonds),
        "Aromatic Bond Fraction": round(
            len(aromatic_bonds) / total_bonds,
            4
        ),
        "Conjugated Atom Fraction": round(
            len(conjugated_atom_indices) / heavy_atom_count,
            4
        ),
        "Heteroaromatic Ring Count": heteroaromatic_rings,
    }


def _descriptor_result_from_mol(molecule_id, smiles, mol):
    atom_counts = count_selected_atoms(mol)
    target_chalcogen_count = (
        atom_counts["Sulfur Count"]
        + atom_counts["Selenium Count"]
        + atom_counts["Tellurium Count"]
    )

    canonical_smiles = Chem.MolToSmiles(mol, canonical=True)

    try:
        inchi = Chem.MolToInchi(mol)
        inchikey = Chem.InchiToInchiKey(inchi) if inchi else ""
    except Exception:
        inchi = ""
        inchikey = ""

    return {
        "Molecule_ID": molecule_id,
        "Original SMILES": smiles,
        "Canonical SMILES": canonical_smiles,
        "InChI": inchi,
        "InChIKey": inchikey,
        "Valid SMILES": True,
        "Molecular Formula": rdMolDescriptors.CalcMolFormula(mol),
        "Molecular Weight": round(Descriptors.MolWt(mol), 4),
        "Exact Molecular Weight": round(Descriptors.ExactMolWt(mol), 4),
        "LogP": round(Descriptors.MolLogP(mol), 4),
        "TPSA": round(Descriptors.TPSA(mol), 4),
        "H-Bond Donors": Descriptors.NumHDonors(mol),
        "H-Bond Acceptors": Descriptors.NumHAcceptors(mol),
        "Rotatable Bonds": Descriptors.NumRotatableBonds(mol),
        "Ring Count": Descriptors.RingCount(mol),
        "Aromatic Ring Count": Descriptors.NumAromaticRings(mol),
        "Aliphatic Ring Count": Descriptors.NumAliphaticRings(mol),
        "Heavy Atom Count": Descriptors.HeavyAtomCount(mol),
        "Heteroatom Count": Descriptors.NumHeteroatoms(mol),
        "Fraction Csp3": round(Descriptors.FractionCSP3(mol), 4),
        "Formal Charge": Chem.GetFormalCharge(mol),
        "Bertz CT": round(Descriptors.BertzCT(mol), 4),
        "Balaban J": round(Descriptors.BalabanJ(mol), 4),
        **atom_counts,
        "Chalcogen Type": classify_chalcogen_type(atom_counts),
        "Target Chalcogen Count": target_chalcogen_count,
        "Contains S": int(atom_counts["Sulfur Count"] > 0),
        "Contains Se": int(atom_counts["Selenium Count"] > 0),
        "Contains Te": int(atom_counts["Tellurium Count"] > 0),
        **calculate_chalcogen_aware_descriptors(mol),
        "Status": "Valid"
    }


def calculate_single_molecule_descriptors(molecule_id, smiles):
    if pd.isna(smiles):
        return None, {
            "Molecule_ID": molecule_id,
            "SMILES": "",
            "Status": "Missing SMILES"
        }

    smiles = str(smiles).strip()

    if not smiles:
        return None, {
            "Molecule_ID": molecule_id,
            "SMILES": "",
            "Status": "Empty SMILES"
        }

    mol = Chem.MolFromSmiles(smiles)

    if mol is None:
        return None, {
            "Molecule_ID": molecule_id,
            "SMILES": smiles,
            "Status": "Invalid SMILES"
        }

    return _descriptor_result_from_mol(molecule_id, smiles, mol), None


def process_molecular_dataset(input_df):
    validate_input_dataframe(input_df)

    valid_results = []
    invalid_results = []

    for row in input_df.itertuples(index=False):
        valid_result, invalid_result = calculate_single_molecule_descriptors(
            molecule_id=row.Molecule_ID,
            smiles=row.SMILES
        )

        if valid_result is not None:
            valid_results.append(valid_result)

        if invalid_result is not None:
            invalid_results.append(invalid_result)

    valid_df = pd.DataFrame(valid_results)
    invalid_df = pd.DataFrame(invalid_results)

    if not valid_df.empty:
        # Structure-level duplicate detection follows the standardized identity
        # used by the curated database. InChIKey is preferred when available;
        # canonical SMILES is used only as a fallback.
        identity = valid_df["InChIKey"].fillna("").astype(str).str.strip()
        fallback = valid_df["Canonical SMILES"].fillna("").astype(str)
        valid_df["_Structure_Identity"] = identity.where(identity.ne(""), fallback)
        duplicate_mask = valid_df.duplicated(
            subset=["_Structure_Identity"],
            keep=False
        )
        valid_df["Duplicate Flag"] = duplicate_mask.astype(bool)
        valid_df = valid_df.drop(columns=["_Structure_Identity"])

    total_records = len(input_df)

    summary = {
        "Total Records": total_records,
        "Valid Molecules": len(valid_df),
        "Invalid Molecules": len(invalid_df),
        "Duplicate Molecules": int(
            valid_df["Duplicate Flag"].sum()
        ) if not valid_df.empty else 0,
        "Success Rate (%)": round(
            100 * len(valid_df) / total_records,
            2
        ) if total_records else 0
    }

    return valid_df, invalid_df, summary


def fingerprint_to_array(fingerprint):
    array = np.zeros(fingerprint.GetNumBits(), dtype=np.uint8)
    DataStructs.ConvertToNumpyArray(fingerprint, array)
    return array


def fingerprint_to_list(fingerprint):
    return fingerprint_to_array(fingerprint).tolist()


def calculate_fingerprint_vector(mol):
    morgan_fp = morgan_generator.GetFingerprint(mol)
    maccs_fp = MACCSkeys.GenMACCSKeys(mol)
    morgan_bits = fingerprint_to_array(morgan_fp)
    maccs_bits = fingerprint_to_array(maccs_fp)
    return np.concatenate((morgan_bits, maccs_bits))


def fingerprint_column_names(vector_length):
    morgan_count = MORGAN_N_BITS
    maccs_count = vector_length - morgan_count
    return (
        [f"Morgan_{index}" for index in range(morgan_count)]
        + [f"MACCS_{index}" for index in range(maccs_count)]
    )


def calculate_fingerprint_columns(mol):
    vector = calculate_fingerprint_vector(mol)
    columns = fingerprint_column_names(len(vector))
    return dict(zip(columns, vector.tolist()))


def _build_fingerprint_dataframe(metadata_rows, fingerprint_vectors):
    if not metadata_rows:
        return pd.DataFrame()

    matrix = np.vstack(fingerprint_vectors)
    bit_columns = fingerprint_column_names(matrix.shape[1])
    bit_df = pd.DataFrame(matrix, columns=bit_columns, dtype=np.uint8)
    metadata_df = pd.DataFrame(metadata_rows).reset_index(drop=True)
    return pd.concat([metadata_df, bit_df], axis=1)


def create_fingerprint_dataset(input_df):
    validate_input_dataframe(input_df)

    metadata_rows = []
    fingerprint_vectors = []
    invalid_results = []

    for row in input_df.itertuples(index=False):
        molecule_id = row.Molecule_ID
        smiles = row.SMILES

        if pd.isna(smiles):
            invalid_results.append({
                "Molecule_ID": molecule_id,
                "SMILES": "",
                "Status": "Missing SMILES"
            })
            continue

        smiles = str(smiles).strip()
        mol = Chem.MolFromSmiles(smiles)

        if mol is None:
            invalid_results.append({
                "Molecule_ID": molecule_id,
                "SMILES": smiles,
                "Status": "Invalid SMILES"
            })
            continue

        metadata_rows.append({
            "Molecule_ID": molecule_id,
            "Original SMILES": smiles,
            "Canonical SMILES": Chem.MolToSmiles(mol, canonical=True)
        })
        fingerprint_vectors.append(calculate_fingerprint_vector(mol))

    fingerprint_df = _build_fingerprint_dataframe(
        metadata_rows,
        fingerprint_vectors,
    )
    return fingerprint_df, pd.DataFrame(invalid_results)


def process_molecular_dataset_with_fingerprints(input_df):
    """Calculate descriptors and fingerprints from a single RDKit parse per row."""
    validate_input_dataframe(input_df)

    valid_results = []
    invalid_results = []
    fingerprint_metadata = []
    fingerprint_vectors = []

    for _, row in input_df.iterrows():
        molecule_id = row["Molecule_ID"]
        smiles = row["SMILES"]

        if pd.isna(smiles):
            invalid_results.append({
                "Molecule_ID": molecule_id,
                "SMILES": "",
                "Status": "Missing SMILES"
            })
            continue

        smiles = str(smiles).strip()
        if not smiles:
            invalid_results.append({
                "Molecule_ID": molecule_id,
                "SMILES": "",
                "Status": "Empty SMILES"
            })
            continue

        mol = Chem.MolFromSmiles(smiles)
        if mol is None:
            invalid_results.append({
                "Molecule_ID": molecule_id,
                "SMILES": smiles,
                "Status": "Invalid SMILES"
            })
            continue

        descriptor_result = _descriptor_result_from_mol(
            molecule_id=molecule_id,
            smiles=smiles,
            mol=mol,
        )
        valid_results.append(descriptor_result)

        fingerprint_metadata.append({
            "Molecule_ID": molecule_id,
            "Original SMILES": smiles,
            "Canonical SMILES": descriptor_result["Canonical SMILES"],
        })
        fingerprint_vectors.append(calculate_fingerprint_vector(mol))

    valid_df = pd.DataFrame(valid_results)
    invalid_df = pd.DataFrame(invalid_results)
    fingerprint_df = _build_fingerprint_dataframe(
        fingerprint_metadata,
        fingerprint_vectors,
    )

    if not valid_df.empty:
        identity = valid_df["InChIKey"].fillna("").astype(str).str.strip()
        fallback = valid_df["Canonical SMILES"].fillna("").astype(str)
        valid_df["_Structure_Identity"] = identity.where(identity.ne(""), fallback)
        duplicate_mask = valid_df.duplicated(
            subset=["_Structure_Identity"],
            keep=False,
        )
        valid_df["Duplicate Flag"] = duplicate_mask.astype(bool)
        valid_df = valid_df.drop(columns=["_Structure_Identity"])

    total_records = len(input_df)
    summary = {
        "Total Records": total_records,
        "Valid Molecules": len(valid_df),
        "Invalid Molecules": len(invalid_df),
        "Duplicate Molecules": int(
            valid_df["Duplicate Flag"].sum()
        ) if not valid_df.empty else 0,
        "Success Rate (%)": round(
            100 * len(valid_df) / total_records,
            2,
        ) if total_records else 0,
    }

    return valid_df, invalid_df, fingerprint_df, summary



def run_molecular_descriptor_platform(
    input_df,
    project_name="molecular_project"
):
    validate_input_dataframe(input_df)

    safe_project_name = sanitize_project_name(project_name)
    output_dir = tempfile.mkdtemp(prefix="chalmoldb_")

    clean_input_df = input_df.copy()

    clean_input_df["Molecule_ID"] = (
        clean_input_df["Molecule_ID"]
        .astype(str)
        .str.strip()
    )

    valid_df, invalid_df, fingerprint_df, summary = (
        process_molecular_dataset_with_fingerprints(clean_input_df)
    )
    summary_df = pd.DataFrame([summary])

    descriptor_file = os.path.join(output_dir, f"{safe_project_name}_descriptors.csv")
    fingerprint_file = os.path.join(output_dir, f"{safe_project_name}_fingerprints.csv")
    invalid_file = os.path.join(output_dir, f"{safe_project_name}_invalid_smiles.csv")
    summary_file = os.path.join(output_dir, f"{safe_project_name}_summary.csv")
    zip_file = os.path.join(output_dir, f"{safe_project_name}_complete_outputs.zip")

    valid_df.to_csv(descriptor_file, index=False)
    fingerprint_df.to_csv(fingerprint_file, index=False)
    invalid_df.to_csv(invalid_file, index=False)
    summary_df.to_csv(summary_file, index=False)

    output_files = [
        descriptor_file,
        fingerprint_file,
        invalid_file,
        summary_file
    ]

    with zipfile.ZipFile(
        zip_file,
        mode="w",
        compression=zipfile.ZIP_DEFLATED
    ) as zip_output:
        for file_name in output_files:
            if os.path.exists(file_name):
                zip_output.write(file_name, arcname=os.path.basename(file_name))

    return {
        "valid_df": valid_df,
        "invalid_df": invalid_df,
        "fingerprint_df": fingerprint_df,
        "summary_df": summary_df,
        "zip_file": zip_file,
        "output_dir": output_dir,
        "project_name": safe_project_name,
    }
