import os
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

    missing_columns = [
        column
        for column in required_columns
        if column not in input_df.columns
    ]

    if missing_columns:
        raise ValueError(
            "Eksik sütunlar: " + ", ".join(missing_columns)
        )


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
        "Chalcogen Fraction": round(
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

    atom_counts = count_selected_atoms(mol)
    heavy_chalcogen_count = (
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

    result = {
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
        "Heavy Chalcogen Count": heavy_chalcogen_count,
        "Total Chalcogen Count": heavy_chalcogen_count,
        "Contains S": int(atom_counts["Sulfur Count"] > 0),
        "Contains Se": int(atom_counts["Selenium Count"] > 0),
        "Contains Te": int(atom_counts["Tellurium Count"] > 0),
        **calculate_chalcogen_aware_descriptors(mol),
        "Status": "Valid"
    }

    return result, None


def process_molecular_dataset(input_df):
    validate_input_dataframe(input_df)

    valid_results = []
    invalid_results = []

    for _, row in input_df.iterrows():
        valid_result, invalid_result = calculate_single_molecule_descriptors(
            molecule_id=row["Molecule_ID"],
            smiles=row["SMILES"]
        )

        if valid_result is not None:
            valid_results.append(valid_result)

        if invalid_result is not None:
            invalid_results.append(invalid_result)

    valid_df = pd.DataFrame(valid_results)
    invalid_df = pd.DataFrame(invalid_results)

    if not valid_df.empty:
        duplicate_mask = valid_df.duplicated(
            subset=["Canonical SMILES"],
            keep=False
        )
        valid_df["Duplicate Flag"] = duplicate_mask.astype(bool)

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


def fingerprint_to_list(fingerprint):
    array = np.zeros(fingerprint.GetNumBits(), dtype=int)
    DataStructs.ConvertToNumpyArray(fingerprint, array)
    return array.tolist()


def calculate_fingerprint_columns(mol):
    morgan_fp = morgan_generator.GetFingerprint(mol)
    maccs_fp = MACCSkeys.GenMACCSKeys(mol)

    morgan_bits = fingerprint_to_list(morgan_fp)
    maccs_bits = fingerprint_to_list(maccs_fp)

    fingerprint_columns = {}

    for index, value in enumerate(morgan_bits):
        fingerprint_columns[f"Morgan_{index}"] = value

    for index, value in enumerate(maccs_bits):
        fingerprint_columns[f"MACCS_{index}"] = value

    return fingerprint_columns


def create_fingerprint_dataset(input_df):
    validate_input_dataframe(input_df)

    fingerprint_results = []
    invalid_results = []

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
        mol = Chem.MolFromSmiles(smiles)

        if mol is None:
            invalid_results.append({
                "Molecule_ID": molecule_id,
                "SMILES": smiles,
                "Status": "Invalid SMILES"
            })
            continue

        result = {
            "Molecule_ID": molecule_id,
            "Original SMILES": smiles,
            "Canonical SMILES": Chem.MolToSmiles(mol, canonical=True)
        }

        result.update(calculate_fingerprint_columns(mol))
        fingerprint_results.append(result)

    return pd.DataFrame(fingerprint_results), pd.DataFrame(invalid_results)


def run_molecular_descriptor_platform(
    input_df,
    project_name="molecular_project"
):
    validate_input_dataframe(input_df)

    clean_input_df = input_df.copy()

    clean_input_df["Molecule_ID"] = (
        clean_input_df["Molecule_ID"]
        .astype(str)
        .str.strip()
    )

    valid_df, invalid_df, summary = process_molecular_dataset(clean_input_df)
    fingerprint_df, _ = create_fingerprint_dataset(clean_input_df)
    summary_df = pd.DataFrame([summary])

    descriptor_file = f"{project_name}_descriptors.csv"
    fingerprint_file = f"{project_name}_fingerprints.csv"
    invalid_file = f"{project_name}_invalid_smiles.csv"
    summary_file = f"{project_name}_summary.csv"
    zip_file = f"{project_name}_complete_outputs.zip"

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
                zip_output.write(file_name, arcname=file_name)

    return {
        "valid_df": valid_df,
        "invalid_df": invalid_df,
        "fingerprint_df": fingerprint_df,
        "summary_df": summary_df,
        "zip_file": zip_file
    }
