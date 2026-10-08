"""Source annotations for the external collection (pure pandas; idempotent).

Taken from the source publications:
  Ozkilinc & Kayi, J. Mol. Model. 2019, 25, 167 (benzochalcogenadiazole polymers; Methods and Table 4)
  Kayi, Sen & Ozkilinc, J. Mol. Model. 2024, 30, 179 (chalcogendiazoloquinoxaline oligomers and polymers)

Changes:
  * Notes: "Source family: Hakan_2019 / Hakan_2024_quinoxaline" replaced by the bibliographic source.
  * Solvent_or_Conditions: "yes" -> "PCM (acetonitrile)", "no" -> "gas phase" (external records).
  * Experimental_Eg_Min_eV / _Max_eV: numeric columns; Experimental_Eg_eV keeps the value(s) as reported.
  * Experimental_Eg_Source: which literature reference of the 2019 article each experimental value comes from.
  * Curation_Note: why the B3LYP/6-31G(d) tellurium records carry no value; level-of-theory provenance notes.
  * HOMO/LUMO withheld for 7 oligomer records whose published HOMO/LUMO pair is inconsistent with the published Eg
    (HOMO_LUMO_CORRECTIONS); Eg is kept.
  * Explicit values instead of empty cells: Chalcogen_Type "None (control)" for control records, and
    Solvent_or_Conditions "not stated in source" for the development collection.
  * Development labels in the form "D58,A1" written as "D58A1" like the rest of the collection (source label kept
    in Curation_Note).

Run on its own (python scripts/annotate_external_sources.py) or from build_database_v1_1.py.
"""

import re
from pathlib import Path

import pandas as pd

ROOT = Path(__file__).resolve().parents[1]
MASTER = ROOT / "data" / "chalcogen_database_v1_master.csv.gz"

SRC_2019 = "Ozkilinc & Kayi, J. Mol. Model. 2019, 25, 167, Table 4"
SRC_2024 = "Kayi, Sen & Ozkilinc, J. Mol. Model. 2024, 30, 179, Table 4"
NOTE_REPLACEMENTS = {
    "Source family: Hakan_2019": f"Source: {SRC_2019}",
    "Source family: Hakan_2024_quinoxaline": f"Source: {SRC_2024}",
}

# Table 4 of the 2019 article: polymer experimental band gaps and their footnotes
# (a [26], b [27], c [28], f [31, 32], g [30], h [29]; reference numbers of the 2019 article).
EXPERIMENTAL_2019 = {
    "OOO": "1.66 eV, ref. [27]",
    "OSO": "1.61 eV, ref. [28]",
    "OSeO": "1.43 eV, ref. [28]",
    "SOS": "1.62 eV, ref. [27]",
    "SSS": "1.80 eV, refs. [31, 32]; 1.50 eV, ref. [30]",
    "SSeS": "1.46 eV, ref. [30]; 0.83 eV, refs. [31, 32]",
    "SeOSe": "1.59 eV, ref. [27]",
}
# Table 4 of the 2019 article: monomer experimental band gaps (footnotes a [26], c [28], h [29]).
EXPERIMENTAL_2019_MONOMER = {
    "OOO": "2.46 eV, ref. [26]",
    "OSO": "2.10 eV, ref. [28]",
    "OSeO": "2.44 eV, ref. [28]",
    "SOS": "2.47 eV, ref. [26]",
    "SSS": "2.43 eV, ref. [26]",
    "SSeS": "2.29 eV, ref. [29]",
    "SeOSe": "2.41 eV, ref. [26]",
    "SeSSe": "2.33 eV, ref. [26]",
    "SeSeSe": "2.19 eV, ref. [29]",
}
# Internal inconsistencies of Table 4 found during curation. In all three cases the published Eg values reproduce
# the published polymer band gaps by the article's 1/n extrapolation (e.g. OSeO, B3LYP/LANL2DZ: 1.37 eV from Eg,
# 1.54 eV from LUMO - HOMO), so Eg is kept and the HOMO/LUMO pair, which cannot belong to that Eg, is withheld.
HOMO_LUMO_CORRECTIONS = {
    ("OSeO", "LANL2DZ", "no", "*"): (
        "HOMO and LUMO withheld: in the source table LUMO - HOMO exceeds Eg by about 0.2 eV for the OSeO "
        "oligomers at this level, and only the published Eg values reproduce the published polymer gap "
        "(1.37 eV) by 1/n extrapolation. Eg is kept as published."),
    ("SeTeSe", "LANL2DZ", "yes", "6"): (
        "HOMO and LUMO withheld: in the source table they are identical to those of the SeTeSe pentamer at the same "
        "level (LUMO - HOMO = 1.31 eV vs Eg = 1.27 eV); the published Eg follows the series trend and reproduces the "
        "published polymer gap. Eg is kept as published."),
    ("SSeS", "LANL2DZ", "no", "6"): (
        "HOMO and LUMO withheld: in the source table they are identical to those of the SSS hexamer and break the "
        "trend of the SSeS series; the published Eg follows the series trend and reproduces the published polymer "
        "gap. Eg is kept as published."),
}
# Table 4 footnotes d, e, i: computed values that the 2019 article took from earlier work of the same group.
EARLIER_WORK_2019 = {
    ("OSO", "6-31G(d)", "no"): "ref. [35]",
    ("OSO", "LANL2DZ", "no"): "ref. [35]",
    ("OSeO", "LANL2DZ", "no"): "ref. [36]",
    ("OSeO", "LANL2DZ", "yes"): "ref. [36]",
    ("SeSeSe", "LANL2DZ", "no"): "ref. [37]",
    ("TeTeTe", "LANL2DZ", "no"): "ref. [37]",
}
NO_TE_BASIS = ("No value in the source: B3LYP/6-31G(d) was applied only to systems without tellurium "
               "(the 6-31G(d) basis set is not defined for Te); the record documents the level-of-theory grid of "
               f"{SRC_2019}.")
SINGLE_POINT = ("B3LYP/6-31G(d) single-point energies on B3LYP/LANL2DZ geometries "
                f"({SRC_2019}, Methods).")
PCM_NOTE = ("Geometries re-optimized and energies computed with PCM (acetonitrile) at B3LYP/LANL2DZ "
            f"({SRC_2019}, Methods).")


def _append(note, text):
    note = "" if pd.isna(note) else str(note)
    if text in note:
        return note
    return f"{note} | {text}" if note else text


def _parse_values(text):
    if pd.isna(text):
        return []
    return [float(v) for v in str(text).replace(",", ";").split(";") if v.strip()]


def annotate(df):
    df = df.copy()
    ext = df.Record_ID.astype(str).str.startswith("HAKAN_")
    a = ext & df.Family.eq("benzochalcogenadiazole")

    for col in ("Notes", "Solvent_or_Conditions", "Curation_Note"):
        df[col] = df[col].astype("object")
    for old, new in NOTE_REPLACEMENTS.items():
        df.loc[ext, "Notes"] = df.loc[ext, "Notes"].astype(str).str.replace(old, new, regex=False)

    raw_solvent = df["Solvent_or_Conditions"].copy()
    df.loc[ext & raw_solvent.eq("yes"), "Solvent_or_Conditions"] = "PCM (acetonitrile)"
    df.loc[ext & raw_solvent.eq("no"), "Solvent_or_Conditions"] = "gas phase"
    solvent_key = df["Solvent_or_Conditions"].map({"PCM (acetonitrile)": "yes", "gas phase": "no"})

    # numeric experimental columns, placed after Experimental_Eg_eV
    values = df["Experimental_Eg_eV"].map(_parse_values)
    new_cols = {
        "Experimental_Eg_Min_eV": values.map(lambda v: min(v) if v else pd.NA),
        "Experimental_Eg_Max_eV": values.map(lambda v: max(v) if v else pd.NA),
        "Experimental_Eg_Source": pd.Series(pd.NA, index=df.index, dtype="object"),
    }
    for name in reversed(list(new_cols)):
        if name in df.columns:
            df = df.drop(columns=name)
        df.insert(df.columns.get_loc("Experimental_Eg_eV") + 1, name, new_cols[name])
    for col in ("Experimental_Eg_Min_eV", "Experimental_Eg_Max_eV"):
        df[col] = pd.to_numeric(df[col])

    for idx in df.index[a]:
        code, basis, solv = df.at[idx, "System_Code"], df.at[idx, "Basis_Set"], solvent_key.at[idx]
        n = str(df.at[idx, "Oligomer_n"])
        if df.at[idx, "Experimental_Eg_eV"] is not pd.NA and pd.notna(df.at[idx, "Experimental_Eg_eV"]):
            table = EXPERIMENTAL_2019 if n == "polymer" else EXPERIMENTAL_2019_MONOMER
            assert n in ("polymer", "1") and code in table, (df.at[idx, "Record_ID"], code, n)
            cited = sorted(float(v) for v in re.findall(r"([0-9.]+) eV", table[code]))
            assert cited == sorted(_parse_values(df.at[idx, "Experimental_Eg_eV"])), (df.at[idx, "Record_ID"], cited)
            df.at[idx, "Experimental_Eg_Source"] = f"{table[code]} as cited in {SRC_2019}"
        if basis == "6-31G(d)":
            if pd.isna(df.at[idx, "Eg_eV"]):
                assert "Te" in code, df.at[idx, "Record_ID"]
                df.at[idx, "Curation_Note"] = _append(df.at[idx, "Curation_Note"], NO_TE_BASIS)
            else:
                df.at[idx, "Curation_Note"] = _append(df.at[idx, "Curation_Note"], SINGLE_POINT)
        if solv == "yes":
            df.at[idx, "Curation_Note"] = _append(df.at[idx, "Curation_Note"], PCM_NOTE)
        earlier = EARLIER_WORK_2019.get((code, basis, solv)) if n == "polymer" else None
        if earlier:
            df.at[idx, "Curation_Note"] = _append(
                df.at[idx, "Curation_Note"],
                f"Value taken by the source article from earlier work of the same group ({earlier} of {SRC_2019}).")
        if n != "polymer":
            for key in ((code, basis, solv, "*"), (code, basis, solv, n)):
                if key in HOMO_LUMO_CORRECTIONS:
                    df.at[idx, "HOMO_eV"] = pd.NA
                    df.at[idx, "LUMO_eV"] = pd.NA
                    df.at[idx, "Curation_Note"] = _append(df.at[idx, "Curation_Note"], HOMO_LUMO_CORRECTIONS[key])
    controls = df.Chalcogen_Type.isna() & df.Scope_Flag.ne("Core_SSeTe")
    df["Chalcogen_Type"] = df["Chalcogen_Type"].astype("object")
    df.loc[controls, "Chalcogen_Type"] = "None (control)"
    dev = df.Record_ID.astype(str).str.startswith("EROL_")
    df.loc[dev & df.Solvent_or_Conditions.isna(), "Solvent_or_Conditions"] = "not stated in source"
    # one label format for the development collection: "D58,A1" (52 source labels) -> "D58A1"
    df["Molecule_Name"] = df["Molecule_Name"].astype("object")
    for idx in df.index[dev & df.Molecule_Name.astype(str).str.fullmatch(r"D\d+,A\d+")]:
        source_label = str(df.at[idx, "Molecule_Name"])
        df.at[idx, "Molecule_Name"] = source_label.replace(",", "")
        df.at[idx, "Curation_Note"] = _append(
            df.at[idx, "Curation_Note"], f"Source label {source_label} written as {df.at[idx, 'Molecule_Name']}.")
    return df


def main():
    df = pd.read_csv(MASTER, low_memory=False)
    out = annotate(df)
    out.to_csv(MASTER, index=False, compression={"method": "gzip", "mtime": 0})
    ext = out[out.Record_ID.str.startswith("HAKAN_")]
    print("solvent:", ext.Solvent_or_Conditions.value_counts().to_dict())
    print("experimental numeric:", int(out.Experimental_Eg_Min_eV.notna().sum()),
          "multi-valued:", int((out.Experimental_Eg_Min_eV < out.Experimental_Eg_Max_eV).sum()))
    print("records with curation note:", int(out.Curation_Note.fillna("").ne("").sum()))


if __name__ == "__main__":
    main()
