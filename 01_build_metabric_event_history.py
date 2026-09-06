#!/usr/bin/env python3
"""Build event-history and post-five-year CTS5 tables from METABRIC clinical data."""

from __future__ import annotations

import argparse
import hashlib
import json
import warnings
from pathlib import Path

import numpy as np
import pandas as pd

LANDMARK = 60.0

REQUIRED = {"patientId", "RFS_MONTHS", "RFS_STATUS", "OS_MONTHS", "VITAL_STATUS",
            "ER_IHC", "HER2_SNP6", "AGE_AT_DIAGNOSIS", "TUMOR_SIZE", "GRADE",
            "LYMPH_NODES_EXAMINED_POSITIVE", "NPI", "HORMONE_THERAPY",
            "CHEMOTHERAPY", "RADIO_THERAPY", "INFERRED_MENOPAUSAL_STATE",}


def cli() -> argparse.Namespace:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--patient", required=True, type=Path)
    parser.add_argument("--sample", required=True, type=Path)
    parser.add_argument("--output-dir", required=True, type=Path)
    parser.add_argument("--overwrite", action="store_true")
    return parser.parse_args()


def sha256(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as handle:
        for block in iter(lambda: handle.read(1024 * 1024), b""):
            digest.update(block)
    return digest.hexdigest()


def require(data: pd.DataFrame, columns: set[str], label: str) -> None:
    missing = sorted(columns - set(data.columns))
    if missing:
        raise ValueError(f"{label} missing columns: {', '.join(missing)}")


def one_row_per_patient(data: pd.DataFrame, label: str) -> None:
    invalid = data["patientId"].isna().any() or data["patientId"].duplicated().any()
    if invalid:
        raise ValueError(f"{label} must contain one nonmissing row per patientId")


def rfs_event(value: object) -> float:
    if pd.isna(value):
        return np.nan
    value = str(value).strip().upper()
    if value.startswith("0:") or value in {"0", "NOT RECURRED"}:
        return 0.0
    if value.startswith("1:") or value in {"1", "RECURRED"}:
        return 1.0
    return np.nan


def er_positive(value: object) -> bool:
    if pd.isna(value):
        return False
    value = str(value).strip().upper()
    return "POS" in value or value in {"1", "TRUE", "YES"}


def her2_negative(value: object) -> bool:
    if pd.isna(value):
        return False
    value = str(value).strip().upper()
    return not ("GAIN" in value or "AMP" in value) and (
        "NEUTRAL" in value or "LOSS" in value)


def yes(value: object) -> bool:
    return not pd.isna(value) and str(value).strip().upper() in {
        "YES", "Y", "TRUE", "1", "POSITIVE", "RECEIVED"}


def node_category(value: object) -> float:
    if pd.isna(value):
        return np.nan
    value = float(value)
    if value < 0 or not value.is_integer():
        return np.nan
    if value == 0:
        return 0
    if value == 1:
        return 1
    if value <= 3:
        return 2
    if value <= 9:
        return 3
    return 4


def risk_group(score: object) -> object:
    if pd.isna(score):
        return np.nan
    if score < 3.13:
        return "LOW"
    if score <= 3.86:
        return "INTERMEDIATE"
    return "HIGH"


def load_merge(patient_file: Path, sample_file: Path) -> pd.DataFrame:
    patient = pd.read_csv(patient_file, low_memory=False)
    sample = pd.read_csv(sample_file, low_memory=False)

    require(patient, {"patientId"}, "patient table")
    require(sample, {"patientId"}, "sample table")
    one_row_per_patient(patient, "patient table")
    one_row_per_patient(sample, "sample table")

    sample_columns = [column
                      for column in sample.columns
                      if column == "patientId" or column not in patient.columns]

    data = patient.merge(sample[sample_columns], on="patientId", how="left",
                         validate="one_to_one",)
    require(data, REQUIRED, "merged table")
    return data


def event_history(data: pd.DataFrame) -> pd.DataFrame:
    data = data.copy()
    data["RFS_MONTHS_NUM"] = pd.to_numeric(data["RFS_MONTHS"], errors="coerce")
    data["OS_MONTHS_NUM"] = pd.to_numeric(data["OS_MONTHS"], errors="coerce")
    data["RFS_EVENT"] = data["RFS_STATUS"].apply(rfs_event)
    data["ER_POS_IHC"] = data["ER_IHC"].apply(er_positive)
    data["HER2_NEG_SNP6"] = data["HER2_SNP6"].apply(her2_negative)
    data["ERPOS_HER2NEG_FROZEN"] = data["ER_POS_IHC"] & data["HER2_NEG_SNP6"]

    valid_rfs = data["RFS_MONTHS_NUM"].notna() & data["RFS_EVENT"].notna()
    out = data.loc[data["ERPOS_HER2NEG_FROZEN"] & valid_rfs].copy()

    if out["RFS_MONTHS_NUM"].lt(0).any():
        raise ValueError("Negative RFS follow-up found")

    event = out["RFS_EVENT"].eq(1)
    censor = out["RFS_EVENT"].eq(0)
    time = out["RFS_MONTHS_NUM"]

    out["EVENT_AT_OR_BEFORE_5Y"] = (event & time.le(LANDMARK)).astype(int)
    out["EVENT_AFTER_5Y"] = (event & time.gt(LANDMARK)).astype(int)
    out["CENSORED_AT_OR_BEFORE_5Y"] = (censor & time.le(LANDMARK)).astype(int)
    out["REACHED_5Y_LANDMARK"] = time.gt(LANDMARK).astype(int)

    out["CONTROL_10Y"] = (censor & time.ge(120)).astype(int)
    out["CONTROL_15Y"] = (censor & time.ge(180)).astype(int)
    out["CONTROL_20Y"] = (censor & time.ge(240)).astype(int)

    out["EVENT_WINDOW"] = np.select(
        [event & time.le(60), event & time.gt(60) & time.le(120),
         event & time.gt(120) & time.le(180), event & time.gt(180),],
        ["EARLY_0_5Y", "LATE_5_10Y", "VERY_LATE_10_15Y", "EXTREME_LATE_GT15Y",],
        default="NO_RFS_EVENT",)

    nesting_ok = (out["CONTROL_20Y"].le(out["CONTROL_15Y"]) 
                  & out["CONTROL_15Y"].le(out["CONTROL_10Y"])).all()

    if not nesting_ok:
        raise AssertionError("Control nesting C20 subset C15 subset C10 failed")

    return out


def landmark_cts5(master: pd.DataFrame) -> pd.DataFrame:
    out = master.loc[master["REACHED_5Y_LANDMARK"].eq(1)].copy()
    out["POST5_TIME_MONTHS"] = out["RFS_MONTHS_NUM"] - LANDMARK

    event = out["RFS_EVENT"].eq(1)
    competing = (~event
                 & out["VITAL_STATUS"].astype(str).str.contains("OTHER",
                                                                case=False,
                                                                na=False,))

    out["POST5_EVENT_TYPE"] = np.select([event, competing], [1, 2], default=0,).astype(int)

    difference = out["OS_MONTHS_NUM"] - out["RFS_MONTHS_NUM"]

    if (competing & difference.abs().gt(1)).any():
        warnings.warn("Some competing-death OS and RFS times differ by >1 month")

    mapping = {"CTS5_AGE_YEARS": "AGE_AT_DIAGNOSIS", "CTS5_TUMOR_SIZE_MM_RAW": "TUMOR_SIZE",
               "CTS5_GRADE": "GRADE", "CTS5_POSITIVE_NODES_RAW": "LYMPH_NODES_EXAMINED_POSITIVE",}

    for target, source in mapping.items():
        out[target] = pd.to_numeric(out[source], errors="coerce")

    observed = out[list(mapping)].dropna()

    invalid = (~np.isfinite(observed).all(axis=1)
               | observed["CTS5_AGE_YEARS"].le(0)
               | observed["CTS5_TUMOR_SIZE_MM_RAW"].le(0)
               | ~observed["CTS5_GRADE"].isin([1, 2, 3])
               | observed["CTS5_POSITIVE_NODES_RAW"].lt(0)
               | ~np.isclose(observed["CTS5_POSITIVE_NODES_RAW"] % 1, 0))

    if invalid.any():
        raise ValueError(f"Invalid CTS5 values in {int(invalid.sum())} complete rows")

    out["CTS5_TUMOR_SIZE_MM_CAPPED30"] = (out["CTS5_TUMOR_SIZE_MM_RAW"].clip(upper=30))
    out["CTS5_NODE_CATEGORY"] = (out["CTS5_POSITIVE_NODES_RAW"].apply(node_category))

    inputs = ["CTS5_AGE_YEARS", "CTS5_TUMOR_SIZE_MM_RAW", "CTS5_GRADE", "CTS5_POSITIVE_NODES_RAW",
              "CTS5_NODE_CATEGORY",]

    out["CTS5_COMPLETE_CASE"] = (out[inputs].notna().all(axis=1)).astype(int)

    size = out["CTS5_TUMOR_SIZE_MM_CAPPED30"]
    inside = (0.093 * size - 0.001 * size.pow(2) + 0.375 * out["CTS5_GRADE"] + 0.017 * out["CTS5_AGE_YEARS"])

    out["CTS5_SCORE"] = np.where(out["CTS5_COMPLETE_CASE"].eq(1),
                                 0.438 * out["CTS5_NODE_CATEGORY"] + 0.988 * inside, np.nan,)

    out["CTS5_RISK_GROUP_LITERATURE"] = (out["CTS5_SCORE"].apply(risk_group))

    out["CTS5_STRICT_SCOPE_PROXY"] = (out["CTS5_COMPLETE_CASE"].eq(1) 
                                      & out["INFERRED_MENOPAUSAL_STATE"].astype(str)
                                      .str.upper()
                                      .str.startswith("POST")
                                      & out["HORMONE_THERAPY"].apply(yes)).astype(int)

    return out


def counts(master: pd.DataFrame, landmark: pd.DataFrame,) -> dict[str, int]:
    window = master["EVENT_WINDOW"]
    group = landmark["CTS5_RISK_GROUP_LITERATURE"]

    return {
        "full_cohort": len(master),
        "early_events": int(master["EVENT_AT_OR_BEFORE_5Y"].sum()),
        "censored_le5y": int(master["CENSORED_AT_OR_BEFORE_5Y"].sum()),
        "landmark5y": len(landmark),
        "late_events": int(master["EVENT_AFTER_5Y"].sum()),
        "late_5_10y": int(window.eq("LATE_5_10Y").sum()),
        "very_late_10_15y": int(window.eq("VERY_LATE_10_15Y").sum()),
        "extreme_late_gt15y": int(window.eq("EXTREME_LATE_GT15Y").sum()),
        "control_10y": int(master["CONTROL_10Y"].sum()),
        "control_15y": int(master["CONTROL_15Y"].sum()),
        "control_20y": int(master["CONTROL_20Y"].sum()),
        "competing_deaths": int(landmark["POST5_EVENT_TYPE"].eq(2).sum()),
        "other_censor": int(landmark["POST5_EVENT_TYPE"].eq(0).sum()),
        "cts5_complete": int(landmark["CTS5_COMPLETE_CASE"].sum()),
        "cts5_low": int(group.eq("LOW").sum()),
        "cts5_intermediate": int(group.eq("INTERMEDIATE").sum()),
        "cts5_high": int(group.eq("HIGH").sum()),
        "strict_scope": int(landmark["CTS5_STRICT_SCOPE_PROXY"].sum()),}


def main() -> None:
    args = cli()

    for path in (args.patient, args.sample):
        if not path.is_file():
            raise FileNotFoundError(path)

    args.output_dir.mkdir(parents=True, exist_ok=True)

    outputs = {
        "event_history":
            args.output_dir / "metabric_event_history.csv",
        "cts5_landmark":
            args.output_dir / "metabric_cts5_landmark5y.csv",
        "summary":
            args.output_dir / "cohort_summary.csv",
        "manifest":
            args.output_dir / "manifest.json",}

    if not args.overwrite and any(
        path.exists() for path in outputs.values()):
        raise FileExistsError(
            "Outputs exist; use --overwrite or a fresh directory")

    master = event_history(load_merge(args.patient, args.sample))
    landmark = landmark_cts5(master)
    observed = counts(master, landmark)

    master.to_csv(outputs["event_history"], index=False)
    landmark.to_csv(outputs["cts5_landmark"], index=False)

    pd.DataFrame(observed.items(), columns=["metric", "n"],
                ).to_csv(outputs["summary"], index=False,)

    manifest = {
        "inputs": {
            "patient": {
                "file": args.patient.name,
                "sha256": sha256(args.patient),},
            "sample": {
                "file": args.sample.name,
                "sha256": sha256(args.sample),},},
        "builder_sha256": sha256(Path(__file__).resolve()),
        "landmark_months": LANDMARK,
        "outputs": {
            key: {"file": path.name, "sha256": sha256(path),}
            for key, path in outputs.items()
            if key != "manifest"},}

    outputs["manifest"].write_text(json.dumps(manifest, indent=2), encoding="utf-8",)

    print(f"Event-history tables written to: {args.output_dir}")
    print(pd.DataFrame(observed.items(), columns=["metric", "n"],).to_string(index=False))

if __name__ == "__main__":
    main()
