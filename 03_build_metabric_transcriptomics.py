#!/usr/bin/env python3
"""Build the METABRIC transcriptomic dataset used in downstream models."""

from __future__ import annotations

import argparse
from pathlib import Path

import pandas as pd

def cli() -> argparse.Namespace:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--expression", required=True, type=Path)
    parser.add_argument("--event-history", required=True, type=Path)
    parser.add_argument("--output-dir", required=True, type=Path)
    return parser.parse_args()

def combine_symbols(values: pd.Series) -> str:
    symbols = {str(value).strip() for value in values if pd.notna(value) and str(value).strip()}
    return "|".join(sorted(symbols))

def build_gene_matrix(expression: pd.DataFrame) -> tuple[pd.DataFrame, pd.DataFrame]:
    id_columns = {"Hugo_Symbol", "Entrez_Gene_Id"}
    sample_columns = [column for column in expression.columns if column not in id_columns]
    expression["Entrez_Gene_Id"] = pd.to_numeric(expression["Entrez_Gene_Id"], errors="coerce",)
    expression = expression.dropna(subset=["Entrez_Gene_Id"]).copy()
    expression["Entrez_Gene_Id"] = (expression["Entrez_Gene_Id"].astype("int64"))
    values = expression[sample_columns].apply(pd.to_numeric, errors="coerce",)
    values.insert(0, "Entrez_Gene_Id", expression["Entrez_Gene_Id"])
    matrix = (values.groupby("Entrez_Gene_Id", sort=True).median())
    annotation = (expression.groupby("Entrez_Gene_Id", sort=True).agg(SOURCE_HUGO_SYMBOLS=("Hugo_Symbol", 
                                                                                           combine_symbols),
                                                                      N_SOURCE_ROWS=("Hugo_Symbol", "size"),)
                  .reset_index())
    return matrix, annotation

def select_samples(matrix: pd.DataFrame, event_history: pd.DataFrame,) -> pd.DataFrame:
  sample_ids = (event_history["sampleId"].dropna().astype(str).tolist())
  selected = matrix[sample_ids].T
  selected.index.name = "sampleId"
  selected.columns = [f"ENTREZ_{gene_id}" for gene_id in selected.columns]
  return selected.reset_index()


def sample_information(event_history: pd.DataFrame) -> pd.DataFrame:
    columns = ["sampleId", "patientId", "COHORT", "RFS_MONTHS_NUM", "RFS_EVENT", "EVENT_AT_OR_BEFORE_5Y",
               "EVENT_AFTER_5Y", "EARLY_EVENT_LE60", "LANDMARK60_ELIGIBLE", "REACHED_5Y_LANDMARK",
               "LATE_RFS_EVENT", "LATE_RFS_MONTHS_FROM_5Y", "CONTROL_10Y", "CONTROL_15Y", "CONTROL_20Y",
               "EVENT_WINDOW",]
    available = [column for column in columns if column in event_history.columns]
    return event_history[available].copy()


def main() -> None:
    args = cli()
    args.output_dir.mkdir(parents=True, exist_ok=True)
    expression = pd.read_csv(args.expression, sep="\t", low_memory=False,)
    event_history = pd.read_csv(args.event_history, low_memory=False,)
    gene_matrix, annotation = build_gene_matrix(expression)
    model_matrix = select_samples(gene_matrix, event_history)
    samples = sample_information(event_history)
    matrix_file = (args.output_dir / "metabric_expression_entrez_median.csv.gz")
    annotation_file = (args.output_dir / "metabric_gene_annotation.csv")
    samples_file = (args.output_dir / "metabric_transcriptomic_samples.csv")
    summary_file = (args.output_dir / "metabric_transcriptomics_summary.csv")
    model_matrix.to_csv(matrix_file, index=False, compression="gzip", float_format="%.10g",)
    annotation.to_csv(annotation_file, index=False)
    samples.to_csv(samples_file, index=False)
    summary = pd.DataFrame({"metric": ["source_expression_rows", "canonical_entrez_genes", "selected_samples",
                                       "missing_expression_cells",],
                            "value": [len(expression), len(annotation), len(model_matrix),
                                      int(model_matrix.iloc[:, 1:].isna().sum().sum()),],})
    summary.to_csv(summary_file, index=False)
    print(f"{len(model_matrix)} samples and "f"{len(annotation)} genes written to {args.output_dir}")

if __name__ == "__main__":
    main()
