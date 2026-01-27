import os

import numpy as np
import pandas as pd


def rpkm(counts: np.ndarray, lengths: np.ndarray) -> np.ndarray:
    """Calculate reads per kilobase transcript per million reads."""
    total_reads = np.sum(counts, axis=0)
    denominator = total_reads[np.newaxis, :] * lengths[:, np.newaxis]
    return np.divide(1e9 * counts, denominator, out=np.zeros_like(counts, dtype=float), where=denominator != 0)


def _log(message: str) -> None:
    if hasattr(snakemake, "log") and snakemake.log:
        os.makedirs(os.path.dirname(snakemake.log[0]), exist_ok=True)
        with open(snakemake.log[0], "a", encoding="utf-8") as log_handle:
            log_handle.write(f"{message}\n")


def main() -> None:
    input_path = snakemake.input[0]
    output_path = snakemake.output[0]
    os.makedirs(os.path.dirname(output_path), exist_ok=True)

    _log(f"Loading counts from {input_path}")
    counts_df = pd.read_csv(input_path, sep="\t", comment="#")

    if "Length" not in counts_df.columns:
        raise ValueError(f"Missing 'Length' column in {input_path}")

    length_index = counts_df.columns.get_loc("Length")
    sample_columns = list(counts_df.columns[length_index + 1 :])

    if not sample_columns:
        raise ValueError(f"No sample columns found in {input_path}")

    counts_matrix = counts_df[sample_columns].to_numpy(dtype=float)
    lengths = counts_df["Length"].to_numpy(dtype=float)

    _log("Calculating RPKM values")
    rpkm_matrix = rpkm(counts_matrix, lengths)

    for col_index, column in enumerate(sample_columns):
        counts_df[column] = rpkm_matrix[:, col_index].astype(float)

    _log(f"Writing RPKM output to {output_path}")
    counts_df.to_csv(output_path, sep="\t", index=False)
    _log("RPKM calculation complete")


if __name__ == "__main__":
    main()
