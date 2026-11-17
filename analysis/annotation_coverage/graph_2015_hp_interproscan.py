#!/usr/bin/env python3
"""Reproduce the 2,015-hypothetical-protein InterProScan analysis.

Inputs
------
1. hypothetical_proteins_manifest_with_counts.csv
   Required columns: gene, seq_length, n_isolates
2. hypothetical_proteins_combined.faa.tsv
   Standard 15-column InterProScan TSV (no header row)

The script deliberately uses two definitions of a "match":

* all_signature_matches: every InterProScan TSV row for a protein.
* unique_interpro_entries: distinct, nonblank IPR accessions for a protein.

Proteins that never appear in the InterProScan TSV are retained with zero
matches by left-joining the counts onto the complete manifest.

Example
-------
python3 graph_2015_hp_interproscan.py \
    hypothetical_proteins_manifest_with_counts.csv \
    hypothetical_proteins_combined.faa.tsv \
    --outdir hp_analysis_outputs

Dependencies
------------
python3 -m pip install pandas numpy scipy matplotlib seaborn
"""

from __future__ import annotations

import argparse
import json
from pathlib import Path

import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
import seaborn as sns
from scipy import stats


IPR_COLUMNS = [
    "protein_id",
    "md5",
    "tsv_length",
    "member_database",
    "signature_accession",
    "signature_description",
    "start",
    "stop",
    "score",
    "status",
    "run_date",
    "interpro_accession",
    "interpro_description",
    "go_terms",
    "pathways",
]


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(
        description="Graph InterProScan matches across all 2,015 hypothetical proteins."
    )
    parser.add_argument(
        "manifest",
        type=Path,
        help="CSV containing gene, seq_length, and n_isolates.",
    )
    parser.add_argument(
        "interpro_tsv",
        type=Path,
        help="Headerless 15-column InterProScan TSV.",
    )
    parser.add_argument(
        "--outdir",
        type=Path,
        default=Path("hp_analysis_outputs"),
        help="Output directory (default: hp_analysis_outputs).",
    )
    parser.add_argument(
        "--expected-proteins",
        type=int,
        default=2015,
        help="Expected number of manifest proteins (default: 2015).",
    )
    parser.add_argument(
        "--bootstrap-resamples",
        type=int,
        default=3000,
        help="Bootstrap samples for each Spearman 95%% CI (default: 3000).",
    )
    parser.add_argument(
        "--seed",
        type=int,
        default=20260810,
        help="Random seed for reproducible bootstrap CIs (default: 20260810).",
    )
    return parser.parse_args()


def bootstrap_spearman(
    x: np.ndarray,
    y: np.ndarray,
    *,
    seed: int,
    n_boot: int,
) -> list[float]:
    """Return the percentile-bootstrap 95% CI for Spearman's rho."""
    rng = np.random.default_rng(seed)
    n = len(x)
    values = np.empty(n_boot, dtype=float)
    for i in range(n_boot):
        idx = rng.integers(0, n, n)
        values[i] = stats.spearmanr(x[idx], y[idx]).statistic
    return np.quantile(values, [0.025, 0.975]).tolist()


def association_stats(
    x: pd.Series,
    y: pd.Series,
    *,
    seed: int,
    n_boot: int,
) -> dict[str, object]:
    """Quantify monotonic association and a descriptive log-log linear fit."""
    x_arr = x.to_numpy(dtype=float)
    y_arr = y.to_numpy(dtype=float)

    spearman = stats.spearmanr(x_arr, y_arr)
    log_fit = stats.linregress(np.log1p(x_arr), np.log1p(y_arr))

    q1_cut = float(x.quantile(0.25))
    q3_cut = float(x.quantile(0.75))
    low = y[x <= q1_cut]
    high = y[x >= q3_cut]

    return {
        "n": int(len(x)),
        "spearman_rho": float(spearman.statistic),
        "spearman_p": float(spearman.pvalue),
        "spearman_ci95": bootstrap_spearman(
            x_arr, y_arr, seed=seed, n_boot=n_boot
        ),
        "log_log_slope": float(log_fit.slope),
        "log_log_r_squared": float(log_fit.rvalue**2),
        "log_log_p": float(log_fit.pvalue),
        "q1_cutoff": q1_cut,
        "q3_cutoff": q3_cut,
        "low_quartile_n": int(len(low)),
        "high_quartile_n": int(len(high)),
        "low_quartile_mean_matches": float(low.mean()),
        "high_quartile_mean_matches": float(high.mean()),
        "low_quartile_median_matches": float(low.median()),
        "high_quartile_median_matches": float(high.median()),
        "low_quartile_any_match_pct": float((low > 0).mean() * 100),
        "high_quartile_any_match_pct": float((high > 0).mean() * 100),
    }


def load_and_validate(
    manifest_path: Path,
    tsv_path: Path,
    expected_proteins: int,
) -> tuple[pd.DataFrame, pd.DataFrame]:
    """Read the two inputs and stop early if their structure is unexpected."""
    manifest = pd.read_csv(manifest_path)
    missing_columns = {"gene", "seq_length", "n_isolates"} - set(manifest.columns)
    if missing_columns:
        raise ValueError(f"Manifest is missing columns: {sorted(missing_columns)}")

    ipr = pd.read_csv(tsv_path, sep="\t", names=IPR_COLUMNS, dtype=str)

    manifest["gene"] = manifest["gene"].astype(str).str.strip()
    manifest["seq_length"] = pd.to_numeric(
        manifest["seq_length"], errors="raise"
    ).astype(int)
    manifest["n_isolates"] = pd.to_numeric(
        manifest["n_isolates"], errors="raise"
    ).astype(int)
    ipr["protein_id"] = ipr["protein_id"].astype(str).str.strip()
    ipr["tsv_length"] = pd.to_numeric(ipr["tsv_length"], errors="raise").astype(int)

    if len(manifest) != expected_proteins:
        raise ValueError(
            f"Expected {expected_proteins:,} manifest rows; found {len(manifest):,}."
        )
    if manifest["gene"].duplicated().any():
        examples = manifest.loc[
            manifest["gene"].duplicated(False), "gene"
        ].tolist()[:10]
        raise ValueError(f"Manifest protein IDs are not unique. Examples: {examples}")

    missing_ids = sorted(set(ipr["protein_id"]) - set(manifest["gene"]))
    if missing_ids:
        raise ValueError(
            "TSV contains protein IDs absent from the manifest. "
            f"Examples: {missing_ids[:10]}"
        )

    return manifest, ipr


def build_analysis(
    manifest: pd.DataFrame,
    ipr: pd.DataFrame,
) -> tuple[pd.DataFrame, pd.DataFrame, pd.DataFrame, pd.DataFrame, pd.DataFrame]:
    """Create one row per protein, including explicit zero-match proteins."""
    # Broad definition: every match row returned by a member database.
    raw_counts = ipr.groupby("protein_id").size().rename("all_signature_matches")

    # Strict definition: unique curated InterPro IDs only.
    ipr_clean = ipr[
        ipr["interpro_accession"].notna()
        & ~ipr["interpro_accession"].isin(["-", ""])
    ].copy()
    integrated_pairs = (
        ipr_clean[["protein_id", "interpro_accession", "interpro_description"]]
        .drop_duplicates(["protein_id", "interpro_accession"])
        .sort_values(["protein_id", "interpro_accession"])
    )
    integrated_counts = (
        integrated_pairs.groupby("protein_id")["interpro_accession"]
        .nunique()
        .rename("unique_interpro_entries")
    )
    database_counts = (
        ipr.groupby("protein_id")["member_database"]
        .nunique()
        .rename("member_databases_with_matches")
    )

    # This left join is crucial: the manifest contains all 2,015 proteins.
    analysis = manifest.rename(columns={"gene": "protein_id"}).copy()
    analysis = analysis.merge(raw_counts, left_on="protein_id", right_index=True, how="left")
    analysis = analysis.merge(
        integrated_counts, left_on="protein_id", right_index=True, how="left"
    )
    analysis = analysis.merge(
        database_counts, left_on="protein_id", right_index=True, how="left"
    )
    count_columns = [
        "all_signature_matches",
        "unique_interpro_entries",
        "member_databases_with_matches",
    ]
    analysis[count_columns] = analysis[count_columns].fillna(0).astype(int)
    analysis["any_signature_match"] = (
        analysis["all_signature_matches"] > 0
    ).astype(int)
    analysis["any_integrated_interpro_match"] = (
        analysis["unique_interpro_entries"] > 0
    ).astype(int)

    tsv_lengths = ipr[["protein_id", "tsv_length"]].drop_duplicates()
    length_check = analysis[["protein_id", "seq_length"]].merge(
        tsv_lengths, on="protein_id", how="inner"
    )
    mismatches = length_check[
        length_check["seq_length"] != length_check["tsv_length"]
    ]
    if not mismatches.empty:
        raise ValueError(
            f"Found {len(mismatches)} proteins whose manifest and TSV lengths disagree."
        )

    # Ten equally sized groups ranked by protein length.
    length_bins = (
        analysis.assign(
            length_decile=pd.qcut(
                analysis["seq_length"].rank(method="first"), 10, labels=False
            )
            + 1
        )
        .groupby("length_decile", observed=True)
        .agg(
            proteins=("protein_id", "size"),
            min_length_aa=("seq_length", "min"),
            median_length_aa=("seq_length", "median"),
            max_length_aa=("seq_length", "max"),
            mean_signature_matches=("all_signature_matches", "mean"),
            median_signature_matches=("all_signature_matches", "median"),
            pct_any_signature_match=("any_signature_match", lambda s: 100 * s.mean()),
            mean_interpro_entries=("unique_interpro_entries", "mean"),
            median_interpro_entries=("unique_interpro_entries", "median"),
            pct_any_interpro_entry=(
                "any_integrated_interpro_match", lambda s: 100 * s.mean()
            ),
        )
        .reset_index()
    )

    # Logarithmic isolate-count groups; these are bins, not deciles.
    isolate_edges = [0, 1, 2, 4, 8, 16, 32, 64, 128, 256, 512, 1024, np.inf]
    isolate_labels = [
        "1",
        "2",
        "3-4",
        "5-8",
        "9-16",
        "17-32",
        "33-64",
        "65-128",
        "129-256",
        "257-512",
        "513-1,024",
        "1,025+",
    ]
    isolate_bins = (
        analysis.assign(
            isolate_bin=pd.cut(
                analysis["n_isolates"],
                bins=isolate_edges,
                labels=isolate_labels,
                include_lowest=True,
            )
        )
        .groupby("isolate_bin", observed=True)
        .agg(
            proteins=("protein_id", "size"),
            min_isolates=("n_isolates", "min"),
            median_isolates=("n_isolates", "median"),
            max_isolates=("n_isolates", "max"),
            mean_signature_matches=("all_signature_matches", "mean"),
            pct_any_signature_match=("any_signature_match", lambda s: 100 * s.mean()),
            mean_interpro_entries=("unique_interpro_entries", "mean"),
            pct_any_interpro_entry=(
                "any_integrated_interpro_match", lambda s: 100 * s.mean()
            ),
        )
        .reset_index()
    )

    max_count = int(
        max(
            analysis["all_signature_matches"].max(),
            analysis["unique_interpro_entries"].max(),
        )
    )
    dist_index = pd.Index(range(max_count + 1), name="matches_per_protein")
    distribution = pd.DataFrame(index=dist_index)
    distribution["proteins_all_signature_matches"] = (
        analysis["all_signature_matches"]
        .value_counts()
        .reindex(dist_index, fill_value=0)
        .sort_index()
    )
    distribution["proteins_unique_interpro_entries"] = (
        analysis["unique_interpro_entries"]
        .value_counts()
        .reindex(dist_index, fill_value=0)
        .sort_index()
    )

    return (
        analysis,
        integrated_pairs,
        length_bins,
        isolate_bins,
        distribution.reset_index(),
    )


def make_summary(
    analysis: pd.DataFrame,
    ipr: pd.DataFrame,
    integrated_pairs: pd.DataFrame,
    *,
    seed: int,
    n_boot: int,
) -> dict[str, object]:
    association_inputs = {
        "length_vs_signature_matches": ("seq_length", "all_signature_matches"),
        "length_vs_interpro_entries": ("seq_length", "unique_interpro_entries"),
        "isolates_vs_signature_matches": ("n_isolates", "all_signature_matches"),
        "isolates_vs_interpro_entries": ("n_isolates", "unique_interpro_entries"),
    }
    associations = {}
    for offset, (name, (x_col, y_col)) in enumerate(
        association_inputs.items(), start=1
    ):
        associations[name] = association_stats(
            analysis[x_col],
            analysis[y_col],
            seed=seed + offset,
            n_boot=n_boot,
        )

    return {
        "quality_control": {
            "manifest_proteins": int(len(analysis)),
            "tsv_match_rows": int(len(ipr)),
            "proteins_with_signature_match": int(analysis["any_signature_match"].sum()),
            "proteins_without_signature_match": int(
                (analysis["any_signature_match"] == 0).sum()
            ),
            "proteins_with_integrated_interpro_entry": int(
                analysis["any_integrated_interpro_match"].sum()
            ),
            "proteins_without_integrated_interpro_entry": int(
                (analysis["any_integrated_interpro_match"] == 0).sum()
            ),
            "unique_integrated_protein_entry_pairs": int(len(integrated_pairs)),
        },
        "associations": associations,
        "definitions": {
            "all_signature_matches": (
                "Number of InterProScan TSV rows for the protein, including all "
                "member-database match instances."
            ),
            "unique_interpro_entries": (
                "Number of distinct nonblank InterPro accessions (IPR IDs) for the "
                "protein; repeated mappings to the same IPR ID count once."
            ),
        },
    }


def plot_six_panels(
    analysis: pd.DataFrame,
    length_bins: pd.DataFrame,
    isolate_bins: pd.DataFrame,
    summary: dict[str, object],
    outdir: Path,
) -> None:
    """Make the two distributions and four requested relationship plots."""
    sns.set_theme(style="whitegrid", context="talk")
    blue = "#2F6B8A"
    orange = "#D98032"
    dark = "#20262E"

    fig = plt.figure(figsize=(15, 15), constrained_layout=True)
    grid = fig.add_gridspec(3, 2, height_ratios=[0.95, 1.15, 1.15])

    distributions = [
        (
            fig.add_subplot(grid[0, 0]),
            "all_signature_matches",
            blue,
            "A. All member-database match instances",
        ),
        (
            fig.add_subplot(grid[0, 1]),
            "unique_interpro_entries",
            orange,
            "B. Unique curated InterPro entries",
        ),
    ]
    for ax, metric, color, title in distributions:
        counts = analysis[metric].value_counts().sort_index()
        ax.bar(counts.index, counts.values, color=color, width=0.82)
        ax.set_yscale("log")
        ax.set_xlabel("Matches per hypothetical protein")
        ax.set_ylabel("Number of proteins (log scale)")
        ax.set_title(title, loc="left", weight="bold")
        ax.text(
            0.98,
            0.96,
            f"n = {len(analysis):,}\nzero matches = {(analysis[metric] == 0).sum():,}",
            transform=ax.transAxes,
            ha="right",
            va="top",
            fontsize=12,
            color=dark,
        )

    metrics = [
        ("all_signature_matches", "All signature matches", "Blues"),
        ("unique_interpro_entries", "Unique InterPro entries", "Oranges"),
    ]
    for col, (metric, label, cmap) in enumerate(metrics):
        suffix = "signature_matches" if col == 0 else "interpro_entries"
        ax = fig.add_subplot(grid[1, col])
        hexbin = ax.hexbin(
            analysis["seq_length"],
            analysis[metric],
            gridsize=(42, 24),
            mincnt=1,
            bins="log",
            cmap=cmap,
            linewidths=0.15,
        )
        ax.plot(
            length_bins["median_length_aa"],
            length_bins[f"mean_{suffix}"],
            color=dark,
            marker="o",
            linewidth=2.2,
            label="Length-decile mean",
        )
        result = summary["associations"][f"length_vs_{suffix}"]
        ax.set_xlabel("Sequence length (amino acids)")
        ax.set_ylabel(f"{label} per protein")
        ax.set_title(
            f"{'C' if col == 0 else 'D'}. Matches versus sequence length",
            loc="left",
            weight="bold",
        )
        ax.text(
            0.98,
            0.96,
            f"Spearman rho = {result['spearman_rho']:.2f}\n"
            f"log-log R-squared = {result['log_log_r_squared']:.2f}",
            transform=ax.transAxes,
            ha="right",
            va="top",
            fontsize=12,
            color=dark,
        )
        ax.legend(frameon=False, loc="upper left")
        fig.colorbar(hexbin, ax=ax, label="Protein density (log count)")

    for col, (metric, label, cmap) in enumerate(metrics):
        suffix = "signature_matches" if col == 0 else "interpro_entries"
        ax = fig.add_subplot(grid[2, col])
        hexbin = ax.hexbin(
            analysis["n_isolates"],
            analysis[metric],
            xscale="log",
            gridsize=(42, 24),
            mincnt=1,
            bins="log",
            cmap=cmap,
            linewidths=0.15,
        )
        ax.plot(
            isolate_bins["median_isolates"],
            isolate_bins[f"mean_{suffix}"],
            color=dark,
            marker="o",
            linewidth=2.2,
            label="Isolate-bin mean",
        )
        result = summary["associations"][f"isolates_vs_{suffix}"]
        ax.set_xscale("log")
        ax.set_xlabel("Isolates containing the protein group (log scale)")
        ax.set_ylabel(f"{label} per protein")
        ax.set_title(
            f"{'E' if col == 0 else 'F'}. Matches versus number of sequences",
            loc="left",
            weight="bold",
        )
        ax.text(
            0.98,
            0.96,
            f"Spearman rho = {result['spearman_rho']:.2f}\n"
            f"log-log R-squared = {result['log_log_r_squared']:.2f}",
            transform=ax.transAxes,
            ha="right",
            va="top",
            fontsize=12,
            color=dark,
        )
        ax.legend(frameon=False, loc="upper left")
        fig.colorbar(hexbin, ax=ax, label="Protein density (log count)")

    fig.suptitle(
        "InterProScan matches across 2,015 S. pyogenes hypothetical proteins",
        fontsize=22,
        weight="bold",
    )
    fig.savefig(outdir / "hp_match_visualization.png", dpi=220, bbox_inches="tight")
    fig.savefig(outdir / "hp_match_visualization.pdf", bbox_inches="tight")
    plt.close(fig)


def main() -> None:
    args = parse_args()
    args.outdir.mkdir(parents=True, exist_ok=True)

    manifest, ipr = load_and_validate(
        args.manifest, args.interpro_tsv, args.expected_proteins
    )
    analysis, integrated_pairs, length_bins, isolate_bins, distribution = (
        build_analysis(manifest, ipr)
    )
    summary = make_summary(
        analysis,
        ipr,
        integrated_pairs,
        seed=args.seed,
        n_boot=args.bootstrap_resamples,
    )

    analysis.to_csv(args.outdir / "hp_match_analysis.csv", index=False)
    integrated_pairs.to_csv(
        args.outdir / "integrated_interpro_pairs.csv", index=False
    )
    length_bins.to_csv(args.outdir / "length_decile_summary.csv", index=False)
    isolate_bins.to_csv(args.outdir / "isolate_bin_summary.csv", index=False)
    distribution.to_csv(args.outdir / "match_count_distribution.csv", index=False)
    (args.outdir / "analysis_summary.json").write_text(
        json.dumps(summary, indent=2), encoding="utf-8"
    )

    plot_six_panels(analysis, length_bins, isolate_bins, summary, args.outdir)

    qc = summary["quality_control"]
    print(f"Analyzed {qc['manifest_proteins']:,} proteins.")
    print(f"InterProScan TSV rows: {qc['tsv_match_rows']:,}")
    print(
        "Proteins with at least one broad match: "
        f"{qc['proteins_with_signature_match']:,}"
    )
    print(
        "Proteins with at least one curated IPR entry: "
        f"{qc['proteins_with_integrated_interpro_entry']:,}"
    )
    print(f"Outputs written to: {args.outdir.resolve()}")


if __name__ == "__main__":
    main()
