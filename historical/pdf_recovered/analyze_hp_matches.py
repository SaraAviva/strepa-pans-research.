from __future__ import annotations

import json
from pathlib import Path

import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
import seaborn as sns
from scipy import stats


ROOT = Path("/workspace/scratch/84c78e3af748")
UPLOAD = ROOT / "upload"
OUT = ROOT / "analysis_outputs"
OUT.mkdir(exist_ok=True)

MANIFEST_PATH = UPLOAD / "hypothetical_proteins_manifest_with_counts.csv"
TSV_PATH = UPLOAD / "hypothetical_proteins_combined.faa.tsv"

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


def bootstrap_spearman(x: np.ndarray, y: np.ndarray, seed: int = 20260810, n_boot: int = 3000):
    rng = np.random.default_rng(seed)
    n = len(x)
    values = np.empty(n_boot, dtype=float)
    for i in range(n_boot):
        idx = rng.integers(0, n, n)
        values[i] = stats.spearmanr(x[idx], y[idx]).statistic
    return np.quantile(values, [0.025, 0.975]).tolist()


def association_stats(x: pd.Series, y: pd.Series, seed_offset: int = 0):
    x_arr = x.to_numpy(dtype=float)
    y_arr = y.to_numpy(dtype=float)
    sp = stats.spearmanr(x_arr, y_arr)
    log_fit = stats.linregress(np.log1p(x_arr), np.log1p(y_arr))

    q1_cut = float(x.quantile(0.25))
    q3_cut = float(x.quantile(0.75))
    low = y[x <= q1_cut]
    high = y[x >= q3_cut]

    return {
        "n": int(len(x)),
        "spearman_rho": float(sp.statistic),
        "spearman_p": float(sp.pvalue),
        "spearman_ci95": bootstrap_spearman(x_arr, y_arr, seed=20260810 + seed_offset),
        "log_log_slope": float(log_fit.slope),
        "log_log_r_squared": float(log_fit.rvalue**2),
        "log_log_p": float(log_fit.pvalue),
        "q1_cutoff": q1_cut,
        "q3_cutoff": q3_cut,
        "short_or_low_n": int(len(low)),
        "long_or_high_n": int(len(high)),
        "short_or_low_mean_matches": float(low.mean()),
        "long_or_high_mean_matches": float(high.mean()),
        "short_or_low_median_matches": float(low.median()),
        "long_or_high_median_matches": float(high.median()),
        "short_or_low_any_match_pct": float((low > 0).mean() * 100),
        "long_or_high_any_match_pct": float((high > 0).mean() * 100),
    }


manifest = pd.read_csv(MANIFEST_PATH)
ipr = pd.read_csv(TSV_PATH, sep="\t", names=IPR_COLUMNS, dtype=str)

manifest["gene"] = manifest["gene"].astype(str).str.strip()
manifest["seq_length"] = pd.to_numeric(manifest["seq_length"], errors="raise").astype(int)
manifest["n_isolates"] = pd.to_numeric(manifest["n_isolates"], errors="raise").astype(int)
ipr["protein_id"] = ipr["protein_id"].astype(str).str.strip()
ipr["tsv_length"] = pd.to_numeric(ipr["tsv_length"], errors="raise").astype(int)

if len(manifest) != 2015:
    raise ValueError(f"Expected 2,015 manifest rows; found {len(manifest):,}")
if manifest["gene"].duplicated().any():
    dupes = manifest.loc[manifest["gene"].duplicated(False), "gene"].tolist()[:10]
    raise ValueError(f"Manifest protein IDs are not unique. Examples: {dupes}")

missing_ids = sorted(set(ipr["protein_id"]) - set(manifest["gene"]))
if missing_ids:
    raise ValueError(f"TSV contains protein IDs absent from manifest. Examples: {missing_ids[:10]}")

raw_counts = ipr.groupby("protein_id").size().rename("all_signature_matches")

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

analysis = manifest.rename(columns={"gene": "protein_id"}).copy()
analysis = analysis.merge(raw_counts, left_on="protein_id", right_index=True, how="left")
analysis = analysis.merge(integrated_counts, left_on="protein_id", right_index=True, how="left")
analysis = analysis.merge(database_counts, left_on="protein_id", right_index=True, how="left")
for col in ["all_signature_matches", "unique_interpro_entries", "member_databases_with_matches"]:
    analysis[col] = analysis[col].fillna(0).astype(int)

analysis["any_signature_match"] = (analysis["all_signature_matches"] > 0).astype(int)
analysis["any_integrated_interpro_match"] = (analysis["unique_interpro_entries"] > 0).astype(int)
analysis["log10_n_isolates"] = np.log10(analysis["n_isolates"])

tsv_lengths = ipr[["protein_id", "tsv_length"]].drop_duplicates()
length_check = analysis[["protein_id", "seq_length"]].merge(tsv_lengths, on="protein_id", how="inner")
length_mismatches = length_check[length_check["seq_length"] != length_check["tsv_length"]]

analysis["length_quartile"] = pd.qcut(
    analysis["seq_length"], 4, labels=["Q1 shortest", "Q2", "Q3", "Q4 longest"], duplicates="drop"
).astype(str)
analysis["isolate_quartile"] = pd.qcut(
    analysis["n_isolates"].rank(method="first"),
    4,
    labels=["Q1 fewest", "Q2", "Q3", "Q4 most"],
).astype(str)

length_bins = (
    analysis.assign(length_decile=pd.qcut(analysis["seq_length"].rank(method="first"), 10, labels=False) + 1)
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
        pct_any_interpro_entry=("any_integrated_interpro_match", lambda s: 100 * s.mean()),
    )
    .reset_index()
)

isolate_edges = [0, 1, 2, 4, 8, 16, 32, 64, 128, 256, 512, 1024, np.inf]
isolate_labels = ["1", "2", "3–4", "5–8", "9–16", "17–32", "33–64", "65–128", "129–256", "257–512", "513–1,024", "1,025+"]
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
        pct_any_interpro_entry=("any_integrated_interpro_match", lambda s: 100 * s.mean()),
    )
    .reset_index()
)

max_count = int(max(analysis["all_signature_matches"].max(), analysis["unique_interpro_entries"].max()))
dist_index = pd.Index(range(max_count + 1), name="matches_per_protein")
match_distribution = pd.DataFrame(index=dist_index)
match_distribution["proteins_all_signature_matches"] = (
    analysis["all_signature_matches"].value_counts().reindex(dist_index, fill_value=0).sort_index()
)
match_distribution["proteins_unique_interpro_entries"] = (
    analysis["unique_interpro_entries"].value_counts().reindex(dist_index, fill_value=0).sort_index()
)
match_distribution = match_distribution.reset_index()

summary = {
    "quality_control": {
        "manifest_proteins": int(len(manifest)),
        "tsv_match_rows": int(len(ipr)),
        "proteins_with_signature_match": int(analysis["any_signature_match"].sum()),
        "proteins_without_signature_match": int((analysis["any_signature_match"] == 0).sum()),
        "proteins_with_integrated_interpro_entry": int(analysis["any_integrated_interpro_match"].sum()),
        "proteins_without_integrated_interpro_entry": int((analysis["any_integrated_interpro_match"] == 0).sum()),
        "unique_integrated_protein_entry_pairs": int(len(integrated_pairs)),
        "tsv_ids_missing_from_manifest": int(len(missing_ids)),
        "length_mismatches_for_matched_proteins": int(len(length_mismatches)),
    },
    "descriptive": {
        "length_min": int(analysis["seq_length"].min()),
        "length_median": float(analysis["seq_length"].median()),
        "length_max": int(analysis["seq_length"].max()),
        "isolates_min": int(analysis["n_isolates"].min()),
        "isolates_median": float(analysis["n_isolates"].median()),
        "isolates_max": int(analysis["n_isolates"].max()),
        "signature_matches_mean": float(analysis["all_signature_matches"].mean()),
        "signature_matches_median": float(analysis["all_signature_matches"].median()),
        "signature_matches_max": int(analysis["all_signature_matches"].max()),
        "interpro_entries_mean": float(analysis["unique_interpro_entries"].mean()),
        "interpro_entries_median": float(analysis["unique_interpro_entries"].median()),
        "interpro_entries_max": int(analysis["unique_interpro_entries"].max()),
    },
    "associations": {
        "length_vs_signature_matches": association_stats(analysis["seq_length"], analysis["all_signature_matches"], 1),
        "length_vs_interpro_entries": association_stats(analysis["seq_length"], analysis["unique_interpro_entries"], 2),
        "isolates_vs_signature_matches": association_stats(analysis["n_isolates"], analysis["all_signature_matches"], 3),
        "isolates_vs_interpro_entries": association_stats(analysis["n_isolates"], analysis["unique_interpro_entries"], 4),
    },
    "definitions": {
        "all_signature_matches": "Number of InterProScan TSV rows for the protein, including all member-database match instances.",
        "unique_interpro_entries": "Number of distinct non-dash InterPro accessions (IPR IDs) for the protein; multiple member databases mapping to the same IPR ID count once.",
    },
}

analysis.to_csv(OUT / "hp_match_analysis.csv", index=False)
ipr.to_csv(OUT / "interproscan_raw_with_headers.csv", index=False)
integrated_pairs.to_csv(OUT / "integrated_interpro_pairs.csv", index=False)
length_bins.to_csv(OUT / "length_decile_summary.csv", index=False)
isolate_bins.to_csv(OUT / "isolate_decile_summary.csv", index=False)
match_distribution.to_csv(OUT / "match_count_distribution.csv", index=False)
(OUT / "analysis_summary.json").write_text(json.dumps(summary, indent=2), encoding="utf-8")

workbook_payload = {
    "manifest": json.loads(manifest.to_json(orient="records")),
    "raw_interpro": json.loads(ipr.to_json(orient="records")),
    "integrated_pairs": json.loads(integrated_pairs.to_json(orient="records")),
    "analysis": json.loads(analysis.to_json(orient="records")),
    "length_bins": json.loads(length_bins.to_json(orient="records")),
    "isolate_bins": json.loads(isolate_bins.to_json(orient="records")),
    "match_distribution": json.loads(match_distribution.to_json(orient="records")),
    "summary": summary,
}
(OUT / "workbook_payload.json").write_text(json.dumps(workbook_payload), encoding="utf-8")


# Publication-ready figure: distributions and both requested relationships.
sns.set_theme(style="whitegrid", context="talk")
blue = "#2F6B8A"
orange = "#D98032"
dark = "#20262E"

fig = plt.figure(figsize=(15, 15), constrained_layout=True)
gs = fig.add_gridspec(3, 2, height_ratios=[0.95, 1.15, 1.15])

for ax, metric, color, title in [
    (fig.add_subplot(gs[0, 0]), "all_signature_matches", blue, "A. All member-database match instances"),
    (fig.add_subplot(gs[0, 1]), "unique_interpro_entries", orange, "B. Unique curated InterPro entries"),
]:
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

for col, (metric, color, label) in enumerate([
    ("all_signature_matches", blue, "All signature matches"),
    ("unique_interpro_entries", orange, "Unique InterPro entries"),
]):
    ax = fig.add_subplot(gs[1, col])
    hb = ax.hexbin(
        analysis["seq_length"],
        analysis[metric],
        gridsize=(42, 24),
        mincnt=1,
        bins="log",
        cmap="Blues" if col == 0 else "Oranges",
        linewidths=0.15,
    )
    ax.plot(length_bins["median_length_aa"], length_bins[f"mean_{'signature_matches' if col == 0 else 'interpro_entries'}"], color=dark, marker="o", linewidth=2.2, label="Decile mean")
    stat = summary["associations"][f"length_vs_{'signature_matches' if col == 0 else 'interpro_entries'}"]
    ax.set_xlabel("Sequence length (amino acids)")
    ax.set_ylabel(f"{label} per protein")
    ax.set_title(f"{'C' if col == 0 else 'D'}. Matches versus sequence length", loc="left", weight="bold")
    ax.text(0.98, 0.96, f"Spearman ρ = {stat['spearman_rho']:.2f}\nlog–log R² = {stat['log_log_r_squared']:.2f}", transform=ax.transAxes, ha="right", va="top", fontsize=12, color=dark)
    ax.legend(frameon=False, loc="upper left")
    fig.colorbar(hb, ax=ax, label="Protein density (log count)")

for col, (metric, color, label) in enumerate([
    ("all_signature_matches", blue, "All signature matches"),
    ("unique_interpro_entries", orange, "Unique InterPro entries"),
]):
    ax = fig.add_subplot(gs[2, col])
    hb = ax.hexbin(
        analysis["n_isolates"],
        analysis[metric],
        xscale="log",
        gridsize=(42, 24),
        mincnt=1,
        bins="log",
        cmap="Blues" if col == 0 else "Oranges",
        linewidths=0.15,
    )
    ax.plot(isolate_bins["median_isolates"], isolate_bins[f"mean_{'signature_matches' if col == 0 else 'interpro_entries'}"], color=dark, marker="o", linewidth=2.2, label="Decile mean")
    stat = summary["associations"][f"isolates_vs_{'signature_matches' if col == 0 else 'interpro_entries'}"]
    ax.set_xscale("log")
    ax.set_xlabel("Isolates containing the protein group (log scale)")
    ax.set_ylabel(f"{label} per protein")
    ax.set_title(f"{'E' if col == 0 else 'F'}. Matches versus number of sequences", loc="left", weight="bold")
    ax.text(0.98, 0.96, f"Spearman ρ = {stat['spearman_rho']:.2f}\nlog–log R² = {stat['log_log_r_squared']:.2f}", transform=ax.transAxes, ha="right", va="top", fontsize=12, color=dark)
    ax.legend(frameon=False, loc="upper left")
    fig.colorbar(hb, ax=ax, label="Protein density (log count)")

fig.suptitle(
    "InterProScan matches across 2,015 S. pyogenes hypothetical proteins",
    fontsize=22,
    weight="bold",
)
fig.savefig(OUT / "hp_match_visualization.png", dpi=220, bbox_inches="tight")
fig.savefig(OUT / "hp_match_visualization.pdf", bbox_inches="tight")
plt.close(fig)

print(json.dumps(summary, indent=2))
