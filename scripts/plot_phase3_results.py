#!/usr/bin/env python3
"""Create the frozen Phase 3 multi-panel main and supplementary figures."""

from __future__ import annotations

import itertools
import shutil
from pathlib import Path

import matplotlib as mpl
import matplotlib.pyplot as plt
from matplotlib.lines import Line2D
import numpy as np
import pandas as pd
import seaborn as sns


ROOT = Path(__file__).resolve().parents[1]
RESULTS = ROOT / "results" / "phase3_robustness" / "v1"
FIG_DIR = ROOT / "manuscript" / "figures"
SUPP_FIG_DIR = ROOT / "manuscript" / "supplement" / "figures"
SOURCE_DIR = ROOT / "manuscript" / "figure_source_data" / "phase3"

PROGRAM_ORDER = [
    "SIG_ATNF_7GENE",
    "SIG_ATNF_ARJIS5",
    "SIG_ATNF_OSM1",
    "SIG_CROSS_UC4",
    "SIG_UST4_COEF",
    "SIG_UST_HUB4",
    "SIG_UST_MH2",
    "SIG_VDZ4",
    "SIG_VDZ_ANDO7",
]
PROGRAM_LABEL = {
    "SIG_ATNF_7GENE": "Anti-TNF7",
    "SIG_ATNF_ARJIS5": "Arijs5",
    "SIG_ATNF_OSM1": "OSM",
    "SIG_CROSS_UC4": "Cross-UC4",
    "SIG_UST4_COEF": "UST4-coef",
    "SIG_UST_HUB4": "UST-Hub4",
    "SIG_UST_MH2": "UST-MH2",
    "SIG_VDZ4": "VDZ4",
    "SIG_VDZ_ANDO7": "VDZ-Ando7",
}
TASK_ORDER = [
    "T_PURSUIT_CR_W6",
    "T_PROGECT_MH_W6",
    "T_EMTAB7604_ER",
    "T_UNITI2_CR_W8",
    "T_UNIFI_MH_W8",
    "T_UNIFI_REM_W8",
]
TASK_LABEL = {
    "T_PURSUIT_CR_W6": "PURSUIT\nresponse",
    "T_PROGECT_MH_W6": "PROgECT\nhealing",
    "T_EMTAB7604_ER": "E-MTAB-7604\nremission",
    "T_UNITI2_CR_W8": "UNITI-2\nresponse",
    "T_UNIFI_MH_W8": "UNIFI\nhealing",
    "T_UNIFI_REM_W8": "UNIFI\nremission",
}
ROLE_LABEL = {
    "SAME_THERAPY_SAME_ENDPOINT_EXTERNAL": "Strict external",
    "WITHIN_THERAPY_ENDPOINT_OR_DISEASE_TRANSPORT": "Within-therapy",
    "CROSS_THERAPY_TRANSPORT_STRESS_TEST": "Cross-therapy",
}
ROLE_COLOR = {
    "SAME_THERAPY_SAME_ENDPOINT_EXTERNAL": "#183A5A",
    "WITHIN_THERAPY_ENDPOINT_OR_DISEASE_TRANSPORT": "#2E75B6",
    "CROSS_THERAPY_TRANSPORT_STRESS_TEST": "#D9822B",
}
BLUE = "#2E75B6"
NAVY = "#183A5A"
ORANGE = "#D9822B"
RED = "#C44E52"
GREEN = "#4C956C"
GRAY = "#6B7280"
LIGHT = "#DCE6F2"


def setup_style() -> None:
    sns.set_theme(style="whitegrid", context="paper")
    mpl.rcParams.update(
        {
            "font.family": "Arial",
            "font.size": 8.5,
            "axes.labelsize": 9,
            "xtick.labelsize": 7.5,
            "ytick.labelsize": 7.5,
            "legend.fontsize": 7.2,
            "axes.linewidth": 0.8,
            "pdf.fonttype": 42,
            "ps.fonttype": 42,
        }
    )


def panel_letter(ax: plt.Axes, letter: str) -> None:
    ax.text(-0.12, 1.08, letter, transform=ax.transAxes, fontsize=13, fontweight="bold", va="top")


def save_figure(fig: plt.Figure, stem: str) -> None:
    FIG_DIR.mkdir(parents=True, exist_ok=True)
    png = FIG_DIR / f"{stem}.png"
    pdf = FIG_DIR / f"{stem}.pdf"
    fig.savefig(png, dpi=300, bbox_inches="tight", facecolor="white")
    fig.savefig(pdf, bbox_inches="tight", facecolor="white")
    plt.close(fig)


def save_supp_figure(fig: plt.Figure, stem: str) -> None:
    SUPP_FIG_DIR.mkdir(parents=True, exist_ok=True)
    fig.savefig(SUPP_FIG_DIR / f"{stem}.png", dpi=300, bbox_inches="tight", facecolor="white")
    fig.savefig(SUPP_FIG_DIR / f"{stem}.pdf", bbox_inches="tight", facecolor="white")
    plt.close(fig)


def heatmap_matrix(rows: pd.DataFrame, value: str) -> pd.DataFrame:
    matrix = rows.pivot(index="signature_id", columns="task_id", values=value)
    return matrix.reindex(index=PROGRAM_ORDER, columns=TASK_ORDER)


def plot_figure_1(evaluations: pd.DataFrame, tasks: pd.DataFrame, inventory: pd.DataFrame) -> None:
    fig, axes = plt.subplots(2, 2, figsize=(12.0, 8.2), gridspec_kw={"height_ratios": [0.9, 1.1]})
    ax = axes[0, 0]
    panel_letter(ax, "A")
    ax.axis("off")
    boxes = [
        (0.02, "Phase 0", "22 records\n18 clinical units\n836-subject\nlower bound", "#DCE9F5"),
        (0.28, "Frozen grid", "9 programs\n6 clinical units\n54 planned cells", "#E4F0D8"),
        (0.54, "Executed", "5 clinical units\n636 subjects\n6 tasks; 46 cells", "#FFF0C9"),
        (0.80, "Evidence strata", "3 strict\n14 within-therapy\n29 cross-therapy", "#F9DFCF"),
    ]
    for x, heading, body, color in boxes:
        rect = mpl.patches.FancyBboxPatch((x, 0.22), 0.18, 0.58, boxstyle="round,pad=0.012", fc=color, ec="#44546A", lw=1.1)
        ax.add_patch(rect)
        ax.text(x + 0.09, 0.67, heading, ha="center", va="center", fontweight="bold", fontsize=9)
        ax.text(x + 0.09, 0.43, body, ha="center", va="center", fontsize=7.0, linespacing=1.35)
    for x in [0.21, 0.47, 0.73]:
        ax.annotate("", xy=(x + 0.05, 0.51), xytext=(x, 0.51), arrowprops=dict(arrowstyle="-|>", lw=1.2, color="#44546A"))
    ax.text(0.5, 0.06, "Repository reuploads, nested subsets, repeated tissues/time points, and same-trial accessions were blocked.", ha="center", color=GRAY, fontsize=7.2)
    ax.set_xlim(0, 1)
    ax.set_ylim(0, 1)

    ax = axes[0, 1]
    panel_letter(ax, "B")
    task_rows = tasks.set_index("task_id").loc[TASK_ORDER].copy()
    y = np.arange(len(task_rows))
    ax.barh(y, task_rows["n_responder"], color=BLUE, label="Responder")
    ax.barh(y, task_rows["n_nonresponder"], left=task_rows["n_responder"], color=LIGHT, edgecolor="#7A8DA8", label="Nonresponder")
    ax.set_yticks(y, [TASK_LABEL[t].replace("\n", " ") for t in TASK_ORDER])
    ax.invert_yaxis()
    ax.set_xlabel("Subjects in endpoint-specific task")
    for index, total in enumerate((task_rows["n_responder"] + task_rows["n_nonresponder"]).astype(int)):
        ax.text(total + 4, index, str(total), va="center", fontsize=7)
    ax.legend(frameon=False, ncol=2, loc="lower center", bbox_to_anchor=(0.5, 1.01))
    sns.despine(ax=ax)

    ax = axes[1, 0]
    panel_letter(ax, "C")
    role_codes = {
        "SAME_THERAPY_SAME_ENDPOINT_EXTERNAL": 1,
        "WITHIN_THERAPY_ENDPOINT_OR_DISEASE_TRANSPORT": 2,
        "CROSS_THERAPY_TRANSPORT_STRESS_TEST": 3,
    }
    matrix = pd.DataFrame(np.nan, index=PROGRAM_ORDER, columns=TASK_ORDER)
    for row in evaluations.itertuples():
        matrix.loc[row.signature_id, row.task_id] = role_codes[row.evaluation_role]
    cmap = mpl.colors.ListedColormap([NAVY, BLUE, ORANGE])
    sns.heatmap(matrix, ax=ax, cmap=cmap, vmin=1, vmax=3, cbar=False, linewidths=0.6, linecolor="white", mask=matrix.isna())
    ax.set_xticklabels([TASK_LABEL[t] for t in TASK_ORDER], rotation=35, ha="right")
    ax.set_yticklabels([PROGRAM_LABEL[p] for p in PROGRAM_ORDER], rotation=0)
    ax.set_xlabel("")
    ax.set_ylabel("")
    handles = [mpl.patches.Patch(color=color, label=ROLE_LABEL[role]) for role, color in ROLE_COLOR.items()]
    ax.legend(handles=handles, frameon=False, ncol=3, loc="upper center", bbox_to_anchor=(0.5, -0.22))

    ax = axes[1, 1]
    panel_letter(ax, "D")
    program = inventory.set_index("signature_id").loc[PROGRAM_ORDER].copy()
    cell_counts = evaluations.groupby("signature_id").size().reindex(PROGRAM_ORDER).fillna(0)
    source_colors = []
    for therapy in program["therapy"]:
        if "ANTI_TNF" in therapy and "ANTI_INTEGRIN" not in therapy:
            source_colors.append(BLUE)
        elif "USTEKINUMAB" in therapy:
            source_colors.append(GREEN)
        elif "ANTI_INTEGRIN" in therapy and "ANTI_TNF" not in therapy:
            source_colors.append(ORANGE)
        else:
            source_colors.append("#8E6C8A")
    y = np.arange(len(PROGRAM_ORDER))
    ax.barh(y, program["n_genes"].astype(int), color=source_colors, alpha=0.85)
    for idx, (genes, count) in enumerate(zip(program["n_genes"].astype(int), cell_counts.astype(int))):
        ax.text(genes + 0.1, idx, f"{count} cells", va="center", fontsize=7)
    ax.set_yticks(y, [PROGRAM_LABEL[p] for p in PROGRAM_ORDER])
    ax.invert_yaxis()
    ax.set_xlabel("Published program size (genes)")
    legend = [
        mpl.patches.Patch(color=BLUE, label="Anti-TNF"),
        mpl.patches.Patch(color=GREEN, label="Ustekinumab"),
        mpl.patches.Patch(color=ORANGE, label="Vedolizumab"),
        mpl.patches.Patch(color="#8E6C8A", label="Cross-biologic"),
    ]
    ax.set_xlim(0, 8.4)
    ax.legend(handles=legend, frameon=False, ncol=4, loc="lower center", bbox_to_anchor=(0.5, 1.01),
              fontsize=6.6, columnspacing=0.8, handlelength=1.4)
    sns.despine(ax=ax)
    fig.tight_layout(h_pad=2.2, w_pad=2.0)
    save_figure(fig, "figure_1_multicohort_design_and_evidence_geometry")


def plot_figure_2(evaluations: pd.DataFrame) -> None:
    strict = evaluations[evaluations["evaluation_role"] == "SAME_THERAPY_SAME_ENDPOINT_EXTERNAL"].copy()
    strict["label"] = strict["signature_id"].map(PROGRAM_LABEL) + " → " + strict["task_id"].map(lambda x: TASK_LABEL[x].replace("\n", " "))
    strict = strict.sort_values("auroc")
    scores = pd.read_csv(RESULTS / "phase3_strict_subject_scores.csv")
    scores["cell"] = scores["signature_id"].map(PROGRAM_LABEL) + " → " + scores["task_id"].map(lambda x: TASK_LABEL[x].replace("\n", " "))
    rocs = pd.read_csv(RESULTS / "phase3_strict_roc_curves.csv")
    nulls = pd.read_csv(RESULTS / "phase3_strict_permutation_nulls.csv")

    fig, axes = plt.subplots(2, 2, figsize=(11.5, 8.5))
    ax = axes[0, 0]
    panel_letter(ax, "A")
    y = np.arange(len(strict))
    ax.errorbar(
        strict["auroc"], y,
        xerr=[strict["auroc"] - strict["auroc_ci_lower"], strict["auroc_ci_upper"] - strict["auroc"]],
        fmt="o", color=NAVY, ecolor="#7AA6C2", capsize=3, ms=6,
    )
    ax.axvline(0.5, ls="--", color=GRAY, lw=1)
    ax.set_yticks(y, strict["label"])
    ax.set_xlim(0.32, 0.83)
    ax.set_xlabel("AUROC (95% subject-bootstrap interval)")
    for idx, row in enumerate(strict.itertuples()):
        ax.text(0.805, idx, f"Holm P={row.holm_adjusted_p_value:.3f}", ha="right", va="center", fontsize=7)
    sns.despine(ax=ax)

    ax = axes[0, 1]
    panel_letter(ax, "B")
    order = strict["label"].tolist()
    sns.violinplot(data=scores, x="cell", y="score", hue="outcome_label", order=order, split=True, inner=None,
                   palette={0: "#B9C7DA", 1: BLUE}, linewidth=0.7, cut=0, ax=ax)
    sns.boxplot(data=scores, x="cell", y="score", hue="outcome_label", order=order, dodge=True, width=0.35,
                palette={0: "white", 1: "white"}, showfliers=False, linewidth=0.8, ax=ax, legend=False)
    ax.set_xticks(np.arange(len(order)), [label.split(" → ")[0] for label in order], rotation=25, ha="right")
    ax.set_xlabel("")
    ax.set_ylabel("Frozen uncalibrated score")
    handles, _ = ax.get_legend_handles_labels()
    ax.legend(handles[:2], ["Nonresponder", "Responder"], frameon=False, loc="upper left")
    sns.despine(ax=ax)

    ax = axes[1, 0]
    panel_letter(ax, "C")
    for (signature_id, task_id), group in rocs.groupby(["signature_id", "task_id"]):
        auc = float(strict.loc[(strict["signature_id"] == signature_id) & (strict["task_id"] == task_id), "auroc"].iloc[0])
        ax.plot(group["false_positive_rate"], group["true_positive_rate"], lw=1.8,
                label=f"{PROGRAM_LABEL[signature_id]} ({auc:.3f})")
    ax.plot([0, 1], [0, 1], ls="--", color=GRAY, lw=1)
    ax.set_aspect("equal", adjustable="box")
    ax.set_xlabel("False-positive rate")
    ax.set_ylabel("True-positive rate")
    ax.legend(frameon=False, title="AUROC", loc="lower right")
    sns.despine(ax=ax)

    ax = axes[1, 1]
    panel_letter(ax, "D")
    colors = [NAVY, BLUE, ORANGE]
    for offset, ((signature_id, task_id), group) in enumerate(nulls.groupby(["signature_id", "task_id"])):
        values = group["permuted_auroc"].to_numpy()
        hist, bins = np.histogram(values, bins=np.linspace(0.2, 0.8, 45), density=True)
        centers = (bins[:-1] + bins[1:]) / 2
        hist = hist / max(hist) * 0.55
        ax.fill_between(centers, offset, offset + hist, color=colors[offset], alpha=0.45)
        observed = float(strict.loc[(strict["signature_id"] == signature_id) & (strict["task_id"] == task_id), "auroc"].iloc[0])
        ax.vlines(observed, offset, offset + 0.62, color=colors[offset], lw=2)
        ax.text(0.205, offset + 0.18, PROGRAM_LABEL[signature_id], fontsize=7.2, va="center")
    ax.axvline(0.5, ls="--", color=GRAY, lw=1)
    ax.set_yticks([])
    ax.set_xlabel("Permutation-null AUROC; vertical line = observed")
    ax.set_ylim(-0.05, 2.75)
    sns.despine(ax=ax, left=True)
    fig.tight_layout(h_pad=2.1, w_pad=2.0)
    save_figure(fig, "figure_2_strict_external_evidence")


def plot_figure_3(evaluations: pd.DataFrame) -> None:
    fig, axes = plt.subplots(2, 2, figsize=(12.0, 9.0))
    matrices = [
        ("A", "auroc", "AUROC", "vlag", 0.35, 0.75, ".2f"),
        ("B", "ci_width", "95% bootstrap interval width", "mako_r", None, None, ".2f"),
        ("C", "prevalence_normalized_ap", "Average precision / prevalence", "crest", None, None, ".2f"),
    ]
    evaluations = evaluations.copy()
    evaluations["ci_width"] = evaluations["auroc_ci_upper"] - evaluations["auroc_ci_lower"]
    for ax, (letter, value, label, cmap, vmin, vmax, fmt) in zip([axes[0, 0], axes[0, 1], axes[1, 0]], matrices):
        panel_letter(ax, letter)
        matrix = heatmap_matrix(evaluations, value)
        sns.heatmap(matrix, ax=ax, cmap=cmap, vmin=vmin, vmax=vmax, center=0.5 if value == "auroc" else None,
                    annot=True, fmt=fmt, annot_kws={"fontsize": 6.2}, linewidths=0.4, linecolor="white",
                    mask=matrix.isna(), cbar_kws={"label": label, "shrink": 0.72})
        ax.set_xticklabels([TASK_LABEL[t] for t in TASK_ORDER], rotation=35, ha="right")
        ax.set_yticklabels([PROGRAM_LABEL[p] for p in PROGRAM_ORDER], rotation=0)
        ax.set_xlabel("")
        ax.set_ylabel("")

    ax = axes[1, 1]
    panel_letter(ax, "D")
    counts = (
        evaluations.assign(role=evaluations["evaluation_role"].map(ROLE_LABEL))
        .groupby(["role", "transport_classification"]).size().unstack(fill_value=0)
        .reindex(["Strict external", "Within-therapy", "Cross-therapy"])
    )
    supported = counts.get("SUPPORTED_DIRECTION", pd.Series(0, index=counts.index))
    uncertain = counts.get("UNCERTAIN_OR_NONTRANSPORTABLE", pd.Series(0, index=counts.index))
    y = np.arange(len(counts))
    ax.barh(y, supported, color=BLUE, label="Directionally supported")
    ax.barh(y, uncertain, left=supported, color=LIGHT, edgecolor="#7A8DA8", label="Uncertain/nontransportable")
    for idx, (s, u) in enumerate(zip(supported, uncertain)):
        ax.text(s / 2, idx, str(int(s)), ha="center", va="center", color="white", fontweight="bold")
        ax.text(s + u / 2, idx, str(int(u)), ha="center", va="center", color="#222", fontweight="bold")
    ax.set_yticks(y, counts.index)
    ax.invert_yaxis()
    ax.set_xlabel("Executable program-task cells")
    ax.legend(frameon=False, loc="lower right")
    sns.despine(ax=ax)
    fig.tight_layout(h_pad=2.2, w_pad=2.0)
    save_figure(fig, "figure_3_complete_transportability_atlas")


def plot_figure_4(cells: pd.DataFrame) -> None:
    np.random.seed(20260905)  # Horizontal jitter only; matches corrected publication figure.
    groups = pd.read_csv(RESULTS / "phase3_context_alignment_groups.csv")
    contrasts = pd.read_csv(RESULTS / "phase3_within_program_context_contrasts.csv")
    therapy = contrasts[contrasts["context_component"] == "same_therapy"].copy()
    contingency = pd.read_csv(RESULTS / "phase3_unifi_endpoint_contingency.csv")
    unifi = pd.read_csv(RESULTS / "phase3_unifi_paired_program_aurocs.csv")
    fig, axes = plt.subplots(2, 2, figsize=(10.8, 8.5))

    ax = axes[0, 0]
    panel_letter(ax, "A")
    sns.boxplot(data=cells, x="alignment_score", y="frozen_auroc", color=LIGHT, width=0.55, showfliers=False, ax=ax)
    sns.stripplot(data=cells, x="alignment_score", y="frozen_auroc", hue="evaluation_role", palette=ROLE_COLOR,
                  jitter=0.17, size=4.3, alpha=0.85, ax=ax)
    ax.axhline(0.5, ls="--", color=GRAY, lw=1)
    for row in groups.itertuples():
        ax.text(row.alignment_score, 0.785, f"{row.n_supported}/{row.n_cells}", ha="center", fontsize=7)
    ax.text(0.02, 0.98, "supported / total", transform=ax.transAxes, va="top", fontsize=7, color=GRAY)
    ax.set_xlabel("Context-alignment score (therapy + endpoint + disease)")
    ax.set_ylabel("AUROC")
    handles = [Line2D([], [], color=ROLE_COLOR[role], marker="o", linestyle="None", markersize=4, label=ROLE_LABEL[role]) for role in ROLE_LABEL]
    ax.legend(handles=handles, frameon=False, loc="lower right")
    sns.despine(ax=ax)

    ax = axes[0, 1]
    panel_letter(ax, "B")
    therapy = therapy.sort_values("aligned_minus_mismatched")
    program_colors = dict(zip(PROGRAM_ORDER, sns.color_palette("tab10", n_colors=len(PROGRAM_ORDER))))
    for row in therapy.itertuples():
        color = program_colors[row.signature_id]
        ax.plot([0, 1], [row.mismatched_median_auroc, row.aligned_median_auroc], color=color, lw=1.3,
                marker="o", ms=4.5, label=PROGRAM_LABEL[row.signature_id])
    ax.axhline(0.5, ls="--", color=GRAY, lw=1)
    ax.set_xticks([0, 1], ["Cross-therapy", "Same-therapy"])
    ax.set_xlim(-0.15, 1.15)
    ax.set_ylabel("Within-program median AUROC")
    ax.legend(frameon=False, ncol=4, loc="lower center", bbox_to_anchor=(0.5, 1.01), fontsize=6.2,
              columnspacing=0.8, handlelength=1.2)
    sns.despine(ax=ax)

    ax = axes[1, 0]
    panel_letter(ax, "C")
    matrix = contingency.pivot(index="mucosal_healing", columns="clinical_remission", values="n_subjects").reindex(index=[1, 0], columns=[0, 1])
    sns.heatmap(matrix, annot=True, fmt="d", cmap="Blues", cbar=False, linewidths=1.2, linecolor="white", ax=ax,
                annot_kws={"fontsize": 13, "fontweight": "bold"})
    ax.set_xticklabels(["No", "Yes"])
    ax.set_yticklabels(["Yes", "No"], rotation=0)
    ax.set_xlabel("Clinical remission at week 8")
    ax.set_ylabel("Mucosal healing at week 8")
    ax.text(0.02, -0.27, "358 shared subjects; agreement 94.4%; Cohen κ=0.772", transform=ax.transAxes, fontsize=7.5)

    ax = axes[1, 1]
    panel_letter(ax, "D")
    for row in unifi.itertuples():
        color = program_colors[row.signature_id]
        ax.scatter(row.mucosal_healing_auroc, row.clinical_remission_auroc, s=36, color=color,
                   label=PROGRAM_LABEL[row.signature_id])
    ax.plot([0.4, 0.75], [0.4, 0.75], ls="--", color=GRAY, lw=1)
    ax.axvline(0.5, ls=":", color="#A0A0A0", lw=0.8)
    ax.axhline(0.5, ls=":", color="#A0A0A0", lw=0.8)
    ax.set_xlim(0.40, 0.76)
    ax.set_ylim(0.40, 0.76)
    ax.set_xlabel("UNIFI mucosal-healing AUROC")
    ax.set_ylabel("UNIFI clinical-remission AUROC")
    ax.legend(frameon=False, ncol=2, loc="lower right", fontsize=6.1, columnspacing=0.8, handletextpad=0.3)
    sns.despine(ax=ax)
    fig.tight_layout(h_pad=2.4, w_pad=2.2)
    save_figure(fig, "figure_4_context_alignment_and_endpoint_drift")


def symmetric_pair_matrix(rows: pd.DataFrame, value: str, diagonal: float = 1.0) -> pd.DataFrame:
    matrix = pd.DataFrame(np.nan, index=PROGRAM_ORDER, columns=PROGRAM_ORDER)
    np.fill_diagonal(matrix.values, diagonal)
    for row in rows.itertuples():
        first = row.signature_1
        second = row.signature_2
        val = getattr(row, value)
        matrix.loc[first, second] = val
        matrix.loc[second, first] = val
    return matrix


def plot_figure_5() -> None:
    weights = pd.read_csv(RESULTS / "phase3_program_gene_weight_matrix.csv").set_index("signature_id").reindex(PROGRAM_ORDER)
    used = [column for column in weights.columns if (weights[column] != 0).any()]
    weights = weights[used]
    gene_pairs = pd.read_csv(RESULTS / "phase3_program_pair_gene_similarity.csv")
    score_pairs = pd.read_csv(RESULTS / "phase3_program_pair_median_score_correlations.csv")
    perf_pairs = pd.read_csv(RESULTS / "phase3_program_performance_profile_correlations.csv")
    fig, axes = plt.subplots(2, 2, figsize=(13.0, 9.2), gridspec_kw={"height_ratios": [1.05, 1]})

    ax = axes[0, 0]
    panel_letter(ax, "A")
    sns.heatmap(weights, cmap="vlag", center=0, vmin=-1, vmax=1, mask=weights.eq(0), cbar_kws={"label": "Normalized signed weight", "shrink": 0.7},
                linewidths=0.25, linecolor="#EEEEEE", ax=ax)
    ax.set_yticklabels([PROGRAM_LABEL[p] for p in PROGRAM_ORDER], rotation=0)
    ax.set_xticklabels(ax.get_xticklabels(), rotation=70, ha="right", fontsize=5.8)
    ax.set_xlabel("")
    ax.set_ylabel("")

    ax = axes[0, 1]
    panel_letter(ax, "B")
    jaccard = symmetric_pair_matrix(gene_pairs, "jaccard_similarity")
    sns.heatmap(jaccard, cmap="Blues", vmin=0, vmax=0.25, annot=True, fmt=".2f", annot_kws={"fontsize": 6},
                cbar_kws={"label": "Gene-set Jaccard similarity", "shrink": 0.7}, ax=ax)
    ax.set_xticklabels([PROGRAM_LABEL[p] for p in PROGRAM_ORDER], rotation=45, ha="right")
    ax.set_yticklabels([PROGRAM_LABEL[p] for p in PROGRAM_ORDER], rotation=0)
    ax.set_xlabel("")
    ax.set_ylabel("")

    ax = axes[1, 0]
    panel_letter(ax, "C")
    score_matrix = symmetric_pair_matrix(score_pairs, "median_spearman_rho")
    sns.heatmap(score_matrix, cmap="vlag", center=0, vmin=-1, vmax=1, annot=True, fmt=".2f", annot_kws={"fontsize": 6},
                mask=score_matrix.isna(), cbar_kws={"label": "Median subject-score Spearman ρ", "shrink": 0.7}, ax=ax)
    ax.set_xticklabels([PROGRAM_LABEL[p] for p in PROGRAM_ORDER], rotation=45, ha="right")
    ax.set_yticklabels([PROGRAM_LABEL[p] for p in PROGRAM_ORDER], rotation=0)
    ax.set_xlabel("")
    ax.set_ylabel("")

    ax = axes[1, 1]
    panel_letter(ax, "D")
    perf_matrix = symmetric_pair_matrix(perf_pairs, "performance_profile_spearman_rho")
    sns.heatmap(perf_matrix, cmap="vlag", center=0, vmin=-1, vmax=1, annot=True, fmt=".2f", annot_kws={"fontsize": 6},
                mask=perf_matrix.isna(), cbar_kws={"label": "Cross-task AUROC-profile Spearman ρ", "shrink": 0.7}, ax=ax)
    ax.set_xticklabels([PROGRAM_LABEL[p] for p in PROGRAM_ORDER], rotation=45, ha="right")
    ax.set_yticklabels([PROGRAM_LABEL[p] for p in PROGRAM_ORDER], rotation=0)
    ax.set_xlabel("")
    ax.set_ylabel("")
    fig.tight_layout(h_pad=2.0, w_pad=2.3)
    save_figure(fig, "figure_5_program_architecture_and_redundancy")


def plot_figure_6(cells: pd.DataFrame) -> None:
    fig, axes = plt.subplots(3, 2, figsize=(10.8, 12.2))

    for ax, letter, alt, label in [
        (axes[0, 0], "A", "mad_auroc", "Median/MAD AUROC"),
        (axes[0, 1], "B", "rank_auroc", "Centered-rank AUROC"),
    ]:
        panel_letter(ax, letter)
        for role, group in cells.groupby("evaluation_role"):
            ax.scatter(group["frozen_auroc"], group[alt], s=28, alpha=0.8, color=ROLE_COLOR[role], label=ROLE_LABEL[role])
        ax.plot([0.35, 0.78], [0.35, 0.78], ls="--", color=GRAY, lw=1)
        ax.axvline(0.5, ls=":", color="#AAAAAA", lw=0.8)
        ax.axhline(0.5, ls=":", color="#AAAAAA", lw=0.8)
        ax.set_xlim(0.35, 0.78)
        ax.set_ylim(0.35, 0.78)
        ax.set_xlabel("Frozen population-z AUROC")
        ax.set_ylabel(label)
        if letter == "A":
            ax.legend(frameon=False, loc="lower right")
        sns.despine(ax=ax)

    ax = axes[1, 0]
    panel_letter(ax, "C")
    for role, group in cells.groupby("evaluation_role"):
        ax.scatter(group["frozen_auroc"], group["bootstrap_probability_auc_gt_0_5"], color=ROLE_COLOR[role], s=30, alpha=0.82)
    ax.axvline(0.5, ls="--", color=GRAY, lw=1)
    ax.axhline(0.95, ls=":", color=GRAY, lw=1)
    ax.set_xlabel("Frozen AUROC")
    ax.set_ylabel("Bootstrap probability AUROC > 0.50")
    ax.set_ylim(-0.03, 1.03)
    sns.despine(ax=ax)

    ax = axes[1, 1]
    panel_letter(ax, "D")
    supported = cells[(cells["frozen_classification"] == "SUPPORTED_DIRECTION") & cells["loo_gene_min_auroc"].notna()].copy()
    supported["label"] = supported["signature_id"].map(PROGRAM_LABEL) + " | " + supported["task_id"].map(lambda x: TASK_LABEL[x].replace("\n", " "))
    supported = supported.sort_values("frozen_auroc")
    y = np.arange(len(supported))
    ax.hlines(y, supported["loo_gene_min_auroc"], supported["loo_gene_max_auroc"], color="#8CA6C1", lw=2)
    ax.scatter(supported["frozen_auroc"], y, color=NAVY, s=24, zorder=3, label="Full program")
    ax.axvline(0.5, ls="--", color=GRAY, lw=1)
    ax.set_yticks(y, supported["label"], fontsize=5.8)
    ax.set_xlabel("Leave-one-gene-out AUROC range")
    ax.legend(frameon=False, loc="lower right")
    sns.despine(ax=ax)

    ax = axes[2, 0]
    panel_letter(ax, "E")
    for role, group in cells.groupby("evaluation_role"):
        ax.scatter(group["n_subjects"], group["jackknife_max_abs_delta"], color=ROLE_COLOR[role], s=30, alpha=0.82, label=ROLE_LABEL[role])
    ax.set_xlabel("Subjects in task")
    ax.set_ylabel("Largest absolute leave-one-subject AUROC change")
    ax.set_xscale("log")
    ax.legend(frameon=False, loc="upper right")
    sns.despine(ax=ax)

    ax = axes[2, 1]
    panel_letter(ax, "F")
    for role, group in cells.groupby("evaluation_role"):
        ax.scatter(group["frozen_auroc"], group["published_orientation_percentile"], color=ROLE_COLOR[role], s=30, alpha=0.82)
    ax.axvline(0.5, ls="--", color=GRAY, lw=1)
    ax.axhline(50, ls=":", color=GRAY, lw=1)
    ax.axhline(75, ls=":", color="#A0A0A0", lw=0.8)
    ax.set_xlabel("Frozen AUROC")
    ax.set_ylabel("Published orientation percentile\nwithin same-gene sign space")
    ax.set_ylim(0, 100)
    sns.despine(ax=ax)
    fig.tight_layout(h_pad=2.2, w_pad=2.2)
    save_figure(fig, "figure_6_multilevel_robustness_diagnostics")


def plot_supplementary(evaluations: pd.DataFrame, cells: pd.DataFrame) -> None:
    # S1: complete cell-level forest.
    ordered = evaluations.sort_values(["evaluation_role", "auroc"]).reset_index(drop=True)
    fig, ax = plt.subplots(figsize=(9.2, 13.5))
    y = np.arange(len(ordered))
    for role, group in ordered.groupby("evaluation_role"):
        idx = group.index.to_numpy()
        ax.errorbar(group["auroc"], idx,
                    xerr=[group["auroc"] - group["auroc_ci_lower"], group["auroc_ci_upper"] - group["auroc"]],
                    fmt="o", color=ROLE_COLOR[role], ecolor=ROLE_COLOR[role], alpha=0.8, ms=3.3, capsize=1.6,
                    label=ROLE_LABEL[role])
    ax.axvline(0.5, ls="--", color=GRAY, lw=1)
    labels = [f"{PROGRAM_LABEL[r.signature_id]} | {TASK_LABEL[r.task_id].replace(chr(10), ' ')}" for r in ordered.itertuples()]
    ax.set_yticks(y, labels, fontsize=5.8)
    ax.invert_yaxis()
    ax.set_xlabel("AUROC (95% subject-bootstrap interval)")
    ax.legend(frameon=False, ncol=3, loc="lower center", bbox_to_anchor=(0.5, 1.005))
    sns.despine(ax=ax)
    fig.tight_layout()
    save_supp_figure(fig, "figure_s1_all_cell_forest")

    # S2: complete leave-one-gene-out delta heatmap.
    loo = pd.read_csv(RESULTS / "phase3_leave_one_gene_out.csv")
    loo["cell"] = loo["signature_id"].map(PROGRAM_LABEL) + " | " + loo["task_id"].map(lambda x: TASK_LABEL[x].replace("\n", " "))
    matrix = loo.pivot(index="cell", columns="omitted_gene", values="delta_auroc")
    matrix = matrix.loc[matrix.abs().max(axis=1).sort_values(ascending=False).index]
    fig, ax = plt.subplots(figsize=(12.5, 11.5))
    sns.heatmap(matrix, cmap="vlag", center=0, vmin=-0.20, vmax=0.20, mask=matrix.isna(),
                cbar_kws={"label": "AUROC change after omitting gene", "shrink": 0.65}, ax=ax)
    ax.set_xlabel("Omitted gene")
    ax.set_ylabel("")
    ax.set_xticklabels(ax.get_xticklabels(), rotation=70, ha="right", fontsize=6)
    ax.set_yticklabels(ax.get_yticklabels(), fontsize=5.6)
    fig.tight_layout()
    save_supp_figure(fig, "figure_s2_leave_one_gene_out_heatmap")

    # S3: task-specific score correlation matrices.
    rows = pd.read_csv(RESULTS / "phase3_task_score_correlations.csv")
    fig, axes = plt.subplots(2, 3, figsize=(13.2, 8.4))
    for ax, task_id in zip(axes.flat, TASK_ORDER):
        task_rows = rows[rows["task_id"] == task_id]
        matrix = symmetric_pair_matrix(task_rows, "spearman_rho")
        keep = [p for p in PROGRAM_ORDER if matrix.loc[p].notna().any()]
        matrix = matrix.loc[keep, keep]
        sns.heatmap(matrix, cmap="vlag", center=0, vmin=-1, vmax=1, cbar=False, square=True, ax=ax)
        ax.set_title(TASK_LABEL[task_id].replace("\n", " "), fontsize=8)
        ax.set_xticklabels([PROGRAM_LABEL[p] for p in keep], rotation=55, ha="right", fontsize=5.3)
        ax.set_yticklabels([PROGRAM_LABEL[p] for p in keep], rotation=0, fontsize=5.3)
    fig.tight_layout()
    save_supp_figure(fig, "figure_s3_task_specific_score_correlations")

    # S4: complete robustness matrices.
    fig, axes = plt.subplots(2, 2, figsize=(12.0, 8.7))
    fields = [
        ("published_orientation_percentile", "Published-orientation percentile", "viridis", 0, 100),
        ("bootstrap_probability_auc_gt_0_5", "Bootstrap probability AUROC > 0.50", "Blues", 0, 1),
        ("scaling_max_abs_delta", "Maximum AUROC change across scalings", "rocket_r", 0, None),
        ("loo_gene_direction_retention_fraction", "Leave-one-gene direction-retention fraction", "Greens", 0, 1),
    ]
    for ax, (field, label, cmap, vmin, vmax) in zip(axes.flat, fields):
        matrix = cells.pivot(index="signature_id", columns="task_id", values=field).reindex(index=PROGRAM_ORDER, columns=TASK_ORDER)
        sns.heatmap(matrix, cmap=cmap, vmin=vmin, vmax=vmax, mask=matrix.isna(), annot=True, fmt=".2f",
                    annot_kws={"fontsize": 5.8}, cbar_kws={"label": label, "shrink": 0.7}, ax=ax)
        ax.set_xticklabels([TASK_LABEL[t] for t in TASK_ORDER], rotation=35, ha="right")
        ax.set_yticklabels([PROGRAM_LABEL[p] for p in PROGRAM_ORDER], rotation=0)
        ax.set_xlabel("")
        ax.set_ylabel("")
    fig.tight_layout(h_pad=2.0, w_pad=2.0)
    save_supp_figure(fig, "figure_s4_complete_robustness_matrices")


def copy_source_data() -> None:
    SOURCE_DIR.mkdir(parents=True, exist_ok=True)
    sources = [
        ROOT / "results" / "performance" / "phase1_v1" / "study_specific_evaluations.csv",
        ROOT / "results" / "performance" / "phase1_v1" / "subject_level_scores.csv",
        ROOT / "metadata" / "phase1" / "phase1_task_summary.csv",
        ROOT / "literature" / "signature_inventory.csv",
    ] + sorted(RESULTS.glob("phase3_*.csv")) + sorted(RESULTS.glob("phase3_*.json"))
    for source in sources:
        shutil.copy2(source, SOURCE_DIR / source.name)


def main() -> None:
    setup_style()
    evaluations = pd.read_csv(ROOT / "results" / "performance" / "phase1_v1" / "study_specific_evaluations.csv")
    evaluations = evaluations[evaluations["execution_status"] == "EXECUTED"].copy()
    numeric = ["auroc", "auroc_ci_lower", "auroc_ci_upper", "prevalence_normalized_ap", "holm_adjusted_p_value"]
    for column in numeric:
        evaluations[column] = pd.to_numeric(evaluations[column], errors="coerce")
    tasks = pd.read_csv(ROOT / "metadata" / "phase1" / "phase1_task_summary.csv")
    tasks = tasks[tasks["status"] == "EXECUTABLE"].copy()
    inventory = pd.read_csv(ROOT / "literature" / "signature_inventory.csv")
    inventory = inventory[inventory["main_grid_eligible"].astype(str).str.upper() == "TRUE"].copy()
    cells = pd.read_csv(RESULTS / "phase3_cell_robustness.csv")
    for column in ["same_therapy", "same_endpoint", "same_disease", "all_scalings_same_direction"]:
        if column in cells:
            cells[column] = cells[column].astype(str).str.lower().eq("true")

    plot_figure_1(evaluations, tasks, inventory)
    plot_figure_2(evaluations)
    plot_figure_3(evaluations)
    plot_figure_4(cells)
    plot_figure_5()
    plot_figure_6(cells)
    plot_supplementary(evaluations, cells)
    copy_source_data()
    print(f"Created 6 main and 4 supplementary figures in {FIG_DIR} and {SUPP_FIG_DIR}")


if __name__ == "__main__":
    main()
