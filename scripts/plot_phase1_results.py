#!/usr/bin/env python3
"""Create compact study-specific Phase 1 audit figures."""

from __future__ import annotations

import csv
from pathlib import Path

import matplotlib.pyplot as plt
import numpy as np


ROOT = Path(__file__).resolve().parents[1]
SOURCE = ROOT / "results" / "performance" / "phase1_v1"
OUT = ROOT / "results" / "figures" / "phase1_v1"


def read_csv(path: Path) -> list[dict]:
    with path.open(newline="", encoding="utf-8-sig") as handle:
        return list(csv.DictReader(handle))


def strict_forest(rows: list[dict]) -> None:
    labels = [f"{r['signature_id']} → {r['task_id']}" for r in rows]
    auc = np.array([float(r["auroc"]) for r in rows])
    low = np.array([float(r["auroc_ci_lower"]) for r in rows])
    high = np.array([float(r["auroc_ci_upper"]) for r in rows])
    y = np.arange(len(rows))[::-1]
    fig, ax = plt.subplots(figsize=(9.2, 4.2))
    ax.axvline(0.5, color="#666666", linestyle="--", linewidth=1)
    ax.errorbar(auc, y, xerr=[auc - low, high - auc], fmt="o", color="#245b78", ecolor="#78a6bd", capsize=4)
    ax.set_yticks(y, labels)
    ax.set_xlim(0.25, 0.9)
    ax.set_xlabel("Subject-level AUROC (95% stratified-bootstrap interval)")
    ax.set_title("Strict same-therapy, same-endpoint external evaluations")
    for yi, r in zip(y, rows):
        ax.text(0.895, yi, f"Holm p={float(r['holm_adjusted_p_value']):.3f}", ha="right", va="center", fontsize=8)
    ax.spines[["top", "right"]].set_visible(False)
    fig.tight_layout()
    fig.savefig(OUT / "strict_external_forest.png", dpi=220)
    plt.close(fig)


def transport_heatmap(rows: list[dict]) -> None:
    signatures = sorted({r["signature_id"] for r in rows})
    tasks = [
        "T_PURSUIT_CR_W6", "T_PROGECT_MH_W6", "T_EMTAB7604_ER",
        "T_UNITI2_CR_W8", "T_UNIFI_MH_W8", "T_UNIFI_REM_W8",
    ]
    values = np.full((len(signatures), len(tasks)), np.nan)
    classes = {}
    for r in rows:
        i = signatures.index(r["signature_id"])
        j = tasks.index(r["task_id"])
        values[i, j] = float(r["auroc"])
        classes[(i, j)] = r["transport_classification"]
    fig, ax = plt.subplots(figsize=(10.8, 7.0))
    image = ax.imshow(values, vmin=0.35, vmax=0.75, cmap="RdBu", aspect="auto")
    ax.set_xticks(range(len(tasks)), tasks, rotation=35, ha="right")
    ax.set_yticks(range(len(signatures)), signatures)
    for i in range(values.shape[0]):
        for j in range(values.shape[1]):
            if np.isfinite(values[i, j]):
                marker = "●" if classes[(i, j)] == "SUPPORTED_DIRECTION" else ""
                ax.text(j, i, f"{values[i, j]:.2f}{marker}", ha="center", va="center", fontsize=8,
                        color="white" if values[i, j] < 0.45 or values[i, j] > 0.68 else "black")
            else:
                ax.text(j, i, "—", ha="center", va="center", color="#777777")
    ax.set_title("Frozen published programs across executable independent tasks\n● 95% bootstrap interval lies above 0.50; blank cells were not eligible")
    cbar = fig.colorbar(image, ax=ax, shrink=0.8)
    cbar.set_label("AUROC; higher score predicts responder")
    fig.tight_layout()
    fig.savefig(OUT / "transportability_matrix.png", dpi=220)
    plt.close(fig)


def main() -> None:
    OUT.mkdir(parents=True, exist_ok=True)
    strict_forest(read_csv(SOURCE / "strict_external_evaluations.csv"))
    transport_heatmap(read_csv(SOURCE / "study_specific_evaluations.csv"))


if __name__ == "__main__":
    main()
