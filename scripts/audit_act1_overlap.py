#!/usr/bin/env python3
"""Audit ACT1 baseline sample reuse between GSE12251 and GSE23597."""

from __future__ import annotations

import csv
import gzip
import hashlib
import re
from pathlib import Path


ROOT = Path(__file__).resolve().parents[1]
CEL_DIR = ROOT / "data" / "raw" / "representative_expression" / "act1_overlap_cel"


def digest(path: Path, decompress: bool = False) -> str:
    h = hashlib.sha256()
    opener = gzip.open if decompress else open
    with opener(path, "rb") as handle:
        for block in iter(lambda: handle.read(1024 * 1024), b""):
            h.update(block)
    return h.hexdigest()


def load(accession: str) -> list[dict[str, str]]:
    path = ROOT / "data" / "raw" / "geo_metadata" / f"{accession}_samples.csv"
    with path.open(encoding="utf-8-sig", newline="") as handle:
        return list(csv.DictReader(handle))


def subject(row: dict[str, str], accession: str) -> str:
    text = row.get("Sample_title", "") if accession == "GSE12251" else row.get("Sample_characteristics_ch1", "")
    match = re.search(r"(?:^|subject:\s*)(P\d+)", text)
    return match.group(1) if match else ""


def main() -> None:
    records = []
    by_subject: dict[tuple[str, str], list[dict[str, str]]] = {}
    for accession in ("GSE12251", "GSE23597"):
        for row in load(accession):
            if accession == "GSE23597" and "time: W0" not in row.get("Sample_characteristics_ch1", ""):
                continue
            sid = row["sample_id"]
            path = CEL_DIR / f"{accession}_{sid}.CEL.gz"
            if not path.exists():
                continue
            record = {
                "dataset_id": accession,
                "sample_id": sid,
                "subject_id": subject(row, accession),
                "sample_title": row.get("Sample_title", ""),
                "compressed_sha256": digest(path),
                "uncompressed_cel_sha256": digest(path, decompress=True),
                "local_file": path.relative_to(ROOT).as_posix(),
                "local_size": str(path.stat().st_size),
            }
            records.append(record)
            by_subject.setdefault((accession, record["subject_id"]), []).append(record)

    out = ROOT / "results" / "lineage" / "act1_baseline_file_hashes.csv"
    out.parent.mkdir(parents=True, exist_ok=True)
    with out.open("w", encoding="utf-8-sig", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=list(records[0]))
        writer.writeheader()
        writer.writerows(records)

    comparisons = []
    subjects = sorted({r["subject_id"] for r in records if r["dataset_id"] == "GSE12251"})
    for sid in subjects:
        left = by_subject.get(("GSE12251", sid), [])
        right = by_subject.get(("GSE23597", sid), [])
        left_hashes = {r["uncompressed_cel_sha256"] for r in left}
        right_hashes = {r["uncompressed_cel_sha256"] for r in right}
        comparisons.append(
            {
                "subject_id": sid,
                "n_gse12251_baseline_files": len(left),
                "n_gse23597_baseline_files": len(right),
                "n_exact_uncompressed_cel_hash_matches": len(left_hashes & right_hashes),
                "all_gse12251_files_found_in_gse23597": left_hashes <= right_hashes,
                "lineage_conclusion": "SAME_BIOLOGICAL_SAMPLE_REUPLOAD" if left_hashes and left_hashes <= right_hashes else "SUBJECT_MATCH_HASH_MISMATCH_OR_MISSING",
            }
        )

    out2 = ROOT / "results" / "lineage" / "gse12251_gse23597_subject_comparison.csv"
    with out2.open("w", encoding="utf-8-sig", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=list(comparisons[0]))
        writer.writeheader()
        writer.writerows(comparisons)


if __name__ == "__main__":
    main()

