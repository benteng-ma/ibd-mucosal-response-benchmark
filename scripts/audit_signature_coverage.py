#!/usr/bin/env python3
"""Audit whether published signature genes are represented on each platform.

This is a Phase 0 metadata/platform audit. It does not read outcome-stratified
expression values, compute signature scores, fit models, or estimate performance.
"""

from __future__ import annotations

import csv
import gzip
import json
import re
import zipfile
from pathlib import Path


ROOT = Path(__file__).resolve().parents[1]
GEO = ROOT / "data" / "raw" / "geo_metadata"
REP = ROOT / "data" / "raw" / "representative_expression"
OUT = ROOT / "results" / "coverage"


def read_csv(path: Path) -> list[dict]:
    with path.open(newline="", encoding="utf-8-sig") as handle:
        return list(csv.DictReader(handle))


def write_csv(path: Path, rows: list[dict]) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    with path.open("w", newline="", encoding="utf-8-sig") as handle:
        writer = csv.DictWriter(handle, fieldnames=list(rows[0]))
        writer.writeheader()
        writer.writerows(rows)


def split_symbols(value: str) -> set[str]:
    values = re.split(r"\s*///\s*|\s*[;,|]\s*", value or "")
    return {v.strip().upper() for v in values if re.fullmatch(r"[A-Za-z][A-Za-z0-9_.-]*", v.strip())}


def platform_catalog(gse: str, wanted_platform: str) -> set[str]:
    path = GEO / f"{gse}_family.soft.gz"
    symbols: set[str] = set()
    in_platform = False
    in_table = False
    header: list[str] = []
    symbol_idx: int | None = None
    assignment_idx: int | None = None
    with gzip.open(path, "rt", encoding="utf-8", errors="replace") as handle:
        for line in handle:
            if line.startswith("^PLATFORM"):
                in_platform = line.rstrip().endswith(wanted_platform)
                in_table = False
                continue
            if in_platform and line.startswith("!platform_table_begin"):
                in_table = True
                header = []
                continue
            if in_platform and line.startswith("!platform_table_end"):
                break
            if not (in_platform and in_table):
                continue
            cells = line.rstrip("\n").split("\t")
            if not header:
                header = cells
                lowered = [x.casefold() for x in header]
                symbol_idx = lowered.index("gene symbol") if "gene symbol" in lowered else None
                assignment_idx = lowered.index("gene_assignment") if "gene_assignment" in lowered else None
                continue
            if symbol_idx is not None and symbol_idx < len(cells):
                symbols.update(split_symbols(cells[symbol_idx]))
            elif assignment_idx is not None and assignment_idx < len(cells):
                for assignment in cells[assignment_idx].split(" /// "):
                    parts = [x.strip() for x in assignment.split(" // ")]
                    if len(parts) > 1:
                        symbols.update(split_symbols(parts[1]))
    return symbols


def gse234736_catalog() -> set[str]:
    path = REP / "GSE234736_Supplementary_Table_2.csv.gz"
    if not path.exists():
        return set()
    direct_symbols: set[str] = set()
    try:
        with gzip.open(path, "rt", encoding="utf-8-sig", errors="replace") as handle:
            reader = csv.reader(handle)
            next(reader, None)
            for row in reader:
                if row:
                    if row[0].upper().startswith("ENSG"):
                        continue
                    direct_symbols.update(split_symbols(row[0]))
    except (EOFError, OSError):
        # A partial representative download is never used to claim coverage.
        return set()
    # The complete public table is Ensembl keyed; metadata labels are not genes.
    return set()


def emtab7604_catalog() -> set[str]:
    """Return direct symbols only; the provided count files are Ensembl IDs."""
    path = REP / "E-MTAB-7604.processed.1.zip"
    if not path.exists():
        return set()
    with zipfile.ZipFile(path) as archive:
        member = next((n for n in archive.namelist() if n.endswith(".count")), None)
        if not member:
            return set()
        with archive.open(member) as handle:
            identifiers = [
                line.decode("utf-8", "replace").split("\t", 1)[0].strip()
                for line in handle
            ]
    return set() if any(x.upper().startswith("ENSG") for x in identifiers[:100]) else {
        x.upper() for x in identifiers if re.fullmatch(r"[A-Za-z][A-Za-z0-9_.-]*", x)
    }


def main() -> None:
    sources = {
        "GPL570": "GSE16879",
        "GPL6244": "GSE73661",
        "GPL13158": "GSE92415",
        "GPL23159": "GSE107865",
        "GPL17996": "GSE52746",
        "GPL32416": "GSE207465",
    }
    catalogs: dict[str, set[str]] = {}
    catalog_summary: list[dict] = []
    for platform, gse in sources.items():
        genes = platform_catalog(gse, platform)
        catalogs[platform] = genes
        catalog_summary.append({"platform": platform, "source_accession": gse, "n_direct_gene_symbols": len(genes), "status": "PARSED" if genes else "NO_DIRECT_SYMBOLS"})
    catalogs["GPL16791"] = gse234736_catalog()
    catalog_summary.append({"platform": "GPL16791", "source_accession": "GSE234736 Supplementary Table 2", "n_direct_gene_symbols": len(catalogs["GPL16791"]), "status": "ENSEMBL_TO_SYMBOL_CROSSWALK_REQUIRED"})
    catalogs["Illumina HiSeq 4000"] = emtab7604_catalog()
    catalog_summary.append({"platform": "Illumina HiSeq 4000", "source_accession": "E-MTAB-7604 processed counts", "n_direct_gene_symbols": len(catalogs["Illumina HiSeq 4000"]), "status": "PARSED_DIRECT_SYMBOL_COUNTS" if catalogs["Illumina HiSeq 4000"] else "ENSEMBL_TO_SYMBOL_CROSSWALK_REQUIRED"})

    datasets = read_csv(ROOT / "metadata" / "dataset_manifest.csv")
    signatures = read_csv(ROOT / "literature" / "signature_inventory.csv")
    gene_rows = read_csv(ROOT / "literature" / "signature_gene_table.csv")
    sig_genes: dict[str, list[str]] = {}
    for row in gene_rows:
        sig_genes.setdefault(row["signature_id"], []).append(row["gene_symbol"].upper())

    rows: list[dict] = []
    for sig in signatures:
        genes = sorted(set(sig_genes.get(sig["signature_id"], [])))
        for dataset in datasets:
            platform = dataset["platform"]
            available = catalogs.get(platform)
            if available is None:
                mapped: list[str] = []
                status = "NOT_COMPUTABLE_NO_DIRECT_PLATFORM_CATALOG"
                fraction = ""
            elif not available:
                mapped = []
                status = "NOT_COMPUTABLE_CROSSWALK_REQUIRED"
                fraction = ""
            else:
                mapped = [g for g in genes if g in available]
                status = "DIRECT_SYMBOL_AUDIT"
                fraction = f"{len(mapped) / len(genes):.3f}" if genes else ""
            missing = [g for g in genes if g not in set(mapped)]
            pass_primary = bool(fraction) and float(fraction) >= 0.8 and sig["direction_available"].upper() == "TRUE"
            rows.append({
                "signature_id": sig["signature_id"],
                "dataset_id": dataset["dataset_id"],
                "platform": platform,
                "n_signature_genes": len(genes),
                "n_mapped_genes": len(mapped),
                "mapped_genes": ";".join(mapped),
                "missing_genes": ";".join(missing),
                "coverage_fraction": fraction,
                "direction_complete": sig["direction_available"].upper(),
                "mapping_status": status,
                "primary_coverage_pass": str(pass_primary).upper(),
                "notes": "Coverage only; no expression values or outcomes were analyzed.",
            })
    write_csv(ROOT / "metadata" / "signature_coverage_matrix.csv", rows)
    write_csv(OUT / "platform_gene_catalog_summary.csv", catalog_summary)
    with (OUT / "coverage_summary.json").open("w", encoding="utf-8") as handle:
        json.dump({
            "n_signature_dataset_pairs": len(rows),
            "n_directly_audited_pairs": sum(r["mapping_status"] == "DIRECT_SYMBOL_AUDIT" for r in rows),
            "n_primary_coverage_pass_pairs": sum(r["primary_coverage_pass"] == "TRUE" for r in rows),
            "coverage_threshold": 0.8,
            "phase0_only": True,
        }, handle, indent=2)


if __name__ == "__main__":
    main()
