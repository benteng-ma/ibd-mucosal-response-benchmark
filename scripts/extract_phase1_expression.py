#!/usr/bin/env python3
"""Extract only frozen signature genes from each cohort, independently.

No outcomes are read here. Microarray values are repository-normalized GEO
values; multiple directly annotated probes are collapsed by median. RNA-seq
counts are converted to log2(CPM + 0.5) within each sample. Cohorts are never
combined or batch corrected.
"""

from __future__ import annotations

import csv
import gzip
import hashlib
import json
import math
import re
import statistics
import zipfile
from collections import defaultdict
from pathlib import Path


ROOT = Path(__file__).resolve().parents[1]
GEO = ROOT / "data" / "raw" / "geo_metadata"
REP = ROOT / "data" / "raw" / "representative_expression"
OUT = ROOT / "data" / "processed" / "phase1_v1"
LOG = ROOT / "results" / "preprocessing"


def read_csv(path: Path) -> list[dict]:
    with path.open(newline="", encoding="utf-8-sig") as handle:
        return list(csv.DictReader(handle))


def write_csv(path: Path, rows: list[dict], fields: list[str]) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    with path.open("w", newline="", encoding="utf-8-sig") as handle:
        writer = csv.DictWriter(handle, fieldnames=fields)
        writer.writeheader()
        writer.writerows(rows)


def sha256(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as handle:
        for block in iter(lambda: handle.read(1024 * 1024), b""):
            digest.update(block)
    return digest.hexdigest()


def split_symbols(value: str) -> set[str]:
    return {
        token.strip().upper()
        for token in re.split(r"\s*///\s*|\s*[;,|]\s*", value or "")
        if re.fullmatch(r"[A-Za-z][A-Za-z0-9_.-]*", token.strip())
    }


def platform_probe_map(soft: Path, platform: str, wanted: set[str]) -> tuple[dict[str, set[str]], dict[str, set[str]]]:
    probe_to_genes: dict[str, set[str]] = {}
    gene_to_probes: dict[str, set[str]] = defaultdict(set)
    in_platform = False
    in_table = False
    header: list[str] = []
    with gzip.open(soft, "rt", encoding="utf-8", errors="replace") as handle:
        for line in handle:
            if line.startswith("^PLATFORM"):
                in_platform = line.rstrip().endswith(platform)
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
                id_idx = lowered.index("id")
                symbol_idx = lowered.index("gene symbol") if "gene symbol" in lowered else None
                assignment_idx = lowered.index("gene_assignment") if "gene_assignment" in lowered else None
                continue
            symbols: set[str] = set()
            if symbol_idx is not None and symbol_idx < len(cells):
                symbols = split_symbols(cells[symbol_idx])
            elif assignment_idx is not None and assignment_idx < len(cells):
                for assignment in cells[assignment_idx].split(" /// "):
                    parts = [x.strip() for x in assignment.split(" // ")]
                    if len(parts) > 1:
                        symbols.update(split_symbols(parts[1]))
            selected = symbols & wanted
            if selected:
                probe = cells[id_idx]
                probe_to_genes[probe] = selected
                for gene in selected:
                    gene_to_probes[gene].add(probe)
    return probe_to_genes, gene_to_probes


def extract_geo(gse: str, dataset_id: str, platform: str, wanted: set[str]) -> dict:
    soft = GEO / f"{gse}_family.soft.gz"
    probe_to_genes, gene_to_probes = platform_probe_map(soft, platform, wanted)
    sample_values: dict[str, dict[str, list[float]]] = {}
    current_sample = ""
    in_table = False
    table_header = False
    n_samples = 0
    with gzip.open(soft, "rt", encoding="utf-8", errors="replace") as handle:
        for line in handle:
            if line.startswith("^SAMPLE"):
                current_sample = line.split("=", 1)[1].strip()
                sample_values[current_sample] = defaultdict(list)
                in_table = False
                n_samples += 1
                if n_samples % 50 == 0:
                    print(f"{gse}: streamed {n_samples} samples", flush=True)
                continue
            if line.startswith("!sample_table_begin"):
                in_table = True
                table_header = True
                continue
            if line.startswith("!sample_table_end"):
                in_table = False
                continue
            if not (in_table and current_sample):
                continue
            if table_header:
                table_header = False
                continue
            cells = line.rstrip("\n").split("\t")
            if len(cells) < 2 or cells[0] not in probe_to_genes:
                continue
            try:
                value = float(cells[1])
            except ValueError:
                continue
            for gene in probe_to_genes[cells[0]]:
                sample_values[current_sample][gene].append(value)

    rows: list[dict] = []
    for sample_id, gene_values in sample_values.items():
        row = {"sample_id": sample_id}
        for gene in sorted(wanted):
            values = gene_values.get(gene, [])
            row[gene] = f"{statistics.median(values):.8g}" if values else ""
        rows.append(row)
    expression_path = OUT / f"{dataset_id}_signature_genes.csv"
    write_csv(expression_path, rows, ["sample_id", *sorted(wanted)])

    mapping = [
        {
            "dataset_id": dataset_id,
            "platform": platform,
            "gene_symbol": gene,
            "n_direct_probes": len(gene_to_probes.get(gene, set())),
            "probe_ids": ";".join(sorted(gene_to_probes.get(gene, set()))),
            "mapping_status": "MAPPED" if gene_to_probes.get(gene) else "MISSING",
        }
        for gene in sorted(wanted)
    ]
    write_csv(LOG / f"{dataset_id}_gene_mapping.csv", mapping, list(mapping[0]))
    return {
        "dataset_id": dataset_id,
        "source_file": str(soft.relative_to(ROOT)).replace("\\", "/"),
        "source_sha256": sha256(soft),
        "platform": platform,
        "n_samples_extracted": len(rows),
        "n_wanted_genes": len(wanted),
        "n_mapped_genes": sum(bool(gene_to_probes.get(g)) for g in wanted),
        "expression_file": str(expression_path.relative_to(ROOT)).replace("\\", "/"),
        "probe_collapse": "median",
        "normalization": "repository-normalized values; no cross-cohort transformation",
    }


def extract_emtab7604(wanted: set[str]) -> dict:
    archive_path = REP / "E-MTAB-7604.processed.1.zip"
    rows: list[dict] = []
    mapping_counts = {gene: 0 for gene in wanted}
    with zipfile.ZipFile(archive_path) as archive:
        members = sorted(n for n in archive.namelist() if n.endswith(".count"))
        for member in members:
            target: dict[str, float] = {}
            total = 0.0
            with archive.open(member) as handle:
                for raw in handle:
                    cells = raw.decode("utf-8", "replace").rstrip("\n").split("\t")
                    if len(cells) < 2:
                        continue
                    try:
                        count = float(cells[1])
                    except ValueError:
                        continue
                    total += count
                    gene = cells[0].upper()
                    if gene in wanted:
                        target[gene] = count
            row = {"sample_id": Path(member).name}
            for gene in sorted(wanted):
                if gene in target and total > 0:
                    row[gene] = f"{math.log2(target[gene] / total * 1_000_000 + 0.5):.8g}"
                    mapping_counts[gene] += 1
                else:
                    row[gene] = ""
            rows.append(row)
    dataset_id = "E-MTAB-7604_ATNF"
    expression_path = OUT / f"{dataset_id}_signature_genes.csv"
    write_csv(expression_path, rows, ["sample_id", *sorted(wanted)])
    mapping = [
        {
            "dataset_id": dataset_id,
            "platform": "Illumina HiSeq 4000",
            "gene_symbol": gene,
            "n_direct_probes": 1 if mapping_counts[gene] else 0,
            "probe_ids": gene if mapping_counts[gene] else "",
            "mapping_status": "MAPPED" if mapping_counts[gene] else "MISSING",
        }
        for gene in sorted(wanted)
    ]
    write_csv(LOG / f"{dataset_id}_gene_mapping.csv", mapping, list(mapping[0]))
    return {
        "dataset_id": dataset_id,
        "source_file": str(archive_path.relative_to(ROOT)).replace("\\", "/"),
        "source_sha256": sha256(archive_path),
        "platform": "Illumina HiSeq 4000",
        "n_samples_extracted": len(rows),
        "n_wanted_genes": len(wanted),
        "n_mapped_genes": sum(v > 0 for v in mapping_counts.values()),
        "expression_file": str(expression_path.relative_to(ROOT)).replace("\\", "/"),
        "probe_collapse": "not applicable; direct gene symbols",
        "normalization": "within-sample log2(CPM + 0.5); no cross-cohort transformation",
    }


def main() -> None:
    OUT.mkdir(parents=True, exist_ok=True)
    LOG.mkdir(parents=True, exist_ok=True)
    inventory = read_csv(ROOT / "literature" / "signature_inventory.csv")
    eligible = {r["signature_id"] for r in inventory if r["main_grid_eligible"].upper() == "TRUE"}
    genes = {
        r["gene_symbol"].upper()
        for r in read_csv(ROOT / "literature" / "signature_gene_table.csv")
        if r["signature_id"] in eligible
    }
    datasets = [
        ("GSE92415", "GSE92415_PURSUIT", "GPL13158"),
        ("GSE212849", "GSE212849_PROGECT", "GPL570"),
        ("GSE112366", "GSE112366_UNITI2_ILEUM", "GPL13158"),
        ("GSE206285", "GSE206285_UNIFI", "GPL13158"),
    ]
    provenance = [extract_geo(*item, genes) for item in datasets]
    provenance.append(extract_emtab7604(genes))
    with (LOG / "phase1_expression_provenance.json").open("w", encoding="utf-8") as handle:
        json.dump({"version": "phase1_v1", "datasets": provenance}, handle, indent=2)
    print(json.dumps({"datasets": len(provenance), "wanted_genes": len(genes)}, indent=2), flush=True)


if __name__ == "__main__":
    main()
