#!/usr/bin/env python3
"""Download small GEO family SOFT metadata files for Phase 0 only.

No expression values are analysed. Files are streamed, size-checked, hashed, and
recorded in a machine-readable download manifest.
"""

from __future__ import annotations

import argparse
import csv
import datetime as dt
import gzip
import hashlib
import json
import re
import urllib.error
import urllib.request
from pathlib import Path


USER_AGENT = "ibd-treatment-response-phase0-audit/0.1 (metadata audit)"


def geo_bucket(accession: str) -> str:
    return accession[:-3] + "nnn"


def sha256(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as handle:
        for block in iter(lambda: handle.read(1024 * 1024), b""):
            digest.update(block)
    return digest.hexdigest()


def request(url: str, method: str = "GET"):
    req = urllib.request.Request(url, method=method, headers={"User-Agent": USER_AGENT})
    return urllib.request.urlopen(req, timeout=120)


def parse_family_soft(path: Path) -> dict:
    series: dict[str, list[str]] = {}
    samples: list[dict[str, object]] = []
    current: dict[str, object] | None = None
    section = None

    with gzip.open(path, "rt", encoding="utf-8", errors="replace") as handle:
        for raw in handle:
            line = raw.rstrip("\n\r")
            if line.startswith("^SERIES = "):
                section = "series"
                continue
            if line.startswith("^SAMPLE = "):
                if current:
                    samples.append(current)
                section = "sample"
                current = {"sample_id": line.split("=", 1)[1].strip()}
                continue
            if line.startswith("^"):
                if current:
                    samples.append(current)
                    current = None
                section = None
                continue
            if not line.startswith("!") or " = " not in line:
                continue
            key, value = line[1:].split(" = ", 1)
            if section == "series":
                series.setdefault(key, []).append(value)
            elif section == "sample" and current is not None:
                existing = current.get(key)
                if existing is None:
                    current[key] = value
                elif isinstance(existing, list):
                    existing.append(value)
                else:
                    current[key] = [existing, value]
    if current:
        samples.append(current)

    return {"series": series, "samples": samples}


def flatten(value: object) -> str:
    if isinstance(value, list):
        return " || ".join(str(x) for x in value)
    return "" if value is None else str(value)


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("accessions", nargs="+")
    parser.add_argument("--root", default=".")
    parser.add_argument("--max-gb", type=float, default=2.0)
    args = parser.parse_args()

    root = Path(args.root).resolve()
    out_dir = root / "data" / "raw" / "geo_metadata"
    out_dir.mkdir(parents=True, exist_ok=True)
    rows = []

    for accession in args.accessions:
        accession = accession.upper().strip()
        url = (
            "https://ftp.ncbi.nlm.nih.gov/geo/series/"
            f"{geo_bucket(accession)}/{accession}/soft/{accession}_family.soft.gz"
        )
        local = out_dir / f"{accession}_family.soft.gz"
        status = "NOT_EXECUTED"
        error = ""
        remote_size = None
        try:
            with request(url, "HEAD") as response:
                remote_size = int(response.headers.get("Content-Length", "0") or 0)
            if remote_size and remote_size > args.max_gb * 1024**3:
                status = "SKIPPED_GT_SIZE_LIMIT"
            else:
                if not local.exists() or (remote_size and local.stat().st_size != remote_size):
                    with request(url) as response, local.open("wb") as output:
                        while True:
                            block = response.read(1024 * 1024)
                            if not block:
                                break
                            output.write(block)
                parsed = parse_family_soft(local)
                payload = out_dir / f"{accession}_parsed.json"
                payload.write_text(json.dumps(parsed, indent=2, ensure_ascii=False), encoding="utf-8")

                sample_rows = parsed["samples"]
                keys = sorted({key for row in sample_rows for key in row})
                with (out_dir / f"{accession}_samples.csv").open("w", newline="", encoding="utf-8-sig") as handle:
                    writer = csv.DictWriter(handle, fieldnames=keys)
                    writer.writeheader()
                    for row in sample_rows:
                        writer.writerow({key: flatten(row.get(key)) for key in keys})
                status = "DOWNLOADED_PARSED"
        except (urllib.error.URLError, TimeoutError, OSError, ValueError) as exc:
            status = "NETWORK_BLOCKED_OR_PARSE_ERROR"
            error = repr(exc)

        rows.append(
            {
                "dataset_id": accession,
                "source_url": url,
                "remote_file": f"{accession}_family.soft.gz",
                "local_file": local.relative_to(root).as_posix(),
                "remote_size": remote_size or "",
                "local_size": local.stat().st_size if local.exists() else "",
                "sha256": sha256(local) if local.exists() else "",
                "file_type": "GEO_FAMILY_SOFT_GZ",
                "raw_or_processed": "repository_metadata",
                "download_date": dt.date.today().isoformat(),
                "readable": status == "DOWNLOADED_PARSED",
                "n_rows": len(parse_family_soft(local)["samples"]) if local.exists() and status == "DOWNLOADED_PARSED" else "",
                "n_columns": "variable_soft_fields",
                "used_phase0": True,
                "status": status,
                "notes": error,
            }
        )

    manifest = out_dir / "geo_file_manifest.csv"
    with manifest.open("w", newline="", encoding="utf-8-sig") as handle:
        writer = csv.DictWriter(handle, fieldnames=list(rows[0]))
        writer.writeheader()
        writer.writerows(rows)


if __name__ == "__main__":
    main()

