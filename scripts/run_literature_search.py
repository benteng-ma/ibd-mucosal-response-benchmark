#!/usr/bin/env python3
"""Run the locked Phase 0 novelty queries and record database hit counts.

This script retrieves counts only. It does not rank papers or perform any
expression analysis. Search results are frozen to the declared cutoff date.
"""

from __future__ import annotations

import csv
import json
import time
import urllib.parse
import urllib.request
from datetime import date
from pathlib import Path


ROOT = Path(__file__).resolve().parents[1]
OUT = ROOT / "literature" / "search_log.csv"
CUTOFF = "2026-08-31"

QUERIES = [
    "inflammatory bowel disease transcriptomic treatment response external validation",
    "IBD mucosal gene signature biologic response benchmark",
    "anti-TNF transcriptomic response signature external validation UC CD",
    "vedolizumab transcriptomic response signature external validation IBD",
    "ustekinumab transcriptomic response signature external validation IBD",
    "IBD cross-therapy biomarker transcriptomic response",
    "IBD therapy-specific gene signature biologic response",
    "IBD treatment response endpoint heterogeneity transcriptomic",
    "IBD mucosal healing transcriptomic predictor external cohort",
    "GSE12251 GSE23597 overlap",
    "GSE14580 GSE16879 overlap",
    "IBD treatment response benchmark data lineage",
    "IBD biologic response leave-one-study-out transcriptomic",
    "IBD treatment-response signature transportability",
    "golimumab mucosal gene signature PROgECT PURSUIT transcriptomic",
    "olamkicept treatment response transcriptomic IBD mucosa",
    "BIOSTOP anti-TNF withdrawal relapse mucosal transcriptomic",
    "IL11 fibroblast anti-TNF nonresponse IBD transcriptomic",
    "vedolizumab interferon T cell five cohort response transcriptomic",
]


def get_json(url: str, attempts: int = 3) -> dict:
    req = urllib.request.Request(
        url,
        headers={"User-Agent": "ibd-phase0-audit/1.0 (research metadata audit)"},
    )
    error = None
    for attempt in range(attempts):
        try:
            with urllib.request.urlopen(req, timeout=30) as handle:
                return json.load(handle)
        except Exception as exc:  # network failures are logged, never concealed
            error = exc
            time.sleep(1.5 * (attempt + 1))
    raise RuntimeError(f"request failed after {attempts} attempts: {error}")


def pubmed_count(query: str) -> int:
    params = urllib.parse.urlencode(
        {
            "db": "pubmed",
            "term": f"({query}) AND (\"1900\"[Date - Publication] : \"{CUTOFF}\"[Date - Publication])",
            "retmode": "json",
            "retmax": 0,
        }
    )
    data = get_json("https://eutils.ncbi.nlm.nih.gov/entrez/eutils/esearch.fcgi?" + params)
    return int(data["esearchresult"]["count"])


def europe_pmc_count(query: str) -> int:
    params = urllib.parse.urlencode(
        {
            "query": f'({query}) AND FIRST_PDATE:[1900-01-01 TO {CUTOFF}]',
            "format": "json",
            "pageSize": 1,
        }
    )
    data = get_json("https://www.ebi.ac.uk/europepmc/webservices/rest/search?" + params)
    return int(data["hitCount"])


def crossref_count(query: str) -> int:
    params = urllib.parse.urlencode(
        {
            "query.bibliographic": query,
            "filter": f"until-pub-date:{CUTOFF}",
            "rows": 0,
            "mailto": "phase0-audit@example.invalid",
        }
    )
    data = get_json("https://api.crossref.org/works?" + params)
    return int(data["message"]["total-results"])


def main() -> None:
    rows = []
    engines = [
        ("PubMed", pubmed_count),
        ("Europe PMC", europe_pmc_count),
        ("Crossref", crossref_count),
    ]
    for number, query in enumerate(QUERIES, 1):
        for database, counter in engines:
            status = "COMPLETED"
            notes = "count-only API search; screening decisions recorded in literature_inventory.csv"
            count = ""
            try:
                count = counter(query)
            except Exception as exc:
                status = "FAILED_AFTER_RETRY"
                notes = str(exc)
            rows.append(
                {
                    "search_id": f"S{number:02d}-{database.replace(' ', '').upper()}",
                    "database": database,
                    "query": query,
                    "search_date": str(date.today()),
                    "cutoff_date": CUTOFF,
                    "n_hits": count,
                    "n_after_dedup": "NOT_APPLICABLE_PER_DATABASE",
                    "n_full_text_reviewed": "SEE_LITERATURE_INVENTORY",
                    "n_included": "SEE_LITERATURE_INVENTORY",
                    "status": status,
                    "notes": notes,
                }
            )
            time.sleep(0.15)
    OUT.parent.mkdir(parents=True, exist_ok=True)
    with OUT.open("w", newline="", encoding="utf-8-sig") as handle:
        writer = csv.DictWriter(handle, fieldnames=list(rows[0]))
        writer.writeheader()
        writer.writerows(rows)


if __name__ == "__main__":
    main()
