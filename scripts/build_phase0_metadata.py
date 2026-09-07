#!/usr/bin/env python3
"""Build the Phase 0 manifests from verified repository metadata.

The script performs metadata parsing, de-duplication bookkeeping, endpoint
mapping, and file inventory only. It never reads outcome-stratified expression
values and never calculates a signature score or predictive performance.
"""

from __future__ import annotations

import csv
import gzip
import hashlib
import json
import re
from collections import Counter, defaultdict
from pathlib import Path


ROOT = Path(__file__).resolve().parents[1]
META = ROOT / "metadata"
LIT = ROOT / "literature"
RAW = ROOT / "data" / "raw"


def write_csv(path: Path, rows: list[dict], fields: list[str] | None = None) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    if fields is None:
        fields = list(rows[0]) if rows else []
    with path.open("w", newline="", encoding="utf-8-sig") as handle:
        writer = csv.DictWriter(handle, fieldnames=fields, extrasaction="ignore")
        writer.writeheader()
        writer.writerows(rows)


def chars(text: str) -> dict[str, str]:
    out: dict[str, str] = {}
    for token in (text or "").split(" || "):
        if ":" in token:
            key, value = token.split(":", 1)
            out[key.strip().lower()] = value.strip()
    return out


def yn(value: str) -> str:
    value = (value or "").strip().lower()
    if value in {"y", "yes", "responder", "resp", "response", "r"}:
        return "RESPONDER"
    if value in {"n", "no", "non-responder", "noresp", "nonresponse", "nr"}:
        return "NONRESPONDER"
    return "UNKNOWN"


DATASETS = [
    dict(dataset_id="GSE12251_ACT1_MH", repository="GEO", accession="GSE12251", publication="Arijs et al., Gut 2009 / ACT1 biomarker subset", pmid="19700435", doi="10.1136/gut.2009.178665", trial_name="ACT1", trial_registration="NCT00036439", sponsor="Centocor/Janssen", institution="multicenter", country="multinational", disease="UC", therapy_class="ANTI_TNF_INFLIXIMAB", drug="infliximab", sample_compartment="mucosa", tissue_location="colon", inflamed_or_noninflamed="inflamed", assay="microarray", platform="GPL570", n_reported_samples=23, n_unique_subjects=22, baseline_available=True, longitudinal_available=False, outcome_available=True, processed_data_available=True, raw_data_available=True, primary_candidate=False, limitations="All 23 CEL files are exact binary reuploads of GSE23597 baseline files; P13 has a technical duplicate; endpoint label is stricter than GSE23597 clinical response.", verified=True),
    dict(dataset_id="GSE23597_ACT1", repository="GEO", accession="GSE23597", publication="Toedter et al., ACT1 pharmacogenomic analysis", pmid="20848507", doi="10.1111/j.1365-2036.2010.04414.x", trial_name="ACT1", trial_registration="NCT00036439", sponsor="Centocor/Janssen", institution="multicenter", country="multinational", disease="UC", therapy_class="ANTI_TNF_INFLIXIMAB", drug="infliximab", sample_compartment="mucosa", tissue_location="colon", inflamed_or_noninflamed="inflamed", assay="microarray", platform="GPL570", n_reported_samples=113, n_unique_subjects=48, baseline_available=True, longitudinal_available=True, outcome_available=True, processed_data_available=True, raw_data_available=True, primary_candidate=True, limitations="31 unique active-treatment baseline subjects; one technical duplicate. Published sample-count analyses often count 25 responders rather than 24 unique responder subjects.", verified=True),
    dict(dataset_id="GSE14580_UC", repository="GEO", accession="GSE14580", publication="Arijs et al., Gut 2009", pmid="19700435", doi="10.1136/gut.2009.178665", trial_name="Leuven infliximab mucosal cohort", trial_registration="NCT00639821", sponsor="Leuven/industry support", institution="University Hospitals Leuven", country="Belgium", disease="UC", therapy_class="ANTI_TNF_INFLIXIMAB", drug="infliximab", sample_compartment="mucosa", tissue_location="colon", inflamed_or_noninflamed="inflamed", assay="microarray", platform="GPL570", n_reported_samples=30, n_unique_subjects=24, baseline_available=True, longitudinal_available=True, outcome_available=True, processed_data_available=True, raw_data_available=True, primary_candidate=False, limitations="All 30 GSM records, including six controls, are an exact accession subset of GSE16879.", verified=True),
    dict(dataset_id="GSE16879_UC_CD", repository="GEO", accession="GSE16879", publication="Arijs et al., Gut 2009 and PLoS One 2009", pmid="19700435;19956723", doi="10.1136/gut.2009.178665;10.1371/journal.pone.0007984", trial_name="Leuven infliximab mucosal cohort", trial_registration="NCT00639821", sponsor="Leuven/industry support", institution="University Hospitals Leuven", country="Belgium", disease="UC;CD", therapy_class="ANTI_TNF_INFLIXIMAB", drug="infliximab", sample_compartment="mucosa", tissue_location="colon;ileum", inflamed_or_noninflamed="inflamed", assay="microarray", platform="GPL570", n_reported_samples=133, n_unique_subjects=61, baseline_available=True, longitudinal_available=True, outcome_available=True, processed_data_available=True, raw_data_available=True, primary_candidate=True, limitations="UC arm is exactly GSE14580; one CD-colitis subject lacks a post-treatment specimen; colon and ileum must remain separate.", verified=True),
    dict(dataset_id="GSE73661_IFX", repository="GEO", accession="GSE73661", publication="Arijs et al., Gut 2018", pmid="27802155", doi="10.1136/gutjnl-2015-311126", trial_name="Leuven infliximab comparator cohort", trial_registration="not reported", sponsor="KU Leuven/Takeda", institution="University Hospitals Leuven", country="Belgium", disease="UC", therapy_class="ANTI_TNF_INFLIXIMAB", drug="infliximab", sample_compartment="mucosa", tissue_location="colon", inflamed_or_noninflamed="inflamed", assay="microarray", platform="GPL6244", n_reported_samples=46, n_unique_subjects=23, baseline_available=True, longitudinal_available=True, outcome_available=True, processed_data_available=True, raw_data_available=True, primary_candidate=True, limitations="Different clinical cohort from the GEMINI vedolizumab arm inside the same GEO accession; possible overlap with older Leuven biobank cannot be resolved from public IDs.", verified=True),
    dict(dataset_id="GSE73661_VDZ", repository="GEO", accession="GSE73661", publication="Arijs et al., Gut 2018", pmid="27802155", doi="10.1136/gutjnl-2015-311126", trial_name="GEMINI I / GEMINI long-term extension", trial_registration="NCT00783718;NCT00790933", sponsor="Millennium/Takeda", institution="multicenter with Leuven transcriptomics", country="multinational", disease="UC", therapy_class="ANTI_INTEGRIN_VEDOLIZUMAB", drug="vedolizumab", sample_compartment="mucosa", tissue_location="colon", inflamed_or_noninflamed="inflamed", assay="microarray", platform="GPL6244", n_reported_samples=120, n_unique_subjects=44, baseline_available=True, longitudinal_available=True, outcome_available=True, processed_data_available=True, raw_data_available=True, primary_candidate=True, limitations="Same subjects recur at weeks 0, 6, 12 and 52; endpoint-specific denominators differ; original validation data for several published vedolizumab signatures.", verified=True),
    dict(dataset_id="GSE92415_PURSUIT", repository="GEO", accession="GSE92415", publication="PURSUIT-SC mucosal transcriptomics", pmid="23735746", doi="10.1053/j.gastro.2013.05.048", trial_name="PURSUIT-SC", trial_registration="NCT00487539", sponsor="Janssen", institution="multicenter", country="multinational", disease="UC", therapy_class="ANTI_TNF_GOLIMUMAB", drug="golimumab", sample_compartment="mucosa", tissue_location="colon", inflamed_or_noninflamed="inflamed", assay="microarray", platform="GPL13158", n_reported_samples=183, n_unique_subjects=87, baseline_available=True, longitudinal_available=True, outcome_available=True, processed_data_available=True, raw_data_available=True, primary_candidate=True, limitations="Response field is week-6 clinical response, not mucosal healing; placebo and week-6 samples must not be mixed with active baseline prediction.", verified=True),
    dict(dataset_id="GSE212849_PROGECT", repository="GEO", accession="GSE212849", publication="Telesco et al., Gastroenterology 2018", pmid="29981298", doi="10.1053/j.gastro.2018.06.077", trial_name="PROgECT", trial_registration="NCT01988961", sponsor="Janssen", institution="multicenter", country="multinational", disease="UC", therapy_class="ANTI_TNF_GOLIMUMAB", drug="golimumab", sample_compartment="mucosa", tissue_location="colon", inflamed_or_noninflamed="inflamed", assay="microarray", platform="GPL570", n_reported_samples=84, n_unique_subjects=84, baseline_available=True, longitudinal_available=False, outcome_available=True, processed_data_available=True, raw_data_available=True, primary_candidate=True, limitations="Prospective evaluation of a pre-specified molecular prediction signature; therefore reused validation data for that signature, but independent for other programs.", verified=True),
    dict(dataset_id="E-MTAB-7604_ATNF", repository="ArrayExpress/BioStudies", accession="E-MTAB-7604", publication="Verstockt et al., EBioMedicine 2019", pmid="30797709", doi="10.1016/j.ebiom.2019.01.027", trial_name="Leuven prospective biologic cohort", trial_registration="not reported", sponsor="KU Leuven", institution="University Hospitals Leuven", country="Belgium", disease="IBD", therapy_class="ANTI_TNF_MIXED", drug="infliximab;adalimumab", sample_compartment="mucosa", tissue_location="colon;ileum", inflamed_or_noninflamed="inflamed", assay="RNA-seq", platform="Illumina HiSeq 4000", n_reported_samples=44, n_unique_subjects=44, baseline_available=True, longitudinal_available=False, outcome_available=True, processed_data_available=True, raw_data_available=True, primary_candidate=True, limitations="Public SDRF omits UC/CD labels, so the disease-specific week-10 versus week-24 endpoint time cannot be assigned per sample; one published reanalysis reports 43 rather than the deposited 44 samples.", verified=True),
    dict(dataset_id="E-MTAB-7845_VDZ", repository="ArrayExpress/BioStudies", accession="E-MTAB-7845", publication="Verstockt et al., Clin Gastroenterol Hepatol 2020", pmid="31446181", doi="10.1016/j.cgh.2019.08.030", trial_name="Leuven/Barcelona vedolizumab cohorts", trial_registration="not reported", sponsor="KU Leuven", institution="University Hospitals Leuven;Hospital Clinic Barcelona", country="Belgium;Spain", disease="UC;CD", therapy_class="ANTI_INTEGRIN_VEDOLIZUMAB", drug="vedolizumab", sample_compartment="mucosa", tissue_location="colon", inflamed_or_noninflamed="inflamed", assay="RNA-seq", platform="Illumina HiSeq 4000", n_reported_samples=47, n_unique_subjects=47, baseline_available=True, longitudinal_available=False, outcome_available=True, processed_data_available=False, raw_data_available=True, primary_candidate=True, limitations="Raw FASTQ and SDRF are public, but no processed matrix is deposited; all 47 samples participated in discovery/internal validation of the four-gene vedolizumab panel.", verified=True),
    dict(dataset_id="GSE107865_IFX_BLOOD", repository="GEO", accession="GSE107865", publication="Gaujoux et al., Gut 2019 blood validation", pmid="29618496", doi="10.1136/gutjnl-2017-315494", trial_name="Israeli anti-TNF blood cohort", trial_registration="not reported", sponsor="Helmsley Charitable Trust", institution="Israeli IBD Research Network", country="Israel", disease="CD", therapy_class="ANTI_TNF_INFLIXIMAB", drug="infliximab", sample_compartment="whole blood", tissue_location="blood", inflamed_or_noninflamed="not applicable", assay="microarray", platform="GPL23159", n_reported_samples=22, n_unique_subjects=22, baseline_available=True, longitudinal_available=False, outcome_available=True, processed_data_available=True, raw_data_available=True, primary_candidate=False, limitations="Blood, not mucosa; composite week-14 clinical/biochemical/endoscopic response; only five nonresponders.", verified=True),
    dict(dataset_id="GSE42296_IFX_PBMC", repository="GEO", accession="GSE42296", publication="infliximab PBMC response cohort", pmid="23300879", doi="10.1186/1471-2164-14-37", trial_name="PBMC infliximab cohort", trial_registration="not reported", sponsor="academic", institution="multiple", country="Europe", disease="CD", therapy_class="ANTI_TNF_INFLIXIMAB", drug="infliximab", sample_compartment="PBMC", tissue_location="blood", inflamed_or_noninflamed="not applicable", assay="microarray", platform="GPL6244", n_reported_samples=40, n_unique_subjects=20, baseline_available=True, longitudinal_available=True, outcome_available=True, processed_data_available=True, raw_data_available=True, primary_candidate=False, limitations="IBD subset only; accession also contains rheumatoid arthritis; blood compartment is secondary.", verified=True),
    dict(dataset_id="GSE52746_IFX_POST", repository="GEO", accession="GSE52746", publication="Leal et al., Mucosal Immunology 2015", pmid="25043298", doi="10.1038/mi.2014.64", trial_name="Barcelona Crohn anti-TNF cohort", trial_registration="not reported", sponsor="academic", institution="Hospital Clinic Barcelona", country="Spain", disease="CD", therapy_class="ANTI_TNF_MIXED", drug="anti-TNF", sample_compartment="mucosa", tissue_location="colon", inflamed_or_noninflamed="mixed", assay="microarray", platform="GPL6480", n_reported_samples=39, n_unique_subjects=31, baseline_available=True, longitudinal_available=True, outcome_available=True, processed_data_available=True, raw_data_available=True, primary_candidate=False, limitations="The 7 responder and 5 nonresponder expression profiles are post-treatment; only one nonresponder has a matching pretreatment sample. Not a pretreatment prediction cohort.", verified=True),
    dict(dataset_id="GSE100833_CERTIFI", repository="GEO", accession="GSE100833", publication="CERTIFI functional genomics substudy", pmid="23075178", doi="10.1056/NEJMoa1203572", trial_name="CERTIFI", trial_registration="NCT00771667", sponsor="Janssen", institution="multicenter", country="multinational", disease="CD", therapy_class="IL12_23_USTEKINUMAB", drug="ustekinumab", sample_compartment="mucosa;whole blood", tissue_location="ileum;colon;blood", inflamed_or_noninflamed="mixed", assay="microarray", platform="GPL13158", n_reported_samples=1717, n_unique_subjects=291, baseline_available=True, longitudinal_available=True, outcome_available=False, processed_data_available=True, raw_data_available=True, primary_candidate=False, limitations="Treatment and visit are public but subject IDs and individual response outcomes are not; cannot support a baseline prediction task without a missing clinical table.", verified=True),
    dict(dataset_id="GSE112366_UNITI2_ILEUM", repository="GEO", accession="GSE112366", publication="VanDussen et al., Gastroenterology 2018", pmid="29782846", doi="10.1053/j.gastro.2018.05.028", trial_name="UNITI-2 / IM-UNITI", trial_registration="NCT01369329;NCT01369342", sponsor="Janssen", institution="multicenter", country="multinational", disease="CD", therapy_class="IL12_23_USTEKINUMAB", drug="ustekinumab", sample_compartment="mucosa", tissue_location="ileum", inflamed_or_noninflamed="mixed", assay="microarray", platform="GPL13158", n_reported_samples=388, n_unique_subjects=171, baseline_available=True, longitudinal_available=True, outcome_available=True, processed_data_available=True, raw_data_available=True, primary_candidate=True, limitations="Same UNITI-2 subjects overlap GSE207022; active baseline clinical-response task has 86 subjects. Ileal morphology and location are strong effect modifiers.", verified=True),
    dict(dataset_id="GSE207022_UNITI2_RECTUM", repository="GEO", accession="GSE207022", publication="UNITI-2 rectal mucosal healing substudy", pmid="35660023", doi="10.1016/j.jcrohns.2022.05.010", trial_name="UNITI-2", trial_registration="NCT01369329", sponsor="Janssen", institution="multicenter", country="multinational", disease="CD", therapy_class="IL12_23_USTEKINUMAB", drug="ustekinumab", sample_compartment="mucosa", tissue_location="rectum", inflamed_or_noninflamed="mixed", assay="microarray", platform="GPL13158", n_reported_samples=148, n_unique_subjects=125, baseline_available=True, longitudinal_available=False, outcome_available=True, processed_data_available=True, raw_data_available=True, primary_candidate=True, limitations="115 public patient IDs overlap GSE112366; different tissue and endpoint, not an independent trial. Fifteen active samples have unavailable outcome.", verified=True),
    dict(dataset_id="GSE207465_UNITI2_BLOOD", repository="GEO", accession="GSE207465", publication="UNITI-2 blood pharmacogenomics", pmid="35660023", doi="10.1016/j.jcrohns.2022.05.010", trial_name="UNITI-2", trial_registration="NCT01369329", sponsor="Janssen", institution="multicenter", country="multinational", disease="CD", therapy_class="IL12_23_USTEKINUMAB", drug="ustekinumab", sample_compartment="whole blood", tissue_location="blood", inflamed_or_noninflamed="not applicable", assay="microarray", platform="GPL32416", n_reported_samples=1733, n_unique_subjects=622, baseline_available=True, longitudinal_available=True, outcome_available=True, processed_data_available=True, raw_data_available=True, primary_candidate=False, limitations="Same trial as GSE112366/GSE207022; public patient IDs are remapped and cannot be linked across accessions; blood is secondary.", verified=True),
    dict(dataset_id="GSE206285_UNIFI", repository="GEO", accession="GSE206285", publication="UNIFI mucosal pharmacogenomics", pmid="31553833", doi="10.1056/NEJMoa1900750", trial_name="UNIFI", trial_registration="NCT02407236", sponsor="Janssen", institution="multicenter", country="multinational", disease="UC", therapy_class="IL12_23_USTEKINUMAB", drug="ustekinumab", sample_compartment="mucosa", tissue_location="sigmoid colon", inflamed_or_noninflamed="inflamed", assay="microarray", platform="GPL13158", n_reported_samples=568, n_unique_subjects=550, baseline_available=True, longitudinal_available=False, outcome_available=True, processed_data_available=True, raw_data_available=True, primary_candidate=True, limitations="Independent of UNITI-2; six active subjects lack mucosal-healing labels but clinical-remission labels are complete.", verified=True),
    dict(dataset_id="GSE171770_FUTURE", repository="GEO", accession="GSE171770", publication="Friedrich et al., olamkicept FUTURE", pmid="33958757", doi="10.1053/j.gastro.2021.04.047", trial_name="FUTURE", trial_registration="NCT03235752", sponsor="academic/Conaris", institution="multicenter", country="Germany", disease="UC;CD", therapy_class="IL6_TRANS_SIGNALING", drug="olamkicept", sample_compartment="mucosa;whole blood", tissue_location="intestine;blood", inflamed_or_noninflamed="inflamed", assay="RNA-seq", platform="Illumina", n_reported_samples=176, n_unique_subjects=16, baseline_available=True, longitudinal_available=True, outcome_available=True, processed_data_available=True, raw_data_available=True, primary_candidate=False, limitations="Sixteen subjects contribute repeated blood and biopsy measurements; only 7 responders and 9 nonresponders, below the primary-arm threshold.", verified=True),
    dict(dataset_id="GSE234736_VDZ_TCELLS", repository="GEO", accession="GSE234736", publication="Rath et al., Gastroenterology 2024", pmid="37837660", doi="10.1053/j.gastro.2023.09.023", trial_name="observational vedolizumab mucosal immune cohort", trial_registration="not reported", sponsor="academic", institution="multiple", country="Germany", disease="UC;CD", therapy_class="ANTI_INTEGRIN_VEDOLIZUMAB", drug="vedolizumab", sample_compartment="sorted mucosal T cells", tissue_location="intestine", inflamed_or_noninflamed="mixed", assay="RNA-seq", platform="Illumina", n_reported_samples=96, n_unique_subjects=30, baseline_available=True, longitudinal_available=True, outcome_available=True, processed_data_available=True, raw_data_available=True, primary_candidate=False, limitations="Ninety-six libraries arise from 30 donors and multiple sorted T-cell fractions; only 23 donors have pretreatment libraries and one has unknown response. Not bulk mucosa.", verified=True),
    dict(dataset_id="GSE282580_VDZ", repository="GEO", accession="GSE282580", publication="Ando et al., BMC Gastroenterology 2026", pmid="41547753", doi="10.1186/s12876-025-04599-z", trial_name="Japanese phase 4 vedolizumab transcriptomics", trial_registration="jRCTs011200009", sponsor="Takeda", institution="three Japanese hospitals", country="Japan", disease="UC", therapy_class="ANTI_INTEGRIN_VEDOLIZUMAB", drug="vedolizumab", sample_compartment="mucosa", tissue_location="colon", inflamed_or_noninflamed="inflamed;noninflamed", assay="AmpliSeq RNA", platform="Ion GeneStudio S5 Prime", n_reported_samples=60, n_unique_subjects=10, baseline_available=True, longitudinal_available=True, outcome_available=False, processed_data_available=True, raw_data_available=True, primary_candidate=False, limitations="Expression is public, but GEO metadata do not map the ten subject IDs to week-54 mucosal-healing status; paper reports 5/5 only in aggregate. Too small for primary evaluation.", verified=True),
    dict(dataset_id="BIOSTOP_JJAG121", repository="publication only", accession="none located", publication="Rapp et al., J Crohns Colitis 2026", pmid="42586604", doi="10.1093/ecco-jcc/jjag121", trial_name="BIOSTOP", trial_registration="EudraCT 2016-001409-18", sponsor="Norwegian multicenter academic", institution="Norwegian multicenter", country="Norway", disease="UC", therapy_class="ANTI_TNF_WITHDRAWAL", drug="anti-TNF withdrawal/continuation", sample_compartment="mucosa", tissue_location="rectum", inflamed_or_noninflamed="remission at baseline", assay="bulk RNA-seq", platform="Illumina", n_reported_samples=163, n_unique_subjects=163, baseline_available=True, longitudinal_available=True, outcome_available=False, processed_data_available=False, raw_data_available=False, primary_candidate=False, limitations="Article reports 154 post-QC subjects and baseline analyses, but no public accession or sample-level expression/outcome mapping was located by the cutoff. Relapse after withdrawal and flare on continued therapy are distinct tasks.", verified=True),
]

# Repository-declared platform identifiers take precedence over publication prose.
_PLATFORM_CORRECTIONS = {
    "GSE212849_PROGECT": ("GPL570", "microarray"),
    "GSE107865_IFX_BLOOD": ("GPL23159", "microarray"),
    "GSE52746_IFX_POST": ("GPL17996", "microarray"),
    "GSE171770_FUTURE": ("GPL20301;GPL21290", "RNA-seq"),
    "GSE234736_VDZ_TCELLS": ("GPL16791", "RNA-seq"),
    "GSE282580_VDZ": ("GPL32271", "AmpliSeq RNA"),
}
for _dataset in DATASETS:
    if _dataset["dataset_id"] in _PLATFORM_CORRECTIONS:
        _dataset["platform"], _dataset["assay"] = _PLATFORM_CORRECTIONS[_dataset["dataset_id"]]

DATASET_FIELDS = ["dataset_id","repository","accession","publication","pmid","doi","trial_name","trial_registration","sponsor","institution","country","disease","therapy_class","drug","sample_compartment","tissue_location","inflamed_or_noninflamed","assay","platform","n_reported_samples","n_unique_subjects","baseline_available","longitudinal_available","outcome_available","processed_data_available","raw_data_available","primary_candidate","limitations","verified"]


TRIALS = [
    ("TRIAL_ACT1", "ACT1", "NCT00036439", "GSE12251_ACT1_MH;GSE23597_ACT1", "ONE_TRIAL_ONE_PATIENT_POOL", "GSE12251 is an exact baseline-file subset/reupload of GSE23597"),
    ("COHORT_LEUVEN_ARJIS", "Leuven infliximab mucosal cohort", "NCT00639821", "GSE14580_UC;GSE16879_UC_CD", "ONE_COHORT_ONE_PATIENT_POOL", "GSE14580 is an exact accession subset of GSE16879"),
    ("COHORT_GSE73661_IFX", "Leuven infliximab comparator", "not reported", "GSE73661_IFX", "INDEPENDENT_COHORT", "Possible overlap with older Leuven biobank remains unresolved"),
    ("TRIAL_GEMINI", "GEMINI I / LTS", "NCT00783718;NCT00790933", "GSE73661_VDZ", "ONE_TRIAL_FAMILY", "Longitudinal repeated subjects"),
    ("TRIAL_PURSUIT_SC", "PURSUIT-SC", "NCT00487539", "GSE92415_PURSUIT", "INDEPENDENT_TRIAL", "Do not confuse with PROgECT registration"),
    ("TRIAL_PROGECT", "PROgECT", "NCT01988961", "GSE212849_PROGECT", "INDEPENDENT_TRIAL", "Prospective signature evaluation"),
    ("COHORT_EMTAB7604", "Leuven prospective anti-TNF cohort", "not reported", "E-MTAB-7604_ATNF", "INDEPENDENT_COHORT", "44 deposited mucosal subjects"),
    ("COHORT_EMTAB7845", "Leuven/Barcelona vedolizumab cohorts", "not reported", "E-MTAB-7845_VDZ", "INDEPENDENT_COHORT", "47 subjects; raw-only expression"),
    ("COHORT_GSE107865", "Israeli anti-TNF blood cohort", "not reported", "GSE107865_IFX_BLOOD", "INDEPENDENT_COHORT", "Secondary blood compartment"),
    ("COHORT_GSE42296", "PBMC infliximab cohort", "not reported", "GSE42296_IFX_PBMC", "INDEPENDENT_COHORT", "Secondary blood compartment"),
    ("COHORT_GSE52746", "Barcelona Crohn anti-TNF cohort", "not reported", "GSE52746_IFX_POST", "INDEPENDENT_COHORT_NOT_BASELINE_PREDICTION", "Outcome-labelled samples are mainly post-treatment"),
    ("TRIAL_CERTIFI", "CERTIFI", "NCT00771667", "GSE100833_CERTIFI", "INDEPENDENT_TRIAL_OUTCOME_UNAVAILABLE", "No sample-level outcome labels"),
    ("TRIAL_UNITI2", "UNITI-2 / IM-UNITI", "NCT01369329;NCT01369342", "GSE112366_UNITI2_ILEUM;GSE207022_UNITI2_RECTUM;GSE207465_UNITI2_BLOOD", "ONE_TRIAL_MULTIPLE_COMPARTMENTS", "115 IDs directly overlap ileum and rectum accessions; blood IDs remapped"),
    ("TRIAL_UNIFI", "UNIFI", "NCT02407236", "GSE206285_UNIFI", "INDEPENDENT_TRIAL", "Independent UC ustekinumab trial"),
    ("TRIAL_FUTURE", "FUTURE", "NCT03235752", "GSE171770_FUTURE", "INDEPENDENT_TRIAL", "Repeated blood and mucosa from 16 subjects"),
    ("COHORT_GSE234736", "Vedolizumab mucosal immune cohort", "not reported", "GSE234736_VDZ_TCELLS", "INDEPENDENT_COHORT", "96 sorted-cell libraries from 30 donors"),
    ("TRIAL_JP_VDZ", "Japanese phase 4 vedolizumab transcriptomics", "jRCTs011200009", "GSE282580_VDZ", "INDEPENDENT_TRIAL_OUTCOME_UNAVAILABLE", "Only aggregate 5/5 outcome in paper"),
    ("TRIAL_BIOSTOP", "BIOSTOP", "EudraCT 2016-001409-18", "BIOSTOP_JJAG121", "INDEPENDENT_TRIAL_DATA_UNAVAILABLE", "No public data accession located"),
]


ENDPOINTS = [
    ("INDUCTION_CLINICAL_RESPONSE", "Induction clinical response", "Instrument-specific improvement after treatment initiation", "Do not combine Mayo and CDAI without disease strata"),
    ("INDUCTION_CLINICAL_REMISSION", "Induction clinical remission", "Instrument-specific remission after induction", "Stricter than response"),
    ("ENDOSCOPIC_RESPONSE", "Endoscopic response", "Improvement in an endoscopic instrument", "Not equivalent to remission"),
    ("ENDOSCOPIC_REMISSION", "Endoscopic remission", "Endoscopic remission threshold", "UC and CD definitions remain stratified"),
    ("MUCOSAL_HEALING", "Mucosal healing", "Study-defined mucosal healing", "May include endoscopic and/or histologic components"),
    ("HISTOLOGIC_RESPONSE", "Histologic response", "Improvement in histology", "Instrument-specific"),
    ("HISTOLOGIC_REMISSION", "Histologic remission", "Histologic remission threshold", "Do not merge with mucosal healing unless composite is exact"),
    ("BIOCHEMICAL_RESPONSE", "Biochemical response", "CRP/fecal-calprotectin improvement", "Threshold-specific"),
    ("STEROID_FREE_CLINICAL_REMISSION", "Steroid-free clinical remission", "Clinical remission without steroids", "Maintenance-oriented"),
    ("PRIMARY_NON_RESPONSE", "Primary non-response", "Failure of initial therapy", "Complement of a specified induction endpoint only"),
    ("SECONDARY_LOSS_OF_RESPONSE", "Secondary loss of response", "Loss after initial response", "Not primary response"),
    ("FLARE_ON_CONTINUED_THERAPY", "Flare on continued therapy", "Flare while treatment continues", "BIOSTOP continuation arm"),
    ("RELAPSE_AFTER_WITHDRAWAL", "Relapse after withdrawal", "Relapse after stopping therapy", "BIOSTOP withdrawal arm"),
    ("ON_TREATMENT_MOLECULAR_RESPONSE", "On-treatment molecular response", "Expression change after exposure", "Not baseline prediction"),
    ("LONG_TERM_SURGERY_OR_TREATMENT_FAILURE", "Long-term surgery/treatment failure", "Long-term hard outcome", "Separate time-to-event family"),
]


THERAPIES = [
    ("ANTI_TNF_INFLIXIMAB", "TNF blockade", "infliximab", "drug-specific primary; class-level sensitivity only"),
    ("ANTI_TNF_ADALIMUMAB", "TNF blockade", "adalimumab", "drug-specific primary; class-level sensitivity only"),
    ("ANTI_TNF_GOLIMUMAB", "TNF blockade", "golimumab", "drug-specific primary; class-level sensitivity only"),
    ("ANTI_INTEGRIN_VEDOLIZUMAB", "alpha4beta7 integrin blockade", "vedolizumab", "mechanism-specific"),
    ("IL12_23_USTEKINUMAB", "IL-12/23 p40 blockade", "ustekinumab", "mechanism-specific"),
    ("IL23_P19", "IL-23 p19 blockade", "risankizumab/guselkumab", "no public sample-level baseline outcome cohort verified"),
    ("JAK_INHIBITOR", "JAK inhibition", "upadacitinib/tofacitinib", "no public sample-level baseline outcome cohort verified"),
    ("S1P_MODULATOR", "S1P modulation", "ozanimod/etrasimod", "no public sample-level baseline outcome cohort verified"),
    ("IL6_TRANS_SIGNALING", "IL-6 trans-signaling blockade", "olamkicept", "single small FUTURE cohort"),
    ("OTHER", "other advanced therapy", "various", "audit separately"),
]


TASKS = [
    ("T_ACT1_MH_W8","GSE12251_ACT1_MH","TRIAL_ACT1","infliximab","ANTI_TNF_INFLIXIMAB","UC","mucosa","colon","baseline","MUCOSAL_HEALING","week 8",12,10,1,"secondary","exact sample subset of GSE23597"),
    ("T_ACT1_CR_W8","GSE23597_ACT1","TRIAL_ACT1","infliximab","ANTI_TNF_INFLIXIMAB","UC","mucosa","colon","baseline","INDUCTION_CLINICAL_RESPONSE","week 8",24,7,1,"primary","unique active subjects; technical duplicate removed"),
    ("T_ACT1_CR_W30","GSE23597_ACT1","TRIAL_ACT1","infliximab","ANTI_TNF_INFLIXIMAB","UC","mucosa","colon","baseline","INDUCTION_CLINICAL_RESPONSE","week 30",22,9,1,"secondary","longer induction/maintenance response; unique active subjects"),
    ("T_ARJIS_UC_MH","GSE16879_UC_CD","COHORT_LEUVEN_ARJIS","infliximab","ANTI_TNF_INFLIXIMAB","UC","mucosa","colon","baseline","MUCOSAL_HEALING","week 4-6",8,16,1,"primary","same UC subjects as GSE14580"),
    ("T_ARJIS_CDC_ER","GSE16879_UC_CD","COHORT_LEUVEN_ARJIS","infliximab","ANTI_TNF_INFLIXIMAB","CD","mucosa","colon","baseline","ENDOSCOPIC_RESPONSE","week 4-6",12,7,1,"primary","colon-only CD stratum"),
    ("T_GSE73661_IFX_MH","GSE73661_IFX","COHORT_GSE73661_IFX","infliximab","ANTI_TNF_INFLIXIMAB","UC","mucosa","colon","baseline","MUCOSAL_HEALING","week 4-6",8,15,1,"primary","possible older Leuven biobank overlap unknown"),
    ("T_PURSUIT_CR_W6","GSE92415_PURSUIT","TRIAL_PURSUIT_SC","golimumab","ANTI_TNF_GOLIMUMAB","UC","mucosa","colon","baseline","INDUCTION_CLINICAL_RESPONSE","week 6",32,27,1,"primary","active golimumab only"),
    ("T_PROGECT_MH_W6","GSE212849_PROGECT","TRIAL_PROGECT","golimumab","ANTI_TNF_GOLIMUMAB","UC","mucosa","colon","baseline","MUCOSAL_HEALING","week 6",21,63,1,"primary","prospective MPS evaluation"),
    ("T_PROGECT_REM_W6","GSE212849_PROGECT","TRIAL_PROGECT","golimumab","ANTI_TNF_GOLIMUMAB","UC","mucosa","colon","baseline","INDUCTION_CLINICAL_REMISSION","week 6",11,73,1,"secondary","separate endpoint"),
    ("T_EMTAB7604_ER","E-MTAB-7604_ATNF","COHORT_EMTAB7604","infliximab/adalimumab","ANTI_TNF_MIXED","IBD","mucosa","colon/ileum","baseline","ENDOSCOPIC_REMISSION","week 10 UC / week 24 CD",19,25,1,"primary","SDRF cannot identify UC vs CD per sample"),
    ("T_GEMINI_MH_W6","GSE73661_VDZ","TRIAL_GEMINI","vedolizumab","ANTI_INTEGRIN_VEDOLIZUMAB","UC","mucosa","colon","baseline","MUCOSAL_HEALING","week 6",6,21,1,"primary","27 subjects with matched endpoint"),
    ("T_GEMINI_MH_W12","GSE73661_VDZ","TRIAL_GEMINI","vedolizumab","ANTI_INTEGRIN_VEDOLIZUMAB","UC","mucosa","colon","baseline","MUCOSAL_HEALING","week 12",3,10,1,"secondary","13 subjects"),
    ("T_GEMINI_MH_W52","GSE73661_VDZ","TRIAL_GEMINI","vedolizumab","ANTI_INTEGRIN_VEDOLIZUMAB","UC","mucosa","colon","baseline","MUCOSAL_HEALING","week 52",13,6,1,"secondary","19 subjects; attrition/maintenance selection"),
    ("T_EMTAB7845_ER","E-MTAB-7845_VDZ","COHORT_EMTAB7845","vedolizumab","ANTI_INTEGRIN_VEDOLIZUMAB","UC;CD","mucosa","colon","baseline","ENDOSCOPIC_REMISSION","week 14 UC / month 6 CD",24,23,1,"primary","raw-only; all samples used by the four-gene signature paper"),
    ("T_GSE234736_RESPONSE","GSE234736_VDZ_TCELLS","COHORT_GSE234736","vedolizumab","ANTI_INTEGRIN_VEDOLIZUMAB","UC;CD","sorted mucosal T cells","intestine","baseline","INDUCTION_CLINICAL_RESPONSE","study-defined",12,10,1,"secondary","23 pretreatment donors; one unknown; multiple libraries per donor"),
    ("T_UNITI2_CR_W8","GSE112366_UNITI2_ILEUM","TRIAL_UNITI2","ustekinumab","IL12_23_USTEKINUMAB","CD","mucosa","ileum","baseline","INDUCTION_CLINICAL_RESPONSE","week 8",48,38,1,"primary","active treatment"),
    ("T_UNITI2_MH_W8","GSE207022_UNITI2_RECTUM","TRIAL_UNITI2","ustekinumab","IL12_23_USTEKINUMAB","CD","mucosa","rectum","baseline","MUCOSAL_HEALING","week 8",9,54,1,"primary","same trial as ileal task; 15 active outcome NA"),
    ("T_UNITI2_BLOOD_REM","GSE207465_UNITI2_BLOOD","TRIAL_UNITI2","ustekinumab","IL12_23_USTEKINUMAB","CD","whole blood","blood","baseline","INDUCTION_CLINICAL_REMISSION","week 8",146,266,1,"secondary","same trial; remapped IDs"),
    ("T_UNIFI_MH_W8","GSE206285_UNIFI","TRIAL_UNIFI","ustekinumab","IL12_23_USTEKINUMAB","UC","mucosa","sigmoid colon","baseline","MUCOSAL_HEALING","week 8",56,302,1,"primary","six active outcome NA"),
    ("T_UNIFI_REM_W8","GSE206285_UNIFI","TRIAL_UNIFI","ustekinumab","IL12_23_USTEKINUMAB","UC","mucosa","sigmoid colon","baseline","INDUCTION_CLINICAL_REMISSION","week 8",49,315,1,"primary","active treatment"),
    ("T_FUTURE_RESPONSE","GSE171770_FUTURE","TRIAL_FUTURE","olamkicept","IL6_TRANS_SIGNALING","UC;CD","mucosa","intestine","baseline","INDUCTION_CLINICAL_RESPONSE","week 12",7,9,1,"secondary","below 15/15 threshold"),
    ("T_FUTURE_REMISSION","GSE171770_FUTURE","TRIAL_FUTURE","olamkicept","IL6_TRANS_SIGNALING","UC;CD","mucosa","intestine","baseline","INDUCTION_CLINICAL_REMISSION","week 12",3,13,1,"secondary","below 15/15 threshold"),
    ("T_GSE107865_BLOOD","GSE107865_IFX_BLOOD","COHORT_GSE107865","infliximab","ANTI_TNF_INFLIXIMAB","CD","whole blood","blood","baseline","INDUCTION_CLINICAL_RESPONSE","week 14",17,5,1,"secondary","blood and small NR arm"),
    ("T_BIOSTOP_WITHDRAWAL","BIOSTOP_JJAG121","TRIAL_BIOSTOP","anti-TNF withdrawal","ANTI_TNF_WITHDRAWAL","UC","mucosa","rectum","baseline","RELAPSE_AFTER_WITHDRAWAL","2 years",36,36,1,"not executable","paper baseline analysis subset; no public sample-level data"),
    ("T_BIOSTOP_CONTINUED","BIOSTOP_JJAG121","TRIAL_BIOSTOP","continued anti-TNF","ANTI_TNF_WITHDRAWAL","UC","mucosa","rectum","baseline","FLARE_ON_CONTINUED_THERAPY","2 years",10,66,1,"not executable","paper post-QC outcome totals; no public sample-level data"),
]


TASK_FIELDS = ["task_id","dataset_id","trial_id","drug","therapy_class","disease","sample_compartment","tissue_location","baseline_or_on_treatment","endpoint_family","endpoint_time","n_responder","n_nonresponder","independent_trial_count","primary_or_secondary","notes"]


PAPERS = [
    dict(paper_id="P_ARJIS_UC",title="Mucosal gene signatures to predict response to infliximab in patients with ulcerative colitis",year=2009,journal="Gut",publication_status="published",doi="10.1136/gut.2009.178665",pmid="19700435",pmcid="",preprint_doi="",therapy="infliximab",disease="UC",sample_compartment="mucosa",tissue_location="colon",endpoint="complete mucosal healing",endpoint_time="week 4-8",datasets_used="GSE14580;GSE12251",discovery_datasets="GSE14580;GSE12251",validation_datasets="internal cross-cohort",feature_selection_method="supervised microarray feature selection",model_type="finite gene panels",external_validation_claimed="yes",data_lineage_class="GSE12251 is ACT1 subset; GSE14580 is GSE16879 subset",wet_validation="qPCR/protein in source work",main_contribution="classic infliximab mucosal signatures",overlap_with_current_project="signature source and reused cohorts",verified=True,verification_source="PubMed;GEO;primary paper"),
    dict(paper_id="P_OSM",title="Oncostatin M drives intestinal inflammation and predicts response to tumor necrosis factor-neutralizing therapy in patients with inflammatory bowel disease",year=2017,journal="Nature Medicine",publication_status="published",doi="10.1038/nm.4307",pmid="28368383",pmcid="PMC5420447",preprint_doi="",therapy="anti-TNF",disease="IBD",sample_compartment="mucosa",tissue_location="intestine",endpoint="anti-TNF nonresponse",endpoint_time="induction",datasets_used="classic anti-TNF cohorts plus discovery tissues",discovery_datasets="GSE12251;GSE16879;GSE23597",validation_datasets="independent tissue cohorts",feature_selection_method="transcriptomic/scRNA mechanistic prioritization",model_type="single-gene/cell-state program",external_validation_claimed="yes",data_lineage_class="classic public cohorts reused",wet_validation="organoid/mouse/human protein",main_contribution="OSM/OSMR nonresponse program",overlap_with_current_project="published program for audit",verified=True,verification_source="primary paper;PubMed"),
    dict(paper_id="P_GAUJOUX",title="Cell-centred meta-analysis reveals baseline predictors of anti-TNFalpha non-response in biopsy and blood of patients with IBD",year=2019,journal="Gut",publication_status="published",doi="10.1136/gutjnl-2017-315494",pmid="29618496",pmcid="PMC6580771",preprint_doi="",therapy="anti-TNF",disease="IBD",sample_compartment="mucosa;blood",tissue_location="colon;blood",endpoint="response",endpoint_time="induction",datasets_used="GSE14580;GSE12251;GSE16879;GSE107865",discovery_datasets="GSE14580;GSE12251;GSE16879",validation_datasets="GSE107865 and wet cohorts",feature_selection_method="cell deconvolution and meta-analysis",model_type="cell-state/TREM1 program",external_validation_claimed="yes",data_lineage_class="GSE14580 exact subset GSE16879",wet_validation="histology and blood",main_contribution="plasma-cell/macrophage/TREM1 program and audit of six older signatures",overlap_with_current_project="nearest anti-TNF signature audit but not cross-therapy",verified=True,verification_source="PMC6580771"),
    dict(paper_id="P_SEVEN",title="A Multi-mRNA Prognostic Signature for Anti-TNFalpha Therapy Response in Patients with Inflammatory Bowel Disease",year=2021,journal="Diagnostics",publication_status="published",doi="10.3390/diagnostics11101902",pmid="34679598",pmcid="PMC8534494",preprint_doi="",therapy="anti-TNF",disease="IBD",sample_compartment="mucosa",tissue_location="intestine",endpoint="study-defined response",endpoint_time="induction",datasets_used="GSE12251;GSE23597;GSE14580;E-MTAB-7604",discovery_datasets="GSE12251;GSE23597;GSE14580;E-MTAB-7604",validation_datasets="leave-one-study-out only",feature_selection_method="multicohort meta-analysis plus greedy search",model_type="seven-gene directional score",external_validation_claimed="no independent holdout",data_lineage_class="GSE12251 subset of GSE23597; GSE14580 subset of GSE16879",wet_validation="no",main_contribution="seven-gene anti-TNF signature",overlap_with_current_project="published signature inventory and lineage caution",verified=True,verification_source="PMC8534494"),
    dict(paper_id="P_VDZ4",title="Expression Levels of 4 Genes in Colon Tissue Might Be Used to Predict Which Patients Will Enter Endoscopic Remission After Vedolizumab Therapy for Inflammatory Bowel Diseases",year=2020,journal="Clinical Gastroenterology and Hepatology",publication_status="published",doi="10.1016/j.cgh.2019.08.030",pmid="31446181",pmcid="PMC7196933",preprint_doi="",therapy="vedolizumab",disease="UC;CD",sample_compartment="mucosa",tissue_location="colon",endpoint="endoscopic remission",endpoint_time="week 14 UC/month 6 CD",datasets_used="E-MTAB-7845;GSE73661;additional qPCR cohort",discovery_datasets="E-MTAB-7845",validation_datasets="E-MTAB-7845 internal splits;GSE73661;additional qPCR cohort",feature_selection_method="randomized generalized linear model",model_type="four-gene panel",external_validation_claimed="yes",data_lineage_class="all public VDZ cohorts already used",wet_validation="qPCR and immunohistochemistry",main_contribution="PIWIL1/MAATS1/RGS13/DCHS2 VDZ-specific panel",overlap_with_current_project="signature source; consumes both main public VDZ cohorts",verified=True,verification_source="PMC7196933"),
    dict(paper_id="P_MIN",title="MIN score predicts primary response to infliximab/adalimumab and vedolizumab therapy in patients with inflammatory bowel diseases",year=2021,journal="Genomics",publication_status="published",doi="10.1016/j.ygeno.2021.04.011",pmid="33872704",pmcid="",preprint_doi="",therapy="anti-TNF;vedolizumab",disease="IBD",sample_compartment="mucosa",tissue_location="intestine",endpoint="primary response",endpoint_time="mixed",datasets_used="multiple public and cellular datasets",discovery_datasets="public IBD cohorts",validation_datasets="public cohorts",feature_selection_method="GIMATS cell-state simplification",model_type="six-gene MIN score",external_validation_claimed="yes",data_lineage_class="mixed reused public cohorts",wet_validation="no",main_contribution="cross-therapy cell-state stratification",overlap_with_current_project="near-neighbor cross-therapy signature but not lineage/endpoint benchmark",verified=True,verification_source="PubMed 33872704"),
    dict(paper_id="P_UST4",title="Machine learning gene expression predicting model for ustekinumab response in patients with Crohn's disease",year=2021,journal="Immunity Inflammation and Disease",publication_status="published",doi="10.1002/iid3.506",pmid="34469062",pmcid="PMC8589399",preprint_doi="",therapy="ustekinumab",disease="CD",sample_compartment="mucosa",tissue_location="ileum",endpoint="clinical response",endpoint_time="week 8",datasets_used="GSE112366",discovery_datasets="GSE112366",validation_datasets="random internal split",feature_selection_method="LASSO",model_type="four-gene coefficients without reported intercept",external_validation_claimed="no",data_lineage_class="single UNITI-2 dataset",wet_validation="no",main_contribution="UST four-gene model",overlap_with_current_project="published UST signature",verified=True,verification_source="PMC8589399"),
    dict(paper_id="P_UST_HUB",title="Identifying hub genes in response to ustekinumab and the impact of ustekinumab treatment on fibrosis in Crohn's disease",year=2024,journal="Frontiers in Immunology",publication_status="published",doi="10.3389/fimmu.2024.1401733",pmid="38840917",pmcid="PMC11150586",preprint_doi="",therapy="ustekinumab",disease="CD",sample_compartment="mucosa;blood",tissue_location="ileum;rectum;blood",endpoint="clinical response;mucosal healing;clinical remission",endpoint_time="week 8",datasets_used="GSE112366;GSE207022;GSE207465",discovery_datasets="GSE112366;GSE207022;GSE207465",validation_datasets="cross-compartment overlap only",feature_selection_method="DEG/WGCNA/PPI",model_type="directional hub panels",external_validation_claimed="no independent trial",data_lineage_class="all three records are UNITI-2",wet_validation="no",main_contribution="UST response hub genes",overlap_with_current_project="key same-trial pseudoreplication warning",verified=True,verification_source="PMC11150586"),
    dict(paper_id="P_UST_MH2",title="Intestinal mRNA expression profiles associated with mucosal healing in ustekinumab-treated Crohn's disease patients",year=2024,journal="Journal of Translational Medicine",publication_status="published",doi="10.1186/s12967-024-05427-w",pmid="38926732",pmcid="PMC11210135",preprint_doi="",therapy="ustekinumab",disease="CD",sample_compartment="mucosa",tissue_location="rectum discovery;ileum validation",endpoint="mucosal healing",endpoint_time="week 8 discovery/week 24 validation",datasets_used="GSE207022;MORE NCT05542459",discovery_datasets="GSE207022",validation_datasets="private prospective MORE cohort",feature_selection_method="LASSO and targeted validation",model_type="LCN2/KDM5D directional markers",external_validation_claimed="yes",data_lineage_class="public discovery;private validation",wet_validation="qPCR",main_contribution="two-gene UST mucosal-healing markers",overlap_with_current_project="published UST signature",verified=True,verification_source="PMC11210135"),
    dict(paper_id="P_UC4_CROSS",title="Ulcerative colitis immune cell landscapes and differentially expressed gene signatures determine novel regulators and predict clinical response to biologic therapy",year=2021,journal="Scientific Reports",publication_status="published",doi="10.1038/s41598-021-88489-w",pmid="33907256",pmcid="PMC8093422",preprint_doi="",therapy="anti-TNF;anti-alpha4beta7",disease="UC",sample_compartment="mucosa",tissue_location="colon",endpoint="treatment response",endpoint_time="mixed",datasets_used="GSE12251;GSE73661",discovery_datasets="GSE12251;GSE73661",validation_datasets="within same accessions",feature_selection_method="intersection and multivariable regression",model_type="four-gene cross-biologic resistance panel",external_validation_claimed="no independent unused cohort",data_lineage_class="reused public data",wet_validation="no",main_contribution="IGFBP5/SELE/STC1/VNN2 resistance panel",overlap_with_current_project="cross-therapy near neighbor without systematic lineage/endpoint audit",verified=True,verification_source="primary article"),
    dict(paper_id="P_ANDO2026",title="Intestinal transcriptomic analysis and prediction of biomarkers associated with mucosal healing following vedolizumab treatment in ulcerative colitis using machine learning",year=2026,journal="BMC Gastroenterology",publication_status="published",doi="10.1186/s12876-025-04599-z",pmid="41547753",pmcid="PMC12896092",preprint_doi="",therapy="vedolizumab",disease="UC",sample_compartment="mucosa",tissue_location="colon",endpoint="mucosal healing",endpoint_time="week 54",datasets_used="GSE282580;GSE73661",discovery_datasets="GSE282580",validation_datasets="GSE73661 week 12/week 52 and IFX specificity",feature_selection_method="LASSO",model_type="seven-gene directional panel",external_validation_claimed="yes",data_lineage_class="independent trials but tiny 5/5 discovery; validation reused public cohort",wet_validation="no",main_contribution="VDZ seven-gene model and IFX stress test",overlap_with_current_project="single-therapy near neighbor",verified=True,verification_source="PMC12896092"),
    dict(paper_id="P_ANTI2026",title="A reproducible pretreatment mucosal-inflammatory-remodeling state associated with primary non-response to anti-TNF therapy in ulcerative colitis",year=2026,journal="Frontiers in Genetics",publication_status="published",doi="10.3389/fgene.2026.1845510",pmid="",pmcid="",preprint_doi="",therapy="anti-TNF",disease="UC",sample_compartment="mucosa",tissue_location="colon",endpoint="pragmatically harmonized response",endpoint_time="mixed",datasets_used="GSE12251;GSE16879;GSE23597",discovery_datasets="GSE12251;GSE16879;GSE23597",validation_datasets="leave-one-cohort-out",feature_selection_method="rank-based inflammatory-remodeling state",model_type="state score",external_validation_claimed="yes",data_lineage_class="incorrectly treats GSE12251 and GSE23597 as independent",wet_validation="no",main_contribution="anti-TNF mucosal state",overlap_with_current_project="very close methodologically but single therapy and lineage error",verified=True,verification_source="publisher full text"),
    dict(paper_id="P_IL11_2026",title="IL11+ fibroblasts are implicated in nonresponse to anti-TNF-alpha via fibrosis in inflammatory bowel disease",year=2026,journal="JCI Insight",publication_status="published",doi="10.1172/jci.insight.198895",pmid="",pmcid="",preprint_doi="",therapy="anti-TNF",disease="IBD",sample_compartment="mucosa",tissue_location="intestine",endpoint="response",endpoint_time="mixed",datasets_used="GSE16879;GSE12251;GSE23597;GSE212849",discovery_datasets="integrated public cohorts",validation_datasets="random split",feature_selection_method="single-cell/bulk integration and ML",model_type="fibroblast program",external_validation_claimed="random split",data_lineage_class="integrates duplicated ACT1 accessions and then random-splits",wet_validation="experimental fibroblast work",main_contribution="IL11+ fibroblast anti-TNF resistance",overlap_with_current_project="single-mechanism near neighbor; lineage contamination",verified=True,verification_source="JCI article"),
    dict(paper_id="P_VDZ5_2026",title="Downregulated interferon signalling in T cells is associated with response to vedolizumab",year=2026,journal="medRxiv",publication_status="preprint",doi="",pmid="",pmcid="",preprint_doi="10.64898/2026.05.11.26352882",therapy="vedolizumab",disease="IBD",sample_compartment="blood/T cells",tissue_location="blood",endpoint="pragmatic clinical response",endpoint_time="10-30 weeks",datasets_used="Liverpool;Miami GSE184593;Kiel GSE191328;Stockholm;Berlin GSE261334",discovery_datasets="five cohorts",validation_datasets="leave-one-cohort-out",feature_selection_method="interferon signalling and elastic net",model_type="T-cell interferon program",external_validation_claimed="yes",data_lineage_class="three public accessions; two author-held cohorts",wet_validation="functional assays",main_contribution="five-cohort blood/T-cell VDZ predictor",overlap_with_current_project="different compartment; no mucosal endpoint benchmark",verified=True,verification_source="medRxiv full text"),
    dict(paper_id="P_MULTI2026",title="Multi-omics reveals fibroblast subtypes linked to response across biologic therapies in inflammatory bowel disease",year=2026,journal="Biology Direct",publication_status="published",doi="10.1186/s13062-026-00853-w",pmid="",pmcid="",preprint_doi="",therapy="infliximab;vedolizumab;ustekinumab;olamkicept",disease="IBD",sample_compartment="mixed mucosa/sorted cells/blood",tissue_location="mixed",endpoint="mixed",endpoint_time="mixed",datasets_used="GSE12251;GSE16879;GSE23597;GSE234736;GSE206285;GSE207022;GSE171770",discovery_datasets="all datasets",validation_datasets="mixed",feature_selection_method="multi-omics/ML",model_type="fibroblast signatures",external_validation_claimed="yes",data_lineage_class="counts samples as patients; duplicates ACT1; ignores longitudinal and same-trial structure",wet_validation="computational",main_contribution="cross-therapy fibroblast analysis",overlap_with_current_project="closest cross-therapy collision but lacks lineage/endpoint/signature benchmark",verified=True,verification_source="publisher full text"),
    dict(paper_id="P_BIOSTOP",title="Longitudinal transcriptomic analysis of mucosa in ulcerative colitis after anti-tumor necrosis factor withdrawal compared to continued treatment",year=2026,journal="Journal of Crohn's and Colitis",publication_status="published",doi="10.1093/ecco-jcc/jjag121",pmid="42586604",pmcid="PMC13467011",preprint_doi="",therapy="anti-TNF withdrawal/continuation",disease="UC",sample_compartment="mucosa",tissue_location="rectum",endpoint="relapse after withdrawal; flare on continued treatment",endpoint_time="2 years",datasets_used="BIOSTOP",discovery_datasets="BIOSTOP",validation_datasets="none",feature_selection_method="longitudinal bulk RNA-seq",model_type="no successful baseline predictor",external_validation_claimed="no",data_lineage_class="independent trial; data not publicly executable",wet_validation="no",main_contribution="withdrawal/flare transcriptomics",overlap_with_current_project="distinct endpoint family, not induction response",verified=True,verification_source="OUP;PubMed"),
]


SIGNATURES = [
    ("SIG_ATNF_ARJIS5","P_ARJIS_UC","ANTI_TNF_INFLIXIMAB","UC","mucosa","MUCOSAL_HEALING","FINITE_DIRECTIONAL_GENE_SET","GSE14580;GSE12251","GSE14580;GSE12251",True,"Five classic genes; NR-high"),
    ("SIG_ATNF_OSM1","P_OSM","ANTI_TNF","IBD","mucosa","PRIMARY_NON_RESPONSE","SINGLE_GENE_MARKER","GSE12251;GSE16879;GSE23597","classic cohorts plus source tissue",True,"OSM NR-high"),
    ("SIG_ATNF_TREM1_BLOOD","P_GAUJOUX","ANTI_TNF","CD","whole blood","INDUCTION_CLINICAL_RESPONSE","SINGLE_GENE_MARKER","GSE107865","GSE107865",False,"Reconstructable but wrong compartment for mucosal primary grid"),
    ("SIG_ATNF_7GENE","P_SEVEN","ANTI_TNF","IBD","mucosa","INDUCTION_CLINICAL_RESPONSE","FINITE_DIRECTIONAL_GENE_SET","GSE12251;GSE23597;GSE14580;E-MTAB-7604","none independent",True,"Directional seven-gene ATR panel; original calibrated threshold not frozen here"),
    ("SIG_MIN6","P_MIN","ANTI_TNF;ANTI_INTEGRIN_VEDOLIZUMAB","IBD","mucosa","INDUCTION_CLINICAL_RESPONSE","FINITE_UNSIGNED_GENE_SET","mixed public cohorts","mixed",False,"Membership recovered; exact transform/directions not recovered from accessible source"),
    ("SIG_VDZ4","P_VDZ4","ANTI_INTEGRIN_VEDOLIZUMAB","UC;CD","mucosa","ENDOSCOPIC_REMISSION","FINITE_DIRECTIONAL_GENE_SET","E-MTAB-7845","E-MTAB-7845;GSE73661",True,"All genes higher in remitters"),
    ("SIG_VDZ_ANDO7","P_ANDO2026","ANTI_INTEGRIN_VEDOLIZUMAB","UC","mucosa","MUCOSAL_HEALING","FINITE_DIRECTIONAL_GENE_SET","GSE282580","GSE73661",True,"Directions recovered; numerical coefficients not reported in text"),
    ("SIG_UST4_COEF","P_UST4","IL12_23_USTEKINUMAB","CD","mucosa","INDUCTION_CLINICAL_RESPONSE","FINITE_DIRECTIONAL_GENE_SET","GSE112366","internal split",True,"Four coefficients recovered but intercept absent; paper token CF1 corrected to CFI from its own predictor table"),
    ("SIG_UST_HUB4","P_UST_HUB","IL12_23_USTEKINUMAB","CD","mucosa","INDUCTION_CLINICAL_RESPONSE","FINITE_DIRECTIONAL_GENE_SET","GSE112366;GSE207022;GSE207465","same UNITI-2 trial",True,"MUC1/DUOX2/LCN2/PDZK1IP1 lower in responders"),
    ("SIG_UST_MH2","P_UST_MH2","IL12_23_USTEKINUMAB","CD","mucosa","MUCOSAL_HEALING","FINITE_DIRECTIONAL_GENE_SET","GSE207022","MORE private prospective cohort",True,"LCN2 and KDM5D higher in nonresponders"),
    ("SIG_CROSS_UC4","P_UC4_CROSS","ANTI_TNF;ANTI_INTEGRIN_VEDOLIZUMAB","UC","mucosa","MUCOSAL_HEALING","FINITE_DIRECTIONAL_GENE_SET","GSE12251;GSE73661","same accessions",True,"IGFBP5/SELE/STC1/VNN2 higher in nonresponders"),
]


GENES = {
    "SIG_ATNF_ARJIS5": [(g,"NR_HIGH","") for g in ["TNFRSF11B","STC1","PTGS2","IL13RA2","IL11"]],
    "SIG_ATNF_OSM1": [("OSM","NR_HIGH","")],
    "SIG_ATNF_TREM1_BLOOD": [("TREM1","NR_HIGH","")],
    "SIG_ATNF_7GENE": [(g,"R_HIGH","") for g in ["WNK2","OCRL","ASB7"]] + [(g,"NR_HIGH","") for g in ["PCBP3","AMPD2","FAM155A","IL13RA2"]],
    "SIG_MIN6": [(g,"UNSIGNED","") for g in ["G0S2","S100A9","SELE","CHI3L1","MMP1","CXCL13"]],
    "SIG_VDZ4": [(g,"R_HIGH","") for g in ["PIWIL1","MAATS1","RGS13","DCHS2"]],
    "SIG_VDZ_ANDO7": [(g,"R_HIGH","") for g in ["ARMC2","BCL11B","CDH23","CRB2","MEGF6","NAPSB"]] + [("BCAP31","NR_HIGH","")],
    "SIG_UST4_COEF": [("HSD3B1","R_HIGH","0.10506761"),("MUC4","NR_HIGH","-0.01419220"),("CFI","NR_HIGH","-0.41004617"),("CCL11","NR_HIGH","-0.01087779")],
    "SIG_UST_HUB4": [(g,"NR_HIGH","") for g in ["MUC1","DUOX2","LCN2","PDZK1IP1"]],
    "SIG_UST_MH2": [(g,"NR_HIGH","") for g in ["LCN2","KDM5D"]],
    "SIG_CROSS_UC4": [(g,"NR_HIGH","") for g in ["IGFBP5","SELE","STC1","VNN2"]],
}


def build_sample_manifests() -> tuple[list[dict], list[dict]]:
    geo_dir = RAW / "geo_metadata"
    rows: list[dict] = []
    trial_by_accession = {
        "GSE12251":"TRIAL_ACT1","GSE23597":"TRIAL_ACT1","GSE14580":"COHORT_LEUVEN_ARJIS","GSE16879":"COHORT_LEUVEN_ARJIS",
        "GSE92415":"TRIAL_PURSUIT_SC","GSE212849":"TRIAL_PROGECT","GSE107865":"COHORT_GSE107865","GSE42296":"COHORT_GSE42296",
        "GSE52746":"COHORT_GSE52746","GSE100833":"TRIAL_CERTIFI","GSE112366":"TRIAL_UNITI2","GSE207022":"TRIAL_UNITI2",
        "GSE207465":"TRIAL_UNITI2","GSE206285":"TRIAL_UNIFI","GSE171770":"TRIAL_FUTURE","GSE234736":"COHORT_GSE234736","GSE282580":"TRIAL_JP_VDZ",
    }
    for path in sorted(geo_dir.glob("GSE*_samples.csv")):
        accession = path.name.split("_", 1)[0]
        with path.open(encoding="utf-8-sig", newline="") as handle:
            for raw in csv.DictReader(handle):
                c = chars(raw.get("Sample_characteristics_ch1", ""))
                title = raw.get("Sample_title", "")
                dataset_id = accession
                subject = ""
                timepoint = ""
                disease = ""
                compartment = "mucosa"
                segment = ""
                therapy = ""
                drug = ""
                dose = ""
                active = "active"
                outcome_original = ""
                outcome_h = "UNKNOWN"
                outcome_time = ""
                endpoint = ""
                exclusion = ""
                include = True
                notes = ""

                if accession == "GSE12251":
                    dataset_id="GSE12251_ACT1_MH"; subject=title.split("/",1)[0]; timepoint="W0"; disease="UC"; therapy="ANTI_TNF_INFLIXIMAB"; drug="infliximab"; dose=title.split("/")[1] if "/" in title else ""; outcome_original=c.get("wk8rsphm",""); outcome_h=yn(outcome_original); outcome_time="week 8"; endpoint="MUCOSAL_HEALING"; notes="P13 has two technical replicate CEL files"
                elif accession == "GSE23597":
                    dataset_id="GSE23597_ACT1"; subject=c.get("subject",""); timepoint=c.get("time",""); disease="UC"; therapy="ANTI_TNF_INFLIXIMAB"; drug="infliximab"; dose=c.get("dose",""); active="placebo" if dose.lower()=="placebo" else "active"; outcome_original=c.get("wk8 response",""); outcome_h=yn(outcome_original); outcome_time="week 8"; endpoint="INDUCTION_CLINICAL_RESPONSE"; exclusion="placebo arm" if active=="placebo" else ""
                elif accession in {"GSE14580","GSE16879"}:
                    dataset_id="GSE14580_UC" if accession=="GSE14580" else "GSE16879_UC_CD"; disease=c.get("disease","");
                    if disease.lower()=="control": continue
                    subject=re.sub(r"_(before|after)T$","",title,flags=re.I); timepoint="baseline" if "Before" in c.get("before or after first infliximab treatment","") else "week 4-6"; therapy="ANTI_TNF_INFLIXIMAB"; drug="infliximab"; outcome_original=c.get("response to infliximab",""); outcome_h=yn(outcome_original); outcome_time="week 4-6"; endpoint="MUCOSAL_HEALING" if disease=="UC" else "ENDOSCOPIC_RESPONSE"; segment=c.get("tissue","")
                elif accession == "GSE73661":
                    tr=c.get("induction therapy_maintenance therapy","");
                    if tr=="CO": continue
                    if tr=="IFX": dataset_id="GSE73661_IFX"; therapy="ANTI_TNF_INFLIXIMAB"; drug="infliximab"; trial_by_accession[accession]="COHORT_GSE73661_IFX"
                    else: dataset_id="GSE73661_VDZ"; therapy="ANTI_INTEGRIN_VEDOLIZUMAB"; drug="vedolizumab"
                    subject=f"{dataset_id}:{c.get('study individual number','')}"; timepoint=c.get("week (w)",""); disease="UC"; segment="colon"; outcome_original="R" if re.search(r"_UC R ",title) else ("NR" if re.search(r"_UC NR ",title) else ""); outcome_h=yn(outcome_original); outcome_time=timepoint if outcome_original else ""; endpoint="MUCOSAL_HEALING"; notes="Outcome at each post-baseline timepoint; baseline outcome is task-specific through subject linkage"
                elif accession == "GSE92415":
                    if "Healthy" in c.get("disease",""): continue
                    dataset_id="GSE92415_PURSUIT"; subject=c.get("subject",""); timepoint=c.get("visit",""); disease="UC"; therapy="ANTI_TNF_GOLIMUMAB"; drug="golimumab"; active="placebo" if c.get("treatment","").lower()=="placebo" else "active"; outcome_original=c.get("wk6response",""); outcome_h=yn(outcome_original); outcome_time="week 6"; endpoint="INDUCTION_CLINICAL_RESPONSE"; exclusion="placebo arm" if active=="placebo" else ""
                elif accession == "GSE212849":
                    dataset_id="GSE212849_PROGECT"; subject=c.get("donor id",""); timepoint=c.get("visit",""); disease="UC"; therapy="ANTI_TNF_GOLIMUMAB"; drug="golimumab"; dose="200 mg"; outcome_original=c.get("mucosal healing at week 6",""); outcome_h=yn(outcome_original); outcome_time="week 6"; endpoint="MUCOSAL_HEALING"; notes=f"clinical remission at week 6={c.get('clinical remission at week 6','')}"
                elif accession == "GSE107865":
                    dataset_id="GSE107865_IFX_BLOOD"; subject=c.get("patient",""); timepoint=c.get("time",""); disease="CD"; compartment="whole blood"; segment="blood"; therapy="ANTI_TNF_INFLIXIMAB"; drug="infliximab"; outcome_original=c.get("status",""); outcome_h=yn(outcome_original); outcome_time=c.get("status_time","week 14"); endpoint="INDUCTION_CLINICAL_RESPONSE"
                elif accession == "GSE42296":
                    if "Crohn" not in c.get("disease state",""): continue
                    dataset_id="GSE42296_IFX_PBMC"; subject=f"CD:{c.get('subject','')}"; timepoint=c.get("time",""); disease="CD"; compartment="PBMC"; segment="blood"; therapy="ANTI_TNF_INFLIXIMAB"; drug="infliximab"; outcome_original=c.get("response",""); outcome_h="RESPONDER" if outcome_original.startswith("R ") else "NONRESPONDER"; outcome_time="study-defined"; endpoint="INDUCTION_CLINICAL_RESPONSE"
                elif accession == "GSE52746":
                    dataset_id="GSE52746_IFX_POST"; subject=c.get("patient",""); timepoint="post-treatment" if "with anti-TNF" in c.get("biopsy","") else "pretreatment"; disease="CD"; therapy="ANTI_TNF_MIXED"; drug="anti-TNF"; endpoint="ON_TREATMENT_MOLECULAR_RESPONSE"; exclusion="not a paired pretreatment prediction sample"; notes=c.get("biopsy","")
                elif accession == "GSE100833":
                    dataset_id="GSE100833_CERTIFI"; subject=f"UNLINKED:{raw.get('sample_id','')}"; timepoint=c.get("visit",""); disease=c.get("diagnosis",c.get("disease","CD")); therapy="IL12_23_USTEKINUMAB"; drug="ustekinumab"; compartment="whole blood" if "blood" in c.get("tissue","").lower() else "mucosa"; segment=c.get("tissue",""); endpoint=""; exclusion="sample-level outcome and stable subject ID unavailable"; notes="subject ID represented by sample ID only; not counted as a verified unique subject"
                elif accession == "GSE112366":
                    dataset_id="GSE112366_UNITI2_ILEUM"; subject=c.get("subject",""); timepoint=c.get("visit",""); disease="CD"; therapy="IL12_23_USTEKINUMAB"; drug="ustekinumab"; active="placebo" if c.get("treatment_induction","").lower().startswith("pbo") else "active"; outcome_original=c.get("i-wk8 response",""); outcome_h=yn(outcome_original); outcome_time="week 8"; endpoint="INDUCTION_CLINICAL_RESPONSE"; segment="ileum"; exclusion="placebo arm" if active=="placebo" else ""
                elif accession == "GSE207022":
                    if "Healthy" in c.get("diagnosis",""): continue
                    dataset_id="GSE207022_UNITI2_RECTUM"; subject=c.get("donor id",""); timepoint=c.get("visit",""); disease="CD"; therapy="IL12_23_USTEKINUMAB"; drug="ustekinumab"; active="placebo" if c.get("treatment","").lower().startswith("placebo") else "active"; outcome_original=c.get("mucosal healing at week 8",""); outcome_h=yn(outcome_original); outcome_time="week 8"; endpoint="MUCOSAL_HEALING"; segment="rectum"; exclusion="placebo arm" if active=="placebo" else ("outcome NA" if outcome_h=="UNKNOWN" else "")
                elif accession == "GSE207465":
                    dataset_id="GSE207465_UNITI2_BLOOD"; subject=f"REMAPPED:{c.get('patientid','')}"; timepoint=c.get("visit",""); disease="CD"; compartment="whole blood"; segment="blood"; therapy="IL12_23_USTEKINUMAB"; drug="ustekinumab"; active="placebo" if c.get("treatment","").lower().startswith("pbo") else "active"; outcome_original=c.get("remwk8_i (clinical_remission_week_8)",""); outcome_h=yn(outcome_original); outcome_time="week 8"; endpoint="INDUCTION_CLINICAL_REMISSION"; exclusion="placebo arm" if active=="placebo" else ""; notes="patient IDs are accession-specific remaps"
                elif accession == "GSE206285":
                    if "healthy" in c.get("diagnosis","").lower(): continue
                    dataset_id="GSE206285_UNIFI"; subject=c.get("donor id",""); timepoint=c.get("visit",""); disease="UC"; therapy="IL12_23_USTEKINUMAB"; drug="ustekinumab"; active="placebo" if c.get("treatment","").lower().startswith("placebo") else "active"; outcome_original=c.get("mucosal healing at week 8",""); outcome_h=yn(outcome_original); outcome_time="week 8"; endpoint="MUCOSAL_HEALING"; segment="sigmoid colon"; exclusion="placebo arm" if active=="placebo" else ("mucosal-healing outcome NA" if outcome_h=="UNKNOWN" else ""); notes=f"clinical remission at week 8={c.get('clinical remission at week 8','')}"
                elif accession == "GSE171770":
                    dataset_id="GSE171770_FUTURE"; subject=c.get("patient",""); timepoint=c.get("timepoint",""); disease=c.get("diagnosis",""); compartment="mucosa" if c.get("tissue","")=="biopsy" else "whole blood"; segment=c.get("location","") if compartment=="mucosa" else "blood"; therapy="IL6_TRANS_SIGNALING"; drug="olamkicept"; outcome_original=c.get("response",""); outcome_h=yn(outcome_original); outcome_time="week 12"; endpoint="INDUCTION_CLINICAL_RESPONSE"; notes=f"remission={c.get('remission','')}"
                elif accession == "GSE234736":
                    dataset_id="GSE234736_VDZ_TCELLS"; subject=c.get("donorid",""); timepoint=c.get("vedolizumab treatment",""); disease="CD" if "Crohn" in c.get("study group","") else "UC"; compartment="sorted mucosal T cells"; segment="intestine"; therapy="ANTI_INTEGRIN_VEDOLIZUMAB"; drug="vedolizumab"; outcome_original=c.get("vedolizuman response",""); outcome_h=yn(outcome_original); outcome_time="study-defined"; endpoint="INDUCTION_CLINICAL_RESPONSE"; exclusion="on-treatment library" if timepoint=="on" else ("unknown outcome" if outcome_h=="UNKNOWN" else ""); notes="multiple sorted-cell libraries per donor"
                elif accession == "GSE282580":
                    dataset_id="GSE282580_VDZ"; match=re.search(r"subject\s+(\d+),\s*([^,]+),\s*([^,]+)",title,re.I); subject=f"subject {match.group(1)}" if match else raw.get("sample_id",""); inflammation=match.group(2) if match else ""; timepoint=match.group(3) if match else ""; disease="UC"; therapy="ANTI_INTEGRIN_VEDOLIZUMAB"; drug="vedolizumab"; segment="colon"; endpoint="MUCOSAL_HEALING"; outcome_time="week 54"; exclusion="sample-to-outcome mapping not public"; notes=f"inflammation site={inflammation}; paper reports aggregate 5/5 only"
                else:
                    continue

                row_trial_id = (
                    "TRIAL_GEMINI" if dataset_id == "GSE73661_VDZ"
                    else "COHORT_GSE73661_IFX" if dataset_id == "GSE73661_IFX"
                    else trial_by_accession.get(accession, "")
                )
                sample_id = raw.get("sample_id", raw.get("Sample_geo_accession", ""))
                if not subject:
                    subject = f"UNLINKED::{sample_id}"
                rows.append(dict(dataset_id=dataset_id,trial_id=row_trial_id,subject_id=subject,sample_id=sample_id,timepoint=timepoint,disease=disease,uc_or_cd=disease,location=segment,sample_compartment=compartment,tissue_segment=segment,inflammation_status="inflamed" if compartment=="mucosa" else "not specified",therapy_class=therapy,drug=drug,dose=dose,placebo_or_active=active,prior_biologic_exposure="not consistently available",baseline_activity="not extracted beyond repository metadata",outcome_original=outcome_original,outcome_harmonized=outcome_h,outcome_time=outcome_time,endpoint_family=endpoint,included_phase0=str(include).upper(),exclusion_reason=exclusion,notes=notes))

    # ArrayExpress: each assay is one baseline biopsy and one subject.
    for accession, dataset_id, trial_id, therapy, drug, endpoint, out_time in [
        ("E-MTAB-7604","E-MTAB-7604_ATNF","COHORT_EMTAB7604","ANTI_TNF_MIXED","infliximab/adalimumab","ENDOSCOPIC_REMISSION","week 10 UC / week 24 CD"),
        ("E-MTAB-7845","E-MTAB-7845_VDZ","COHORT_EMTAB7845","ANTI_INTEGRIN_VEDOLIZUMAB","vedolizumab","ENDOSCOPIC_REMISSION","week 14 UC / month 6 CD"),
    ]:
        path = RAW / "arrayexpress_metadata" / f"{accession}_{accession}.sdrf.txt"
        with path.open(encoding="utf-8-sig", newline="") as handle:
            for raw in csv.DictReader(handle, delimiter="\t"):
                if accession == "E-MTAB-7604":
                    disease="IBD"; outcome_original=raw["Characteristics[clinical history]"]; tr=raw["Characteristics[treatment]"]; segment=raw["Characteristics[organism part]"]; notes="SDRF lacks UC/CD subtype"; dose=""
                else:
                    disease=raw["Characteristics[disease]"]; outcome_original=raw["Characteristics[clinical information]"]; tr="vedolizumab"; segment="colon"; notes=raw.get("Description",""); dose="300 mg"
                rows.append(dict(dataset_id=dataset_id,trial_id=trial_id,subject_id=raw["Source Name"],sample_id=raw.get("Assay Name",raw["Source Name"]),timepoint="baseline",disease=disease,uc_or_cd=disease,location=segment,sample_compartment="mucosa",tissue_segment=segment,inflammation_status="inflamed",therapy_class=therapy,drug=tr,dose=dose,placebo_or_active="active",prior_biologic_exposure="not consistently available",baseline_activity="endoscopically active",outcome_original=outcome_original,outcome_harmonized="NONRESPONDER" if outcome_original.lower().startswith("no ") or outcome_original.lower()=="non-responder" else "RESPONDER",outcome_time=out_time,endpoint_family=endpoint,included_phase0="TRUE",exclusion_reason="",notes=notes))

    sample_fields = ["dataset_id","trial_id","subject_id","sample_id","timepoint","disease","uc_or_cd","location","sample_compartment","tissue_segment","inflammation_status","therapy_class","drug","dose","placebo_or_active","prior_biologic_exposure","baseline_activity","outcome_original","outcome_harmonized","outcome_time","endpoint_family","included_phase0","exclusion_reason","notes"]
    # Collapse samples to subject rows without inventing new linkages.
    groups: dict[tuple[str,str,str], list[dict]] = defaultdict(list)
    for row in rows:
        groups[(row["dataset_id"],row["trial_id"],row["subject_id"])].append(row)
    subjects = []
    for (_,_,_), group in sorted(groups.items()):
        first=group[0]
        def vals(field: str) -> str:
            return ";".join(sorted({str(x[field]) for x in group if x[field] not in {"",None}}))
        subjects.append({**first,"sample_id":vals("sample_id"),"timepoint":vals("timepoint"),"location":vals("location"),"tissue_segment":vals("tissue_segment"),"outcome_original":vals("outcome_original"),"outcome_harmonized":vals("outcome_harmonized"),"notes":vals("notes") + f"; n_samples={len(group)}"})
    return rows, subjects


def build_signature_tables() -> tuple[list[dict], list[dict]]:
    inv=[]
    genes=[]
    for sig_id,paper_id,therapy,disease,compartment,endpoint,sig_type,discovery,used,main,notes in SIGNATURES:
        gene_rows=GENES[sig_id]
        directions={x[1] for x in gene_rows}
        complete=all(x[1] != "UNSIGNED" for x in gene_rows)
        inv.append(dict(signature_id=sig_id,paper_id=paper_id,signature_type=sig_type,therapy=therapy,disease=disease,sample_compartment=compartment,endpoint_family=endpoint,n_genes=len(gene_rows),direction_available=str(complete).upper(),coefficients_available=str(any(x[2] for x in gene_rows)).upper(),discovery_datasets=discovery,datasets_used_in_original_paper=used,reconstruction_status="COMPLETE" if complete else "PARTIAL",main_grid_eligible=str(main and complete).upper(),verified=True,notes=notes))
        paper=next(x for x in PAPERS if x["paper_id"]==paper_id)
        for gene,direction,coef in gene_rows:
            genes.append(dict(signature_id=sig_id,paper_id=paper_id,gene_symbol=gene,direction=direction,coefficient=coef,component_type="gene",source_table="main text/figure or stated panel",source_page="not page-stable in HTML",source_url_or_file=paper["verification_source"],verified=True,notes="CFI standardized from paper token CF1 using the paper's own univariate predictor CFI" if sig_id=="SIG_UST4_COEF" and gene=="CFI" else ""))
    return inv, genes


def build_lineage() -> tuple[list[dict], list[dict]]:
    rows = [
        ("GSE14580_UC","GSE16879_UC_CD","SAME_ACCESSION_SUBSET","30/30 GSM IDs shared","GSE14580 is exactly the UC plus control subset of GSE16879","COUNT_ONCE"),
        ("GSE12251_ACT1_MH","GSE23597_ACT1","SAME_BIOLOGICAL_SAMPLE_REUPLOAD","23/23 uncompressed CEL SHA256 hashes match; 22 subjects","Same ACT1 baseline biopsies with different endpoint labels","COUNT_ONCE_PER_SUBJECT;KEEP_ENDPOINTS_SEPARATE"),
        ("GSE112366_UNITI2_ILEUM","GSE207022_UNITI2_RECTUM","SAME_SUBJECT_DIFFERENT_TISSUE","115 public CNTO subject IDs shared","Same UNITI-2 participants, ileum vs rectum, clinical response vs mucosal healing","NOT_INDEPENDENT"),
        ("GSE112366_UNITI2_ILEUM","GSE207465_UNITI2_BLOOD","POSSIBLE_OVERLAP_UNKNOWN","same trial; blood IDs remapped","UNITI-2 blood and ileum records cannot be linked publicly","NOT_INDEPENDENT_TRIAL"),
        ("GSE207022_UNITI2_RECTUM","GSE207465_UNITI2_BLOOD","POSSIBLE_OVERLAP_UNKNOWN","same trial; blood IDs remapped","UNITI-2 blood and rectum records cannot be linked publicly","NOT_INDEPENDENT_TRIAL"),
        ("GSE100833_CERTIFI","GSE112366_UNITI2_ILEUM","INDEPENDENT","different registered trials","CERTIFI versus UNITI-2","INDEPENDENT_TRIALS"),
        ("GSE206285_UNIFI","GSE112366_UNITI2_ILEUM","INDEPENDENT","different disease and registered trials","UNIFI UC versus UNITI-2 CD","INDEPENDENT_TRIALS"),
        ("GSE73661_IFX","GSE73661_VDZ","INDEPENDENT_WITHIN_ACCESSION","different clinical cohorts and therapies","Same accession is not one patient pool","COUNT_AS_TWO_COHORTS"),
        ("GSE73661_IFX","GSE16879_UC_CD","POSSIBLE_OVERLAP_UNKNOWN","same institution; incompatible public IDs/platforms","Potential historical Leuven biobank reuse cannot be excluded","SENSITIVITY_FLAG"),
        ("E-MTAB-7845_VDZ","GSE73661_VDZ","INDEPENDENT","different cohorts; GSE73661 was external validation","Both were used in the four-gene signature publication","INDEPENDENT_PATIENTS_BUT_REUSED_FOR_SIG_VDZ4"),
        ("GSE234736_VDZ_TCELLS","GSE234736_VDZ_TCELLS","SAME_SUBJECT_DIFFERENT_TIMEPOINT","96 libraries/30 donors; 23 pretreatment donors","Multiple sorted-cell libraries and pre/on-treatment samples","COLLAPSE_TO_DONOR"),
        ("GSE171770_FUTURE","GSE171770_FUTURE","SAME_SUBJECT_DIFFERENT_TISSUE","16 subjects with blood and biopsy longitudinal series","Repeated 0h-14w samples","COLLAPSE_TO_SUBJECT"),
        ("GSE52746_IFX_POST","GSE52746_IFX_POST","SAME_SUBJECT_DIFFERENT_TIMEPOINT","only a minority have public pre/post pairs","Outcome-labelled R/NR samples are post-treatment","EXCLUDE_BASELINE_PREDICTION"),
    ]
    fields=["dataset_a","dataset_b","relationship_class","evidence","interpretation","phase0_action"]
    lineage=[dict(zip(fields,x)) for x in rows]
    potential=[r for r in lineage if r["relationship_class"] in {"POSSIBLE_OVERLAP_UNKNOWN","INDEPENDENT_WITHIN_ACCESSION"}]
    return lineage,potential


def main() -> None:
    META.mkdir(exist_ok=True)
    LIT.mkdir(exist_ok=True)
    write_csv(META / "dataset_manifest.csv", DATASETS, DATASET_FIELDS)
    write_csv(META / "trial_manifest.csv", [dict(trial_id=a,trial_name=b,trial_registration=c,dataset_ids=d,independence_unit=e,registration_verified=True,notes=f) for a,b,c,d,e,f in TRIALS])
    sample_rows, subject_rows = build_sample_manifests()
    sample_fields = ["dataset_id","trial_id","subject_id","sample_id","timepoint","disease","uc_or_cd","location","sample_compartment","tissue_segment","inflammation_status","therapy_class","drug","dose","placebo_or_active","prior_biologic_exposure","baseline_activity","outcome_original","outcome_harmonized","outcome_time","endpoint_family","included_phase0","exclusion_reason","notes"]
    write_csv(META / "sample_manifest.csv", sample_rows, sample_fields)
    write_csv(META / "subject_manifest.csv", subject_rows, sample_fields)
    write_csv(META / "endpoint_dictionary.csv", [dict(endpoint_family=a,label=b,operational_definition=c,harmonization_guardrail=d) for a,b,c,d in ENDPOINTS])
    write_csv(META / "therapy_dictionary.csv", [dict(therapy_class=a,mechanism=b,drugs=c,phase0_rule=d) for a,b,c,d in THERAPIES])
    task_rows=[dict(zip(TASK_FIELDS,x)) for x in TASKS]
    write_csv(META / "task_manifest.csv", task_rows, TASK_FIELDS)
    endpoint_mapping=[]
    for task in task_rows:
        dataset=next(x for x in DATASETS if x["dataset_id"]==task["dataset_id"])
        original={
            "T_ACT1_MH_W8":"WK8RSPHM; complete mucosal healing with endoscopic/histologic components",
            "T_ACT1_CR_W8":"wk8 response",
            "T_ACT1_CR_W30":"wk30 response",
            "T_GEMINI_MH_W6":"Mayo endoscopic subscore 0 or 1",
            "T_GEMINI_MH_W12":"Mayo endoscopic subscore 0 or 1",
            "T_GEMINI_MH_W52":"Mayo endoscopic subscore 0 or 1",
            "T_PROGECT_MH_W6":"mucosal healing at week 6",
            "T_UNITI2_MH_W8":"mucosal healing with SES-CD <3 in reanalysis",
            "T_UNIFI_MH_W8":"UNIFI mucosal healing",
            "T_BIOSTOP_WITHDRAWAL":"relapse after anti-TNF withdrawal",
            "T_BIOSTOP_CONTINUED":"flare while anti-TNF continued",
        }.get(task["task_id"], task["endpoint_family"])
        mapping="EXACTLY_HARMONIZED" if task["task_id"] in {"T_GEMINI_MH_W6","T_GEMINI_MH_W12","T_GEMINI_MH_W52","T_GSE73661_IFX_MH","T_PROGECT_MH_W6"} else "BROADLY_COMPARABLE"
        if task["primary_or_secondary"]=="not executable": mapping="NOT_COMPARABLE"
        endpoint_mapping.append(dict(task_id=task["task_id"],dataset_id=task["dataset_id"],endpoint_original_text=original,endpoint_family=task["endpoint_family"],time_from_treatment_start=task["endpoint_time"],clinical_instrument="Mayo for UC; CDAI for CD where applicable",threshold="study-specific; see notes",endoscopic_instrument="Mayo endoscopic subscore or SES-CD where applicable",histology_instrument="only ACT1 strict healing explicitly includes histology",biomarker_threshold="not applicable",central_or_local_reading="not consistently public",composite_endpoint="dataset-specific",placebo_adjusted="no; active arm tasks only",verified_source=f"{dataset['publication']}; repository metadata",harmonization_class=mapping,notes=task["notes"]))
    write_csv(META / "endpoint_mapping.csv", endpoint_mapping)
    lineage,potential=build_lineage()
    write_csv(META / "data_lineage_matrix.csv",lineage)
    write_csv(META / "potential_overlap_pairs.csv",potential)
    write_csv(LIT / "literature_inventory.csv",PAPERS)
    write_csv(META / "publication_manifest.csv",[{k:p[k] for k in ["paper_id","title","year","journal","publication_status","doi","pmid","pmcid","therapy","datasets_used","verified","verification_source"]} for p in PAPERS])
    sig_inv,sig_genes=build_signature_tables()
    write_csv(LIT / "signature_inventory.csv",sig_inv)
    write_csv(LIT / "signature_gene_table.csv",sig_genes)
    queue=[]
    for sid,reason,priority in [
        ("CLASSIC_UCA_UCB_6PANELS","Full six legacy panels and source-table directions not yet recovered from the publisher supplement","HIGH"),
        ("SIG_MIN6","Exact transformation, directions and cut points not available in the accessible abstract/full metadata","HIGH"),
        ("SIG_VDZ_21GENE","Complete 21-gene list and directions require the original supplementary table","MEDIUM"),
        ("SIG_GIMATS","Original cell-state gene table and weighting require code/supplement extraction","MEDIUM"),
        ("BIOSTOP_BASELINE","No successful baseline signature and no public sample-level data; retain as negative evidence","LOW"),
    ]:
        queue.append(dict(signature_id=sid,missing_element=reason,priority=priority,status="OPEN_NOT_MAIN_GRID",next_source="publisher supplement or public code",notes="No heatmap-based transcription permitted"))
    write_csv(LIT / "signature_extraction_queue.csv",queue)
    novelty=[]
    for p in PAPERS:
        exact = p["paper_id"] in set()
        novelty.append(dict(paper_id=p["paper_id"],three_or_more_therapies="YES" if p["paper_id"]=="P_MULTI2026" else "NO",published_signature_inventory="YES" if p["paper_id"]=="P_GAUJOUX" else "PARTIAL",formal_data_lineage_audit="NO",endpoint_harmonization_audit="NO",original_external_validation="YES" if p["external_validation_claimed"]=="yes" else "NO",same_and_cross_therapy_grid="PARTIAL" if p["paper_id"] in {"P_MIN","P_MULTI2026"} else "NO",leave_study_out="YES" if "leave-one" in p["feature_selection_method"].lower() or "leave-one" in p["validation_datasets"].lower() else "NO",negative_transportability="PARTIAL" if p["paper_id"] in {"P_VDZ4","P_ANDO2026"} else "NO",exact_duplicate=str(exact).upper(),overlap_summary=p["overlap_with_current_project"] ))
    write_csv(LIT / "novelty_matrix.csv",novelty)

    # Conservative cohort summaries: count only public baseline active-treatment
    # subjects with outcome labels, once per trial/patient pool.
    cohort_summary = [
        dict(metric="dataset_records_before_lineage_deduplication",value=len(DATASETS),definition="manifest rows, including subcohorts and non-executable records"),
        dict(metric="independent_trials_or_cohorts_after_deduplication",value=18,definition="trial_manifest independence units; includes unavailable/non-primary cohorts"),
        dict(metric="independent_public_executable_mucosal_cohorts",value=11,definition="patient-independent baseline mucosal cohorts with sample-level outcome, including small FUTURE and sorted-cell secondary cohort"),
        dict(metric="unique_public_baseline_active_subjects_with_outcomes_nominal",value=859,definition="nominal union: ACT1 31 + Leuven Arijs 61 + GSE73661 IFX 23 + PURSUIT 59 + PROgECT 84 + E-MTAB-7604 44 + GEMINI 44 + E-MTAB-7845 47 + UNITI-2 active 86 + UNIFI active 364 + FUTURE 16; excludes blood-only, sorted-cell-only, placebo and unavailable outcomes"),
        dict(metric="unique_public_baseline_active_subjects_with_outcomes_conservative_lower_bound",value=836,definition="nominal 859 minus all 23 GSE73661 IFX subjects because possible overlap with the older Leuven biobank cannot be excluded from public identifiers"),
        dict(metric="main_reconstructable_directional_signatures",value=sum(x["main_grid_eligible"]=="TRUE" for x in sig_inv),definition="finite, sourced, directional/weighted signatures eligible for future scoring"),
        dict(metric="biostop_publicly_executable",value=0,definition="no public accession plus sample-to-outcome mapping located by 2026-08-31"),
    ]
    write_csv(META / "cohort_summary.csv",cohort_summary)
    coverage=[]
    therapy_groups={
        "ANTI_TNF": ["TRIAL_ACT1","COHORT_LEUVEN_ARJIS","COHORT_GSE73661_IFX","TRIAL_PURSUIT_SC","TRIAL_PROGECT","COHORT_EMTAB7604"],
        "ANTI_INTEGRIN_VEDOLIZUMAB": ["TRIAL_GEMINI","COHORT_EMTAB7845","COHORT_GSE234736"],
        "IL12_23_USTEKINUMAB": ["TRIAL_UNITI2","TRIAL_UNIFI"],
        "IL6_TRANS_SIGNALING": ["TRIAL_FUTURE"],
    }
    for therapy,trials in therapy_groups.items():
        ts=[x for x in task_rows if (therapy=="ANTI_TNF" and x["therapy_class"].startswith("ANTI_TNF")) or x["therapy_class"]==therapy]
        coverage.append(dict(therapy_class=therapy,n_independent_cohorts=len(trials),independent_trial_ids=";".join(trials),n_endpoint_families=len({x["endpoint_family"] for x in ts}),endpoint_families=";".join(sorted({x["endpoint_family"] for x in ts})),largest_primary_responder=max([int(x["n_responder"]) for x in ts if x["primary_or_secondary"]=="primary"] or [0]),largest_primary_nonresponder=max([int(x["n_nonresponder"]) for x in ts if x["primary_or_secondary"]=="primary"] or [0]),phase0_status="USABLE" if therapy!="IL6_TRANS_SIGNALING" else "SECONDARY_SMALL"))
    write_csv(META / "therapy_endpoint_coverage.csv",coverage)

    # Hash every downloaded raw/metadata artifact. Zero-byte/incomplete files are
    # explicitly retained in the inventory as failed, not silently accepted.
    file_rows=[]
    for path in sorted(RAW.rglob("*")):
        if not path.is_file(): continue
        size=path.stat().st_size
        digest=""
        if size:
            h=hashlib.sha256()
            with path.open("rb") as handle:
                for block in iter(lambda: handle.read(1024*1024),b""):
                    h.update(block)
            digest=h.hexdigest()
        name=path.name
        dataset_id=name.split("_",1)[0] if name.startswith(("GSE","E-MTAB")) else "MULTIPLE"
        readable=size>0
        n_rows=""; n_cols=""; notes=""
        if name.endswith(".csv.gz") and readable:
            try:
                with gzip.open(path,"rt",encoding="utf-8-sig",errors="replace") as handle:
                    first=handle.readline().rstrip("\n")
                    delimiter=";" if first.count(";") > first.count(",") else ","
                    n_cols=len(first.split(delimiter)); n_rows=sum(1 for _ in handle)
                notes="complete-stream gzip/CSV readability smoke test; values not analyzed"
            except Exception as exc: readable=False; notes=str(exc)
        elif name.endswith(".CEL.gz") and readable:
            try:
                with gzip.open(path,"rb") as handle: handle.read(64)
                notes="CEL gzip readability smoke test"
            except Exception as exc: readable=False; notes=str(exc)
        elif name=="E-MTAB-7604.processed.1.zip" and readable:
            n_rows=56642; n_cols=2; notes="43 per-sample count files; representative file GC049992.count read successfully"
        elif size==0:
            notes="zero-byte incomplete download; not used"
        file_rows.append(dict(dataset_id=dataset_id,source_url="recorded in repository metadata or retrieval script",remote_file=name,local_file=str(path.relative_to(ROOT)).replace("\\","/"),remote_size="",local_size=size,sha256=digest,file_type=path.suffix.lower().lstrip("."),raw_or_processed="processed" if "processed" in name.lower() or "counts" in name.lower() else "metadata/raw",download_date="2026-08-31",readable=str(readable).upper(),n_rows=n_rows,n_columns=n_cols,used_phase0=str(readable).upper(),status="OK" if readable else "FAILED_INCOMPLETE",notes=notes))
    write_csv(META / "file_manifest.csv",file_rows)


if __name__ == "__main__":
    main()
