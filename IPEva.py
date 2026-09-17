#!/usr/bin/env python3
"""
IPEva: Imputation Performance Evaluator
Version: 1.0
Auther:Huang Xiaowen
Time:2026.7.14
Truth-based, single-sample SNP imputation performance evaluation.
"""
from __future__ import annotations

import argparse
import gzip
import html
import json
import math
import sys
from collections import Counter, defaultdict
from dataclasses import dataclass, asdict
from datetime import datetime
from pathlib import Path
from typing import Any, Dict, List, Optional, Tuple

VERSION = "1.0"
NA = "NA"
DEFAULT_GPTH = 0.9
DEFAULT_GPSCAN = [round(x / 100, 2) for x in range(90, 100)]
DEFAULT_MAFBIN = [0, 0.001, 0.003, 0.005, 0.007, 0.01, 0.02, 0.05, 0.10, 0.30, 0.50]
DEFAULT_SHETBIN = [0, 0.2, 0.4, 1]


@dataclass
class VariantRecord:
    chrom: str
    pos: int
    ref: str
    alt: str
    gt_code: int


@dataclass
class SiteRecord:
    chrom: str
    pos: int
    ref: str
    alt: str
    true_gt: int
    imputed_gt: int
    gp: float
    correct: bool
    genotype_class: str
    swapped: bool = False
    maf: Optional[float] = None
    maf_bin: Optional[str] = None
    maf_range: Optional[str] = None
    shet: Optional[float] = None
    shet_bin: Optional[str] = None
    shet_range: Optional[str] = None


def open_text(path: str):
    return gzip.open(path, "rt") if path.endswith(".gz") else open(path, "rt")


def is_snp(ref: str, alt: str) -> bool:
    return len(ref) == 1 and len(alt) == 1 and alt not in {".", "*"} and "," not in alt


def parse_gt_code(gt_value: Optional[str]) -> Optional[int]:
    if gt_value is None or gt_value in {".", "./.", ".|."}:
        return None
    gt = gt_value.split(":", 1)[0].replace("|", "/")
    if "/" not in gt:
        return None
    parts = gt.split("/")
    if len(parts) != 2 or any(x == "." for x in parts):
        return None
    try:
        vals = [int(x) for x in parts]
    except ValueError:
        return None
    if any(x not in (0, 1) for x in vals):
        return None
    return sum(vals)


def genotype_class(gt_code: int) -> str:
    return "Heterozygous" if gt_code == 1 else "Homozygous"


def read_vcf_header_samples(path: str) -> List[str]:
    with open_text(path) as fh:
        for line in fh:
            if line.startswith("#CHROM"):
                return line.rstrip("\n").split("\t")[9:]
    raise SystemExit(f"ERROR: no #CHROM header line found in VCF: {path}")


def ensure_single_sample_vcfs(truth_vcf: str, imputed_vcf: str, warnings: List[str]) -> Tuple[str, str]:
    truth_samples = read_vcf_header_samples(truth_vcf)
    imputed_samples = read_vcf_header_samples(imputed_vcf)
    if len(truth_samples) != 1:
        raise SystemExit(
            "ERROR: IPEva only supports single-sample VCF files. "
            f"The truth VCF contains {len(truth_samples)} samples: {', '.join(truth_samples[:10])}. "
            "Please extract one target sample before running IPEva."
        )
    if len(imputed_samples) != 1:
        raise SystemExit(
            "ERROR: IPEva only supports single-sample VCF files. "
            f"The imputed VCF contains {len(imputed_samples)} samples: {', '.join(imputed_samples[:10])}. "
            "Please provide a single-sample imputed VCF."
        )
    if truth_samples[0] != imputed_samples[0]:
        msg = (
            "WARNING: The truth VCF and imputed VCF contain different sample names. "
            f"Truth sample: {truth_samples[0]}; Imputed sample: {imputed_samples[0]}. "
            "IPEva will treat them as the same target sample because both files are single-sample VCFs."
        )
        print(msg, file=sys.stderr)
        warnings.append(msg)
    return truth_samples[0], imputed_samples[0]


def read_truth_vcf(path: str) -> Tuple[Dict[Tuple[str, int], VariantRecord], Dict[str, int]]:
    records: Dict[Tuple[str, int], VariantRecord] = {}
    qc = Counter()
    with open_text(path) as fh:
        for line in fh:
            if line.startswith("#"):
                continue
            fields = line.rstrip("\n").split("\t")
            if len(fields) < 10:
                qc["malformed"] += 1
                continue
            chrom, pos_s, ref, alt = fields[0], fields[1], fields[3], fields[4]
            try:
                pos = int(pos_s)
            except ValueError:
                qc["bad_pos"] += 1
                continue
            if not is_snp(ref, alt):
                qc["non_biallelic_snp"] += 1
                continue
            key = (chrom, pos)
            if key in records:
                qc["duplicate_position"] += 1
                continue
            fmt_keys = fields[8].split(":")
            sample_values = fields[9].split(":")
            fmt_map = dict(zip(fmt_keys, sample_values))
            gt_code = parse_gt_code(fmt_map.get("GT", sample_values[0] if sample_values else None))
            if gt_code is None:
                qc["missing_or_invalid_gt"] += 1
                continue
            records[key] = VariantRecord(chrom=chrom, pos=pos, ref=ref, alt=alt, gt_code=gt_code)
            qc["loaded"] += 1
    return records, dict(qc)


def parse_gp(fmt_keys: List[str], sample_values: List[str]) -> Optional[float]:
    fmt_map = dict(zip(fmt_keys, sample_values))
    gp_raw = fmt_map.get("GP")
    if gp_raw is None or gp_raw in {".", ""}:
        return None
    sep = "," if "," in gp_raw else ";" if ";" in gp_raw else None
    parts = gp_raw.split(sep) if sep else [gp_raw]
    vals: List[float] = []
    try:
        for p in parts:
            if p not in {"", "."}:
                vals.append(float(p))
    except ValueError:
        return None
    return max(vals) if vals else None


def read_imputed_vcf_and_match(path: str, truth_records: Dict[Tuple[str, int], VariantRecord]) -> Tuple[List[SiteRecord], Dict[str, int], List[Dict[str, Any]]]:
    sites: List[SiteRecord] = []
    qc = Counter()
    excluded: List[Dict[str, Any]] = []
    seen_imp = set()
    with open_text(path) as fh:
        for line in fh:
            if line.startswith("#"):
                continue
            fields = line.rstrip("\n").split("\t")
            if len(fields) < 10:
                qc["imputed_malformed"] += 1
                continue
            chrom, pos_s, ref, alt = fields[0], fields[1], fields[3], fields[4]
            try:
                pos = int(pos_s)
            except ValueError:
                qc["imputed_bad_pos"] += 1
                continue
            if not is_snp(ref, alt):
                qc["imputed_non_biallelic_snp"] += 1
                continue
            key = (chrom, pos)
            if key in seen_imp:
                qc["imputed_duplicate_position"] += 1
                excluded.append({"CHROM": chrom, "POS": pos, "reason": "imputed_duplicate_position"})
                continue
            seen_imp.add(key)
            truth = truth_records.get(key)
            if truth is None:
                qc["not_in_truth"] += 1
                continue
            fmt_keys = fields[8].split(":")
            sample_values = fields[9].split(":")
            fmt_map = dict(zip(fmt_keys, sample_values))
            imp_code = parse_gt_code(fmt_map.get("GT", sample_values[0] if sample_values else None))
            if imp_code is None:
                qc["missing_or_invalid_imputed_gt"] += 1
                excluded.append({"CHROM": chrom, "POS": pos, "reason": "missing_or_invalid_imputed_gt"})
                continue
            gp = parse_gp(fmt_keys, sample_values)
            if gp is None:
                qc["missing_or_invalid_gp"] += 1
                excluded.append({"CHROM": chrom, "POS": pos, "reason": "missing_or_invalid_gp"})
                continue
            swapped = False
            if truth.ref == ref and truth.alt == alt:
                adj_imp_code = imp_code
            elif truth.ref == alt and truth.alt == ref:
                adj_imp_code = 2 - imp_code
                swapped = True
                qc["ref_alt_swapped"] += 1
            else:
                qc["ref_alt_mismatch"] += 1
                excluded.append({
                    "CHROM": chrom, "POS": pos, "reason": "ref_alt_mismatch",
                    "truth_ref": truth.ref, "truth_alt": truth.alt,
                    "imputed_ref": ref, "imputed_alt": alt,
                })
                continue
            sites.append(SiteRecord(
                chrom=truth.chrom,
                pos=truth.pos,
                ref=truth.ref,
                alt=truth.alt,
                true_gt=truth.gt_code,
                imputed_gt=adj_imp_code,
                gp=gp,
                correct=(truth.gt_code == adj_imp_code),
                genotype_class=genotype_class(truth.gt_code),
                swapped=swapped,
            ))
            qc["matched_tested"] += 1
    qc["truth_not_in_imputed"] = len(set(truth_records) - seen_imp)
    return sites, dict(qc), excluded


def parse_float_list(value: str, name: str) -> List[float]:
    try:
        vals = [float(x.strip()) for x in value.split(",") if x.strip() != ""]
    except ValueError:
        raise SystemExit(f"ERROR: {name} must be a comma-separated list of numeric values.")
    if not vals:
        raise SystemExit(f"ERROR: {name} cannot be empty.")
    return vals


def validate_gpth(value: float) -> None:
    if not 0 <= value <= 1:
        raise SystemExit("ERROR: --GPth must be a numeric value between 0 and 1.")


def validate_gpscan(vals: List[float]) -> List[float]:
    for v in vals:
        if not 0 <= v <= 1:
            raise SystemExit("ERROR: --GPscan values must be between 0 and 1.")
    return sorted(set(round(v, 6) for v in vals))


def validate_mafbin(vals: List[float]) -> List[float]:
    if len(vals) < 2:
        raise SystemExit("ERROR: --MAFbin must contain at least two complete boundaries.")
    if abs(vals[0] - 0.0) > 1e-12 or abs(vals[-1] - 0.5) > 1e-12:
        raise SystemExit("ERROR: --MAFbin must provide complete boundaries starting with 0 and ending with 0.5. Example: --MAFbin 0,0.01,0.05,0.5")
    if any(v < 0 or v > 0.5 for v in vals):
        raise SystemExit("ERROR: --MAFbin values must be between 0 and 0.5.")
    if any(vals[i] >= vals[i + 1] for i in range(len(vals) - 1)):
        raise SystemExit("ERROR: --MAFbin values must be strictly increasing.")
    return vals


def validate_shetbin(vals: List[float]) -> List[float]:
    if len(vals) < 2:
        raise SystemExit("ERROR: --sHetbin must contain at least two complete boundaries.")
    if abs(vals[0] - 0.0) > 1e-12 or abs(vals[-1] - 1.0) > 1e-12:
        raise SystemExit("ERROR: --sHetbin must provide complete boundaries starting with 0 and ending with 1. Example: --sHetbin 0,0.2,0.4,1")
    if any(v < 0 or v > 1 for v in vals):
        raise SystemExit("ERROR: --sHetbin values must be between 0 and 1.")
    if any(vals[i] >= vals[i + 1] for i in range(len(vals) - 1)):
        raise SystemExit("ERROR: --sHetbin values must be strictly increasing.")
    return vals


def read_annotation(path: str, value_col_default: str) -> Dict[Tuple[Any, ...], float]:
    data: Dict[Tuple[Any, ...], float] = {}
    with open_text(path) as fh:
        first = fh.readline()
        if not first:
            return data
        first_parts = first.strip().split()
        lower = [x.lower() for x in first_parts]
        has_header = any(x in lower for x in ["chrom", "chr", "#chrom", "pos", "position", value_col_default.lower(), "maf", "shet"])
        rows: List[List[str]] = []
        header: List[str] = []
        if has_header:
            header = lower
        else:
            rows.append(first_parts)
        for line in fh:
            if not line.strip() or line.startswith("#"):
                continue
            rows.append(line.strip().split())
    if has_header:
        def find(names, required=True):
            for n in names:
                if n in header:
                    return header.index(n)
            if required:
                raise SystemExit(f"ERROR: Annotation file {path} is missing required column among: {names}")
            return None
        chrom_i = find(["chrom", "chr", "#chrom"])
        pos_i = find(["pos", "position", "bp"])
        ref_i = find(["ref", "reference"], required=False)
        alt_i = find(["alt", "alternate"], required=False)
        val_i = find([value_col_default.lower(), "maf", "shet", "value", "freq", "af"])
    else:
        chrom_i, pos_i = 0, 1
        if len(rows[0]) >= 5:
            ref_i, alt_i, val_i = 2, 3, 4
        elif len(rows[0]) >= 3:
            ref_i, alt_i, val_i = None, None, 2
        else:
            raise SystemExit(f"ERROR: Annotation file {path} must have at least 3 columns: CHROM POS VALUE.")
    max_needed = max(i for i in [chrom_i, pos_i, val_i] if i is not None)
    for parts in rows:
        if len(parts) <= max_needed:
            continue
        try:
            chrom = parts[chrom_i]
            pos = int(parts[pos_i])
            val = float(parts[val_i])
        except ValueError:
            continue
        if ref_i is not None and alt_i is not None and len(parts) > max(ref_i, alt_i):
            key = (chrom, pos, parts[ref_i], parts[alt_i])
        else:
            key = (chrom, pos)
        data[key] = val
    return data


def bin_label(prefix: str, i: int, total: int) -> str:
    if prefix == "sHet" and total == 3:
        return ["Low", "Medium", "High"][i]
    return f"{prefix}_bin_{i + 1}"


def assign_maf(sites: List[SiteRecord], maf_file: Optional[str], bins: List[float]) -> Dict[str, int]:
    qc = Counter()
    if not maf_file:
        return dict(qc)
    ann = read_annotation(maf_file, "MAF")
    for s in sites:
        val = ann.get((s.chrom, s.pos, s.ref, s.alt))
        if val is None:
            val = ann.get((s.chrom, s.pos, s.alt, s.ref))
        if val is None:
            val = ann.get((s.chrom, s.pos))
        if val is None:
            qc["MAF_missing"] += 1
            continue
        if 0.5 < val <= 1:
            val = min(val, 1 - val)
        s.maf = val
        for i, (lo, hi) in enumerate(zip(bins[:-1], bins[1:])):
            if lo < val <= hi:
                s.maf_bin = f"MAF_bin_{i + 1}"
                s.maf_range = f"({lo:g}, {hi:g}]"
                qc["MAF_assigned"] += 1
                break
        if s.maf_bin is None:
            qc["MAF_outside_bins_or_zero"] += 1
    return dict(qc)


def assign_shet(sites: List[SiteRecord], shet_file: Optional[str], bins: List[float]) -> Dict[str, int]:
    qc = Counter()
    if not shet_file:
        return dict(qc)
    ann = read_annotation(shet_file, "sHet")
    total = len(bins) - 1
    for s in sites:
        val = ann.get((s.chrom, s.pos, s.ref, s.alt))
        if val is None:
            val = ann.get((s.chrom, s.pos, s.alt, s.ref))
        if val is None:
            val = ann.get((s.chrom, s.pos))
        if val is None:
            qc["sHet_missing"] += 1
            continue
        s.shet = val
        for i, (lo, hi) in enumerate(zip(bins[:-1], bins[1:])):
            if (i == 0 and lo <= val <= hi) or (i > 0 and lo < val <= hi):
                s.shet_bin = bin_label("sHet", i, total)
                s.shet_range = f"[{lo:g}, {hi:g}]" if i == 0 else f"({lo:g}, {hi:g}]"
                qc["sHet_assigned"] += 1
                break
        if s.shet_bin is None:
            qc["sHet_outside_bins"] += 1
    return dict(qc)


def safe_div(a: float, b: float) -> Optional[float]:
    return None if b == 0 else a / b


def pearson(x: List[float], y: List[float]) -> Optional[float]:
    if len(x) < 2 or len(y) < 2:
        return None
    mx, my = sum(x) / len(x), sum(y) / len(y)
    vx = sum((v - mx) ** 2 for v in x)
    vy = sum((v - my) ** 2 for v in y)
    if vx <= 0 or vy <= 0:
        return None
    cov = sum((a - mx) * (b - my) for a, b in zip(x, y))
    return cov / math.sqrt(vx * vy)


def harmonic_mean(vals: List[Optional[float]]) -> Optional[float]:
    if any(v is None for v in vals):
        return None
    vals2 = [max(0.0, float(v)) for v in vals]
    if any(v == 0 for v in vals2):
        return 0.0
    return len(vals2) / sum(1 / v for v in vals2)


def compute_metrics(sites: List[SiteRecord], gpth: float, analysis: str, group: str) -> Dict[str, Any]:
    n_tested = len(sites)
    retained = [s for s in sites if s.gp >= gpth]
    tp = sum(1 for s in sites if s.gp >= gpth and s.correct)
    fp = sum(1 for s in sites if s.gp >= gpth and not s.correct)
    fn = sum(1 for s in sites if s.gp < gpth and s.correct)
    tn = sum(1 for s in sites if s.gp < gpth and not s.correct)
    ip = safe_div(tp, tp + fp)
    gkr = safe_div(tp + fp, n_tested)
    ca = safe_div(tp + tn, n_tested)
    tpr = safe_div(tp, tp + fn)
    tnr = safe_div(tn, tn + fp)
    fpr = safe_div(fp, fp + tn)
    fnr = safe_div(fn, tp + fn)
    fdr = safe_div(fp, tp + fp)
    cor = pearson([s.true_gt for s in retained], [s.imputed_gt for s in retained])
    cis = harmonic_mean([cor, ip, gkr, ca, tpr])
    return {
        "Analysis": analysis,
        "Group": group,
        "GPth": gpth,
        "N_tested": n_tested,
        "N_retained": tp + fp,
        "N_retained_correct": tp,
        "N_incorrect_retained": fp,
        "COR": cor,
        "IP": ip,
        "GKR": gkr,
        "CA": ca,
        "TP": tp,
        "FP": fp,
        "FN": fn,
        "TN": tn,
        "TPR": tpr,
        "TNR": tnr,
        "FPR": fpr,
        "FNR": fnr,
        "FDR": fdr,
        "CIS": cis,
    }


def group_by_genotype(sites: List[SiteRecord]) -> Dict[str, List[SiteRecord]]:
    return {
        "Homozygous": [s for s in sites if s.genotype_class == "Homozygous"],
        "Heterozygous": [s for s in sites if s.genotype_class == "Heterozygous"],
    }


def group_by_attr(sites: List[SiteRecord], attr: str) -> Dict[str, List[SiteRecord]]:
    out = defaultdict(list)
    for s in sites:
        key = getattr(s, attr)
        if key is not None:
            out[key].append(s)
    return dict(out)


def parse_region(region: str) -> Tuple[str, int, int, str]:
    try:
        chrom, coords = region.split(":")
        start_s, end_s = coords.replace(",", "").split("-")
        start, end = int(start_s), int(end_s)
    except Exception:
        raise SystemExit("ERROR: --region must be formatted as CHR:START-END, e.g. Chr1:1000000-2000000")
    if start >= end:
        raise SystemExit("ERROR: --region start must be smaller than end.")
    return chrom, start, end, f"{chrom}:{start}-{end}"


def read_regions(region: Optional[str], region_bed: Optional[str]) -> List[Tuple[str, int, int, str]]:
    regions: List[Tuple[str, int, int, str]] = []
    if region:
        regions.append(parse_region(region))
    if region_bed:
        with open_text(region_bed) as fh:
            for line in fh:
                if not line.strip() or line.startswith("#"):
                    continue
                parts = line.strip().split()
                if len(parts) < 3:
                    continue
                chrom, start, end = parts[0], int(parts[1]), int(parts[2])
                name = f"{chrom}:{start}-{end}"
                regions.append((chrom, start, end, name))
    return regions


def compute_region_metrics(sites: List[SiteRecord], regions: List[Tuple[str, int, int, str]], gpth: float) -> List[Dict[str, Any]]:
    out = []
    for chrom, start, end, name in regions:
        sub = [s for s in sites if s.chrom == chrom and start <= s.pos <= end]
        m = compute_metrics(sub, gpth, "Region", name)
        m.update({"CHROM": chrom, "Start": start, "End": end, "Region": name})
        out.append(m)
    return out


def compute_whet_windows(sites: List[SiteRecord], window_size: int, window_step: int, min_sites: int, gpth: float) -> List[Dict[str, Any]]:
    by_chrom = defaultdict(list)
    for s in sites:
        by_chrom[s.chrom].append(s)
    out: List[Dict[str, Any]] = []
    for chrom, chrom_sites in sorted(by_chrom.items()):
        chrom_sites = sorted(chrom_sites, key=lambda s: s.pos)
        if not chrom_sites:
            continue
        max_pos = chrom_sites[-1].pos
        start = 0
        while start <= max_pos:
            end = start + window_size
            sub = [s for s in chrom_sites if start <= s.pos < end]
            if len(sub) >= min_sites:
                n_het = sum(1 for s in sub if s.true_gt == 1)
                whet = n_het / len(sub)
                m = compute_metrics(sub, gpth, "wHet", f"{chrom}:{start}-{end}")
                m.update({
                    "CHROM": chrom,
                    "Start": start,
                    "End": end,
                    "N_truth_callable": len(sub),
                    "N_truth_heterozygous": n_het,
                    "wHet": whet,
                })
                out.append(m)
            start += window_step
    return out


PERCENT_METRICS = {"IP", "GKR", "CA", "TPR", "TNR", "FPR", "FNR", "FDR"}
DECIMAL4_METRICS = {"COR", "CIS"}


def fmt_num(v: Any, digits: int = 6) -> str:
    if v is None:
        return NA
    if isinstance(v, str):
        return v
    if isinstance(v, int):
        return str(v)
    if isinstance(v, float):
        if math.isnan(v) or math.isinf(v):
            return NA
        return f"{v:.{digits}g}"
    return str(v)


def fmt_metric(key: str, v: Any, digits: int = 6, html_display: bool = True) -> str:
    """Format values for report tables. Internal calculations stay as 0-1 fractions."""
    if v is None:
        return NA
    if isinstance(v, str):
        return v
    if isinstance(v, bool):
        return str(int(v))
    if isinstance(v, int):
        return str(v)
    if isinstance(v, float):
        if math.isnan(v) or math.isinf(v):
            return NA
        if html_display and key in PERCENT_METRICS:
            return f"{v * 100:.2f}%"
        if html_display and key in DECIMAL4_METRICS:
            return f"{v:.4f}"
        if key == "GPth":
            return f"{v:.4g}"
        return f"{v:.{digits}g}"
    return str(v)


def html_escape(x: Any) -> str:
    return html.escape(str(x))


def html_table(headers: List[str], rows: List[List[Any]]) -> str:
    thead = "".join(f"<th>{html_escape(h)}</th>" for h in headers)
    trs = []
    for row in rows:
        trs.append("<tr>" + "".join(f"<td>{html_escape(x)}</td>" for x in row) + "</tr>")
    return f"<div class='table-wrap'><table><thead><tr>{thead}</tr></thead><tbody>{''.join(trs)}</tbody></table></div>"


def metrics_rows_to_table(rows: List[Dict[str, Any]], keys: List[str]) -> str:
    return html_table(keys, [[fmt_metric(k, r.get(k), 6) for k in keys] for r in rows])


def vertical_metric_table(m: Dict[str, Any], keys: List[str]) -> str:
    return html_table(["Metric", "Value"], [[k, fmt_metric(k, m.get(k), 6)] for k in keys])


def confusion_matrix_html(m: Dict[str, Any], title: str) -> str:
    total = m.get("N_tested", 0) or 0
    def cell(k: str) -> str:
        v = m.get(k, 0) or 0
        pct = (v / total * 100) if total else 0
        return f"<strong>{v}</strong><br><span>{pct:.2f}%</span>"
    return f"""
<div class='matrix-card'>
  <h4>{html_escape(title)}</h4>
  <table class='matrix'>
    <thead><tr><th></th><th>Retained<br>GP ≥ {fmt_num(m.get('GPth'))}</th><th>Filtered<br>GP &lt; {fmt_num(m.get('GPth'))}</th></tr></thead>
    <tbody>
      <tr><th>Correct</th><td class='tp'>{cell('TP')}</td><td class='fn'>{cell('FN')}</td></tr>
      <tr><th>Incorrect</th><td class='fp'>{cell('FP')}</td><td class='tn'>{cell('TN')}</td></tr>
    </tbody>
  </table>
</div>
"""


def make_svg_bar(rows: List[Dict[str, Any]], x_key: str, y_key: str, title: str, y_label: str, percent: bool = True, width: int = 680, height: int = 310) -> str:
    vals = []
    for r in rows:
        y = r.get(y_key)
        if y is None:
            continue
        vals.append((str(r.get(x_key, "")), y * 100 if percent else y))
    if not vals:
        return f"<div class='note'>No valid values for {html_escape(y_key)}.</div>"
    ml, mr, mt, mb = 60, 25, 35, 75
    if percent:
        min_y, max_y = 0.0, 100.0
    else:
        min_y, max_y = 0.0, max(v for _, v in vals) or 1.0
        max_y *= 1.15
        if max_y == min_y:
            max_y = min_y + 1.0
    bar_area = width - ml - mr
    gap = bar_area / max(1, len(vals))
    bw = gap * 0.62
    def sy(y): return height - mb - (y - min_y) / (max_y - min_y) * (height - mt - mb)
    grid = []
    for i in range(5):
        yy = min_y + i * (max_y - min_y) / 4
        yp = sy(yy)
        grid.append(f"<line class='grid' x1='{ml}' x2='{width-mr}' y1='{yp:.2f}' y2='{yp:.2f}'/>")
        grid.append(f"<text class='tick' x='{ml-7}' y='{yp+4:.2f}' text-anchor='end'>{yy:.0f}</text>" if percent else f"<text class='tick' x='{ml-7}' y='{yp+4:.2f}' text-anchor='end'>{yy:.2g}</text>")
    bars = []
    for i, (name, val) in enumerate(vals):
        cx = ml + gap * i + gap / 2
        x = cx - bw / 2
        y = sy(val)
        h = height - mb - y
        bars.append(f"<rect x='{x:.2f}' y='{y:.2f}' width='{bw:.2f}' height='{h:.2f}' rx='4'><title>{html_escape(name)}: {val:.3g}</title></rect>")
        bars.append(f"<text class='tick' x='{cx:.2f}' y='{height-mb+18}' text-anchor='middle'>{html_escape(name)}</text>")
        bars.append(f"<text class='tick' x='{cx:.2f}' y='{max(12, y-5):.2f}' text-anchor='middle'>{val:.2f}</text>")
    return f"""
<div class='chart-card'><h4>{html_escape(title)}</h4>
<svg viewBox='0 0 {width} {height}'>
<style>.axis{{stroke:#333;stroke-width:1.1}} .grid{{stroke:#e7e7e7}} rect{{fill:#8e7cc3}} .tick{{font-size:11px;fill:#444}} .label{{font-size:13px;font-weight:600;fill:#222}}</style>
{''.join(grid)}<line class='axis' x1='{ml}' x2='{width-mr}' y1='{height-mb}' y2='{height-mb}'/><line class='axis' x1='{ml}' x2='{ml}' y1='{mt}' y2='{height-mb}'/>{''.join(bars)}
<text class='label' transform='translate(16,{height/2}) rotate(-90)' text-anchor='middle'>{html_escape(y_label)}</text>
</svg></div>"""



def make_svg_index_line(rows: List[Dict[str, Any]], y_key: str, title: str, y_label: str, percent: bool = True, width: int = 720, height: int = 310) -> str:
    vals = []
    for r in rows:
        y = r.get(y_key)
        x = r.get("bin_index")
        if x is None or y is None:
            continue
        label = str(r.get("Group", x)).replace("MAF_bin_", "")
        tooltip = str(r.get("MAF_range", r.get("Group", label)))
        vals.append((int(x), float(y) * 100 if percent else float(y), label, tooltip))
    if not vals:
        return f"<div class='note'>No valid values for {html_escape(y_key)}.</div>"
    vals = sorted(vals, key=lambda x: x[0])
    ml, mr, mt, mb = 62, 26, 35, 62
    ys = [v for _, v, _, _ in vals]
    if percent:
        min_y, max_y = 0.0, 100.0
    else:
        min_y, max_y = min(ys), max(ys)
        if min_y == max_y:
            min_y -= 1; max_y += 1
        ypad = (max_y - min_y) * 0.08
        min_y -= ypad; max_y += ypad
    min_x = min(x for x, _, _, _ in vals)
    max_x = max(x for x, _, _, _ in vals)
    if min_x == max_x:
        min_x -= 1; max_x += 1
    def sx(x): return ml + (x - min_x) / (max_x - min_x) * (width - ml - mr)
    def sy(y): return height - mb - (y - min_y) / (max_y - min_y) * (height - mt - mb)
    grid = []
    for i in range(5):
        yy = min_y + i * (max_y - min_y) / 4
        yp = sy(yy)
        grid.append(f"<line class='grid' x1='{ml}' x2='{width-mr}' y1='{yp:.2f}' y2='{yp:.2f}'/>")
        grid.append(f"<text class='tick' x='{ml-7}' y='{yp+4:.2f}' text-anchor='end'>{yy:.0f}</text>" if percent else f"<text class='tick' x='{ml-7}' y='{yp+4:.2f}' text-anchor='end'>{yy:.2g}</text>")
    for x, _, label, _ in vals:
        xp = sx(x)
        grid.append(f"<text class='tick' x='{xp:.2f}' y='{height-mb+17}' text-anchor='middle'>{html_escape(label)}</text>")
    path_d = " ".join(("M" if i == 0 else "L") + f"{sx(x):.2f},{sy(y):.2f}" for i, (x, y, _, _) in enumerate(vals))
    path = f"<path class='line' d='{path_d}'/>"
    circles = "".join(f"<circle cx='{sx(x):.2f}' cy='{sy(y):.2f}' r='3.5'><title>MAF bin {html_escape(label)} {html_escape(tooltip)}, {html_escape(y_label)}={y:.4g}</title></circle>" for x, y, label, tooltip in vals)
    return f"""
<div class='chart-card'><h4>{html_escape(title)}</h4>
<svg viewBox='0 0 {width} {height}'>
<style>.axis{{stroke:#333;stroke-width:1.1}} .grid{{stroke:#e7e7e7}} .line{{fill:none;stroke:#4a6fa5;stroke-width:2.2}} circle{{fill:#4a6fa5;opacity:.82}} .tick{{font-size:11px;fill:#444}} .label{{font-size:13px;font-weight:600;fill:#222}}</style>
{''.join(grid)}<line class='axis' x1='{ml}' x2='{width-mr}' y1='{height-mb}' y2='{height-mb}'/><line class='axis' x1='{ml}' x2='{ml}' y1='{mt}' y2='{height-mb}'/>{path}{circles}
<text class='label' x='{width/2}' y='{height-8}' text-anchor='middle'>MAF bin index</text>
<text class='label' transform='translate(16,{height/2}) rotate(-90)' text-anchor='middle'>{html_escape(y_label)}</text>
</svg></div>"""

def make_svg_line(rows: List[Dict[str, Any]], x_key: str, y_key: str, title: str, y_label: str, percent: bool = True, width: int = 720, height: int = 310, zoom_percent: bool = False, x_min_limit: Optional[float] = None, x_max_limit: Optional[float] = None) -> str:
    pts = []
    for r in rows:
        x, y = r.get(x_key), r.get(y_key)
        if x is None or y is None:
            continue
        pts.append((float(x), y * 100 if percent else float(y)))
    if not pts:
        return f"<div class='note'>No valid values for {html_escape(y_key)}.</div>"
    pts = sorted(pts, key=lambda p: p[0])
    return make_xy_svg(
        pts,
        title,
        x_label=x_key,
        y_label=y_label,
        connect=True,
        width=width,
        height=height,
        percent=percent,
        zoom_percent=zoom_percent,
        x_min_limit=x_min_limit,
        x_max_limit=x_max_limit,
    )


def make_svg_scatter(rows: List[Dict[str, Any]], x_key: str, y_key: str, title: str, x_label: str, y_label: str, percent: bool = True, width: int = 720, height: int = 310, zoom_percent: bool = False, x_min_limit: Optional[float] = None, x_max_limit: Optional[float] = None) -> str:
    pts = []
    for r in rows:
        x, y = r.get(x_key), r.get(y_key)
        if x is None or y is None:
            continue
        pts.append((float(x), y * 100 if percent else float(y)))
    if not pts:
        return f"<div class='note'>No valid values for {html_escape(y_key)}.</div>"
    return make_xy_svg(
        pts,
        title,
        x_label=x_label,
        y_label=y_label,
        connect=False,
        width=width,
        height=height,
        percent=percent,
        zoom_percent=zoom_percent,
        x_min_limit=x_min_limit,
        x_max_limit=x_max_limit,
    )


def nice_tick_label(value: float, percent: bool = False) -> str:
    if percent:
        # Keep zoomed percentage axes readable without scientific notation.
        if abs(value - round(value)) < 1e-8:
            return f"{value:.0f}"
        return f"{value:.1f}" if abs(value) >= 10 else f"{value:.2f}".rstrip('0').rstrip('.')
    return f"{value:.3g}"


def make_xy_svg(points: List[Tuple[float, float]], title: str, x_label: str, y_label: str, connect: bool, width: int, height: int, percent: bool = False, zoom_percent: bool = False, x_min_limit: Optional[float] = None, x_max_limit: Optional[float] = None) -> str:
    ml, mr, mt, mb = 62, 26, 35, 52
    xs = [x for x, _ in points]
    ys = [y for _, y in points]
    min_x, max_x = min(xs), max(xs)
    if min_x == max_x:
        xpad = max(abs(min_x) * 0.05, 0.01)
    else:
        xpad = (max_x - min_x) * 0.03
    min_x -= xpad
    max_x += xpad
    if x_min_limit is not None:
        min_x = max(x_min_limit, min_x)
    if x_max_limit is not None:
        max_x = min(x_max_limit, max_x)
    if min_x == max_x:
        min_x -= 0.01
        max_x += 0.01
        if x_min_limit is not None:
            min_x = max(x_min_limit, min_x)
        if x_max_limit is not None:
            max_x = min(x_max_limit, max_x)

    if percent and not zoom_percent:
        min_y, max_y = 0.0, 100.0
    else:
        min_y, max_y = min(ys), max(ys)
        if min_y == max_y:
            ypad = max(abs(min_y) * 0.05, 1.0 if percent else 0.01)
        else:
            ypad = (max_y - min_y) * 0.12
        min_y -= ypad
        max_y += ypad
        if percent:
            min_y = max(0.0, min_y)
            max_y = min(100.0, max_y)
        if min_y == max_y:
            if percent:
                min_y = max(0.0, min_y - 1.0)
                max_y = min(100.0, max_y + 1.0)
            else:
                min_y -= 1.0
                max_y += 1.0
    if min_y == max_y:
        max_y = min_y + 1.0

    def sx(x): return ml + (x - min_x) / (max_x - min_x) * (width - ml - mr)
    def sy(y): return height - mb - (y - min_y) / (max_y - min_y) * (height - mt - mb)
    grid = []
    for i in range(5):
        yy = min_y + i * (max_y - min_y) / 4
        yp = sy(yy)
        grid.append(f"<line class='grid' x1='{ml}' x2='{width-mr}' y1='{yp:.2f}' y2='{yp:.2f}'/>")
        grid.append(f"<text class='tick' x='{ml-7}' y='{yp+4:.2f}' text-anchor='end'>{nice_tick_label(yy, percent)}</text>")
    for i in range(5):
        xx = min_x + i * (max_x - min_x) / 4
        xp = sx(xx)
        grid.append(f"<text class='tick' x='{xp:.2f}' y='{height-mb+17}' text-anchor='middle'>{xx:.3g}</text>")
    path = ""
    if connect:
        path_d = " ".join(("M" if i == 0 else "L") + f"{sx(x):.2f},{sy(y):.2f}" for i, (x, y) in enumerate(points))
        path = f"<path class='line' d='{path_d}'/>"
    circles = "".join(f"<circle cx='{sx(x):.2f}' cy='{sy(y):.2f}' r='3.5'><title>{html_escape(x_label)}={x:.4g}, {html_escape(y_label)}={y:.4g}</title></circle>" for x, y in points)
    return f"""
<div class='chart-card'><h4>{html_escape(title)}</h4>
<svg viewBox='0 0 {width} {height}'>
<style>.axis{{stroke:#333;stroke-width:1.1}} .grid{{stroke:#e7e7e7}} .line{{fill:none;stroke:#4a6fa5;stroke-width:2.2}} circle{{fill:#4a6fa5;opacity:.82}} .tick{{font-size:11px;fill:#444}} .label{{font-size:13px;font-weight:600;fill:#222}}</style>
{''.join(grid)}<line class='axis' x1='{ml}' x2='{width-mr}' y1='{height-mb}' y2='{height-mb}'/><line class='axis' x1='{ml}' x2='{ml}' y1='{mt}' y2='{height-mb}'/>{path}{circles}
<text class='label' x='{width/2}' y='{height-8}' text-anchor='middle'>{html_escape(x_label)}</text>
<text class='label' transform='translate(16,{height/2}) rotate(-90)' text-anchor='middle'>{html_escape(y_label)}</text>
</svg></div>"""


def transpose_scan_table(scan_rows: List[Dict[str, Any]], metrics: List[str]) -> List[List[str]]:
    header = ["Metric"] + [f"GP={r['GPth']:g}" for r in scan_rows]
    table = [header]
    for metric in metrics:
        table.append([metric] + [fmt_metric(metric, r.get(metric), 6) for r in scan_rows])
    return table


def site_level_table(path: Path, sites: List[SiteRecord]) -> None:
    with gzip.open(path, "wt") as fh:
        cols = ["CHROM", "POS", "REF", "ALT", "true_GT", "imputed_GT", "GP", "correct", "genotype_class", "MAF", "MAF_bin", "sHet", "sHet_bin"]
        fh.write("\t".join(cols) + "\n")
        for s in sites:
            vals = [s.chrom, s.pos, s.ref, s.alt, s.true_gt, s.imputed_gt, s.gp, int(s.correct), s.genotype_class, s.maf, s.maf_bin, s.shet, s.shet_bin]
            fh.write("\t".join(fmt_num(v, 6) for v in vals) + "\n")


def write_tsv(path: Path, rows: List[Dict[str, Any]]) -> None:
    if not rows:
        return
    keys = list(rows[0].keys())
    for r in rows[1:]:
        for k in r.keys():
            if k not in keys:
                keys.append(k)
    with open(path, "w") as fh:
        fh.write("\t".join(keys) + "\n")
        for r in rows:
            fh.write("\t".join(fmt_num(r.get(k), 10) for k in keys) + "\n")


def build_html_report(results: Dict[str, Any], outpath: Path) -> None:
    css = """
    body{font-family:Arial,Helvetica,sans-serif;margin:0;background:#f7f8fb;color:#1f2937} main{max-width:1180px;margin:0 auto;padding:26px}
    header{background:#111827;color:white;padding:28px 34px} h1{margin:0;font-size:30px}.subtitle{opacity:.85;margin-top:8px}
    section{background:white;border:1px solid #e5e7eb;border-radius:16px;padding:22px;margin:22px 0;box-shadow:0 1px 4px rgba(0,0,0,.04)}
    h2{margin-top:0;border-bottom:1px solid #e5e7eb;padding-bottom:9px} h3{margin-top:24px}.note{color:#5b6472;font-size:14px}.warn{background:#fff7ed;border:1px solid #fed7aa;padding:10px;border-radius:8px;margin:8px 0}
    .cards{display:grid;grid-template-columns:repeat(3,220px);gap:18px;justify-content:center;align-items:stretch;margin:24px auto}.card{background:#f9fafb;border:1px solid #e5e7eb;border-radius:12px;padding:16px 18px;min-height:82px}.card .k{font-size:13px;color:#6b7280}.card .v{font-size:24px;font-weight:700;margin-top:8px}@media (max-width:820px){.cards{grid-template-columns:1fr}}
    .table-wrap{overflow:auto;margin:12px 0} table{border-collapse:collapse;width:100%;font-size:13px} th,td{border:1px solid #e5e7eb;padding:7px 9px;text-align:right;white-space:nowrap} th:first-child,td:first-child{text-align:left} th{background:#f3f4f6}
    .matrix-wrap{display:grid;grid-template-columns:repeat(auto-fit,minmax(280px,1fr));gap:16px}.matrix-card{border:1px solid #e5e7eb;border-radius:12px;padding:14px;background:#fcfcfd}.matrix td,.matrix th{text-align:center}.matrix .tp{background:#dcfce7}.matrix .fp{background:#fee2e2}.matrix .fn{background:#fef3c7}.matrix .tn{background:#e0f2fe}.matrix span{font-size:12px;color:#4b5563}
    .chart-grid{display:grid;grid-template-columns:repeat(auto-fit,minmax(380px,1fr));gap:16px}.chart-card{border:1px solid #e5e7eb;border-radius:14px;padding:14px;background:#fff} svg{width:100%;height:auto} footer{text-align:center;color:#6b7280;padding:20px}
    """
    parts = ["<!doctype html><html><head><meta charset='utf-8'>", "<title>IPEva (Imputation Performance Evaluator) Report</title>", f"<style>{css}</style></head><body>"]
    meta = results["metadata"]
    parts.append(f"<header><h1>IPEva (Imputation Performance Evaluator) Report</h1><div class='subtitle'>Generated {html_escape(results['generated_at'])}</div></header><main>")

    parts.append("<section><h2>1. Report overview</h2>")
    for w in results.get("warnings", []):
        parts.append(f"<div class='warn'>{html_escape(w)}</div>")
    overview_keys = ["Mode", "Site_analysis", "Interval_analysis", "GPth"]
    if meta.get("Mode") == "all" or "GP" in str(meta.get("Site_analysis", "")):
        overview_keys.append("GPscan")
    overview_keys += ["Truth_sample_name", "Imputed_sample_name", "Truth_VCF", "Imputed_VCF"]
    parts.append(vertical_metric_table(meta, overview_keys))
    parts.append("</section>")

    parts.append("<section><h2>2. SNP matching summary</h2>")
    qc_rows = [[k, v] for k, v in sorted(results.get("matching_qc", {}).items())]
    parts.append(html_table(["Item", "Count"], qc_rows))
    parts.append("</section>")

    if results.get("overall"):
        m = results["overall"]
        parts.append("<section><h2>3. Overall evaluation</h2>")
        card_keys = ["N_tested", "N_retained", "N_retained_correct", "COR", "IP", "GKR", "CA", "TPR", "CIS"]
        parts.append("<div class='cards'>" + "".join(f"<div class='card'><div class='k'>{k}</div><div class='v'>{fmt_metric(k, m.get(k), 5)}</div></div>" for k in card_keys) + "</div>")
        parts.append(vertical_metric_table(m, ["GPth", "N_tested", "N_retained", "N_retained_correct", "N_incorrect_retained", "COR", "IP", "GKR", "CA", "TPR", "TNR", "FPR", "FNR", "FDR", "CIS"]))
        parts.append("<h3>2×2 confusion matrix</h3><div class='matrix-wrap'>")
        parts.append(confusion_matrix_html(m, "Overall"))
        parts.append("</div></section>")

    if results.get("site"):
        parts.append("<section><h2>Site-stratified evaluation</h2>")
        site = results["site"]
        if site.get("Genotype"):
            rows = site["Genotype"]
            parts.append("<h3>Genotype-stratified evaluation</h3>")
            parts.append("<p class='note'>Heterozygous sites have no variance in the truth genotype column because all true genotypes are coded as 1. Therefore COR and CIS are reported as NA for the Heterozygous group.</p>")
            parts.append(metrics_rows_to_table(rows, ["Group", "N_tested", "N_retained", "N_retained_correct", "COR", "IP", "GKR", "CA", "TPR", "TNR", "FPR", "FNR", "FDR", "CIS"]))
            parts.append("<div class='matrix-wrap'>" + "".join(confusion_matrix_html(r, r["Group"]) for r in rows) + "</div>")
            parts.append("<div class='chart-grid'>")
            for metric in ["IP", "GKR", "CA", "TPR"]:
                parts.append(make_svg_bar(rows, "Group", metric, f"Genotype vs {metric}", f"{metric} (%)"))
            parts.append("</div>")
        if site.get("MAF"):
            rows = site["MAF"]
            parts.append("<h3>MAF-stratified evaluation</h3>")
            parts.append(metrics_rows_to_table(rows, ["Group", "MAF_range", "N_tested", "N_retained", "N_retained_correct", "COR", "IP", "GKR", "CA", "TPR", "TNR", "FPR", "FDR", "CIS"]))
            parts.append("<div class='matrix-wrap'>" + "".join(confusion_matrix_html(r, r["Group"]) for r in rows) + "</div>")
            parts.append("<div class='chart-grid'>")
            for metric in ["IP", "GKR", "CA", "TPR", "COR", "CIS"]:
                parts.append(make_svg_index_line(rows, metric, f"MAF bin vs {metric}", f"{metric} (%)" if metric != "COR" else "COR", percent=(metric != "COR")))
            parts.append("</div>")
        if site.get("sHet"):
            rows = site["sHet"]
            parts.append("<h3>sHet-stratified evaluation</h3>")
            parts.append(metrics_rows_to_table(rows, ["Group", "sHet_range", "N_tested", "N_retained", "N_retained_correct", "COR", "IP", "GKR", "CA", "TPR", "TNR", "FPR", "FDR", "CIS"]))
            parts.append("<div class='matrix-wrap'>" + "".join(confusion_matrix_html(r, r["Group"]) for r in rows) + "</div>")
            parts.append("<div class='chart-grid'>")
            for metric in ["IP", "GKR", "CA", "TPR", "COR", "CIS"]:
                parts.append(make_svg_bar(rows, "Group", metric, f"sHet bin vs {metric}", f"{metric} (%)" if metric != "COR" else "COR", percent=(metric != "COR")))
            parts.append("</div>")
        if site.get("GP"):
            rows = site["GP"]
            parts.append("<h3>GP threshold sensitivity analysis</h3>")
            scan_table = transpose_scan_table(rows, ["N_retained", "N_retained_correct", "COR", "IP", "GKR", "CA", "TPR", "TNR", "FPR", "FNR", "FDR", "CIS", "TP", "FP", "FN", "TN"])
            parts.append(html_table(scan_table[0], scan_table[1:]))
            parts.append("<div class='matrix-wrap'>" + "".join(confusion_matrix_html(r, r["Group"]) for r in rows) + "</div>")
            parts.append("<div class='chart-grid'>")
            for metric in ["IP", "GKR", "CA", "TPR", "COR", "CIS"]:
                parts.append(make_svg_line(rows, "GPth", metric, f"GP threshold vs {metric}", f"{metric} (%)" if metric != "COR" else "COR", percent=(metric != "COR"), zoom_percent=True, x_min_limit=0, x_max_limit=1))
            parts.append("</div>")
        parts.append("</section>")

    if results.get("interval"):
        parts.append("<section><h2>Interval-based evaluation</h2>")
        interval = results["interval"]
        if interval.get("Region"):
            rows = interval["Region"]
            parts.append("<h3>User-defined region evaluation</h3>")
            parts.append(metrics_rows_to_table(rows, ["Region", "CHROM", "Start", "End", "N_tested", "N_retained", "N_retained_correct", "COR", "IP", "GKR", "CA", "TPR", "CIS"]))
            parts.append("<div class='matrix-wrap'>" + "".join(confusion_matrix_html(r, r["Region"]) for r in rows) + "</div>")
        if interval.get("wHet"):
            rows = interval["wHet"]
            shown = rows[:200]
            parts.append("<h3>Window heterozygosity evaluation</h3>")
            parts.append("<p class='note'>wHet is calculated as the proportion of truth heterozygous genotypes among comparable SNPs in each window. Scatter plots do not connect points because each window is an independent observation.</p>")
            parts.append(metrics_rows_to_table(shown, ["CHROM", "Start", "End", "N_tested", "N_retained", "N_truth_heterozygous", "wHet", "COR", "IP", "GKR", "CA", "TPR", "CIS"]))
            if len(rows) > 200:
                parts.append("<p class='note'>Only the first 200 windows are displayed in HTML. Use --export-tables to export all windows.</p>")
            parts.append("<div class='chart-grid'>")
            for metric in ["IP", "GKR", "CA", "TPR", "COR", "CIS"]:
                parts.append(make_svg_scatter(rows, "wHet", metric, f"wHet vs {metric}", "wHet", f"{metric} (%)" if metric != "COR" else "COR", percent=(metric != "COR"), zoom_percent=True, x_min_limit=0, x_max_limit=1))
            parts.append("</div>")
        parts.append("</section>")

    parts.append("</main></body></html>")
    outpath.write_text("".join(parts), encoding="utf-8")


def parse_site_analysis(value: str) -> List[str]:
    allowed = {"Genotype", "MAF", "sHet", "GP"}
    items = [x.strip() for x in value.split(",") if x.strip()]
    if not items:
        raise SystemExit("ERROR: --site-analysis cannot be empty.")
    bad = [x for x in items if x not in allowed]
    if bad:
        raise SystemExit(f"ERROR: invalid --site-analysis value(s): {', '.join(bad)}. Allowed values are case-sensitive: Genotype,MAF,sHet,GP")
    return items


def parse_interval_analysis(value: str) -> List[str]:
    allowed = {"Region", "wHet"}
    items = [x.strip() for x in value.split(",") if x.strip()]
    if not items:
        raise SystemExit("ERROR: --interval-analysis cannot be empty.")
    bad = [x for x in items if x not in allowed]
    if bad:
        raise SystemExit(f"ERROR: invalid --interval-analysis value(s): {', '.join(bad)}. Allowed values are case-sensitive: Region,wHet")
    return items


def export_tables(outdir: Path, results: Dict[str, Any], sites: List[SiteRecord], excluded: List[Dict[str, Any]], export_site: bool) -> None:
    tables = outdir / "tables"
    tables.mkdir(parents=True, exist_ok=True)
    if results.get("overall"):
        write_tsv(tables / "overall.tsv", [results["overall"]])
        write_tsv(tables / "overall_confusion_matrix.tsv", [{k: results["overall"].get(k) for k in ["Analysis", "Group", "GPth", "TP", "FP", "FN", "TN", "N_tested"]}])
    for section_name in ["site", "interval"]:
        if results.get(section_name):
            for name, rows in results[section_name].items():
                if rows:
                    write_tsv(tables / f"{section_name}_{name}.tsv", rows)
                    conf_rows = [{k: r.get(k) for k in ["Analysis", "Group", "GPth", "TP", "FP", "FN", "TN", "N_tested"] if k in r} for r in rows]
                    write_tsv(tables / f"{section_name}_{name}_confusion_matrix.tsv", conf_rows)
    if excluded:
        write_tsv(tables / "excluded_sites.tsv", excluded)
    if export_site:
        site_level_table(tables / "site_level.tsv.gz", sites)


def main() -> None:
    parser = argparse.ArgumentParser(description="IPEva: single-sample truth-based SNP genotype imputation performance evaluator.")
    parser.add_argument("--truth", required=True, help="Single-sample truth VCF/VCF.GZ from hcWGS genotype calls.")
    parser.add_argument("--imputed", required=True, help="Single-sample imputed VCF/VCF.GZ containing GT and GP fields.")
    parser.add_argument("--mode", default="overall", choices=["overall", "site", "interval", "all"], help="Evaluation mode. Default: overall.")
    parser.add_argument("--site-analysis", default="Genotype", help="Case-sensitive comma-separated site analyses: Genotype,MAF,sHet,GP. Default: Genotype.")
    parser.add_argument("--interval-analysis", default="wHet", help="Case-sensitive comma-separated interval analyses: Region,wHet. Default: wHet.")
    parser.add_argument("--GPth", type=float, default=DEFAULT_GPTH, help="GP filtering threshold for main analyses. Default: 0.9.")
    parser.add_argument("--GPscan", default=",".join(f"{x:.2f}" for x in DEFAULT_GPSCAN), help="Comma-separated GP thresholds for site-analysis GP. Default: 0.90,0.91,...,0.99.")
    parser.add_argument("--MAFfile", help="MAF annotation table with CHROM POS [REF ALT] MAF.")
    parser.add_argument("--MAFbin", default=",".join(map(str, DEFAULT_MAFBIN)), help="Complete MAF bin boundaries, starting with 0 and ending with 0.5.")
    parser.add_argument("--sHetfile", help="sHet annotation table with CHROM POS [REF ALT] sHet.")
    parser.add_argument("--sHetbin", default=",".join(map(str, DEFAULT_SHETBIN)), help="Complete sHet bin boundaries, starting with 0 and ending with 1.")
    parser.add_argument("--region", help="Single user-defined region, e.g. Chr1:1000000-2000000.")
    parser.add_argument("--region-bed", help="BED file with user-defined regions.")
    parser.add_argument("--window-size", type=int, default=50000, help="Window size for wHet analysis. Default: 50000.")
    parser.add_argument("--window-step", type=int, default=50000, help="Window step for wHet analysis. Default: 50000.")
    parser.add_argument("--min-window-sites", type=int, default=100, help="Minimum comparable SNPs required in a wHet window. Default: 100.")
    parser.add_argument("--outdir", default="IPEva_result", help="Output directory. Default: IPEva_result.")
    parser.add_argument("--export-tables", action="store_true", help="Export TSV summary tables in addition to the HTML report.")
    parser.add_argument("--export-site-table", action="store_true", help="Export compressed site-level table. Implies --export-tables.")
    args = parser.parse_args()

    validate_gpth(args.GPth)
    gpscan = validate_gpscan(parse_float_list(args.GPscan, "--GPscan"))
    mafbin = validate_mafbin(parse_float_list(args.MAFbin, "--MAFbin"))
    shetbin = validate_shetbin(parse_float_list(args.sHetbin, "--sHetbin"))
    if args.window_size <= 0 or args.window_step <= 0 or args.min_window_sites <= 0:
        raise SystemExit("ERROR: --window-size, --window-step and --min-window-sites must be positive integers.")

    if args.mode == "all":
        missing = []
        if not args.MAFfile:
            missing.append("--MAFfile")
        if not args.sHetfile:
            missing.append("--sHetfile")
        if missing:
            raise SystemExit("ERROR: --mode all requires " + " and ".join(missing) + ".")

    site_requested = parse_site_analysis(args.site_analysis)
    interval_requested = parse_interval_analysis(args.interval_analysis)
    if args.mode == "all":
        site_requested = ["Genotype", "MAF", "sHet", "GP"]
        interval_requested = ["wHet"] + (["Region"] if (args.region or args.region_bed) else [])

    if args.mode == "site":
        if "MAF" in site_requested and not args.MAFfile:
            raise SystemExit("ERROR: --site-analysis MAF requires --MAFfile.")
        if "sHet" in site_requested and not args.sHetfile:
            raise SystemExit("ERROR: --site-analysis sHet requires --sHetfile.")
    if args.mode == "interval" and "Region" in interval_requested and not (args.region or args.region_bed):
        raise SystemExit("ERROR: --interval-analysis Region requires --region or --region-bed.")

    warnings: List[str] = []
    truth_sample, imputed_sample = ensure_single_sample_vcfs(args.truth, args.imputed, warnings)
    truth_records, truth_qc = read_truth_vcf(args.truth)
    sites, match_qc, excluded = read_imputed_vcf_and_match(args.imputed, truth_records)
    if not sites:
        raise SystemExit("ERROR: No comparable SNPs were found between truth and imputed VCFs.")

    maf_qc = assign_maf(sites, args.MAFfile, mafbin) if args.MAFfile else {}
    shet_qc = assign_shet(sites, args.sHetfile, shetbin) if args.sHetfile else {}

    results: Dict[str, Any] = {
        "generated_at": datetime.now().strftime("%Y-%m-%d %H:%M:%S"),
        "warnings": warnings,
        "metadata": {
            "Mode": args.mode,
            "Site_analysis": ",".join(site_requested) if args.mode in {"site", "all"} else NA,
            "Interval_analysis": ",".join(interval_requested) if args.mode in {"interval", "all"} else NA,
            "GPth": args.GPth,
            "GPscan": ",".join(f"{v:g}" for v in gpscan),
            "MAFbin": ",".join(f"{v:g}" for v in mafbin),
            "sHetbin": ",".join(f"{v:g}" for v in shetbin),
            "Truth_sample_name": truth_sample,
            "Imputed_sample_name": imputed_sample,
            "Truth_VCF": args.truth,
            "Imputed_VCF": args.imputed,
            "MAFfile": args.MAFfile or NA,
            "sHetfile": args.sHetfile or NA,
        },
        "matching_qc": {**{f"truth_{k}": v for k, v in truth_qc.items()}, **match_qc, **maf_qc, **shet_qc},
        "overall": None,
        "site": {},
        "interval": {},
    }

    if args.mode in {"overall", "all"}:
        results["overall"] = compute_metrics(sites, args.GPth, "Overall", "All")

    if args.mode in {"site", "all"}:
        site_results: Dict[str, List[Dict[str, Any]]] = {}
        if "Genotype" in site_requested:
            rows = []
            for group, sub in group_by_genotype(sites).items():
                rows.append(compute_metrics(sub, args.GPth, "Genotype", group))
            site_results["Genotype"] = sorted(rows, key=lambda r: 0 if r["Group"] == "Homozygous" else 1)
        if "MAF" in site_requested:
            groups = group_by_attr(sites, "maf_bin")
            rows = []
            for i, (lo, hi) in enumerate(zip(mafbin[:-1], mafbin[1:])):
                group = f"MAF_bin_{i + 1}"
                sub = groups.get(group, [])
                m = compute_metrics(sub, args.GPth, "MAF", group)
                m["MAF_range"] = f"({lo:g}, {hi:g}]"
                m["bin_index"] = i + 1
                rows.append(m)
            site_results["MAF"] = rows
        if "sHet" in site_requested:
            groups = group_by_attr(sites, "shet_bin")
            rows = []
            total = len(shetbin) - 1
            for i, (lo, hi) in enumerate(zip(shetbin[:-1], shetbin[1:])):
                group = bin_label("sHet", i, total)
                sub = groups.get(group, [])
                m = compute_metrics(sub, args.GPth, "sHet", group)
                m["sHet_range"] = f"[{lo:g}, {hi:g}]" if i == 0 else f"({lo:g}, {hi:g}]"
                m["bin_index"] = i + 1
                rows.append(m)
            site_results["sHet"] = rows
        if "GP" in site_requested:
            rows = [compute_metrics(sites, gp, "GP", f"GP>={gp:g}") for gp in gpscan]
            site_results["GP"] = rows
        results["site"] = site_results

    if args.mode in {"interval", "all"}:
        interval_results: Dict[str, List[Dict[str, Any]]] = {}
        if "Region" in interval_requested:
            regions = read_regions(args.region, args.region_bed)
            if regions:
                interval_results["Region"] = compute_region_metrics(sites, regions, args.GPth)
        if "wHet" in interval_requested:
            interval_results["wHet"] = compute_whet_windows(sites, args.window_size, args.window_step, args.min_window_sites, args.GPth)
        results["interval"] = interval_results

    outdir = Path(args.outdir)
    outdir.mkdir(parents=True, exist_ok=True)
    (outdir / "IPEva_parameters.json").write_text(json.dumps(results["metadata"], indent=2, ensure_ascii=False), encoding="utf-8")
    with open(outdir / "IPEva_run.log", "w") as log:
        log.write(f"IPEva version {VERSION}\n")
        for w in warnings:
            log.write(w + "\n")
        log.write(json.dumps(results["matching_qc"], indent=2) + "\n")
    build_html_report(results, outdir / "IPEva_report.html")
    if args.export_site_table:
        args.export_tables = True
    if args.export_tables:
        export_tables(outdir, results, sites, excluded, args.export_site_table)
    print(f"IPEva report generated: {outdir / 'IPEva_report.html'}")


if __name__ == "__main__":
    main()
