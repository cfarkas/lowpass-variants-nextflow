#!/usr/bin/env python3
import argparse
import csv
import gzip
import os
from collections import defaultdict, Counter

TAB = chr(9)
NL = chr(10)
CR = chr(13)
CALLER_ORDER = ["MUTECT2", "FREEBAYES", "BCFTOOLS"]
CALLER_PRIORITY = {"MUTECT2": 3, "FREEBAYES": 2, "BCFTOOLS": 1}


def strip_line(line):
    return line.rstrip(NL).rstrip(CR)


def open_text(path):
    return gzip.open(path, "rt") if str(path).endswith(".gz") else open(path, "rt")


def norm_chrom(chrom):
    c = str(chrom).strip()
    if c.startswith("chr"):
        c = c[3:]
    if c == "M":
        c = "MT"
    return c


def chrom_sort_key(chrom):
    c = norm_chrom(chrom)
    order = {str(i): i for i in range(1, 23)}
    order.update({"X": 23, "Y": 24, "MT": 25})
    return (order.get(c, 999), c)


def parse_info(s):
    d = {}
    if s in ("", "."):
        return d
    for item in str(s).split(";"):
        if not item:
            continue
        if "=" in item:
            k, v = item.split("=", 1)
            d[k] = v
        else:
            d[item] = True
    return d


def parse_format(fmt, vals):
    if fmt in ("", ".") or vals in ("", "."):
        return {}
    keys = str(fmt).split(":")
    values = str(vals).split(":")
    return dict(zip(keys, values))


def numeric_list(v, cast=float):
    if v is None or v in ("", "."):
        return []
    out = []
    for x in str(v).split(","):
        if x in ("", "."):
            continue
        try:
            out.append(cast(float(x)))
        except Exception:
            continue
    return out


def safe_int(v):
    try:
        if v is None or v in ("", "."):
            return None
        return int(float(v))
    except Exception:
        return None


def safe_float(v):
    try:
        if v is None or v in ("", "."):
            return None
        return float(v)
    except Exception:
        return None


def depth_from_fields(fm, info):
    dp = safe_int(fm.get("DP"))
    if dp is not None:
        return dp
    dp = safe_int(info.get("DP"))
    if dp is not None:
        return dp
    ad = numeric_list(fm.get("AD"), int)
    if ad:
        return sum(ad)
    ro = safe_int(fm.get("RO"))
    ao = numeric_list(fm.get("AO"), int)
    if ro is not None and ao:
        return ro + sum(ao)
    return None


def alt_count_from_fields(fm, info, alt_i):
    ad = numeric_list(fm.get("AD"), int)
    if ad:
        idx = alt_i + 1
        if idx < len(ad):
            return ad[idx]
        if len(ad) >= 2:
            return max(ad[1:])
    ao = numeric_list(fm.get("AO"), int)
    if ao:
        if alt_i < len(ao):
            return ao[alt_i]
        return max(ao)
    dv = safe_int(fm.get("DV"))
    if dv is not None:
        return dv
    ao_info = numeric_list(info.get("AO"), int)
    if ao_info:
        if alt_i < len(ao_info):
            return ao_info[alt_i]
        return max(ao_info)
    return None


def af_from_fields(fm, info, alt_reads, dp, alt_i):
    for key in ("AF", "VAF"):
        vals = numeric_list(fm.get(key), float)
        if vals:
            return vals[alt_i] if alt_i < len(vals) else max(vals)
    for key in ("AF", "VAF"):
        vals = numeric_list(info.get(key), float)
        if vals:
            return vals[alt_i] if alt_i < len(vals) else max(vals)
    if alt_reads is not None and dp and dp > 0:
        return float(alt_reads) / float(dp)
    return None


def calculated_af(dp, alt_reads):
    if dp is None or alt_reads is None:
        return None
    try:
        dp = float(dp)
        alt_reads = float(alt_reads)
    except Exception:
        return None
    if dp <= 0:
        return None
    return max(0.0, min(1.0, alt_reads / dp))


def is_ct_or_ga(ref, alt):
    ref = str(ref).upper()
    alt = str(alt).upper()
    return (ref == "C" and alt == "T") or (ref == "G" and alt == "A")


def is_simple_variant(ref, alt):
    if not ref or not alt or ref == "." or alt == "." or alt == "*":
        return False
    if alt.startswith("<") or ref.startswith("<"):
        return False
    valid = set("ACGTN")
    return set(ref.upper()) <= valid and set(alt.upper()) <= valid


def read_ref_contigs(fai):
    contigs = []
    if fai and os.path.exists(fai):
        with open(fai, "rt") as fh:
            for line in fh:
                f = line.rstrip("\n").split("\t")
                if len(f) >= 2:
                    try:
                        contigs.append((f[0], int(f[1])))
                    except Exception:
                        pass
    return contigs


def read_targets(path):
    targets = defaultdict(list)
    if not path or path == "." or not os.path.exists(path) or os.path.getsize(path) == 0:
        return targets
    with open(path, "rt") as fh:
        for line in fh:
            if not line.strip() or line.startswith("#"):
                continue
            f = strip_line(line).split(TAB)
            if len(f) < 3:
                continue
            try:
                start = int(f[1])
                end = int(f[2])
            except Exception:
                continue
            chrom = f[0]
            region_id = f[3] if len(f) > 3 and f[3] else f"{chrom}:{start}-{end}"
            region_type = f[4] if len(f) > 4 and f[4] else "."
            gene = f[5] if len(f) > 5 and f[5] else "."
            targets[norm_chrom(chrom)].append((start, end, gene, region_type, region_id))
    for chrom in targets:
        targets[chrom].sort(key=lambda x: (x[0], x[1], x[4]))
    return targets


def join_unique(values, default="."):
    out = []
    seen = set()
    for v in values:
        v = str(v)
        if not v or v == ".":
            continue
        if v not in seen:
            seen.add(v)
            out.append(v)
    return ",".join(out) if out else default


def target_annotation(targets, chrom, pos1):
    if not targets:
        return {"gene": ".", "region_type": "genome_wide", "region_id": "."}
    c = norm_chrom(chrom)
    pos0 = int(pos1) - 1
    genes = []
    types = []
    ids = []
    for start, end, gene, region_type, region_id in targets.get(c, []):
        if end <= pos0:
            continue
        if start > pos0:
            break
        genes.append(gene)
        types.append(region_type)
        ids.append(region_id)
    return {
        "gene": join_unique(genes),
        "region_type": join_unique(types),
        "region_id": join_unique(ids),
    }


def canonical_caller(caller):
    c = str(caller).lower()
    if c.startswith("mutect2"):
        return "MUTECT2"
    if c == "freebayes":
        return "FREEBAYES"
    if c == "bcftools":
        return "BCFTOOLS"
    return str(caller).upper()


def parse_vcf(vcf, sample, caller, status):
    calls = []
    if not vcf or vcf == "." or not os.path.exists(vcf) or os.path.getsize(vcf) == 0:
        return calls
    canonical = canonical_caller(caller)
    with open_text(vcf) as fh:
        for line in fh:
            if not line.strip() or line.startswith("#"):
                continue
            f = strip_line(line).split(TAB)
            if len(f) < 8:
                continue
            chrom, pos, _vid, ref, alts, qual, filt, info_s = f[:8]
            try:
                pos_i = int(pos)
            except Exception:
                continue
            info = parse_info(info_s)
            fmt = f[8] if len(f) > 8 else "."
            sample_vals = f[9] if len(f) > 9 else "."
            fm = parse_format(fmt, sample_vals)
            dp = depth_from_fields(fm, info)
            for alt_i, alt in enumerate(str(alts).split(",")):
                alt = alt.strip()
                ref_u = str(ref).upper()
                alt_u = str(alt).upper()
                if not is_simple_variant(ref_u, alt_u):
                    continue
                alt_reads = alt_count_from_fields(fm, info, alt_i)
                af_observed = af_from_fields(fm, info, alt_reads, dp, alt_i)
                af_calc = calculated_af(dp, alt_reads)
                if caller == "mutect2_raw_unfiltered" or str(status).startswith("raw_unfiltered"):
                    caller_pass = False
                else:
                    caller_pass = filt in ("PASS", ".")
                calls.append({
                    "sample": sample,
                    "caller": canonical,
                    "caller_raw": caller,
                    "status": status,
                    "chrom": chrom,
                    "chrom_norm": norm_chrom(chrom),
                    "pos": pos_i,
                    "ref": ref_u,
                    "alt": alt_u,
                    "filter": filt,
                    "qual": qual,
                    "dp": dp,
                    "alt_reads": alt_reads,
                    "af_observed": af_observed,
                    "af_calculated": af_calc,
                    "caller_pass": caller_pass,
                })
    return calls


def key_of(call):
    return (call["sample"], call["chrom_norm"], int(call["pos"]), call["ref"], call["alt"])


def call_af(call):
    if call is None:
        return None
    if call.get("af_calculated") is not None:
        return call.get("af_calculated")
    return call.get("af_observed")


def evidence_ok(call, args):
    if call is None or not call.get("caller_pass"):
        return False
    dp = call.get("dp") if call.get("dp") is not None else -1
    alt = call.get("alt_reads") if call.get("alt_reads") is not None else -1
    af = call_af(call)
    af = af if af is not None else -1.0
    return dp >= args.min_dp and alt >= args.min_alt_reads and af >= args.min_af


def best_call(calls, args=None):
    best = None
    best_score = None
    for call in calls:
        af = call_af(call)
        af = af if af is not None else -1.0
        alt = call.get("alt_reads") if call.get("alt_reads") is not None else -1
        dp = call.get("dp") if call.get("dp") is not None else -1
        score = (
            1 if (args is not None and evidence_ok(call, args)) else 0,
            1 if call.get("caller_pass") else 0,
            CALLER_PRIORITY.get(call.get("caller", ""), 0),
            af,
            alt,
            dp,
        )
        if best is None or score > best_score:
            best = call
            best_score = score
    return best


def best_by_caller(calls, args):
    grouped = defaultdict(list)
    for call in calls:
        grouped[call["caller"]].append(call)
    return {caller: best_call(vals, args) for caller, vals in grouped.items()}


def metric_max(calls, field, default):
    vals = [c.get(field) for c in calls if c is not None and c.get(field) is not None]
    return max(vals) if vals else default


def af_max(calls):
    vals = [call_af(c) for c in calls if c is not None and call_af(c) is not None]
    return max(vals) if vals else None


def fmt(value):
    if value is None:
        return "."
    if isinstance(value, float):
        return f"{value:.6f}"
    return str(value)


def info_escape(value):
    s = str(value)
    if s in ("", "None"):
        return "."
    for ch in [";", "=", "\t", "\n", "\r", " "]:
        s = s.replace(ch, "_")
    return s if s else "."


def vcf_sample_format(dp, alt_reads, af):
    dp_i = None
    alt_i = None
    try:
        if dp is not None and str(dp) != ".":
            dp_i = max(0, int(round(float(dp))))
    except Exception:
        dp_i = None
    try:
        if alt_reads is not None and str(alt_reads) != ".":
            alt_i = max(0, int(round(float(alt_reads))))
    except Exception:
        alt_i = None
    if dp_i is not None and alt_i is not None:
        ref_i = max(0, dp_i - alt_i)
        ad = f"{ref_i},{alt_i}"
        if af is None:
            af = calculated_af(dp_i, alt_i)
    else:
        ad = "."
    if af is None:
        af_s = "."
        gt = "0/1"
    else:
        af_s = f"{float(af):.6f}"
        gt = "1/1" if float(af) >= 0.80 else "0/1"
    dp_s = "." if dp_i is None else str(dp_i)
    return f"{gt}:{dp_s}:{ad}:{af_s}"


def decide_variant(ref, alt, calls_by_caller, args):
    calls = list(calls_by_caller.values())
    source_callers = [c for c in CALLER_ORDER if c in calls_by_caller]
    supporting_callers = [c for c in CALLER_ORDER if c in calls_by_caller and evidence_ok(calls_by_caller[c], args)]
    pass_callers = [c for c in CALLER_ORDER if c in calls_by_caller and calls_by_caller[c].get("caller_pass")]
    vote_count = len(supporting_callers)
    max_dp = metric_max(calls, "dp", -1)
    max_alt = metric_max(calls, "alt_reads", -1)
    max_af = af_max(calls)
    max_af_num = max_af if max_af is not None else -1.0
    ffpe_mode = args.sample_mode == "ffpe"
    ctga = is_ct_or_ga(ref, alt) if ffpe_mode else False
    obvious = []
    review = []
    if not pass_callers:
        obvious.append("no_PASS_caller")
    if vote_count < 1:
        obvious.append("no_caller_with_minimum_evidence")
    if max_dp < args.min_dp:
        obvious.append(f"DP<{args.min_dp}")
    if max_alt < args.min_alt_reads:
        obvious.append(f"ALT<{args.min_alt_reads}")
    if max_af_num < args.min_af:
        obvious.append(f"AF<{args.min_af}")
    if vote_count < args.vote_threshold:
        review.append(f"supporting_callers<{args.vote_threshold}")
    ffpe_low_support = False
    if ffpe_mode and ctga and len(supporting_callers) < 3:
        if max_af_num < args.ffpe_ct_vaf_keep or max_alt < args.ffpe_ct_min_alt:
            ffpe_low_support = True
            review.append(f"FFPE_CtoT_GtoA_low_support_AF<{args.ffpe_ct_vaf_keep}_or_ALT<{args.ffpe_ct_min_alt}")
    if obvious:
        decision = "LOW_CONFIDENCE_OR_ARTIFACT"
        reasons = obvious + review
    elif review:
        decision = "MODERATE_REVIEW"
        reasons = review
    else:
        decision = "HIGH_CONFIDENCE"
        reasons = []
    decision_row = {
        "decision": decision,
        "source": ",".join(source_callers) if source_callers else ".",
        "supporting_callers": ",".join(supporting_callers) if supporting_callers else ".",
        "pass_callers": ",".join(pass_callers) if pass_callers else ".",
        "caller_count": len(source_callers),
        "vote_count": vote_count,
        "max_dp": max_dp,
        "max_alt_reads": max_alt,
        "max_af": max_af,
        "pipeline_mode": args.sample_mode,
        "drop_reason": ";".join(reasons) if reasons else ".",
    }
    if ffpe_mode:
        decision_row.update(
            is_CT_or_GA="yes" if ctga else "no",
            ffpe_risk="FFPE_CtoT_GtoA_risk" if ctga else "none",
            ffpe_low_support="yes" if ffpe_low_support else "no",
        )
    return decision_row


def caller_fields(calls_by_caller):
    out = {}
    for caller in CALLER_ORDER:
        prefix = caller.lower()
        call = calls_by_caller.get(caller)
        if call is None:
            out[f"{prefix}_filter"] = "."
            out[f"{prefix}_dp"] = "."
            out[f"{prefix}_alt_reads"] = "."
            out[f"{prefix}_af_observed"] = "."
            out[f"{prefix}_af_calculated"] = "."
            out[f"{prefix}_status"] = "."
        else:
            out[f"{prefix}_filter"] = call.get("filter", ".")
            out[f"{prefix}_dp"] = fmt(call.get("dp"))
            out[f"{prefix}_alt_reads"] = fmt(call.get("alt_reads"))
            out[f"{prefix}_af_observed"] = fmt(call.get("af_observed"))
            out[f"{prefix}_af_calculated"] = fmt(call.get("af_calculated"))
            out[f"{prefix}_status"] = call.get("status", ".")
    return out


def write_tsv(path, fields, rows):
    with open(path, "w", newline="") as out:
        writer = csv.DictWriter(out, delimiter=TAB, fieldnames=fields, lineterminator=NL, extrasaction="ignore")
        writer.writeheader()
        for row in rows:
            writer.writerow({f: row.get(f, ".") for f in fields})


def write_vcf(path, rows, sample, contigs, include_artifacts, sample_mode):
    rows_to_write = []
    seen = set()
    for row in sorted(rows, key=lambda r: (chrom_sort_key(r["chrom"]), int(r["pos"]), r["ref"], r["alt"])):
        if row["decision"] == "LOW_CONFIDENCE_OR_ARTIFACT" and not include_artifacts:
            continue
        key = (norm_chrom(row["chrom"]), int(row["pos"]), row["ref"], row["alt"])
        if key in seen:
            continue
        seen.add(key)
        rows_to_write.append(row)
    with open(path, "w") as out:
        out.write("##fileformat=VCFv4.2" + NL)
        out.write("##source=lowpass-variants-nextflow" + NL)
        out.write(f"##lowpass_variants_mode={sample_mode}" + NL)
        out.write('##FILTER=<ID=PASS,Description="High-confidence final variant for this low-pass workflow">' + NL)
        out.write('##FILTER=<ID=REVIEW,Description="Review-level final variant retained for this low-pass workflow">' + NL)
        out.write('##FILTER=<ID=ARTIFACT,Description="Low-confidence or obvious artifact retained only because --include-obvious-artifacts was used">' + NL)
        for chrom, length in contigs:
            out.write(f"##contig=<ID={chrom},length={length}>" + NL)
        out.write('##INFO=<ID=SOURCE,Number=.,Type=String,Description="Uppercase callers that found the variant: MUTECT2, FREEBAYES, BCFTOOLS">' + NL)
        out.write('##INFO=<ID=SUPPORTING_CALLERS,Number=.,Type=String,Description="Callers passing caller filter and minimum DP/ALT/AF evidence">' + NL)
        out.write('##INFO=<ID=PASS_CALLERS,Number=.,Type=String,Description="Callers with PASS or missing FILTER before evidence thresholds">' + NL)
        out.write('##INFO=<ID=DECISION,Number=1,Type=String,Description="Pipeline decision: HIGH_CONFIDENCE, MODERATE_REVIEW, or LOW_CONFIDENCE_OR_ARTIFACT">' + NL)
        out.write('##INFO=<ID=CALLER_COUNT,Number=1,Type=Integer,Description="Number of callers that found this normalized variant">' + NL)
        out.write('##INFO=<ID=VOTE_COUNT,Number=1,Type=Integer,Description="Number of callers supporting this variant after caller filter and minimum evidence thresholds">' + NL)
        out.write('##INFO=<ID=PIPELINE_MODE,Number=1,Type=String,Description="Selected sample mode: ffpe or fresh">' + NL)
        out.write('##INFO=<ID=REP_CALLER,Number=1,Type=String,Description="Caller selected as representative for FORMAT DP/AD/AF">' + NL)
        out.write('##INFO=<ID=AF_CALCULATED,Number=1,Type=Float,Description="Allele fraction calculated from representative ALT_READS/DP">' + NL)
        out.write('##INFO=<ID=MAX_AF,Number=1,Type=Float,Description="Maximum calculated/observed AF across callers">' + NL)
        out.write('##INFO=<ID=MAX_DP,Number=1,Type=Integer,Description="Maximum depth across callers">' + NL)
        out.write('##INFO=<ID=MAX_ALT_READS,Number=1,Type=Integer,Description="Maximum alternate read support across callers">' + NL)
        if sample_mode == "ffpe":
            out.write('##INFO=<ID=FFPE_RISK,Number=1,Type=String,Description="C>T/G>A FFPE-risk status">' + NL)
            out.write('##INFO=<ID=FFPE_LOW_SUPPORT,Number=1,Type=String,Description="FFPE C>T/G>A support assessment">' + NL)
        out.write('##INFO=<ID=DROP_REASON,Number=.,Type=String,Description="Reasons for review or artifact classification">' + NL)
        out.write('##INFO=<ID=GENE,Number=.,Type=String,Description="Gene annotation when --genes/--bed annotation is available">' + NL)
        out.write('##INFO=<ID=REGION_TYPE,Number=.,Type=String,Description="Target region type when available">' + NL)
        out.write('##INFO=<ID=REGION_ID,Number=.,Type=String,Description="Target region identifier when available">' + NL)
        for caller in CALLER_ORDER:
            p = caller
            out.write(f'##INFO=<ID={p}_FILTER,Number=1,Type=String,Description="Original FILTER value for {caller}">' + NL)
            out.write(f'##INFO=<ID={p}_DP,Number=1,Type=String,Description="Depth parsed for {caller}">' + NL)
            out.write(f'##INFO=<ID={p}_ALT_READS,Number=1,Type=String,Description="Alternate reads parsed for {caller}">' + NL)
            out.write(f'##INFO=<ID={p}_AF,Number=1,Type=String,Description="Calculated AF for {caller}, or observed AF if calculation was unavailable">' + NL)
        out.write('##FORMAT=<ID=GT,Number=1,Type=String,Description="Approximate genotype derived for downstream compatibility; not a clinical germline genotype">' + NL)
        out.write('##FORMAT=<ID=DP,Number=1,Type=Integer,Description="Representative depth used for calculated AF">' + NL)
        out.write('##FORMAT=<ID=AD,Number=R,Type=Integer,Description="Representative REF and ALT read counts estimated as DP-ALT,ALT">' + NL)
        out.write('##FORMAT=<ID=AF,Number=A,Type=Float,Description="Calculated allele fraction as ALT_READS/DP from representative evidence">' + NL)
        out.write(TAB.join(["#CHROM", "POS", "ID", "REF", "ALT", "QUAL", "FILTER", "INFO", "FORMAT", sample]) + NL)
        for row in rows_to_write:
            decision = row["decision"]
            filt = "PASS" if decision == "HIGH_CONFIDENCE" else ("REVIEW" if decision == "MODERATE_REVIEW" else "ARTIFACT")
            info = [
                f"SOURCE={info_escape(row.get('source', '.'))}",
                f"SUPPORTING_CALLERS={info_escape(row.get('supporting_callers', '.'))}",
                f"PASS_CALLERS={info_escape(row.get('pass_callers', '.'))}",
                f"DECISION={info_escape(decision)}",
                f"CALLER_COUNT={info_escape(row.get('caller_count', '.'))}",
                f"VOTE_COUNT={info_escape(row.get('vote_count', '.'))}",
                f"PIPELINE_MODE={info_escape(row.get('pipeline_mode', sample_mode))}",
                f"REP_CALLER={info_escape(row.get('representative_caller', '.'))}",
                f"AF_CALCULATED={info_escape(row.get('representative_af_calculated', '.'))}",
                f"MAX_AF={info_escape(row.get('max_af', '.'))}",
                f"MAX_DP={info_escape(row.get('max_dp', '.'))}",
                f"MAX_ALT_READS={info_escape(row.get('max_alt_reads', '.'))}",
                f"DROP_REASON={info_escape(row.get('drop_reason', '.'))}",
                f"GENE={info_escape(row.get('gene', '.'))}",
                f"REGION_TYPE={info_escape(row.get('region_type', '.'))}",
                f"REGION_ID={info_escape(row.get('region_id', '.'))}",
            ]
            if sample_mode == "ffpe":
                info.extend(
                    [
                        f"FFPE_RISK={info_escape(row.get('ffpe_risk', '.'))}",
                        f"FFPE_LOW_SUPPORT={info_escape(row.get('ffpe_low_support', '.'))}",
                    ]
                )
            for caller in CALLER_ORDER:
                p = caller.lower()
                info.append(f"{caller}_FILTER={info_escape(row.get(p + '_filter', '.'))}")
                info.append(f"{caller}_DP={info_escape(row.get(p + '_dp', '.'))}")
                info.append(f"{caller}_ALT_READS={info_escape(row.get(p + '_alt_reads', '.'))}")
                # prefer calculated AF for per-caller AF in VCF INFO
                caller_af = row.get(p + "_af_calculated", ".")
                if caller_af == ".":
                    caller_af = row.get(p + "_af_observed", ".")
                info.append(f"{caller}_AF={info_escape(caller_af)}")
            sample_fmt = vcf_sample_format(
                row.get("representative_dp"),
                row.get("representative_alt_reads"),
                safe_float(row.get("representative_af_calculated")),
            )
            out.write(TAB.join([
                str(row["chrom"]),
                str(row["pos"]),
                ".",
                str(row["ref"]),
                str(row["alt"]),
                ".",
                filt,
                ";".join(info),
                "GT:DP:AD:AF",
                sample_fmt,
            ]) + NL)


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--sample", required=True)
    ap.add_argument("--sample-mode", choices=("ffpe", "fresh"), required=True)
    ap.add_argument("--mutect-vcf", default="")
    ap.add_argument("--mutect-status", default="not_run")
    ap.add_argument("--freebayes-vcf", default="")
    ap.add_argument("--freebayes-status", default="completed")
    ap.add_argument("--bcftools-vcf", default="")
    ap.add_argument("--bcftools-status", default="completed")
    ap.add_argument("--targets", default="")
    ap.add_argument("--ref-fai", required=True)
    ap.add_argument("--outdir", required=True)
    ap.add_argument("--out-vcf", required=True)
    ap.add_argument("--min-dp", type=int, required=True)
    ap.add_argument("--min-alt-reads", type=int, required=True)
    ap.add_argument("--min-af", type=float, required=True)
    ap.add_argument("--vote-threshold", type=int, required=True)
    ap.add_argument("--ffpe-ct-vaf-keep", type=float)
    ap.add_argument("--ffpe-ct-min-alt", type=int)
    ap.add_argument("--include-obvious-artifacts", action="store_true")
    args = ap.parse_args()

    if args.sample_mode == "ffpe":
        if args.ffpe_ct_vaf_keep is None or args.ffpe_ct_min_alt is None:
            ap.error("FFPE mode requires --ffpe-ct-vaf-keep and --ffpe-ct-min-alt")

    os.makedirs(args.outdir, exist_ok=True)
    out_vcf_parent = os.path.dirname(os.path.abspath(args.out_vcf))
    if out_vcf_parent:
        os.makedirs(out_vcf_parent, exist_ok=True)
    targets = read_targets(args.targets)
    contigs = read_ref_contigs(args.ref_fai)

    all_calls = []
    all_calls.extend(parse_vcf(args.mutect_vcf, args.sample, "mutect2", args.mutect_status) if args.mutect_vcf else [])
    all_calls.extend(parse_vcf(args.freebayes_vcf, args.sample, "freebayes", args.freebayes_status) if args.freebayes_vcf else [])
    all_calls.extend(parse_vcf(args.bcftools_vcf, args.sample, "bcftools", args.bcftools_status) if args.bcftools_vcf else [])

    grouped = defaultdict(list)
    for call in all_calls:
        grouped[key_of(call)].append(call)

    rows = []
    for key, calls in sorted(grouped.items(), key=lambda x: (chrom_sort_key(x[0][1]), x[0][2], x[0][3], x[0][4])):
        _sample, _chrom_norm, pos, ref, alt = key
        calls_by_caller = best_by_caller(calls, args)
        source_call = best_call(list(calls_by_caller.values()), args)
        if source_call is None:
            continue
        decision = decide_variant(ref, alt, calls_by_caller, args)
        rep = best_call(list(calls_by_caller.values()), args)
        rep_af_calc = calculated_af(rep.get("dp"), rep.get("alt_reads")) if rep else None
        if rep_af_calc is None and rep:
            rep_af_calc = rep.get("af_calculated") if rep.get("af_calculated") is not None else rep.get("af_observed")
        annot = target_annotation(targets, source_call["chrom"], pos)
        row = {
            "sample": args.sample,
            "chrom": source_call["chrom"],
            "pos": pos,
            "ref": ref,
            "alt": alt,
            "gene": annot["gene"],
            "region_type": annot["region_type"],
            "region_id": annot["region_id"],
            "decision": decision["decision"],
            "source": decision["source"],
            "supporting_callers": decision["supporting_callers"],
            "pass_callers": decision["pass_callers"],
            "caller_count": decision["caller_count"],
            "vote_count": decision["vote_count"],
            "max_dp": fmt(decision["max_dp"] if decision["max_dp"] >= 0 else None),
            "max_alt_reads": fmt(decision["max_alt_reads"] if decision["max_alt_reads"] >= 0 else None),
            "max_af": fmt(decision["max_af"]),
            "pipeline_mode": decision["pipeline_mode"],
            "drop_reason": decision["drop_reason"],
            "representative_caller": rep.get("caller", ".") if rep else ".",
            "representative_filter": rep.get("filter", ".") if rep else ".",
            "representative_dp": fmt(rep.get("dp") if rep else None),
            "representative_alt_reads": fmt(rep.get("alt_reads") if rep else None),
            "representative_af_observed": fmt(rep.get("af_observed") if rep else None),
            "representative_af_calculated": fmt(rep_af_calc),
        }
        if args.sample_mode == "ffpe":
            row.update(
                is_CT_or_GA=decision["is_CT_or_GA"],
                ffpe_risk=decision["ffpe_risk"],
                ffpe_low_support=decision["ffpe_low_support"],
            )
        row.update(caller_fields(calls_by_caller))
        rows.append(row)

    fields = [
        "sample", "chrom", "pos", "ref", "alt", "gene", "region_type", "region_id",
        "decision", "source", "supporting_callers", "pass_callers", "caller_count", "vote_count",
        "max_dp", "max_alt_reads", "max_af", "pipeline_mode", "drop_reason",
        "representative_caller", "representative_filter", "representative_dp", "representative_alt_reads",
        "representative_af_observed", "representative_af_calculated",
        "mutect2_filter", "mutect2_dp", "mutect2_alt_reads", "mutect2_af_observed", "mutect2_af_calculated", "mutect2_status",
        "freebayes_filter", "freebayes_dp", "freebayes_alt_reads", "freebayes_af_observed", "freebayes_af_calculated", "freebayes_status",
        "bcftools_filter", "bcftools_dp", "bcftools_alt_reads", "bcftools_af_observed", "bcftools_af_calculated", "bcftools_status",
    ]
    if args.sample_mode == "ffpe":
        insert_at = fields.index("drop_reason")
        fields[insert_at:insert_at] = ["is_CT_or_GA", "ffpe_risk", "ffpe_low_support"]
    prefix = os.path.join(args.outdir, args.sample)
    retained = [r for r in rows if args.include_obvious_artifacts or r["decision"] != "LOW_CONFIDENCE_OR_ARTIFACT"]
    write_tsv(prefix + ".final_variants.all_candidates.tsv", fields, rows)
    write_tsv(prefix + ".final_variants.retained.tsv", fields, retained)

    counts = Counter(r["decision"] for r in rows)
    with open(prefix + ".final_variants.counts.tsv", "w", newline="") as out:
        w = csv.writer(out, delimiter=TAB, lineterminator=NL)
        w.writerow(["sample", "metric", "count"])
        for key in ["HIGH_CONFIDENCE", "MODERATE_REVIEW", "LOW_CONFIDENCE_OR_ARTIFACT"]:
            w.writerow([args.sample, key, counts.get(key, 0)])
        w.writerow([args.sample, "total_candidates", len(rows)])
        w.writerow([args.sample, "total_retained_in_final_vcf", len(retained)])

    write_vcf(args.out_vcf, rows, args.sample, contigs, args.include_obvious_artifacts, args.sample_mode)


if __name__ == "__main__":
    main()
