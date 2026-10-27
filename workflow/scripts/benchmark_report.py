#!/usr/bin/env python3
"""
Build an HTML report of runtime, CPU and memory usage from the Snakemake
benchmark files written by Pipeasm (benchmarks/{sample}.<rule>.txt).

Only the jobs that actually ran leave a benchmark file, and which jobs run
depends on the inputs (HiFi, Hi-C, ONT) and options (solo_asm, run_mitohifi,
gxdb, ...). The report is therefore built from whatever is on disk: each file
is matched to its rule through the `benchmark:` directives in
workflow/rules/*.smk, and each rule is assigned to a pipeline stage by the .smk
file that defines it. Nothing is hard-coded about which rules ran.

Usage (from the Pipeasm working directory):
    python workflow/scripts/benchmark_report.py -c config/config.yaml
"""

import argparse
import csv
import datetime
import html
import logging
import os
import re
import statistics
import sys
from collections import defaultdict

import yaml


# ── Constants ──────────────────────────────────────────────────────────────────

METRICS = [
    "s",
    "max_rss",
    "max_vms",
    "max_uss",
    "max_pss",
    "io_in",
    "io_out",
    "mean_load",
    "cpu_time",
]

# Pipeline stages in execution order; the order also fixes each stage's color,
# so a stage keeps its hue whichever stages were run.
STAGES = [
    ("trimming_qc", "Trimming and QC"),
    ("kmer", "K-mer profiling"),
    ("assembly", "Assembly"),
    ("assembly_stats", "Assembly statistics"),
    ("decontamination", "Decontamination"),
    ("hic_scaffolding", "Hi-C mapping and scaffolding"),
    ("scaffolding_stats", "Scaffolding statistics"),
]
STAGE_ORDER = [key for key, _ in STAGES]
STAGE_LABELS = dict(STAGES, other="Other")

# Rule file (workflow/rules/<name>.smk) -> stage.
SMK_TO_STAGE = {
    "hifi_trimming_qc": "trimming_qc",
    "hic_trimming_qc": "trimming_qc",
    "ont_trimming_qc": "trimming_qc",
    "hifi_kmer_evaluation": "kmer",
    "phased_asm": "assembly",
    "solo_asm": "assembly",
    "ont_phased_asm": "assembly",
    "ont_solo_asm": "assembly",
    "phased_mitohifi": "assembly",
    "solo_mitohifi": "assembly",
    "download_mito_reference": "assembly",
    "phased_stats": "assembly_stats",
    "solo_stats": "assembly_stats",
    "phased_completeness": "assembly_stats",
    "solo_completeness": "assembly_stats",
    "phased_merqury": "assembly_stats",
    "solo_merqury": "assembly_stats",
    "phased_asm_plots": "assembly_stats",
    "solo_asm_plots": "assembly_stats",
    "phased_decontamination": "decontamination",
    "solo_decontamination": "decontamination",
    "hic_initital_contacts": "hic_scaffolding",
    "hic_scaffolding": "hic_scaffolding",
    "hic_final_contacts": "hic_scaffolding",
    "scaffolding_stats": "scaffolding_stats",
    "scaffolding_completeness": "scaffolding_stats",
    "scaffolding_merqury": "scaffolding_stats",
    "scaffolding_plots": "scaffolding_stats",
}

# Used for rules in a new .smk file or for benchmark files that match no rule:
# the stage is guessed from the rule name (first match wins).
RULE_TO_STAGE = [
    (r"^(trimming|qc)_", "trimming_qc"),
    (r"^(fastk|genescope|smudgeplot|kat)", "kmer"),
    (r"fcs", "decontamination"),
    (r"(assembly|mitohifi|mito_reference)", "assembly"),
    (r"^(initial_contacts|final_contacts|yahs)_", "hic_scaffolding"),
    (r"^scaffolding_", "scaffolding_stats"),
    (r"^(solo|phased)_", "assembly_stats"),
]

# Rule name tokens -> assembly the job works on.
ASSEMBLY_TOKENS = {
    "hap1": "Hap1",
    "hap2": "Hap2",
    "prim": "Primary",
    "primary": "Primary",
    "alt": "Alternate",
}

# Matches the `sample` wildcard constraint in workflow/rules/common.smk.
WILDCARD_REGEX = {"sample": r"[A-Za-z0-9_]+"}

# Run settings shown in the report (tool parameters are left out).
CONFIG_FIELDS = [
    "sample",
    "genome_size",
    "hifi_reads",
    "hic_r1",
    "hic_r2",
    "ont_reads",
    "solo_asm",
    "diff_species_hic",
    "run_mitohifi",
    "buscodb",
    "taxid",
    "gxdb",
]
PATH_FIELDS = {"hifi_reads", "hic_r1", "hic_r2", "ont_reads", "gxdb"}

# A job using at least this share of its requested memory is flagged.
MEMORY_WARN = 0.9


# ── Parsing ────────────────────────────────────────────────────────────────────


def parse_arguments():
    parser = argparse.ArgumentParser(
        description="Generate an HTML report from Pipeasm benchmark files"
    )
    parser.add_argument(
        "-b", "--benchmarks", default="benchmarks", help="Benchmarks directory"
    )
    parser.add_argument(
        "-c",
        "--config",
        default=None,
        help="Pipeasm config file. Defaults to config/config.yaml when present",
    )
    parser.add_argument(
        "-r", "--rules", default="workflow/rules", help="Directory with the .smk rule files"
    )
    parser.add_argument(
        "-o", "--output", default="workflow/report/benchmark.html", help="Output HTML"
    )
    parser.add_argument(
        "--threads", default="config/threads.yaml", help="Threads per rule (config/threads.yaml)"
    )
    parser.add_argument(
        "--slurm-params",
        default="config/slurm_params.yaml",
        help="SLURM profile template with the default memory (config/slurm_params.yaml)",
    )
    parser.add_argument(
        "--resources",
        default="config/resources",
        help="Directory with the SLURM memory per rule for each genome_size",
    )
    parser.add_argument(
        "--all-samples",
        action="store_true",
        help="Include samples found in the benchmarks directory other than the config sample",
    )
    parser.add_argument(
        "--top", type=int, default=15, help="Number of rules in the peak memory chart"
    )
    return parser.parse_args()


def pattern_to_regex(pattern, bench_dirname):
    """Convert a Snakemake benchmark path pattern into a regex over paths
    relative to the benchmarks directory. Returns (None, 0) for malformed patterns."""
    for prefix in (bench_dirname.rstrip("/") + "/", "benchmarks/"):
        if pattern.startswith(prefix):
            pattern = pattern[len(prefix):]
            break

    regex = []
    seen = set()
    literal_len = 0
    for part in re.split(r"(\{[^{}]*\})", pattern):
        if part.startswith("{") and part.endswith("}"):
            name = part[1:-1].split(",")[0].strip()
            if not name.isidentifier():
                return None, 0
            if name in seen:
                regex.append(f"(?P={name})")
            else:
                seen.add(name)
                regex.append(f"(?P<{name}>{WILDCARD_REGEX.get(name, r'[^/]+?')})")
        else:
            if "{" in part or "}" in part:
                return None, 0
            regex.append(re.escape(part))
            literal_len += len(part)
    return re.compile("^" + "".join(regex) + "$"), literal_len


def load_rule_patterns(rules_dir, bench_dirname):
    """Read `rule NAME:` blocks and their `benchmark:` paths from every .smk file.
    Also returns the rules that have no benchmark directive."""
    patterns, unmeasured = [], []
    if not os.path.isdir(rules_dir):
        logging.warning(f"Rules directory not found: {rules_dir}")
        return patterns, unmeasured

    for fname in sorted(os.listdir(rules_dir)):
        if not fname.endswith(".smk"):
            continue
        smk = fname[: -len(".smk")]
        with open(os.path.join(rules_dir, fname)) as handle:
            lines = handle.read().splitlines()

        rule = None
        rules_in_file, benchmarked = [], set()
        for i, line in enumerate(lines):
            match = re.match(r"^\s*(?:rule|checkpoint)\s+(\w+)\s*:", line)
            if match:
                rule = match.group(1)
                rules_in_file.append(rule)
                continue
            match = re.match(r"^\s*benchmark\s*:\s*(.*)$", line)
            if not (match and rule):
                continue
            text = match.group(1)
            j = i
            while not text.strip() and j + 1 < len(lines):
                j += 1
                text = lines[j]
            quoted = re.search(r"[\"']([^\"']+)[\"']", text)
            if not quoted:
                continue
            regex, literal_len = pattern_to_regex(quoted.group(1), bench_dirname)
            if regex is None:
                logging.warning(
                    f"Skipping malformed benchmark pattern in {fname} (rule {rule}): {quoted.group(1)}"
                )
                continue
            benchmarked.add(rule)
            patterns.append(
                {
                    "rule": rule,
                    "stage": SMK_TO_STAGE.get(smk) or stage_from_rule(rule),
                    "regex": regex,
                    "specificity": literal_len,
                }
            )
        unmeasured += [r for r in rules_in_file if r not in benchmarked]

    # Most specific patterns first so e.g. "{sample}.fastk_hap1" wins over "{sample}.fastk".
    patterns.sort(key=lambda p: -p["specificity"])
    return patterns, unmeasured


def read_benchmark(path):
    """Return the mean of each metric over the rows (repeats) of a benchmark file."""
    with open(path) as handle:
        rows = list(csv.DictReader(handle, delimiter="\t"))
    if not rows or "s" not in rows[0]:
        return None

    values = {}
    for metric in METRICS:
        numbers = []
        for row in rows:
            try:
                numbers.append(float(row.get(metric, "")))
            except (TypeError, ValueError):
                pass
        values[metric] = statistics.mean(numbers) if numbers else None
    if values["s"] is None:
        return None
    return values


def stage_from_rule(rule):
    for pattern, stage in RULE_TO_STAGE:
        if re.search(pattern, rule):
            return stage
    return "other"


def assembly_of(rule):
    """Assembly a job works on, from the rule name (e.g. ..._hap1_r1 -> Hap1)."""
    tokens = rule.lower().split("_")
    for token in tokens:
        if token in ASSEMBLY_TOKENS:
            label = ASSEMBLY_TOKENS[token]
            # Solo rules name the primary (p_ctg) and alternate (a_ctg) contigs hap1/hap2
            if tokens[0] == "solo":
                label = {"Hap1": "Primary", "Hap2": "Alternate"}.get(label, label)
            return label
    return ""


def classify(relpath, patterns):
    """Return (sample, rule, stage, matched) for a benchmark path."""
    for pattern in patterns:
        match = pattern["regex"].match(relpath)
        if match:
            sample = match.groupdict().get("sample") or ""
            return sample, pattern["rule"], pattern["stage"], True

    # Pipeasm names benchmark files {sample}.<rule>.txt
    stem = re.sub(r"\.(benchmark|txt|tsv)$", "", relpath.split("/")[-1])
    sample, _, rule = stem.partition(".")
    if not rule:
        sample, rule = "", stem
    return sample, rule, stage_from_rule(rule), False


def load_config(path):
    with open(path) as handle:
        return yaml.safe_load(handle) or {}


def load_threads(path):
    """{rule: threads} from config/threads.yaml (rules not listed use 1)."""
    threads = {}
    if not os.path.exists(path):
        return threads
    with open(path) as handle:
        for line in handle:
            match = re.match(r"^\s+(\w+):\s*(\d+)\s*(#.*)?$", line)
            if match:
                threads[match.group(1)] = int(match.group(2))
    return threads


def load_requested_memory(slurm_params, resources_dir, genome_size):
    """(default_mb, {rule: mem_mb}) requested in SLURM mode for this genome_size."""
    default = None
    if os.path.exists(slurm_params):
        with open(slurm_params) as handle:
            match = re.search(r"^\s+mem_mb:\s*(\d+)", handle.read(), re.M)
            if match:
                default = int(match.group(1))

    per_rule = {}
    path = os.path.join(resources_dir, f"{genome_size}.yaml")
    if genome_size and os.path.exists(path):
        with open(path) as handle:
            for line in handle:
                match = re.match(r"^\s+(\w+):\s*\{[^}]*\bmem_mb:\s*(\d+)", line)
                if match:
                    per_rule[match.group(1)] = int(match.group(2))
    return default, per_rule


# ── Aggregation ────────────────────────────────────────────────────────────────


def summarize(jobs):
    walls = [j["s"] for j in jobs]
    rss = [j["max_rss"] for j in jobs if j["max_rss"] is not None]
    cpu = [j["cpu_time"] for j in jobs if j["cpu_time"] is not None]
    return {
        "jobs": len(jobs),
        "wall": sum(walls),
        "max": max(walls),
        "cpu": sum(cpu),
        "rss_peak": max(rss) if rss else None,
        "peak_job": max(jobs, key=lambda j: j["max_rss"] or 0),
        "io_in": sum(j["io_in"] or 0 for j in jobs),
        "io_out": sum(j["io_out"] or 0 for j in jobs),
        "cores": (sum(cpu) / sum(walls)) if cpu and sum(walls) else None,
    }


def group_by(jobs, key):
    groups = defaultdict(list)
    for job in jobs:
        groups[key(job)].append(job)
    return groups


def stage_sort_key(stage):
    if stage in STAGE_ORDER:
        return (0, STAGE_ORDER.index(stage), stage)
    return (1, 0, stage)


def natural_key(text):
    return [int(t) if t.isdigit() else t.lower() for t in re.split(r"(\d+)", str(text))]


# ── Formatting ─────────────────────────────────────────────────────────────────


def esc(value):
    return html.escape(str(value), quote=True)


def fmt_duration(seconds):
    if seconds is None:
        return "–"
    seconds = float(seconds)
    if seconds < 60:
        return f"{seconds:.1f} s"
    if seconds < 3600:
        minutes, sec = divmod(int(round(seconds)), 60)
        return f"{minutes}m {sec:02d}s"
    hours, rem = divmod(int(round(seconds)), 3600)
    if hours < 48:
        return f"{hours}h {rem // 60:02d}m"
    return f"{seconds / 86400:.1f} d"


def fmt_hours(seconds):
    hours = seconds / 3600
    if hours < 10:
        return f"{hours:,.2f} h"
    return f"{hours:,.1f} h"


def fmt_mem(mb):
    """Decimal units, like the mem_mb values in config/resources (200000 MB = 200 GB)."""
    if mb is None:
        return "–"
    if mb < 1000:
        return f"{mb:,.0f} MB"
    return f"{mb / 1000:,.1f} GB"


def fmt_num(value, digits=2):
    if value is None:
        return "–"
    return f"{value:,.{digits}f}"


def stage_label(stage):
    return STAGE_LABELS.get(stage, stage.replace("_", " ").capitalize())


def stage_var(stage):
    if stage in STAGE_ORDER:
        return f"var(--c{STAGE_ORDER.index(stage) + 1})"
    return "var(--c-other)"


def td(value, sort=None, cls="num"):
    sort_attr = f' data-sort="{esc(sort)}"' if sort is not None else ""
    cls_attr = f' class="{cls}"' if cls else ""
    return f"<td{cls_attr}{sort_attr}>{value}</td>"


def swatch_label(stage):
    return f'<span class="swatch" style="background:{stage_var(stage)}"></span>{esc(stage_label(stage))}'


# ── HTML building blocks ───────────────────────────────────────────────────────


def build_kpis(jobs, samples):
    total = summarize(jobs)
    longest = max(jobs, key=lambda j: j["s"])
    peak = total["peak_job"]

    def where(job):
        return job["rule"] if len(samples) == 1 else f"{job['rule']} · {job['sample']}"

    if len(samples) == 1:
        sample_tile = ("Sample", samples[0], "with benchmark data")
    else:
        sample_tile = ("Samples", f"{len(samples)}", "with benchmark data")

    tiles = [
        sample_tile,
        ("Jobs", f"{len(jobs):,}", f"{len({j['rule'] for j in jobs})} distinct rules"),
        ("Cumulative runtime", fmt_hours(total["wall"]), "sum of job wall-clock times"),
        ("CPU time", fmt_hours(total["cpu"]), f"avg. {fmt_num(total['cores'])} cores per running job"),
        ("Peak memory", fmt_mem(peak["max_rss"]), where(peak)),
        ("Longest job", fmt_duration(longest["s"]), where(longest)),
        ("Data read / written", f"{fmt_mem(total['io_in'])} / {fmt_mem(total['io_out'])}", "sum over all jobs"),
    ]
    out = ['<section class="kpis">']
    for label, value, note in tiles:
        out.append(
            f'<div class="kpi"><div class="kpi-label">{esc(label)}</div>'
            f'<div class="kpi-value" title="{esc(value)}">{esc(value)}</div>'
            f'<div class="kpi-note" title="{esc(note)}">{esc(note)}</div></div>'
        )
    out.append("</section>")
    return "\n".join(out)


def build_config_section(config):
    if not config:
        return ""
    items = []
    for field in CONFIG_FIELDS:
        value = config.get(field)
        text = "" if value is None else str(value).strip()
        if field in PATH_FIELDS and text:
            shown = os.path.basename(text.rstrip("/")) or text
        else:
            shown = text
        if not text:
            items.append(f'<div class="kv"><dt>{esc(field)}</dt><dd class="muted">not set</dd></div>')
        else:
            items.append(f'<div class="kv"><dt>{esc(field)}</dt><dd title="{esc(text)}">{esc(shown)}</dd></div>')
    return f"""
<section>
  <h2>Run configuration</h2>
  <p class="sub">Inputs and options that decide which rules run. Hover a value to see the full path.</p>
  <dl class="kv-grid">{''.join(items)}</dl>
</section>"""


def build_stage_section(jobs, stages):
    by_stage = group_by(jobs, lambda j: j["stage"])
    total_wall = sum(j["s"] for j in jobs)
    max_wall = max(sum(j["s"] for j in g) for g in by_stage.values())

    bars, rows = [], []
    for stage in stages:
        stats = summarize(by_stage[stage])
        share = stats["wall"] / total_wall * 100 if total_wall else 0
        width = stats["wall"] / max_wall * 100 if max_wall else 0
        tip = f"{stage_label(stage)}|Runtime {fmt_hours(stats['wall'])} ({share:.1f}%)|{stats['jobs']:,} jobs · peak memory {fmt_mem(stats['rss_peak'])}"
        bars.append(
            f'<div class="hbar-row"><div class="hbar-label">{esc(stage_label(stage))}</div>'
            f'<div class="hbar-track"><div class="hbar" tabindex="0" style="width:{width:.2f}%;background:{stage_var(stage)}" data-tip="{esc(tip)}"></div>'
            f'<span class="hbar-value">{esc(fmt_hours(stats["wall"]))} · {share:.1f}%</span></div></div>'
        )
        rows.append(
            "<tr>"
            f"<td>{swatch_label(stage)}</td>"
            + td(f"{stats['jobs']:,}", stats["jobs"])
            + td(fmt_hours(stats["wall"]), stats["wall"])
            + td(f"{share:.1f}%", share)
            + td(fmt_duration(stats["max"]), stats["max"])
            + td(fmt_hours(stats["cpu"]), stats["cpu"])
            + td(fmt_mem(stats["rss_peak"]), stats["rss_peak"] or 0)
            + "</tr>"
        )

    return f"""
<section>
  <h2>Runtime by stage</h2>
  <p class="sub">Cumulative wall-clock time of all jobs in each stage. Jobs run in parallel, so this measures compute spent, not elapsed time.</p>
  <div class="chart hbars">{''.join(bars)}</div>
  <div class="table-wrap"><table class="sortable">
    <thead><tr><th>Stage</th><th class="num">Jobs</th><th class="num">Runtime</th><th class="num">Share</th>
    <th class="num">Longest job</th><th class="num">CPU time</th><th class="num">Peak memory</th></tr></thead>
    <tbody>{''.join(rows)}</tbody>
  </table></div>
</section>"""


def build_rule_section(jobs, samples, threads, requested):
    by_rule = group_by(jobs, lambda j: j["rule"])
    total_wall = sum(j["s"] for j in jobs)
    ordered = sorted(by_rule.items(), key=lambda kv: -sum(j["s"] for j in kv[1]))
    multi = len(samples) > 1

    rows = []
    for rule, group in ordered:
        stats = summarize(group)
        share = stats["wall"] / total_wall * 100 if total_wall else 0
        req = requested(rule)
        used = stats["rss_peak"] / req if req and stats["rss_peak"] is not None else None
        used_cls = "num warn" if used is not None and used >= MEMORY_WARN else "num"
        cells = [
            f'<td class="mono">{esc(rule)}</td>',
            f'<td class="nowrap">{swatch_label(group[0]["stage"])}</td>',
            f'<td>{esc(group[0]["assembly"]) or "–"}</td>',
        ]
        if multi:
            cells.append(td(f"{stats['jobs']:,}", stats["jobs"]))
        cells += [
            td(fmt_duration(stats["wall"]), stats["wall"]),
            td(f"{share:.1f}%", share),
            td(fmt_hours(stats["cpu"]), stats["cpu"]),
            td(fmt_num(stats["cores"]), stats["cores"] or 0),
            td(f"{threads.get(rule, 1)}", threads.get(rule, 1)),
            td(fmt_mem(stats["rss_peak"]), stats["rss_peak"] or 0),
            td(fmt_mem(req) if req else "–", req or 0),
            td(f"{used * 100:.0f}%" if used is not None else "–", used if used is not None else -1, used_cls),
        ]
        rows.append("<tr>" + "".join(cells) + "</tr>")

    jobs_head = '<th class="num">Jobs</th>' if multi else ""
    return f"""
<section>
  <h2>Resources by rule</h2>
  <p class="sub">One row per Snakemake rule, sorted by runtime. Avg. cores is CPU time divided by wall-clock time; compare it with the threads
  set in <span class="mono">config/threads.yaml</span>. Requested is the memory asked for in SLURM mode for this <span class="mono">genome_size</span>
  (<span class="mono">config/resources/</span>); jobs using {MEMORY_WARN:.0%} or more of it are highlighted. Click a column header to sort.</p>
  <div class="table-wrap"><table class="sortable">
    <thead><tr><th>Rule</th><th>Stage</th><th>Assembly</th>{jobs_head}<th class="num">Runtime</th><th class="num">Share</th>
    <th class="num">CPU time</th><th class="num">Avg. cores</th><th class="num">Threads</th>
    <th class="num">Peak memory</th><th class="num">Requested</th><th class="num">Used</th></tr></thead>
    <tbody>{''.join(rows)}</tbody>
  </table></div>
</section>"""


def build_memory_section(jobs, requested, top):
    by_rule = group_by(jobs, lambda j: j["rule"])
    peaks = {r: summarize(g)["rss_peak"] or 0 for r, g in by_rule.items()}
    ordered = sorted(peaks, key=lambda r: -peaks[r])[:top]
    if not ordered or not peaks[ordered[0]]:
        return ""
    max_peak = peaks[ordered[0]]

    bars = []
    for rule in ordered:
        stage = by_rule[rule][0]["stage"]
        req = requested(rule)
        width = peaks[rule] / max_peak * 100
        note = f" of {fmt_mem(req)}" if req else ""
        tip = f"{rule}|Peak memory {fmt_mem(peaks[rule])}{note}|{stage_label(stage)}"
        bars.append(
            f'<div class="hbar-row"><div class="hbar-label mono" title="{esc(rule)}">{esc(rule)}</div>'
            f'<div class="hbar-track"><div class="hbar" tabindex="0" style="width:{width:.2f}%;background:{stage_var(stage)}" data-tip="{esc(tip)}"></div>'
            f'<span class="hbar-value">{esc(fmt_mem(peaks[rule]))}{esc(note)}</span></div></div>'
        )

    return f"""
<section>
  <h2>Peak memory by rule</h2>
  <p class="sub">The {len(ordered)} rules with the highest peak memory (max_rss), with the memory requested in SLURM mode when known.</p>
  <div class="chart hbars wide-labels">{''.join(bars)}</div>
</section>"""


def build_sample_section(jobs, stages, samples):
    """Only shown when the benchmarks directory holds more than one sample."""
    if len(samples) < 2:
        return ""
    by_sample = group_by(jobs, lambda j: j["sample"])
    walls = {s: sum(j["s"] for j in by_sample[s]) for s in samples}
    max_wall = max(walls.values()) or 1

    legend = "".join(
        f'<span class="legend-item">{swatch_label(s)}</span>' for s in stages
    )
    bars, rows = [], []
    for sample in samples:
        group = by_sample[sample]
        by_stage = group_by(group, lambda j: j["stage"])
        segments = []
        for stage in stages:
            if stage not in by_stage:
                continue
            wall = sum(j["s"] for j in by_stage[stage])
            tip = f"{sample}|{stage_label(stage)}: {fmt_hours(wall)} ({wall / walls[sample] * 100:.1f}%)|{len(by_stage[stage]):,} jobs"
            segments.append(
                f'<span class="seg" tabindex="0" style="width:{wall / max_wall * 100:.3f}%;background:{stage_var(stage)}" data-tip="{esc(tip)}"></span>'
            )
        bars.append(
            f'<div class="hbar-row"><div class="hbar-label mono">{esc(sample)}</div>'
            f'<div class="hbar-track"><div class="stack">{"".join(segments)}</div>'
            f'<span class="hbar-value">{esc(fmt_hours(walls[sample]))}</span></div></div>'
        )
        stats = summarize(group)
        rows.append(
            "<tr>"
            f'<td class="mono">{esc(sample)}</td>'
            + td(f"{stats['jobs']:,}", stats["jobs"])
            + td(fmt_hours(stats["wall"]), stats["wall"])
            + td(fmt_hours(stats["cpu"]), stats["cpu"])
            + td(fmt_mem(stats["rss_peak"]), stats["rss_peak"] or 0)
            + "</tr>"
        )

    return f"""
<section>
  <h2>Samples</h2>
  <p class="sub">Cumulative runtime per sample, split by stage. Hover a segment for details.</p>
  <div class="legend">{legend}</div>
  <div class="chart hbars stacked">{''.join(bars)}</div>
  <div class="table-wrap"><table class="sortable">
    <thead><tr><th>Sample</th><th class="num">Jobs</th><th class="num">Runtime</th><th class="num">CPU time</th><th class="num">Peak memory</th></tr></thead>
    <tbody>{''.join(rows)}</tbody>
  </table></div>
</section>"""


CSS = """
:root {
  color-scheme: light;
  --page: #f9f9f7; --surface: #fcfcfb; --ink: #0b0b0b; --ink-2: #52514e; --muted: #898781;
  --grid: #e1e0d9; --axis: #c3c2b7; --border: rgba(11,11,11,0.10); --hover: rgba(11,11,11,0.035);
  --accent: #2a78d6; --warn: #b54708;
  --c1: #2a78d6; --c2: #eb6834; --c3: #1baf7a; --c4: #eda100;
  --c5: #e87ba4; --c6: #008300; --c7: #4a3aa7; --c8: #e34948; --c-other: #898781;
}
@media (prefers-color-scheme: dark) {
  :root:not([data-theme="light"]) {
    color-scheme: dark;
    --page: #0d0d0d; --surface: #1a1a19; --ink: #ffffff; --ink-2: #c3c2b7; --muted: #898781;
    --grid: #2c2c2a; --axis: #383835; --border: rgba(255,255,255,0.10); --hover: rgba(255,255,255,0.04);
    --accent: #3987e5; --warn: #f79009;
    --c1: #3987e5; --c2: #d95926; --c3: #199e70; --c4: #c98500;
    --c5: #d55181; --c6: #008300; --c7: #9085e9; --c8: #e66767;
  }
}
:root[data-theme="dark"] {
  color-scheme: dark;
  --page: #0d0d0d; --surface: #1a1a19; --ink: #ffffff; --ink-2: #c3c2b7; --muted: #898781;
  --grid: #2c2c2a; --axis: #383835; --border: rgba(255,255,255,0.10); --hover: rgba(255,255,255,0.04);
  --accent: #3987e5; --warn: #f79009;
  --c1: #3987e5; --c2: #d95926; --c3: #199e70; --c4: #c98500;
  --c5: #d55181; --c6: #008300; --c7: #9085e9; --c8: #e66767;
}
* { box-sizing: border-box; }
html { -webkit-text-size-adjust: 100%; }
body {
  margin: 0; background: var(--page); color: var(--ink);
  font: 14px/1.5 system-ui, -apple-system, "Segoe UI", Roboto, sans-serif;
}
main { max-width: 1200px; margin: 0 auto; padding: 32px 24px 64px; }
header.top { display: flex; justify-content: space-between; align-items: flex-end; gap: 24px; flex-wrap: wrap;
  padding-bottom: 20px; border-bottom: 1px solid var(--grid); margin-bottom: 24px; }
.eyebrow { font-size: 12px; letter-spacing: .08em; text-transform: uppercase; color: var(--muted); margin: 0 0 4px; }
h1 { font-size: 26px; font-weight: 600; margin: 0; letter-spacing: -0.01em; }
h2 { font-size: 17px; font-weight: 600; margin: 0 0 4px; }
.meta { color: var(--ink-2); font-size: 13px; text-align: right; min-width: 0; max-width: 100%; }
.meta div { white-space: nowrap; overflow: hidden; text-overflow: ellipsis; max-width: 460px; }
section { min-width: 0; background: var(--surface); border: 1px solid var(--border); border-radius: 10px; padding: 20px 22px; margin-bottom: 20px; }
.sub { color: var(--ink-2); margin: 0 0 16px; font-size: 13px; }
.muted { color: var(--muted); }
.mono { font-family: ui-monospace, SFMono-Regular, Menlo, Consolas, monospace; font-size: 12.5px; }
.nowrap { white-space: nowrap; }

.kpis { display: grid; grid-template-columns: repeat(auto-fit, minmax(min(150px, 100%), 1fr)); gap: 0; padding: 0; overflow: hidden; }
.kpi { padding: 16px 18px; border-right: 1px solid var(--grid); border-bottom: 1px solid var(--grid); margin: 0 -1px -1px 0; min-width: 0; }
.kpi-label { font-size: 12px; color: var(--ink-2); }
.kpi-value { font-size: 22px; font-weight: 600; margin: 2px 0; white-space: nowrap; overflow: hidden; text-overflow: ellipsis; }
.kpi-note { font-size: 12px; color: var(--muted); overflow-wrap: anywhere; }

.kv-grid { display: grid; grid-template-columns: repeat(auto-fill, minmax(min(260px, 100%), 1fr)); gap: 0 24px; margin: 0; }
.kv { display: flex; justify-content: space-between; gap: 12px; padding: 6px 0; border-bottom: 1px solid var(--grid); min-width: 0; }
.kv dt { color: var(--ink-2); font-family: ui-monospace, SFMono-Regular, Menlo, Consolas, monospace; font-size: 12px; flex-shrink: 0; }
.kv dd { margin: 0; text-align: right; white-space: nowrap; overflow: hidden; text-overflow: ellipsis; min-width: 0; }

.chart { margin: 4px 0 20px; }
.hbar-row { display: grid; grid-template-columns: 190px 1fr; align-items: center; gap: 12px; min-height: 26px; }
.wide-labels .hbar-row { grid-template-columns: 300px 1fr; }
.stacked .hbar-row { min-height: 22px; }
.hbar-label { color: var(--ink-2); font-size: 13px; white-space: nowrap; overflow: hidden; text-overflow: ellipsis; text-align: right; }
.hbar-track { display: flex; align-items: center; gap: 8px; min-width: 0; border-left: 1px solid var(--axis); padding: 3px 0; }
.hbar { height: 16px; min-width: 2px; border-radius: 0 4px 4px 0; flex-shrink: 0; }
.stack { display: flex; gap: 2px; flex: 0 1 auto; width: calc(100% - 72px); }
.seg { height: 14px; min-width: 1px; display: block; }
.seg:last-child { border-radius: 0 4px 4px 0; }
.hbar-value { font-size: 12px; color: var(--ink-2); white-space: nowrap; font-variant-numeric: tabular-nums; }
.hbars .hbar-track > .hbar { max-width: calc(100% - 150px); }
.hbar:hover, .seg:hover, .hbar:focus, .seg:focus { filter: brightness(1.12); outline: none; }

.legend { display: flex; flex-wrap: wrap; gap: 6px 16px; margin: 0 0 10px; font-size: 12px; color: var(--ink-2); }
.legend-item { display: inline-flex; align-items: center; }
.swatch { display: inline-block; width: 10px; height: 10px; border-radius: 2px; margin-right: 6px; vertical-align: -1px; flex-shrink: 0; }

.table-wrap { overflow-x: auto; border: 1px solid var(--grid); border-radius: 8px; }
table { border-collapse: collapse; width: 100%; font-size: 13px; }
th, td { padding: 7px 10px; text-align: left; border-bottom: 1px solid var(--grid); vertical-align: middle; }
td { max-width: 340px; overflow: hidden; text-overflow: ellipsis; white-space: nowrap; }
tbody tr:last-child td { border-bottom: 0; }
th { font-weight: 600; font-size: 12px; color: var(--ink-2); background: var(--page); position: sticky; top: 0; white-space: nowrap; }
.sortable th { cursor: pointer; user-select: none; }
.sortable th:hover { color: var(--ink); }
th[aria-sort=ascending]::after { content: " ↑"; }
th[aria-sort=descending]::after { content: " ↓"; }
.num { text-align: right; font-variant-numeric: tabular-nums; }
tbody tr:hover { background: var(--hover); }
.warn { color: var(--warn); font-weight: 600; }
.notes { color: var(--ink-2); font-size: 12.5px; margin: 0; padding-left: 18px; }
.notes li { margin: 3px 0; }
footer { color: var(--muted); font-size: 12px; text-align: center; margin-top: 28px; }

#tip { position: fixed; pointer-events: none; z-index: 10; background: var(--surface); color: var(--ink);
  border: 1px solid var(--border); border-radius: 6px; padding: 8px 10px; font-size: 12px; line-height: 1.45;
  box-shadow: 0 4px 16px rgba(0,0,0,.12); opacity: 0; transition: opacity .08s; max-width: 320px; }
#tip b { display: block; font-weight: 600; margin-bottom: 2px; }
#tip span { display: block; color: var(--ink-2); }

@media (max-width: 640px) {
  main { padding: 20px 16px 48px; }
  section { padding: 16px; }
  .meta { text-align: left; }
  .hbar-row, .wide-labels .hbar-row { grid-template-columns: 110px 1fr; gap: 8px; }
}
@media print {
  body { background: #fff; }
  section { break-inside: avoid; border-color: #ddd; }
  #tip { display: none; }
}
"""

JS = """
(function () {
  // Tooltip for chart marks
  var tip = document.getElementById('tip');
  function show(el, x, y) {
    var parts = el.getAttribute('data-tip').split('|');
    tip.innerHTML = '';
    parts.forEach(function (p, i) {
      var n = document.createElement(i === 0 ? 'b' : 'span');
      n.textContent = p; tip.appendChild(n);
    });
    tip.style.opacity = 1;
    var r = tip.getBoundingClientRect();
    var left = Math.min(x + 14, window.innerWidth - r.width - 8);
    var top = y + 14 + r.height > window.innerHeight ? y - r.height - 10 : y + 14;
    tip.style.left = Math.max(8, left) + 'px'; tip.style.top = top + 'px';
  }
  document.querySelectorAll('[data-tip]').forEach(function (el) {
    el.addEventListener('mousemove', function (e) { show(el, e.clientX, e.clientY); });
    el.addEventListener('mouseleave', function () { tip.style.opacity = 0; });
    el.addEventListener('focus', function () { var r = el.getBoundingClientRect(); show(el, r.left + r.width / 2, r.bottom); });
    el.addEventListener('blur', function () { tip.style.opacity = 0; });
  });

  // Sortable tables (uses data-sort when present)
  document.querySelectorAll('table.sortable').forEach(function (table) {
    table.querySelectorAll('thead th').forEach(function (th, col) {
      th.addEventListener('click', function () {
        var asc = th.getAttribute('aria-sort') !== 'ascending';
        table.querySelectorAll('thead th').forEach(function (h) { h.removeAttribute('aria-sort'); });
        th.setAttribute('aria-sort', asc ? 'ascending' : 'descending');
        var body = table.tBodies[0];
        var rows = Array.prototype.slice.call(body.rows);
        function key(row) {
          var cell = row.cells[col];
          if (!cell) return '';
          var v = cell.getAttribute('data-sort');
          if (v === null) v = cell.textContent.trim();
          var n = parseFloat(v);
          return isNaN(n) || !/^-?[\\d.]+(e[-+]?\\d+)?$/i.test(v) ? v.toLowerCase() : n;
        }
        rows.sort(function (a, b) {
          var x = key(a), y = key(b);
          if (typeof x === 'number' && typeof y === 'number') return asc ? x - y : y - x;
          return asc ? String(x).localeCompare(String(y), undefined, {numeric: true})
                     : String(y).localeCompare(String(x), undefined, {numeric: true});
        });
        rows.forEach(function (r) { body.appendChild(r); });
      });
    });
  });
})();
"""


def build_html(context):
    sections = "\n".join(s for s in context["sections"] if s)
    notes = ""
    if context["notes"]:
        items = "".join(f"<li>{n}</li>" for n in context["notes"])
        notes = f'<section><h2>Notes</h2><ul class="notes">{items}</ul></section>'

    meta = [f"Generated {esc(context['generated'])}", f"Benchmarks: <span class=\"mono\">{esc(context['bench_dir'])}</span>"]
    if context["config_path"]:
        meta.append(f"Config: <span class=\"mono\" title=\"{esc(context['config_path'])}\">{esc(context['config_path'])}</span>")

    return f"""<!doctype html>
<html lang="en">
<head>
<meta charset="utf-8">
<meta name="viewport" content="width=device-width, initial-scale=1">
<title>Pipeasm Benchmark Report</title>
<style>{CSS}</style>
</head>
<body>
<main>
  <header class="top">
    <div>
      <p class="eyebrow">Pipeasm</p>
      <h1>Benchmark report{esc(context['title_suffix'])}</h1>
    </div>
    <div class="meta">{''.join(f'<div>{m}</div>' for m in meta)}</div>
  </header>
  {sections}
  {notes}
  <footer>Generated by workflow/scripts/benchmark_report.py from Snakemake benchmark files.</footer>
</main>
<div id="tip" role="tooltip"></div>
<script>{JS}</script>
</body>
</html>
"""


# ── Main ───────────────────────────────────────────────────────────────────────


def main():
    args = parse_arguments()
    logging.basicConfig(
        level=logging.INFO,
        stream=sys.stdout,
        format="%(asctime)s [%(levelname)s] %(message)s",
        datefmt="%Y-%m-%d %H:%M:%S",
    )

    bench_dir = os.path.normpath(args.benchmarks)
    if not os.path.isdir(bench_dir):
        logging.error(f"Benchmarks directory not found: {bench_dir}")
        sys.exit(1)

    config_path = args.config
    if config_path is None and os.path.exists("config/config.yaml"):
        config_path = "config/config.yaml"
    config = {}
    if config_path:
        if not os.path.exists(config_path):
            logging.error(f"Config file not found: {config_path}")
            sys.exit(1)
        config = load_config(config_path)
        logging.info(f"Loaded config {config_path}")
    config_sample = str(config.get("sample") or "").strip()
    genome_size = str(config.get("genome_size") or "").strip()

    patterns, unmeasured = load_rule_patterns(args.rules, os.path.basename(bench_dir))
    logging.info(f"Parsed {len(patterns)} benchmark patterns from {args.rules}")

    threads = load_threads(args.threads)
    default_mem, rule_mem = load_requested_memory(args.slurm_params, args.resources, genome_size)

    def requested(rule):
        return rule_mem.get(rule, default_mem)

    jobs, unreadable, unmatched = [], [], []
    for root, _, files in os.walk(bench_dir):
        for fname in files:
            path = os.path.join(root, fname)
            relpath = os.path.relpath(path, bench_dir).replace(os.sep, "/")
            try:
                values = read_benchmark(path)
            except (OSError, UnicodeDecodeError, csv.Error):
                values = None
            if values is None:
                unreadable.append(relpath)
                continue
            sample, rule, stage, matched = classify(relpath, patterns)
            if not matched:
                unmatched.append(relpath)
            jobs.append(
                {**values, "sample": sample, "rule": rule, "stage": stage, "assembly": assembly_of(rule)}
            )

    notes = []
    if config_sample and not args.all_samples:
        others = sorted({j["sample"] for j in jobs} - {config_sample}, key=natural_key)
        if others:
            notes.append(
                f"Benchmarks of {len(others)} other sample(s) found in the same directory were left out "
                f"(use <span class=\"mono\">--all-samples</span> to include them): "
                + ", ".join(f'<span class="mono">{esc(s)}</span>' for s in others[:20])
                + (" …" if len(others) > 20 else "")
            )
        jobs = [j for j in jobs if j["sample"] == config_sample]

    if not jobs:
        logging.error("No benchmark data found" + (f" for sample {config_sample}" if config_sample else ""))
        sys.exit(1)

    samples = sorted({j["sample"] for j in jobs}, key=natural_key)
    stages = sorted({j["stage"] for j in jobs}, key=stage_sort_key)

    if unmatched:
        notes.append(
            f"{len(unmatched)} benchmark file(s) did not match any rule in "
            f"<span class=\"mono\">{esc(args.rules)}</span>; their rule names were taken from the file names: "
            + ", ".join(f'<span class="mono">{esc(p)}</span>' for p in unmatched[:10])
            + (" …" if len(unmatched) > 10 else "")
        )
    if unreadable:
        notes.append(
            f"{len(unreadable)} file(s) could not be read as Snakemake benchmarks and were skipped: "
            + ", ".join(f'<span class="mono">{esc(p)}</span>' for p in unreadable[:10])
            + (" …" if len(unreadable) > 10 else "")
        )
    if not genome_size or not rule_mem:
        notes.append(
            "The requested memory could not be read (missing <span class=\"mono\">genome_size</span> in the config "
            f"or <span class=\"mono\">{esc(args.resources)}/&lt;genome_size&gt;.yaml</span>)"
            + ("; only the default memory is shown." if default_mem else ".")
        )
    notes.append(
        "Only the jobs that ran have benchmark files, so rules skipped by the chosen inputs and options "
        "(e.g., no Hi-C, <span class=\"mono\">solo_asm: \"No\"</span>, empty <span class=\"mono\">gxdb</span>) do not appear."
    )
    if unmeasured:
        notes.append(
            "Rules without a <span class=\"mono\">benchmark:</span> directive are not measured: "
            + ", ".join(f'<span class="mono">{esc(r)}</span>' for r in sorted(set(unmeasured)))
            + "."
        )
    notes.append(
        "Runtimes are summed over jobs; because Snakemake runs jobs in parallel, they reflect compute "
        "spent rather than elapsed time. Memory is the maximum resident set size (max_rss) of a single job. "
        "Requested memory applies to SLURM mode only."
    )

    context = {
        "generated": datetime.datetime.now().strftime("%Y-%m-%d %H:%M"),
        "bench_dir": bench_dir,
        "config_path": config_path,
        "title_suffix": f" · {samples[0]}" if len(samples) == 1 else "",
        "notes": notes,
        "sections": [
            build_kpis(jobs, samples),
            build_config_section(config),
            build_stage_section(jobs, stages),
            build_rule_section(jobs, samples, threads, requested),
            build_memory_section(jobs, requested, args.top),
            build_sample_section(jobs, stages, samples),
        ],
    }

    os.makedirs(os.path.dirname(os.path.abspath(args.output)), exist_ok=True)
    with open(args.output, "w") as handle:
        handle.write(build_html(context))

    logging.info(
        f"Report written to {args.output} ({len(jobs)} jobs, {len(samples)} sample(s), {len(stages)} stages)"
    )


if __name__ == "__main__":
    main()
