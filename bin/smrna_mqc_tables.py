#!/usr/bin/env python3
"""Reshape the miRDeep2 collapsed-read table into MultiQC custom-content sections.

Input is the long-format table emitted by the awk step in modules/mirdeep2:

    library  readlength  base1  miRBaseMatch  #distinctReads  #reads

`library` is the three-character config code that mapper.pl prefixes onto every
collapsed read, so it is remapped back to the sample-sheet label where possible.

Emits *_mqc.yaml rather than *_mqc.tsv: the YAML form states sample/x/value
explicitly, so it does not depend on MultiQC's row-vs-column inference for the
TSV custom-content parser.
"""

import argparse
import csv
import glob
import json
import os
import sys
from collections import defaultdict

BASES = ["A", "C", "G", "T", "N"]

# The canonical mature-miRNA window. Used for the summary "% in miRNA range"
# metric only; the full distribution is always plotted.
MIRNA_MIN, MIRNA_MAX = 21, 23


def parse_args():
    p = argparse.ArgumentParser(description=__doc__,
                                formatter_class=argparse.RawDescriptionHelpFormatter)
    p.add_argument("table",
                   help="awk output: (quant_)mirmap_firstbase_readlengthcounts.txt")
    p.add_argument("--sample-sheet",
                   help="pipeline sample sheet, to map config codes back to labels")
    p.add_argument("--outdir", default=".", help="where to write the *_mqc.yaml files")
    p.add_argument("--prefix", default="smrna", help="id prefix for the MultiQC sections")
    p.add_argument("--fastp-dir",
                   help="directory of <sample>.fastp.json files. Adds the raw and m10 "
                        "read counts, which the collapsed table cannot supply.")
    p.add_argument("--mirbase", action="store_true",
                   help="QUANT ran, so the miRBaseMatch column is meaningful. Emits the "
                        "miRBase panel and metric even when nothing matched, so a run "
                        "that mapped nothing reads as 0%% rather than as a missing "
                        "section. Without it the flag is hard-coded 0 and both are "
                        "omitted, because the annotation was never performed.")
    return p.parse_args()


def read_sample_map(path):
    """config code -> sample label, from the pipeline sample sheet."""
    if not path or not os.path.exists(path):
        return {}
    mapping = {}
    with open(path, newline="") as fh:
        for row in csv.DictReader(fh):
            code = (row.get("config") or "").strip()
            label = (row.get("label") or "").strip()
            if code and label:
                mapping[code] = label
    return mapping


def read_fastp(dirname):
    """sample -> (raw, m10) from fastp JSONs, named <sample>.fastp.json."""
    out = {}
    for path in sorted(glob.glob(os.path.join(dirname or "", "*.fastp.json"))):
        sample = os.path.basename(path)[: -len(".fastp.json")]
        try:
            doc = json.load(open(path))
            # passed_filter_reads, not summary.after_filtering.total_reads: this is
            # the field MultiQC's fastp module shows as Reads After Filtering, so
            # m10 matches that column by construction rather than by coincidence.
            out[sample] = (doc["summary"]["before_filtering"]["total_reads"],
                           doc["filtering_result"]["passed_filter_reads"])
        except (ValueError, KeyError) as exc:
            sys.stderr.write(f"{path}: skipping, could not read summary ({exc})\n")
    return out


def read_table(path):
    """Yield (library, readlength, base1, mirbase_match, distinct, reads)."""
    with open(path) as fh:
        first = fh.readline()
        if first and not first.startswith("library"):
            fh.seek(0)  # no header, e.g. a pre-concatenated table
        for lineno, line in enumerate(fh, start=2):
            line = line.rstrip("\n")
            if not line.strip():
                continue
            fields = line.split("\t")
            if len(fields) < 6:
                sys.stderr.write(f"{path}:{lineno}: skipping short line\n")
                continue
            try:
                yield (fields[0],
                       int(fields[1]),
                       fields[2].upper(),
                       int(fields[3]),
                       int(fields[4]),
                       int(fields[5]))
            except ValueError:
                sys.stderr.write(f"{path}:{lineno}: skipping unparseable line\n")


def yaml_scalar(value):
    if isinstance(value, bool):
        return "true" if value else "false"
    if isinstance(value, (int, float)):
        return repr(value) if isinstance(value, float) else str(value)
    return '"{}"'.format(str(value).replace('"', '\\"'))


def yaml_key(value):
    """Numeric keys stay unquoted so MultiQC reads them as x-axis positions."""
    if isinstance(value, int):
        return str(value)
    return '"{}"'.format(str(value).replace('"', '\\"'))


def write_section(path, meta, data):
    """Write one custom-content section: scalar/nested meta, then data."""
    with open(path, "w") as out:
        for key, value in meta.items():
            if isinstance(value, dict):
                out.write(f"{key}:\n")
                for subkey, subvalue in value.items():
                    out.write(f"    {subkey}: {yaml_scalar(subvalue)}\n")
            else:
                out.write(f"{key}: {yaml_scalar(value)}\n")
        out.write("data:\n")
        for sample in sorted(data):
            out.write(f"    {yaml_key(sample)}:\n")
            for x in data[sample]:
                out.write(f"        {yaml_key(x)}: {data[sample][x]}\n")


def write_general_stats(path, headers, data):
    """generalstats sections take their column config as a list under pconfig."""
    with open(path, "w") as out:
        out.write('plot_type: "generalstats"\n')
        out.write("pconfig:\n")
        for name, cfg in headers:
            out.write(f"    - {name}:\n")
            for key, value in cfg.items():
                out.write(f"        {key}: {yaml_scalar(value)}\n")
        out.write("data:\n")
        for sample in sorted(data):
            out.write(f"    {yaml_key(sample)}:\n")
            for name, _ in headers:
                if name in data[sample]:
                    out.write(f"        {name}: {data[sample][name]}\n")


def pct(numerator, denominator):
    return round(100.0 * numerator / denominator, 2) if denominator else 0.0


def main():
    args = parse_args()
    sample_map = read_sample_map(args.sample_sheet)
    fastp = read_fastp(args.fastp_dir) if args.fastp_dir else {}

    reads_by_len = defaultdict(lambda: defaultdict(int))
    reads_by_base = defaultdict(lambda: defaultdict(int))
    matched_by_len = defaultdict(lambda: defaultdict(int))

    total_reads = defaultdict(int)
    total_distinct = defaultdict(int)
    total_matched = defaultdict(int)
    mirna_range_reads = defaultdict(int)

    lengths = set()
    any_match = False

    for library, length, base, match, distinct, reads in read_table(args.table):
        sample = sample_map.get(library, library)
        lengths.add(length)

        reads_by_len[sample][length] += reads
        reads_by_base[sample][base if base in BASES else "N"] += reads

        total_reads[sample] += reads
        total_distinct[sample] += distinct
        if MIRNA_MIN <= length <= MIRNA_MAX:
            mirna_range_reads[sample] += reads

        if match:
            any_match = True
            matched_by_len[sample][length] += reads
            total_matched[sample] += reads

    if not total_reads:
        sys.exit(f"error: no usable rows parsed from {args.table}")

    samples = sorted(total_reads)
    # Zero-fill the length axis so lines are continuous rather than gappy.
    length_axis = list(range(min(lengths), max(lengths) + 1))
    os.makedirs(args.outdir, exist_ok=True)
    prefix = args.prefix

    def out(name):
        return os.path.join(args.outdir, f"{prefix}_{name}_mqc.yaml")

    write_section(
        out("length_reads"),
        {
            "id": f"{prefix}_length_reads",
            "section_name": "smRNA read length distribution",
            "description": ("Total reads per length, after collapsing. Mature miRNAs "
                            "peak at 21-23 nt; a broad 28-34 nt shoulder usually means "
                            "degradation products or rRNA/tRNA fragments."),
            "plot_type": "linegraph",
            "pconfig": {
                "id": f"{prefix}_length_reads_plot",
                "title": "smRNA: read length distribution",
                "xlab": "Read length (nt)",
                "ylab": "Reads",
            },
        },
        {s: {l: reads_by_len[s].get(l, 0) for l in length_axis} for s in samples},
    )

    write_section(
        out("first_base"),
        {
            "id": f"{prefix}_first_base",
            "section_name": "smRNA 5' nucleotide bias",
            "description": ("First base of each read, weighted by read count. Mature "
                            "miRNAs and piRNAs are strongly 5'-U; a flat profile "
                            "suggests degradation or untrimmed random-mer adapters."),
            "plot_type": "bargraph",
            "pconfig": {
                "id": f"{prefix}_first_base_plot",
                "title": "smRNA: 5' nucleotide bias",
                "ylab": "% of reads",
                "cpswitch": False,
            },
        },
        {s: {b: pct(reads_by_base[s].get(b, 0), total_reads[s]) for b in BASES}
         for s in samples},
    )

    if args.mirbase:
        write_section(
            out("mirbase_by_length"),
            {
                "id": f"{prefix}_mirbase_by_length",
                "section_name": "smRNA miRBase-mappable fraction by length",
                "description": ("Percentage of reads at each length that aligned to a "
                                "miRBase hairpin. Should be high across 21-23 nt and "
                                "fall away outside it."),
                "plot_type": "linegraph",
                "pconfig": {
                    "id": f"{prefix}_mirbase_by_length_plot",
                    "title": "smRNA: miRBase-mappable fraction by length",
                    "xlab": "Read length (nt)",
                    "ylab": "% of reads at this length",
                    "ymax": 100,
                    "ymin": 0,
                },
            },
            {s: {l: pct(matched_by_len[s].get(l, 0), reads_by_len[s].get(l, 0))
                 for l in length_axis}
             for s in samples},
        )

    # shared_key read_count hands formatting to MultiQC's read_count_multiplier,
    # the same mechanism behind fastp's "Reads After Filtering" column, so every
    # count here renders in M with identical precision.
    headers = []
    if fastp:
        headers += [
            (f"{prefix}_raw", {
                "title": "Raw reads",
                "description": "Reads into fastp (before_filtering.total_reads)",
                "scale": "Greys",
                "shared_key": "read_count",
            }),
            (f"{prefix}_m10", {
                "title": "m10",
                "description": "Reads surviving fastp trimming "
                               "(filtering_result.passed_filter_reads); the same field "
                               "MultiQC's fastp module shows as Reads After Filtering",
                "scale": "Blues",
                "shared_key": "read_count",
            }),
        ]
    headers += [
        (f"{prefix}_total_reads", {
            "title": "m18",
            "description": "Collapsed reads assigned to this library in the miRDeep2 "
                           "table. Named m18 by convention, but the floor is whatever "
                           "reached mapper.pl: pass -l 18 to enforce 18 nt, otherwise it "
                           "is fastp's --length_required and m18 will track m10",
            "scale": "Blues",
            "shared_key": "read_count",
        }),
        (f"{prefix}_distinct", {
            "title": "Distinct seqs",
            "description": "Unique collapsed sequences",
            "scale": "Purples",
            "shared_key": "read_count",
        }),
        (f"{prefix}_pct_mirna_len", {
            "title": "% 21-23 nt",
            "description": f"Reads in the {MIRNA_MIN}-{MIRNA_MAX} nt mature-miRNA window",
            "suffix": "%",
            "max": 100,
            "min": 0,
            "scale": "RdYlGn",
        }),
        (f"{prefix}_pct_5p_u", {
            "title": "% 5'-U",
            "description": "Reads beginning with U (T), the mature-miRNA signature",
            "suffix": "%",
            "max": 100,
            "min": 0,
            "scale": "RdYlGn",
        }),
    ]
    if args.mirbase:
        headers.append((f"{prefix}_pct_mirbase", {
            "title": "% miRBase",
            "description": "Reads aligning to a miRBase hairpin",
            "suffix": "%",
            "max": 100,
            "min": 0,
            "scale": "RdYlGn",
        }))

    stats = {}
    for s in samples:
        row = {
            f"{prefix}_total_reads": total_reads[s],
            f"{prefix}_distinct": total_distinct[s],
            f"{prefix}_pct_mirna_len": pct(mirna_range_reads[s], total_reads[s]),
            f"{prefix}_pct_5p_u": pct(reads_by_base[s].get("T", 0), total_reads[s]),
        }
        if args.mirbase:
            row[f"{prefix}_pct_mirbase"] = pct(total_matched[s], total_reads[s])
        if s in fastp:
            raw, m10 = fastp[s]
            row[f"{prefix}_raw"] = raw
            row[f"{prefix}_m10"] = m10
        stats[s] = row

    write_general_stats(out("stats"), headers, stats)

    sys.stderr.write(
        "smrna_mqc_tables: {} librar{}, lengths {}-{}, miRBase flag {}\n".format(
            len(samples), "y" if len(samples) == 1 else "ies",
            min(lengths), max(lengths),
            ("annotated, {} matched".format("some" if any_match else "none")
             if args.mirbase else "not annotated (MAPPER-only run)")))
    if fastp:
        missing = [s for s in samples if s not in fastp]
        sys.stderr.write("smrna_mqc_tables: fastp counts for {}/{} samples{}\n".format(
            len(samples) - len(missing), len(samples),
            "; missing " + ", ".join(missing) if missing else ""))


if __name__ == "__main__":
    main()
