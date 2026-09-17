#!/usr/bin/env python3
"""Fixtures for the MultiQC smoke test: a collapsed-read table and fastp JSONs."""
import json
import os

here = os.path.join(os.path.dirname(os.path.abspath(__file__)), "fixtures")
os.makedirs(os.path.join(here, "fastp"), exist_ok=True)

rows = ["library\treadlength\tbase1\tmiRBaseMatch\t#distinctReads\t#reads"]
for lib, scale in (("101", 24_000_000), ("102", 30_000_000)):
    for length in range(18, 52):
        for base in "ACGT":
            weight = 6 if 21 <= length <= 23 and base == "T" else (3 if 21 <= length <= 23 else 1)
            reads = scale * weight // 900
            match = 1 if 21 <= length <= 23 else 0
            rows.append(f"{lib}\t{length}\t{base}\t{match}\t{reads // 40}\t{reads}")
open(os.path.join(here, "table.txt"), "w").write("\n".join(rows) + "\n")

open(os.path.join(here, "sample-sheet.csv"), "w").write(
    "label,fastq1,config\nSS1,SS1_R1.fastq.gz,101\nSS2,SS2_R1.fastq.gz,102\n")

for sample, raw in (("SS1", 24_609_228), ("SS2", 30_591_338)):
    passed = int(raw * 0.991)
    json.dump({"summary": {"before_filtering": {"total_reads": raw},
                           "after_filtering": {"total_reads": passed}},
               "filtering_result": {"passed_filter_reads": passed},
               "command": "fastp -i SS_R1.fastq.gz"},
              open(os.path.join(here, "fastp", f"{sample}.fastp.json"), "w"))

open(os.path.join(here, "software_versions_mqc.yml"), "w").write(
    'id: pipeline_software_versions\nsection_name: Software Versions\n'
    'plot_type: html\ndata: "<table><tr><td>miRDeep2</td><td>2.0.1.3</td></tr></table>"\n')
print("fixtures written to", here)
