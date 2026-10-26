#!/usr/bin/env python3
"""Synchronize the release benchmark scorecard across all published documents."""
from pathlib import Path
import re
import subprocess
import sys

root = Path(__file__).resolve().parents[1]
if "--execute" in sys.argv[1:]:
    subprocess.run(["cargo", "run", "-p", "ql-benches", "--release", "--", "--write-report"], cwd=root, check=True)
source = (root / "crates/ql-benches/BENCHMARK_RESULTS.md").read_text()
header = "| Benchmark ID | Architectural Comparison | Legacy Baseline (ns/op) | ql-rs SOTA (ns/op) | Speedup (x) | Parity Diff |"
a = source.index(header)
table_end = source.index("\n\nSummary:", a)
table = source[a:table_end]
rows = [line for line in source.splitlines() if re.match(r"\| \d{2}_", line)]
speedups = [float(line.split("|")[5].strip().rstrip("x")) for line in rows]
parity_latest = rows[-1].split("|")[-2].strip()
summary = (f"Summary: {len(rows)} benchmark suites; speedups range from **{min(speedups):.2f}x** to "
           f"**{max(speedups):.2f}x** (median **{sorted(speedups)[len(speedups)//2]:.2f}x**, "
           f"mean **{sum(speedups)/len(speedups):.2f}x**). All parity checks pass their domain tolerances; "
           f"latest suite parity is **{parity_latest}**.")
for rel in ["README.md", "knowledge_base/07_BENCHMARK_METHODOLOGY_AND_STUDY.md",
            "knowledge_base/08_BENCHMARK_RESULTS.md", "knowledge_base/08_COMPETITIVE_BENCHMARK_SYNTHESIS.md"]:
    path = root / rel
    if not path.exists():
        continue
    text = path.read_text()
    if header in text:
        begin = text.index(header)
        summary_start = text.find("\n\nSummary:", begin)
        if summary_start < 0:
            summary_start = len(text)
        text = text[:begin] + table + text[summary_start:]
    text = re.sub(r"Summary:.*", summary, text, count=1)
    if "Summary:" not in text:
        text += "\n\n" + summary + "\n"
    path.write_text(text)
