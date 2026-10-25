#!/usr/bin/env python3
"""Synchronize benchmark scorecard counts and suite 37 measurements."""
from pathlib import Path
import subprocess
import sys
root = Path(__file__).resolve().parents[1]
args = sys.argv[1:]
if "--execute" in args:
    subprocess.run(["cargo", "run", "-p", "ql-benches", "--release", "--", "--write-report"], cwd=root, check=True)
# Keep narrative count markers synchronized without touching historical phase reviews.
for rel in ["README.md", "knowledge_base/00_INDEX_AND_VERDICT.md", "knowledge_base/05_BLUEPRINT_MODERN_RUST_QUANTLIB.md", "knowledge_base/07_BENCHMARK_METHODOLOGY_AND_STUDY.md", "knowledge_base/08_BENCHMARK_RESULTS.md", "knowledge_base/08_COMPETITIVE_BENCHMARK_SYNTHESIS.md", "knowledge_base/09_BOOK_LIBRARY_AND_PHASE21_ROADMAP.md"]:
    p = root / rel
    if not p.exists():
        continue
    s = p.read_text()
    s = s.replace("36 benchmark suites", "37 benchmark suites").replace("36 suites", "37 suites").replace("36 benchmark", "37 benchmark")
    s = s.replace("34 workstation tabs", "35 workstation tabs").replace("34 tabs", "35 tabs").replace("34-tab", "35-tab")
    s = s.replace("36 parity checks", "37 parity checks")
    p.write_text(s)
