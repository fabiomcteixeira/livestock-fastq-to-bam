#!/usr/bin/env python3
"""Summarize observed Nextflow tasks without inferring missing-task causes."""
import csv
import re
import sys
from pathlib import Path

manifest_path, trace_path, output_dir = map(Path, sys.argv[1:4])
output_dir.mkdir(parents=True, exist_ok=True)

def read_tsv(path):
    with path.open(newline="") as handle:
        return list(csv.DictReader(handle, delimiter="\t"))

def write_tsv(name, fields, rows):
    with (output_dir / name).open("w", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=fields, delimiter="\t")
        writer.writeheader()
        writer.writerows(rows)

samples = read_tsv(manifest_path)
trace_available = trace_path.is_file() and trace_path.stat().st_size > 0
tasks = read_tsv(trace_path) if trace_available else []

events = []
lookup = {}

for task in tasks:
    name = task.get("name", "")
    match = re.match(r"^(.*?) \((.*)\)$", name)
    process = (match.group(1) if match else name).split(":")[-1]
    tag = match.group(2) if match else ""
    workdir = task.get("workdir", "")
    stderr = (
        str(Path(workdir) / ".command.err")
        if workdir and workdir != "-" else ""
    )
    event = {
        "process": process,
        "tag": tag,
        "status": task.get("status", ""),
        "exit": task.get("exit", ""),
        "attempt": task.get("attempt", ""),
        "workdir": workdir,
        "stderr_path": stderr,
    }
    events.append(event)
    lookup.setdefault((process, tag), []).append(event)

fields = ["process", "tag", "status", "exit", "attempt",
          "workdir", "stderr_path"]
write_tsv("task_status.tsv", fields, events)
write_tsv(
    "failed_tasks.tsv", fields,
    [event for event in events if event["status"] == "FAILED"]
)

def last_status(process, tag):
    observed = lookup.get((process, tag), [])
    if observed:
        return observed[-1]["status"]
    return "NOT_OBSERVED" if trace_available else "TRACE_UNAVAILABLE"

run_rows = []
for sample in samples:
    animal = sample["animal_id"]
    run = sample["run"]
    run_rows.append({
        "animal_id": animal,
        "run": run,
        "FASTQC_RAW": last_status("FASTQC_RAW", run),
        "NORMALIZE_LB": last_status("NORMALIZE_LB", animal + "/" + run),
        "TRIM": last_status("TRIM", animal + "/" + run),
        "FASTQC_TRIM": last_status("FASTQC_TRIM", run),
        "MAP": last_status("MAP", animal + "/" + run),
    })

write_tsv(
    "run_status.tsv",
    ["animal_id", "run", "FASTQC_RAW", "NORMALIZE_LB",
     "TRIM", "FASTQC_TRIM", "MAP"],
    run_rows
)

animal_rows = []
for animal in sorted({sample["animal_id"] for sample in samples}):
    runs = [row for row in run_rows if row["animal_id"] == animal]
    mapped = sum(row["MAP"] in ("COMPLETED", "CACHED") for row in runs)
    animal_rows.append({
        "animal_id": animal,
        "expected_runs": len(runs),
        "mapped_runs": mapped,
        "all_runs_mapped": "YES" if mapped == len(runs) else "NO",
        "MERGE_RUNS": last_status("MERGE_RUNS", animal),
        "DEDUP": last_status("DEDUP", animal),
    })

write_tsv(
    "animal_status.tsv",
    ["animal_id", "expected_runs", "mapped_runs", "all_runs_mapped",
     "MERGE_RUNS", "DEDUP"],
    animal_rows
)

if not trace_available:
    print("WARNING: trace unavailable; task outcomes are unknown.",
          file=sys.stderr)
print("Execution summaries: {}".format(output_dir))
