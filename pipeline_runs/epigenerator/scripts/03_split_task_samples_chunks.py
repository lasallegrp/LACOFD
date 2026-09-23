#!/usr/bin/env python3
r"""
Step 3: split task_samples.yaml into chunks of ~N samples for running CpG_Me2 Part 1
(trim + align) as several parallel Snakemake instances, one chunk per instance.

All lanes of a sample stay in the same chunk. A sample is identified the same way CpG_Me2
merges lanes: the text before the first "_" of each YAML entry
(LA001_S1_L004_JLLW_Nova1479P_Williams_L2 -> LA001).

Usage (run from the run directory, e.g. pipeline_runs/epigenerator/EPI_AZ_01):
python3 ../scripts/03_split_task_samples_chunks.py --chunk-size 10

Writes (in --outdir, default chunks/):
  chunk_01.yaml, chunk_02.yaml, ...   same format as task_samples.yaml, subset of samples
  chunks.tsv                          chunk, number of samples, number of lanes, sample IDs

Existing chunk files are never overwritten (running Snakemake instances read them) unless
--overwrite is given.
"""

import argparse
import os
import re
import sys
from collections import defaultdict
from typing import Dict, List

SAFE = re.compile(r'^[A-Za-z0-9/][A-Za-z0-9_./-]*$')


def scalar(v: str) -> str:
    return v if SAFE.match(v) and not re.fullmatch(r'[0-9.eE+-]+', v) else "'" + v.replace("'", "''") + "'"


def read_task_yaml(path: str) -> Dict:
    """Read task_samples.yaml. Uses PyYAML if available, else a parser for this simple format."""
    try:
        import yaml
        with open(path) as fh:
            return yaml.safe_load(fh)
    except ImportError:
        pass
    data, key = {"samples": []}, None
    with open(path) as fh:
        for raw in fh:
            line = raw.rstrip("\n")
            if not line.strip() or line.strip() == "---" or line.lstrip().startswith("#"):
                continue
            if line.lstrip().startswith("- "):
                if key != "samples":
                    sys.exit(f"❌ Unexpected list item outside 'samples:' in {path}: {line}")
                data["samples"].append(line.lstrip()[2:].strip().strip("'\""))
            else:
                k, _, v = line.partition(":")
                key = k.strip()
                if v.strip():
                    data[key] = v.strip().strip("'\"")
    return data


def write_task_yaml(path: str, sequence_data: str, genome: str, samples: List[str]) -> None:
    lines = ["---", f"sequence_data: {scalar(sequence_data)}", f"genome: {scalar(genome)}", "samples:"]
    lines += [f"- {scalar(s)}" for s in samples]
    with open(path, "w") as fh:
        fh.write("\n".join(lines) + "\n")


def natural_key(s: str):
    return [int(t) if t.isdigit() else t for t in re.split(r'(\d+)', s)]


def main() -> None:
    ap = argparse.ArgumentParser(description="Split task_samples.yaml into per-sample chunks for parallel CpG_Me2 Part 1 runs.")
    ap.add_argument("--yaml", default="task_samples.yaml", help="Full sample list (default: task_samples.yaml)")
    ap.add_argument("--chunk-size", type=int, default=10, help="Samples (not lanes) per chunk (default: 10)")
    ap.add_argument("--outdir", default="chunks", help="Output directory (default: chunks/)")
    ap.add_argument("--overwrite", action="store_true", help="Replace existing chunk files")
    args = ap.parse_args()

    if args.chunk_size < 1:
        sys.exit("--chunk-size must be >= 1")
    cfg = read_task_yaml(args.yaml)
    for k in ("sequence_data", "genome", "samples"):
        if not cfg.get(k):
            sys.exit(f"❌ '{k}' missing or empty in {args.yaml}")

    lanes_by_id = defaultdict(list)
    for lane in cfg["samples"]:
        lanes_by_id[str(lane).split("_")[0]].append(str(lane))
    ids = sorted(lanes_by_id, key=natural_key)

    # Lane counts per sample: flag samples that differ from the most common count
    counts = defaultdict(list)
    for sid in ids:
        counts[len(lanes_by_id[sid])].append(sid)
    typical = max(counts, key=lambda c: len(counts[c]))
    odd = {c: s for c, s in counts.items() if c != typical}
    if odd:
        print(f"⚠️  Most samples have {typical} lanes, but:")
        for c, s in sorted(odd.items()):
            print(f"   {len(s)} sample(s) have {c}: {', '.join(s[:10])}{' ...' if len(s) > 10 else ''}")

    chunks = [ids[i:i + args.chunk_size] for i in range(0, len(ids), args.chunk_size)]
    os.makedirs(args.outdir, exist_ok=True)
    width = max(2, len(str(len(chunks))))
    paths = [os.path.join(args.outdir, f"chunk_{i + 1:0{width}d}.yaml") for i in range(len(chunks))]
    existing = [p for p in paths if os.path.exists(p)]
    if existing and not args.overwrite:
        sys.exit(f"❌ {len(existing)} chunk file(s) already exist in {args.outdir}/ (e.g. {existing[0]}). "
                 "Running Snakemake instances may be using them. Use --overwrite only if none are running.")

    tsv = ["chunk\tn_samples\tn_lanes\tsample_ids"]
    for path, chunk_ids in zip(paths, chunks):
        lanes = [lane for sid in chunk_ids for lane in sorted(lanes_by_id[sid])]
        write_task_yaml(path, cfg["sequence_data"], cfg["genome"], lanes)
        name = os.path.splitext(os.path.basename(path))[0]
        tsv.append(f"{name}\t{len(chunk_ids)}\t{len(lanes)}\t{','.join(chunk_ids)}")
    with open(os.path.join(args.outdir, "chunks.tsv"), "w") as fh:
        fh.write("\n".join(tsv) + "\n")

    print(f"Samples: {len(ids)}   Lanes: {len(cfg['samples'])}   Chunks: {len(chunks)} "
          f"(up to {args.chunk_size} samples each)")
    print(f"Wrote {args.outdir}/chunk_*.yaml and {args.outdir}/chunks.tsv")


if __name__ == "__main__":
    main()