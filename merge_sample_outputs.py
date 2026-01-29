import argparse
import csv
import json
import os
from typing import Dict, List, Optional


def read_tsv_header(path: str) -> List[str]:
    with open(path, "r", newline="") as f:
        r = csv.reader(f, delimiter="\t")
        return next(r)


def append_tsv_rows(src_path: str, dst_fh, *, header: List[str], wrote_header: bool) -> bool:
    with open(src_path, "r", newline="") as f:
        r = csv.reader(f, delimiter="\t")
        src_header = next(r)
        if src_header != header:
            raise ValueError(f"Header mismatch in {src_path}")
        w = csv.writer(dst_fh, delimiter="\t")
        if not wrote_header:
            w.writerow(header)
            wrote_header = True
        for row in r:
            w.writerow(row)
    return wrote_header


def main() -> None:
    p = argparse.ArgumentParser()
    p.add_argument("--type", required=True, choices=["numt", "nupt"])
    p.add_argument("--out-dir", required=True)
    p.add_argument("--samples", default="", help="Comma-separated sample IDs to merge (default: all under out-dir)")
    args = p.parse_args()

    type_dir = os.path.join(args.out_dir, args.type)
    if not os.path.isdir(type_dir):
        raise FileNotFoundError(type_dir)

    samples: Optional[List[str]] = None
    if args.samples.strip():
        samples = [x.strip() for x in args.samples.split(",") if x.strip()]
    else:
        samples = sorted([x for x in os.listdir(type_dir) if os.path.isdir(os.path.join(type_dir, x))])

    region_paths: List[str] = []
    qc_paths: List[str] = []
    for s in samples:
        region_path = os.path.join(type_dir, s, "region_validation.tsv")
        qc_path = os.path.join(type_dir, s, "qc.json")
        if os.path.isfile(region_path):
            region_paths.append(region_path)
        if os.path.isfile(qc_path):
            qc_paths.append(qc_path)

    if not region_paths:
        raise FileNotFoundError(f"No region_validation.tsv found under {type_dir}")

    header = read_tsv_header(region_paths[0])
    combined_region_path = os.path.join(type_dir, "all_samples_region_validation.tsv")
    with open(combined_region_path, "w", newline="") as out_f:
        wrote = False
        for path in region_paths:
            wrote = append_tsv_rows(path, out_f, header=header, wrote_header=wrote)

    qc_all: List[Dict[str, object]] = []
    for path in qc_paths:
        with open(path, "r") as f:
            qc_all.append(json.load(f))
    qc_all_path = os.path.join(type_dir, "qc_all_samples.json")
    with open(qc_all_path, "w") as f:
        json.dump(qc_all, f, indent=2, sort_keys=True)


if __name__ == "__main__":
    main()

