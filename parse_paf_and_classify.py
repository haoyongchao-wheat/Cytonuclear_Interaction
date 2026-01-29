import argparse
import csv
import glob
import json
import os
from dataclasses import dataclass
from typing import Dict, Iterable, List, Optional, Set, Tuple


@dataclass(frozen=True)
class RegionMeta:
    sample: str
    region_id: str
    category: str
    orig_chrom: str
    orig_start: int
    orig_end: int
    orig_length: int
    extr_chrom: str
    extr_start: int
    extr_end: int
    extr_length: int
    left_flank: int
    right_flank: int
    total_flank: int


@dataclass
class AlignmentRecord:
    qname: str
    qlen: int
    qstart: int
    qend: int
    strand: str
    tname: str
    tlen: int
    tstart: int
    tend: int
    nmatch: int
    alen: int
    mapq: int
    tp: Optional[str]


def parse_region_id_from_tname(tname: str) -> str:
    return tname.split("::", 1)[0]


def parse_sample_from_region_id(region_id: str) -> str:
    return region_id.split("|", 1)[0]


def overlap_len(a0: int, a1: int, b0: int, b1: int) -> int:
    lo = max(a0, b0)
    hi = min(a1, b1)
    return max(0, hi - lo)


def load_position_mapping(csv_path: str) -> Dict[Tuple[str, str], RegionMeta]:
    mapping: Dict[Tuple[str, str], RegionMeta] = {}
    with open(csv_path, "r", newline="") as f:
        reader = csv.DictReader(f)
        required = {
            "sample",
            "region_id",
            "category",
            "orig_chrom",
            "orig_start",
            "orig_end",
            "orig_length",
            "extr_chrom",
            "extr_start",
            "extr_end",
            "extr_length",
            "left_flank",
            "right_flank",
            "total_flank",
        }
        missing = required - set(reader.fieldnames or [])
        if missing:
            raise ValueError(f"Missing columns in {csv_path}: {sorted(missing)}")
        for row in reader:
            sample = row["sample"]
            region_id = row["region_id"]
            mapping[(sample, region_id)] = RegionMeta(
                sample=sample,
                region_id=region_id,
                category=row["category"],
                orig_chrom=row["orig_chrom"],
                orig_start=int(row["orig_start"]),
                orig_end=int(row["orig_end"]),
                orig_length=int(row["orig_length"]),
                extr_chrom=row["extr_chrom"],
                extr_start=int(row["extr_start"]),
                extr_end=int(row["extr_end"]),
                extr_length=int(row["extr_length"]),
                left_flank=int(row["left_flank"]),
                right_flank=int(row["right_flank"]),
                total_flank=int(row["total_flank"]),
            )
    return mapping


def parse_tp_tag(fields: List[str]) -> Optional[str]:
    for x in fields[12:]:
        if x.startswith("tp:A:") and len(x) >= 6:
            return x[5:]
    return None


def iter_paf_records(paf_path: str) -> Iterable[AlignmentRecord]:
    with open(paf_path, "r") as f:
        for line in f:
            if not line or line.startswith("#"):
                continue
            parts = line.rstrip("\n").split("\t")
            if len(parts) < 12:
                continue
            tp = parse_tp_tag(parts)
            yield AlignmentRecord(
                qname=parts[0],
                qlen=int(parts[1]),
                qstart=int(parts[2]),
                qend=int(parts[3]),
                strand=parts[4],
                tname=parts[5],
                tlen=int(parts[6]),
                tstart=int(parts[7]),
                tend=int(parts[8]),
                nmatch=int(parts[9]),
                alen=int(parts[10]),
                mapq=int(parts[11]),
                tp=tp,
            )


def choose_better_alignment(a: AlignmentRecord, b: AlignmentRecord) -> AlignmentRecord:
    a_span = a.tend - a.tstart
    b_span = b.tend - b.tstart
    if b_span != a_span:
        return b if b_span > a_span else a
    if b.alen != a.alen:
        return b if b.alen > a.alen else a
    if b.mapq != a.mapq:
        return b if b.mapq > a.mapq else a
    return a


def compute_evidence(
    aln: AlignmentRecord,
    meta: RegionMeta,
    *,
    edge_window: int,
    min_target_cov: float,
    min_overlap_bp: int,
    insert_min_overlap_for_single_junction: int,
) -> Tuple[float, int, bool, bool]:
    tlen = aln.tlen
    tstart = aln.tstart
    tend = aln.tend
    target_cov = (tend - tstart) / tlen if tlen > 0 else 0.0

    left_flank = max(0, min(meta.left_flank, tlen))
    right_flank = max(0, min(meta.right_flank, tlen))
    insert_start = left_flank
    insert_end = max(insert_start, tlen - right_flank)
    right_start = insert_end

    left_ok = left_flank == 0 or overlap_len(tstart, tend, 0, left_flank) >= min_overlap_bp
    insert_ok = overlap_len(tstart, tend, insert_start, insert_end) >= min_overlap_bp
    right_ok = right_flank == 0 or overlap_len(tstart, tend, right_start, tlen) >= min_overlap_bp

    full_span = (
        target_cov >= min_target_cov
        and tstart <= edge_window
        and tend >= tlen - edge_window
        and left_ok
        and insert_ok
        and right_ok
    )

    insert_overlap = overlap_len(tstart, tend, insert_start, insert_end)
    left_junction_supported = tstart < insert_start and tend > insert_start
    right_junction_supported = tstart < insert_end and tend > insert_end
    single_junction = insert_overlap >= insert_min_overlap_for_single_junction and (
        left_junction_supported or right_junction_supported
    )

    return target_cov, insert_overlap, full_span, single_junction


def safe_mkdir(path: str) -> None:
    os.makedirs(path, exist_ok=True)


def write_tsv(path: str, header: List[str], rows: Iterable[List[str]]) -> None:
    with open(path, "w", newline="") as f:
        w = csv.writer(f, delimiter="\t")
        w.writerow(header)
        for row in rows:
            w.writerow(row)


def classify_one_paf(
    paf_path: str,
    *,
    mapping: Dict[Tuple[str, str], RegionMeta],
    regions_by_sample: Dict[str, List[RegionMeta]],
    type_label: str,
    out_dir: str,
    mapq_min: int,
    include_tp: Set[str],
    edge_window: int,
    min_target_cov: float,
    min_overlap_bp: int,
    insert_min_overlap_for_single_junction: int,
    min_support_reads: int,
    write_read_details: bool,
) -> Tuple[List[List[str]], Dict[str, object]]:
    base = os.path.basename(paf_path)
    sample = base.split("_", 1)[0]

    sample_dir = os.path.join(out_dir, type_label, sample)
    safe_mkdir(sample_dir)

    regions = regions_by_sample.get(sample, [])
    meta_by_region_id = {m.region_id: m for m in regions}

    best_by_pair: Dict[Tuple[str, str], AlignmentRecord] = {}
    qc = {
        "paf_path": paf_path,
        "sample": sample,
        "type": type_label,
        "mapq_min": mapq_min,
        "include_tp": sorted(include_tp),
        "total_records": 0,
        "kept_records": 0,
        "unique_pairs": 0,
        "missing_region_meta": 0,
        "tlen_vs_extr_length_mismatch_count": 0,
        "tlen_vs_extr_length_max_abs_diff": 0,
    }

    for aln in iter_paf_records(paf_path):
        qc["total_records"] = int(qc["total_records"]) + 1
        if aln.mapq < mapq_min:
            continue
        if aln.tp is not None and aln.tp not in include_tp:
            continue
        region_id = parse_region_id_from_tname(aln.tname)
        key = (aln.qname, region_id)
        existing = best_by_pair.get(key)
        if existing is None:
            best_by_pair[key] = aln
        else:
            best_by_pair[key] = choose_better_alignment(existing, aln)
        qc["kept_records"] = int(qc["kept_records"]) + 1

    qc["unique_pairs"] = len(best_by_pair)

    region_full_span_reads: Dict[str, Set[str]] = {}
    region_single_junction_reads: Dict[str, Set[str]] = {}
    region_max_target_cov: Dict[str, float] = {}
    region_max_insert_overlap: Dict[str, int] = {}
    region_n_pairs: Dict[str, int] = {}
    region_tlen_seen: Dict[str, int] = {}

    read_detail_header = [
        "sample",
        "type",
        "region_id",
        "qname",
        "tlen",
        "tstart",
        "tend",
        "target_cov",
        "insert_overlap",
        "full_span_evidence",
        "single_junction_evidence",
        "mapq",
        "tp",
    ]

    read_details_path = os.path.join(sample_dir, "read_details.tsv")
    read_details_fh = None
    read_details_writer = None
    if write_read_details:
        read_details_fh = open(read_details_path, "w", newline="")
        read_details_writer = csv.writer(read_details_fh, delimiter="\t")
        read_details_writer.writerow(read_detail_header)

    for (qname, region_id), aln in best_by_pair.items():
        meta = meta_by_region_id.get(region_id)
        if meta is None:
            qc["missing_region_meta"] = int(qc["missing_region_meta"]) + 1
            continue

        diff = abs(aln.tlen - meta.extr_length)
        if diff != 0:
            qc["tlen_vs_extr_length_mismatch_count"] = int(qc["tlen_vs_extr_length_mismatch_count"]) + 1
            qc["tlen_vs_extr_length_max_abs_diff"] = max(int(qc["tlen_vs_extr_length_max_abs_diff"]), diff)

        target_cov, insert_overlap, full_span, single_junction = compute_evidence(
            aln,
            meta,
            edge_window=edge_window,
            min_target_cov=min_target_cov,
            min_overlap_bp=min_overlap_bp,
            insert_min_overlap_for_single_junction=insert_min_overlap_for_single_junction,
        )

        region_n_pairs[region_id] = region_n_pairs.get(region_id, 0) + 1
        region_max_target_cov[region_id] = max(region_max_target_cov.get(region_id, 0.0), target_cov)
        region_max_insert_overlap[region_id] = max(region_max_insert_overlap.get(region_id, 0), insert_overlap)
        region_tlen_seen.setdefault(region_id, aln.tlen)

        if full_span:
            region_full_span_reads.setdefault(region_id, set()).add(qname)
        if single_junction:
            region_single_junction_reads.setdefault(region_id, set()).add(qname)

        if read_details_writer is not None:
            read_details_writer.writerow(
                [
                    sample,
                    type_label,
                    region_id,
                    qname,
                    str(aln.tlen),
                    str(aln.tstart),
                    str(aln.tend),
                    f"{target_cov:.6f}",
                    str(insert_overlap),
                    "1" if full_span else "0",
                    "1" if single_junction else "0",
                    str(aln.mapq),
                    aln.tp or "",
                ]
            )

    if read_details_fh is not None:
        read_details_fh.close()
    with open(os.path.join(sample_dir, "qc.json"), "w") as f:
        json.dump(qc, f, indent=2, sort_keys=True)

    region_rows: List[List[str]] = []
    region_header = [
        "sample",
        "type",
        "region_id",
        "validated_level",
        "extr_chrom",
        "extr_start",
        "extr_end",
        "extr_length",
        "left_flank",
        "right_flank",
        "insert_length",
        "tlen",
        "n_full_span_reads",
        "n_single_junction_reads",
        "max_target_cov",
        "max_insert_overlap",
        "n_read_region_pairs",
        "orig_chrom",
        "orig_start",
        "orig_end",
        "orig_length",
        "category",
    ]

    for meta in regions:
        region_id = meta.region_id
        extr_len = meta.extr_length
        insert_len = max(0, extr_len - meta.left_flank - meta.right_flank)

        full_span_reads = region_full_span_reads.get(region_id, set())
        single_junction_reads = region_single_junction_reads.get(region_id, set())

        tlen = str(region_tlen_seen.get(region_id, ""))

        validated = "Not_Validated"
        if extr_len <= 15000:
            if len(full_span_reads) >= min_support_reads:
                validated = "Validated_Full_Span"
        else:
            if len(full_span_reads) >= 1:
                validated = "Validated_Full_Span"
            elif len(single_junction_reads) >= min_support_reads:
                validated = "Validated_Single_Junction"

        region_rows.append(
            [
                sample,
                type_label,
                region_id,
                validated,
                meta.extr_chrom,
                str(meta.extr_start),
                str(meta.extr_end),
                str(meta.extr_length),
                str(meta.left_flank),
                str(meta.right_flank),
                str(insert_len),
                tlen,
                str(len(full_span_reads)),
                str(len(single_junction_reads)),
                f"{region_max_target_cov.get(region_id, 0.0):.6f}",
                str(region_max_insert_overlap.get(region_id, 0)),
                str(region_n_pairs.get(region_id, 0)),
                meta.orig_chrom,
                str(meta.orig_start),
                str(meta.orig_end),
                str(meta.orig_length),
                meta.category,
            ]
        )

    write_tsv(os.path.join(sample_dir, "region_validation.tsv"), region_header, region_rows)

    return region_rows, qc


def main() -> None:
    p = argparse.ArgumentParser()
    p.add_argument("--type", required=True, choices=["numt", "nupt"])
    p.add_argument("--paf-dir", required=True)
    p.add_argument("--position-mapping-csv", required=True)
    p.add_argument("--out-dir", required=True)
    p.add_argument("--samples", default="", help="Comma-separated sample IDs to process (default: all)")
    p.add_argument("--mapq-min", type=int, default=20)
    p.add_argument("--write-read-details", action="store_true")
    p.add_argument("--no-combined", action="store_true")
    p.add_argument("--include-tp", default="P", help="Comma-separated tp tags to include, e.g. P or P,S")
    p.add_argument("--edge-window", type=int, default=200)
    p.add_argument("--min-target-cov", type=float, default=0.95)
    p.add_argument("--min-overlap-bp", type=int, default=1)
    p.add_argument("--insert-min-overlap", type=int, default=10000)
    p.add_argument("--min-support-reads", type=int, default=2)
    args = p.parse_args()

    mapping = load_position_mapping(args.position_mapping_csv)
    regions_by_sample: Dict[str, List[RegionMeta]] = {}
    for (sample, _), meta in mapping.items():
        regions_by_sample.setdefault(sample, []).append(meta)
    include_tp = {x.strip() for x in args.include_tp.split(",") if x.strip()}

    safe_mkdir(args.out_dir)
    paf_paths = sorted(glob.glob(os.path.join(args.paf_dir, f"*_{args.type}_mapped.paf")))
    if not paf_paths:
        raise FileNotFoundError(f"No PAF files matched in {args.paf_dir} for type {args.type}")
    if args.samples.strip():
        keep = {x.strip() for x in args.samples.split(",") if x.strip()}
        paf_paths = [p for p in paf_paths if os.path.basename(p).split("_", 1)[0] in keep]
        if not paf_paths:
            raise FileNotFoundError(f"No PAF files matched requested samples: {sorted(keep)}")

    all_region_rows: List[List[str]] = []
    qc_rows: List[Dict[str, object]] = []

    for paf_path in paf_paths:
        region_rows, qc = classify_one_paf(
            paf_path,
            mapping=mapping,
            regions_by_sample=regions_by_sample,
            type_label=args.type,
            out_dir=args.out_dir,
            mapq_min=args.mapq_min,
            include_tp=include_tp,
            edge_window=args.edge_window,
            min_target_cov=args.min_target_cov,
            min_overlap_bp=args.min_overlap_bp,
            insert_min_overlap_for_single_junction=args.insert_min_overlap,
            min_support_reads=args.min_support_reads,
            write_read_details=args.write_read_details,
        )
        all_region_rows.extend(region_rows)
        qc_rows.append(qc)

    combined_dir = os.path.join(args.out_dir, args.type)
    safe_mkdir(combined_dir)
    if not args.no_combined:
        combined_region_path = os.path.join(combined_dir, "all_samples_region_validation.tsv")
        if all_region_rows:
            header = [
                "sample",
                "type",
                "region_id",
                "validated_level",
                "extr_chrom",
                "extr_start",
                "extr_end",
                "extr_length",
                "left_flank",
                "right_flank",
                "insert_length",
                "tlen",
                "n_full_span_reads",
                "n_single_junction_reads",
                "max_target_cov",
                "max_insert_overlap",
                "n_read_region_pairs",
                "orig_chrom",
                "orig_start",
                "orig_end",
                "orig_length",
                "category",
            ]
            write_tsv(combined_region_path, header, all_region_rows)

        with open(os.path.join(combined_dir, "qc_all_samples.json"), "w") as f:
            json.dump(qc_rows, f, indent=2, sort_keys=True)


if __name__ == "__main__":
    main()
