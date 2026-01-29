#!/usr/bin/env bash
set -euo pipefail

CONC="${1:-4}"
LOGDIR="07_paf_structural_validation_v1/out/logs"

SAMPLES=(AMN Abo BJ8 CM42 HD6172 JM22 JM47 KF11 MZM NC4 S4185 XN6028 XY6 YM158 ZM16 ZM22 ZM366)

mkdir -p "$LOGDIR"

printf '%s\n' "${SAMPLES[@]}" | xargs -n 1 -P "$CONC" bash -lc '
  python3 07_paf_structural_validation_v1/parse_paf_and_classify.py \
    --type nupt \
    --paf-dir 06_minimap2_mapping/nupt \
    --position-mapping-csv 04_bed_by_flanking_status_v2/nupt_position_mapping.csv \
    --out-dir 07_paf_structural_validation_v1/out \
    --samples "$1" \
    --no-combined \
    > "07_paf_structural_validation_v1/out/logs/nupt_${1}.log" 2>&1
' _ 

python3 07_paf_structural_validation_v1/merge_sample_outputs.py \
  --type nupt \
  --out-dir 07_paf_structural_validation_v1/out \
  > "$LOGDIR/nupt_MERGE.log" 2>&1

