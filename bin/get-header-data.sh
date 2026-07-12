#!/usr/bin/env bash

# get run info (flowcell, rundate, basecaller model) from a fastq file header
#
# usage: get-header-data.sh <fastq_dir> [<regex>]

set -euo pipefail

FASTQDIR="${1:?usage: get-header-data.sh <fastq_dir> [<regex>]}"
REGEX="${2:-fast(q|q.gz)$}"

FASTQFILE=$(find "$FASTQDIR" -type f 2>/dev/null | grep -E "$REGEX" | sort | head -n 1 || true)

if [[ -z "$FASTQFILE" ]]; then
    echo "NA,NA,NA"
    exit 0
fi


HEADER=$(gzip -cd "$FASTQFILE" | head -n 1)
echo header: "$HEADER"

FLOWCELL=$(printf '%s\n' "$HEADER" | grep -oE 'PU:Z:[^[:space:]]+' | head -n 1 | cut -d: -f3)
if [[ -z "$FLOWCELL" ]]; then
    FLOWCELL=$(printf '%s\n' "$HEADER" | grep -oE 'flow_cell_id=[^[:space:]]+' | head -n 1 | cut -d= -f2)
fi

RUNDATE=$(printf '%s\n' "$HEADER" | grep -oE 'DT:Z:[^[:space:]]+' | head -n 1 | cut -d: -f3 | cut -dT -f1)
if [[ -z "$RUNDATE" ]]; then
    RUNDATE=$(printf '%s\n' "$HEADER" | grep -oE 'start_time=[^[:space:]]+' | head -n 1 | cut -d= -f2 | cut -dT -f1)
fi

BC_MODEL=$(printf '%s\n' "$HEADER" | grep -oE 'RG:Z:[^[:space:]]+' | head -n 1 | cut -d: -f3)
if [[ -n "$BC_MODEL" ]]; then
    BC_MODEL="${BC_MODEL#*_}"
    BC_MODEL="${BC_MODEL%_barcode*}"
else
    BC_MODEL=$(printf '%s\n' "$HEADER" | grep -oE 'model_version_id=[^[:space:]]+' | head -n 1 | cut -d= -f2)
fi

FLOWCELL="${FLOWCELL:-NA}"
RUNDATE="${RUNDATE:-NA}"
BC_MODEL="${BC_MODEL:-NA}"

echo "${FLOWCELL},${RUNDATE},${BC_MODEL}"
