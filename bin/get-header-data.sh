#!/usr/bin/env bash

# get run info (platform, flowcell, rundate, basecaller model) from a fastq/fasta/bam file
# platform is auto-detected as one of: ont, pacbio, illumina, unknown
#
# usage: get-header-data.sh <fastq_dir>

FASTQDIR="${1:?usage: get-header-data.sh <fastq_dir>}"
# Matches the Nextflow pattern: *.{bam,fasta,fastq,fastq.gz,fq,fq.gz}
REGEX='\.(fastq|fq|fasta|bam)(\.gz)?$'

FASTQFILE=$(find -L "$FASTQDIR" -type f 2>/dev/null | grep -E "$REGEX" | sort | head -n 1 || true)

if [[ -z "$FASTQFILE" ]]; then
    echo "Error: no fastq/bam/fasta file found in $FASTQDIR matching regex $REGEX" >&2
    exit 1
fi

# ---------------------------------------------------------------------------
# BAM vs Text-stream unified reader logic
# ---------------------------------------------------------------------------
if [[ "$FASTQFILE" =~ \.bam$ ]]; then
    if ! command -v samtools &> /dev/null; then
        echo "Error: samtools is required to process BAM files but was not found in PATH." >&2
        exit 1
    fi
    # Stitch the file's global @RG headers and the first alignment line together 
    HEADER=$(printf "%s\n%s" "$(samtools view -H "$FASTQFILE" | grep '^@RG')" "$(samtools view "$FASTQFILE" | head -n 1)")
else
    # Standard stream reader for fastq/fasta
    HEADER=$(zcat -f "$FASTQFILE" | head -n 1)
fi

if [[ -z "$HEADER" ]]; then
    echo "Error: could not read header data from $FASTQFILE" >&2
    exit 1
fi

# ---------------------------------------------------------------------------
# platform detection
# ---------------------------------------------------------------------------
detect_platform() {
    local h="$1"

    if [[ "$h" =~ PL:Z:PACBIO ]] || \
       [[ "$h" =~ PL:PACBIO ]] || \
       [[ "$h" =~ ^@?m[0-9]{5,6}_[0-9]{6}_[0-9]{6}(_s[0-9]+)?/ ]] || \
       [[ "$h" =~ /ccs([[:space:]]|$) ]]; then
        echo "pacbio"
        return
    fi

    if [[ "$h" =~ PL:Z:ONT ]] || \
       [[ "$h" =~ PL:ONT ]] || \
       [[ "$h" =~ flow_cell_id= ]] || \
       [[ "$h" =~ runid= ]] || \
       [[ "$h" =~ start_time= ]] || \
       [[ "$h" =~ model_version_id= ]] || \
       [[ "$h" =~ RG:Z:[^[:space:]]*_(dna|rna)_r[0-9] ]] || \
       [[ "$h" =~ st:Z: ]] || \
       [[ "$h" =~ \bt:Z: ]] || \
       [[ "$h" =~ fn:Z: ]]; then
        echo "ont"
        return
    fi

    if [[ "$h" =~ PL:Z:ILLUMINA ]] || \
       [[ "$h" =~ PL:ILLUMINA ]] || \
       [[ "$h" =~ ^@[A-Za-z0-9_-]+:[0-9]+:[A-Za-z0-9_-]+:[0-9]+:[0-9]+:[0-9]+:[0-9]+([[:space:]]|$) ]]; then
        echo "illumina"
        return
    fi

    echo "unknown"
}

PLATFORM=$(detect_platform "$HEADER")

FLOWCELL=""
RUNDATE=""
BC_MODEL=""

case "$PLATFORM" in
    illumina)
        # Check standard SAM header PU tag first, otherwise parse legacy raw illumina fastq read ID
        FLOWCELL=$(printf '%s\n' "$HEADER" | grep -oE 'PU:[^[:space:]]+' | head -n 1 | sed -E 's/PU:(Z:)?//')
        if [[ -z "$FLOWCELL" ]]; then
            FLOWCELL=$(printf '%s\n' "$HEADER" | awk -F: '{print $3}')
        fi
        ;;

    ont)
        # Capture either SAM block 'PU:flowcell' format or per-read tag 'PU:Z:flowcell'
        FLOWCELL=$(printf '%s\n' "$HEADER" | grep -oE 'PU:[^[:space:]]+' | head -n 1 | sed -E 's/PU:(Z:)?//')
        if [[ -z "$FLOWCELL" ]]; then
            FLOWCELL=$(printf '%s\n' "$HEADER" | grep -oE 'flow_cell_id=[^[:space:]]+' | head -n 1 | cut -d= -f2)
        fi
        if [[ -z "$FLOWCELL" ]]; then
            FN_TAG=$(printf '%s\n' "$HEADER" | grep -oE 'fn:Z:[^[:space:]]+' | head -n 1 | cut -d: -f3)
            if [[ -n "$FN_TAG" ]]; then
                FLOWCELL=$(printf '%s\n' "$FN_TAG" | cut -d_ -f1)
            fi
        fi

        # Extract Run Date from @RG DT field or ONT native read properties
        RUNDATE=$(printf '%s\n' "$HEADER" | grep -oE 'DT:[^[:space:]]+' | head -n 1 | sed -E 's/DT:(Z:)?//' | cut -dT -f1)
        if [[ -z "$RUNDATE" ]]; then
            RUNDATE=$(printf '%s\n' "$HEADER" | grep -oE 'start_time=[^[:space:]]+' | head -n 1 | cut -d= -f2 | cut -dT -f1)
        fi
        if [[ -z "$RUNDATE" ]]; then
            RUNDATE=$(printf '%s\n' "$HEADER" | grep -oE 'st:Z:[0-9]{4}-[0-9]{2}-[0-9]{2}' | head -n 1 | cut -d: -f3)
        fi
        if [[ -z "$RUNDATE" ]]; then
            RUNDATE=$(printf '%s\n' "$HEADER" | grep -oE '\bt:Z:[0-9]{4}-[0-9]{2}-[0-9]{2}' | head -n 1 | cut -d: -f3)
        fi

        BC_MODEL=$(printf '%s\n' "$HEADER" | grep -oE 'RG:Z:[^[:space:]]+' | head -n 1 | cut -d: -f3)
        if [[ -n "$BC_MODEL" ]]; then
            BC_MODEL="${BC_MODEL#*_}"
            BC_MODEL=$(printf '%s\n' "$BC_MODEL" | sed -E 's/(@v[0-9.]+)(_.*)?/\1/')
            BC_MODEL="${BC_MODEL%_barcode*}"
        else
            BC_MODEL=$(printf '%s\n' "$HEADER" | grep -oE 'model_version_id=[^[:space:]]+' | head -n 1 | cut -d= -f2)
        fi
        ;;

    pacbio)
        FLOWCELL=$(printf '%s\n' "$HEADER" | grep -oE 'PU:[^[:space:]]+' | head -n 1 | sed -E 's/PU:(Z:)?//')
        BC_MODEL=$(printf '%s\n' "$HEADER" | grep -oE 'BC:Z:[^[:space:]]+' | head -n 1 | cut -d: -f3)

        READID=$(printf '%s\n' "$HEADER" | grep -oE '^@?m[0-9]{5,6}_[0-9]{6}_[0-9]{6}(_s[0-9]+)?' | head -n 1)
        if [[ -n "$READID" ]]; then
            MOVIE="${READID#@}"
            if [[ -z "$FLOWCELL" ]]; then
                FLOWCELL="$MOVIE"
            fi
            if [[ -z "$RUNDATE" ]]; then
                YYMMDD=$(printf '%s\n' "$MOVIE" | cut -d_ -f2)
                if [[ "$YYMMDD" =~ ^[0-9]{6}$ ]]; then
                    RUNDATE="20${YYMMDD:0:2}-${YYMMDD:2:2}-${YYMMDD:4:2}"
                fi
            fi
        fi
        ;;

    *)
        FLOWCELL=$(printf '%s\n' "$HEADER" | grep -oE 'PU:[^[:space:]]+' | head -n 1 | sed -E 's/PU:(Z:)?//')
        if [[ -z "$FLOWCELL" ]]; then
            FLOWCELL=$(printf '%s\n' "$HEADER" | grep -oE 'flow_cell_id=[^[:space:]]+' | head -n 1 | cut -d= -f2)
        fi
        if [[ -z "$FLOWCELL" ]]; then
            FN_TAG=$(printf '%s\n' "$HEADER" | grep -oE 'fn:Z:[^[:space:]]+' | head -n 1 | cut -d: -f3)
            if [[ -n "$FN_TAG" ]]; then
                FLOWCELL=$(printf '%s\n' "$FN_TAG" | cut -d_ -f1)
            fi
        fi

        RUNDATE=$(printf '%s\n' "$HEADER" | grep -oE 'DT:[^[:space:]]+' | head -n 1 | sed -E 's/DT:(Z:)?//' | cut -dT -f1)
        if [[ -z "$RUNDATE" ]]; then
            RUNDATE=$(printf '%s\n' "$HEADER" | grep -oE 'start_time=[^[:space:]]+' | head -n 1 | cut -d= -f2 | cut -dT -f1)
        fi
        if [[ -z "$RUNDATE" ]]; then
            RUNDATE=$(printf '%s\n' "$HEADER" | grep -oE 'st:Z:[0-9]{4}-[0-9]{2}-[0-9]{2}' | head -n 1 | cut -d: -f3)
        fi
        if [[ -z "$RUNDATE" ]]; then
            RUNDATE=$(printf '%s\n' "$HEADER" | grep -oE '\bt:Z:[0-9]{4}-[0-9]{2}-[0-9]{2}' | head -n 1 | cut -d: -f3)
        fi

        BC_MODEL=$(printf '%s\n' "$HEADER" | grep -oE 'RG:Z:[^[:space:]]+' | head -n 1 | cut -d: -f3)
        if [[ -n "$BC_MODEL" ]]; then
            BC_MODEL="${BC_MODEL#*_}"
            BC_MODEL=$(printf '%s\n' "$BC_MODEL" | sed -E 's/(@v[0-9.]+)(_.*)?/\1/')
            BC_MODEL="${BC_MODEL%_barcode*}"
        else
            BC_MODEL=$(printf '%s\n' "$HEADER" | grep -oE 'model_version_id=[^[:space:]]+' | head -n 1 | cut -d= -f2)
        fi
        ;;
esac

FLOWCELL="${FLOWCELL:-NA}"
RUNDATE="${RUNDATE:-NA}"
BC_MODEL="${BC_MODEL:-NA}"

echo "${PLATFORM},${FLOWCELL},${RUNDATE},${BC_MODEL}"
