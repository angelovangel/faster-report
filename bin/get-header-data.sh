#!/usr/bin/env bash

# get run info (platform, flowcell, rundate, basecaller model, modifications, filetype)
# from a fastq/fasta/bam file
# platform is auto-detected as one of: ont, pacbio, illumina, unknown
#
# usage: get-header-data.sh <fastq_dir>
# output: platform,flowcell,rundate,bc_model,mods,filetype

FASTQDIR="${1:?usage: get-header-data.sh <fastq_dir>}"
# Matches the Nextflow pattern: *.{bam,fasta,fastq,fastq.gz,fq,fq.gz}
REGEX='\.(fastq|fq|fasta|bam)(\.gz)?$'

FASTQFILE=$(find -L "$FASTQDIR" -type f 2>/dev/null | grep -E "$REGEX" | sort | head -n 1 || true)

if [[ -z "$FASTQFILE" ]]; then
    echo "Error: no fastq/bam/fasta file found in $FASTQDIR matching regex $REGEX" >&2
    exit 1
fi

# ---------------------------------------------------------------------------
# modification code -> label
# letters: SAM spec single-letter codes; numbers: ChEBI IDs
# ---------------------------------------------------------------------------
label_mods() {
    # stdin: one code per line; stdout: labels joined by "/"
    sed -E '
        s/^m$/5mC/; s/^h$/5hmC/; s/^f$/5fC/; s/^c$/5caC/;
        s/^g$/5hmU/; s/^e$/5fU/; s/^b$/5caU/;
        s/^a$/6mA/; s/^o$/8oxoG/; s/^n$/Xao/;
        s/^21839$/4mC/; s/^27551$/5mC/; s/^76792$/5hmC/; s/^28871$/6mA/;
        s/^80961$/5fC/; s/^76794$/5caC/
    ' | sort -u | paste -sd/ -
}

# stdin: text containing MM:Z: tags; stdout: distinct mod codes, one per line
# (letter codes and ChEBI IDs; multi-letter codes such as C+mh are split into m and h)
mm_codes() {
    grep -o 'MM:Z:[^[:space:]]*' | sed 's/^MM:Z://' | tr ';' '\n' \
        | sed -nE 's/^[ACGTUN][+-]([a-z0-9]+)[.?].*/\1/p' \
        | awk '/^[a-z]+$/ { n = split($0, c, ""); for (i = 1; i <= n; i++) print c[i]; next } { print }' \
        | sort -u
}

# ---------------------------------------------------------------------------
# BAM vs text-stream reader
# ---------------------------------------------------------------------------
MOD_LABELS=""

if [[ "$FASTQFILE" =~ \.bam$ ]]; then
    if ! command -v samtools &> /dev/null; then
        echo "Error: samtools is required to process BAM files but was not found in PATH." >&2
        exit 1
    fi

    BAMHEADER=$(samtools view -H "$FASTQFILE")

    # @RG header lines + first alignment line
    HEADER=$(printf "%s\n%s" "$(printf '%s\n' "$BAMHEADER" | grep '^@RG')" "$(samtools view "$FASTQFILE" | head -n 1)")

    # Preferred: modbase_models in @RG DS, e.g.
    #   ..._sup@v5.2.0_6mA@v1,..._sup@v5.2.0_4mC_5mC@v1  ->  6mA, 4mC, 5mC
    MODS=$(printf '%s\n' "$BAMHEADER" | grep -oE 'modbase_models=[^[:space:]]+' | head -n 1 | cut -d= -f2 \
        | tr ',' '\n' | sed -E 's/@v[0-9.]+$//; s/.*@v[0-9.]+_//' | tr '_' '\n' | sort -u)

    if [[ -n "$MODS" ]]; then
        MOD_LABELS=$(printf '%s\n' "$MODS" | sort -u | paste -sd/ -)
    else
        # Fallback: MM tags from the first 1000 reads (letter codes and ChEBI IDs;
        # multi-letter codes such as C+mh are split into m and h)
        MODS=$(samtools view "$FASTQFILE" | head -n 1000 | mm_codes)
        [[ -n "$MODS" ]] && MOD_LABELS=$(printf '%s\n' "$MODS" | label_mods)
    fi
else
    # gzip -dcf works on both GNU and BSD/macOS (zcat -f does not on macOS)
    HEADER=$(gzip -dcf "$FASTQFILE" 2>/dev/null | head -n 1)

    # Mods from MM tags in the first 1000 read headers (ONT fastq written with --emit-moves/modbase tags)
    MODS=$(gzip -dcf "$FASTQFILE" 2>/dev/null | head -n 4000 | awk 'NR % 4 == 1' | mm_codes)
    [[ -n "$MODS" ]] && MOD_LABELS=$(printf '%s\n' "$MODS" | label_mods)
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
       [[ "$h" =~ (^|[[:space:]])t:Z: ]] || \
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

# ---------------------------------------------------------------------------
# field extractors
# ---------------------------------------------------------------------------
# flowcell from PU tag, flow_cell_id=, or fn:Z: (ONT)
extract_flowcell_generic() {
    local h="$1" fc fn
    fc=$(printf '%s\n' "$h" | grep -oE 'PU:[^[:space:]]+' | head -n 1 | sed -E 's/PU:(Z:)?//')
    if [[ -z "$fc" ]]; then
        fc=$(printf '%s\n' "$h" | grep -oE 'flow_cell_id=[^[:space:]]+' | head -n 1 | cut -d= -f2)
    fi
    if [[ -z "$fc" ]]; then
        fn=$(printf '%s\n' "$h" | grep -oE 'fn:Z:[^[:space:]]+' | head -n 1 | cut -d: -f3)
        [[ -n "$fn" ]] && fc=$(printf '%s\n' "$fn" | cut -d_ -f1)
    fi
    printf '%s' "$fc"
}

# run date from @RG DT, start_time=, st:Z:, or t:Z:
extract_rundate_generic() {
    local h="$1" d
    d=$(printf '%s\n' "$h" | grep -oE 'DT:[^[:space:]]+' | head -n 1 | sed -E 's/DT:(Z:)?//' | cut -dT -f1)
    if [[ -z "$d" ]]; then
        d=$(printf '%s\n' "$h" | grep -oE 'start_time=[^[:space:]]+' | head -n 1 | cut -d= -f2 | cut -dT -f1)
    fi
    if [[ -z "$d" ]]; then
        d=$(printf '%s\n' "$h" | grep -oE 'st:Z:[0-9]{4}-[0-9]{2}-[0-9]{2}' | head -n 1 | cut -d: -f3)
    fi
    if [[ -z "$d" ]]; then
        d=$(printf '%s\n' "$h" | grep -oE '(^|[[:space:]])t:Z:[0-9]{4}-[0-9]{2}-[0-9]{2}' | head -n 1 | sed -E 's/.*t:Z://')
    fi
    printf '%s' "$d"
}

# basecaller model from RG:Z: (strip run id prefix and barcode suffix) or model_version_id=
extract_bcmodel_generic() {
    local h="$1" m
    m=$(printf '%s\n' "$h" | grep -oE 'RG:Z:[^[:space:]]+' | head -n 1 | cut -d: -f3)
    if [[ -n "$m" ]]; then
        m="${m#*_}"
        m=$(printf '%s\n' "$m" | sed -E 's/(@v[0-9.]+)(_.*)?/\1/')
        m="${m%_barcode*}"
    else
        m=$(printf '%s\n' "$h" | grep -oE 'model_version_id=[^[:space:]]+' | head -n 1 | cut -d= -f2)
    fi
    printf '%s' "$m"
}

# last resort: ONT flowcell ID (e.g. PBM21646, FAX12345) in the file path / name
extract_flowcell_from_path() {
    printf '%s\n' "$1" | grep -oE '(^|[^A-Za-z0-9])[A-Z]{3}[0-9]{4,6}([^A-Za-z0-9]|$)' | head -n 1 | grep -oE '[A-Z]{3}[0-9]{4,6}'
}

PLATFORM=$(detect_platform "$HEADER")

FLOWCELL=""
RUNDATE=""
BC_MODEL=""

case "$PLATFORM" in
    illumina)
        FLOWCELL=$(printf '%s\n' "$HEADER" | grep -oE 'PU:[^[:space:]]+' | head -n 1 | sed -E 's/PU:(Z:)?//')
        if [[ -z "$FLOWCELL" ]]; then
            FLOWCELL=$(printf '%s\n' "$HEADER" | awk -F: '{print $3}')
        fi
        ;;

    pacbio)
        FLOWCELL=$(printf '%s\n' "$HEADER" | grep -oE 'PU:[^[:space:]]+' | head -n 1 | sed -E 's/PU:(Z:)?//')
        BC_MODEL=$(printf '%s\n' "$HEADER" | grep -oE 'BC:Z:[^[:space:]]+' | head -n 1 | cut -d: -f3)

        READID=$(printf '%s\n' "$HEADER" | grep -oE '^@?m[0-9]{5,6}_[0-9]{6}_[0-9]{6}(_s[0-9]+)?' | head -n 1)
        if [[ -n "$READID" ]]; then
            MOVIE="${READID#@}"
            [[ -z "$FLOWCELL" ]] && FLOWCELL="$MOVIE"
            YYMMDD=$(printf '%s\n' "$MOVIE" | cut -d_ -f2)
            if [[ "$YYMMDD" =~ ^[0-9]{6}$ ]]; then
                RUNDATE="20${YYMMDD:0:2}-${YYMMDD:2:2}-${YYMMDD:4:2}"
            fi
        fi
        ;;

    ont|*)
        FLOWCELL=$(extract_flowcell_generic "$HEADER")
        RUNDATE=$(extract_rundate_generic "$HEADER")
        BC_MODEL=$(extract_bcmodel_generic "$HEADER")
        if [[ -z "$FLOWCELL" ]]; then
            FLOWCELL=$(extract_flowcell_from_path "$FASTQFILE")
        fi
        ;;
esac

FLOWCELL="${FLOWCELL:-NA}"
RUNDATE="${RUNDATE:-NA}"
BC_MODEL="${BC_MODEL:-NA}"
MOD_LABELS="${MOD_LABELS:--}"

if [[ "$FASTQFILE" =~ \.bam$ ]]; then
    FILETYPE="bam"
else
    FILETYPE="fastq"
fi

echo "${PLATFORM},${FLOWCELL},${RUNDATE},${BC_MODEL},${MOD_LABELS},${FILETYPE}"
