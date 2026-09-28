#!/usr/bin/env bash
set -euo pipefail
export LC_ALL=C

# ============================================================
# SARS-CoV-2 collection benchmark
#
# Compares:
#   1. JARVIS3 default configuration (level 7)
#   2. Optional manually/expert-tuned JARVIS3 configuration
#   3. OptJV3 automatically optimized configuration
#
# Expected in current directory:
#   coronas.fasta
#   JARVIS3
#   OptJV3
#
# Run:
#   chmod +x run.sh
#   ./run.sh
#
# Stronger optimization:
#   POPULATION=100 GENERATIONS=100 THREADS=8 ./run.sh
#
# Optional manual baseline:
#   MANUAL_PARAMS='-lr ... -cm ... -rm ...' ./run.sh
# ============================================================


# ------------------------------------------------------------
# Configuration
# ------------------------------------------------------------

FASTA="${FASTA:-coronas.fasta}"
JARVIS="${JARVIS:-./JARVIS3}"
OPT="${OPT:-./OptJV3}"

THREADS="${THREADS:-8}"
POPULATION="${POPULATION:-32}"
GENERATIONS="${GENERATIONS:-10}"
SEED="${SEED:-12345}"

# Optimize using this many MB from the collection.
SAMPLE_MB="${SAMPLE_MB:-10}"

OUTDIR="${OUTDIR:-benchmark_coronas}"

# Optional expert/manual JARVIS3 parameters.
MANUAL_PARAMS="${MANUAL_PARAMS:-}"

mkdir -p "$OUTDIR"

SEQ="$OUTDIR/coronas.seq"
SAMPLE="$OUTDIR/coronas_${SAMPLE_MB}MB.seq"

RESULTS="$OUTDIR/results.tsv"
FILTER_REPORT="$OUTDIR/filter_report.tsv"

BEST="$OUTDIR/best.txt"
HISTORY="$OUTDIR/history.csv"
OPTLOG="$OUTDIR/optjv3.log"
WORK="$OUTDIR/work"


# ------------------------------------------------------------
# Checks
# ------------------------------------------------------------

[[ -f "$FASTA" ]] || {
    echo "ERROR: $FASTA not found"
    exit 1
}

[[ -x "$JARVIS" ]] || {
    echo "ERROR: $JARVIS not found or not executable"
    exit 1
}

[[ -x "$OPT" ]] || {
    echo "ERROR: $OPT not found or not executable"
    exit 1
}


# ------------------------------------------------------------
# Prepare genomic sequence collection
#
# Important:
#   - FASTA headers are removed.
#   - Sequence is converted to uppercase.
#   - ALL non-ACGT symbols are removed:
#       N, R, Y, K, M, S, W, B, D, H, V, gaps, etc.
# ------------------------------------------------------------

echo "============================================================"
echo "Preparing SARS-CoV-2 genomic collection"
echo "============================================================"

N_GENOMES=$(grep -c '^>' "$FASTA" || true)

# Count all sequence symbols before filtering
# (excluding headers and whitespace).
RAW_SYMBOLS=$(
    awk '
        /^>/ { next }
        {
            line = $0
            gsub(/[[:space:]]/, "", line)
            n += length(line)
        }
        END { print n+0 }
    ' "$FASTA"
)

# Create raw ACGT-only sequence.
awk '
    /^>/ { next }
    {
        line = toupper($0)
        gsub(/[^ACGT]/, "", line)
        printf "%s", line
    }
' "$FASTA" > "$SEQ"

INPUT_BYTES=$(wc -c < "$SEQ")
REMOVED=$(( RAW_SYMBOLS - INPUT_BYTES ))

if (( INPUT_BYTES == 0 )); then
    echo "ERROR: no ACGT sequence data found after filtering"
    exit 1
fi

REMOVED_PERCENT=$(
    awk \
        -v r="$REMOVED" \
        -v n="$RAW_SYMBOLS" \
        'BEGIN {
            if (n > 0)
                printf "%.6f", 100*r/n;
            else
                printf "0";
        }'
)

echo "FASTA:              $FASTA"
echo "Genomes:            $N_GENOMES"
echo "Symbols before QC:  $RAW_SYMBOLS"
echo "ACGT symbols:       $INPUT_BYTES"
echo "Removed non-ACGT:   $REMOVED ($REMOVED_PERCENT %)"

awk -v n="$INPUT_BYTES" \
    'BEGIN {
        printf "Filtered size:      %.2f MB\n", n/1000000
    }'

# Save QC information.
printf \
    "genomes\traw_symbols\tacgt_symbols\tremoved\tremoved_percent\n" \
    > "$FILTER_REPORT"

printf \
    "%s\t%s\t%s\t%s\t%s\n" \
    "$N_GENOMES" \
    "$RAW_SYMBOLS" \
    "$INPUT_BYTES" \
    "$REMOVED" \
    "$REMOVED_PERCENT" \
    >> "$FILTER_REPORT"


# ------------------------------------------------------------
# Sanity check: guarantee that only ACGT remains
# ------------------------------------------------------------

if grep -q '[^ACGT]' "$SEQ"; then
    echo "ERROR: non-ACGT symbols remain after filtering"
    exit 1
fi

echo "ACGT validation:    OK"


# ------------------------------------------------------------
# Representative sample
#
# Take SAMPLE_MB from approximately the middle of the
# concatenated genomic collection.
# ------------------------------------------------------------

SAMPLE_BYTES=$(( SAMPLE_MB * 1000 * 1000 ))

if (( SAMPLE_BYTES >= INPUT_BYTES )); then

    cp "$SEQ" "$SAMPLE"

else

    START_BYTE=$(( (INPUT_BYTES - SAMPLE_BYTES) / 2 ))

    dd \
        if="$SEQ" \
        of="$SAMPLE" \
        iflag=skip_bytes,count_bytes \
        skip="$START_BYTE" \
        count="$SAMPLE_BYTES" \
        status=none
fi

ACTUAL_SAMPLE_BYTES=$(wc -c < "$SAMPLE")

awk -v n="$ACTUAL_SAMPLE_BYTES" \
    'BEGIN {
        printf "Opt. sample:        %.2f MB\n", n/1000000
    }'

echo


# ------------------------------------------------------------
# Helper functions
# ------------------------------------------------------------

calc_bps()
{
    local COMPRESSED_BYTES="$1"

    awk \
        -v c="$COMPRESSED_BYTES" \
        -v n="$INPUT_BYTES" \
        'BEGIN {
            printf "%.6f", (8.0*c)/n
        }'
}


calc_gain()
{
    local BASE="$1"
    local TEST="$2"

    awk \
        -v b="$BASE" \
        -v t="$TEST" \
        'BEGIN {
            if (b > 0)
                printf "%.3f", 100.0*(b-t)/b;
            else
                printf "NA";
        }'
}


run_jarvis()
{
    local NAME="$1"
    shift

    local OUT="$OUTDIR/${NAME}.jc"
    local LOG="$OUTDIR/${NAME}.log"

    rm -f "$OUT" "$LOG"

    echo
    echo "------------------------------------------------------------"
    echo "$NAME"
    echo "------------------------------------------------------------"

    echo "Command:"
    printf "%q " "$JARVIS" -v "$@" -o "$OUT" "$SEQ"
    echo
    echo

    local START_TIME
    local END_TIME
    local ELAPSED
    local SIZE
    local BPS

    START_TIME=$(date +%s.%N)

    "$JARVIS" \
        -v \
        "$@" \
        -o "$OUT" \
        "$SEQ" \
        2>&1 | tee "$LOG"

    END_TIME=$(date +%s.%N)

    if [[ ! -f "$OUT" ]]; then
        echo
        echo "ERROR: JARVIS3 did not produce:"
        echo "  $OUT"
        echo
        echo "See:"
        echo "  $LOG"
        exit 1
    fi

    ELAPSED=$(
        awk \
            -v a="$START_TIME" \
            -v b="$END_TIME" \
            'BEGIN {
                printf "%.3f", b-a
            }'
    )

    SIZE=$(stat -c '%s' "$OUT")
    BPS=$(calc_bps "$SIZE")

    printf \
        "%s\t%s\t%s\t%s\n" \
        "$NAME" \
        "$SIZE" \
        "$BPS" \
        "$ELAPSED" \
        >> "$RESULTS"

    echo
    echo "Result:"
    echo "  compressed bytes = $SIZE"
    echo "  BPS              = $BPS"
    echo "  time             = $ELAPSED s"

    rm -f "$OUT"
}


# ------------------------------------------------------------
# Results table
# ------------------------------------------------------------

printf \
    "method\tcompressed_bytes\tbps\ttime_s\n" \
    > "$RESULTS"


# ============================================================
# 1. JARVIS3 DEFAULT
# ============================================================

echo
echo "============================================================"
echo "1. JARVIS3 default configuration"
echo "============================================================"
echo
echo "JARVIS3 default = level 7"

run_jarvis \
    "JARVIS3_default_L7" \
    -l 7

DEFAULT_BPS=$(
    awk -F'\t' \
        '$1=="JARVIS3_default_L7" {print $3}' \
        "$RESULTS"
)


# ============================================================
# 2. OPTIONAL MANUAL / EXPERT CONFIGURATION
# ============================================================

MANUAL_BPS=""

if [[ -n "$MANUAL_PARAMS" ]]; then

    echo
    echo "============================================================"
    echo "2. Expert/manual JARVIS3 configuration"
    echo "============================================================"

    echo
    echo "Parameters:"
    echo "$MANUAL_PARAMS"

    read -r -a MANUAL_ARRAY <<< "$MANUAL_PARAMS"

    run_jarvis \
        "JARVIS3_manual" \
        "${MANUAL_ARRAY[@]}"

    MANUAL_BPS=$(
        awk -F'\t' \
            '$1=="JARVIS3_manual" {print $3}' \
            "$RESULTS"
    )

else

    echo
    echo "============================================================"
    echo "2. Expert/manual JARVIS3 configuration"
    echo "============================================================"
    echo
    echo "Skipped: MANUAL_PARAMS was not supplied."

fi


# ============================================================
# 3. OptJV3 optimization
# ============================================================

echo
echo "============================================================"
echo "3. OptJV3 optimization"
echo "============================================================"

echo
echo "Input sample:  $SAMPLE"
echo "Sample bases:  $ACTUAL_SAMPLE_BYTES"
echo "Population:    $POPULATION"
echo "Generations:   $GENERATIONS"
echo "Threads:       $THREADS"
echo "Seed:          $SEED"
echo

rm -rf "$WORK"
mkdir -p "$WORK"

rm -f \
    "$BEST" \
    "$HISTORY" \
    "$OPTLOG"

echo "Starting OptJV3..."
echo "Progress will be reported after each completed generation."
echo


# ------------------------------------------------------------
# Run OptJV3
#
# Output is line-buffered so that progress appears immediately.
# Whenever OptJV3 reports:
#
#   => Generation X of Y ...
#
# an additional progress line is printed:
#
#   [PROGRESS] Generation X / Y (Z%)
#
# This does NOT modify the optimization algorithm.
# ------------------------------------------------------------

stdbuf -oL -eL "$OPT" \
    --jarvis "$JARVIS" \
    --input "$SAMPLE" \
    --population "$POPULATION" \
    --generations "$GENERATIONS" \
    --threads "$THREADS" \
    --seed "$SEED" \
    --elite 4 \
    --tournament 3 \
    --crossover 0.85 \
    --mutation 0.12 \
    --toggle 0.08 \
    --objective bps \
    --max-cmodels 4 \
    --min-cmodels 1 \
    --max-rmodels 2 \
    --min-rmodels 1 \
    --global-bounds \
        "hs=8:128,lr=0.0:0.15,seed=1:1000" \
    --cm-bounds \
        "ctx=1:14,den=1:5000,ir=0:2,gamma=0.01:0.99,edits=0:20,eden=1:5000,eir=0:1,egamma=0.01:0.99" \
    --rm-bounds \
        "nr=1:500,ctx=1:14,beta=0.01:0.99,limit=1:20,gamma=0.01:0.99,ir=0:2,weight=0.01:0.99,cache=1:12" \
    --workdir "$WORK" \
    --best-out "$BEST" \
    --history-out "$HISTORY" \
    2>&1 |
awk -v total="$GENERATIONS" '
{
    # Preserve the original OptJV3 output.
    print
    fflush()

    # Expected OptJV3 format:
    # => Generation 5 of 100 best: ...
    if ($2 == "Generation" && $4 == "of") {

        gen = $3 + 0
        tot = $5 + 0

        if (tot <= 0)
            tot = total

        pct = 100.0 * gen / tot

        printf "\n"
        printf "[PROGRESS] Generation %d / %d (%.1f%%)\n", \
               gen, tot, pct
        printf "\n"

        fflush()
    }
}
' | tee "$OPTLOG"


# ------------------------------------------------------------
# Extract optimized parameters
# ------------------------------------------------------------

if [[ ! -s "$BEST" ]]; then
    echo
    echo "ERROR: OptJV3 did not create:"
    echo "  $BEST"
    exit 1
fi

OPT_PARAMS=$(
    sed -n \
        's/^params:[[:space:]]*//p' \
        "$BEST" |
    head -n 1
)

if [[ -z "$OPT_PARAMS" ]]; then
    echo
    echo "ERROR: could not extract params from:"
    echo "  $BEST"
    echo
    echo "Content:"
    cat "$BEST"
    exit 1
fi


echo
echo "============================================================"
echo "Best parameters found by OptJV3"
echo "============================================================"
echo
echo "$OPT_PARAMS"


# ============================================================
# 4. Apply optimized configuration to FULL collection
# ============================================================

echo
echo "============================================================"
echo "4. Testing optimized configuration on full collection"
echo "============================================================"

read -r -a OPT_ARRAY <<< "$OPT_PARAMS"

run_jarvis \
    "OptJV3_optimized" \
    "${OPT_ARRAY[@]}"

OPT_BPS=$(
    awk -F'\t' \
        '$1=="OptJV3_optimized" {print $3}' \
        "$RESULTS"
)


# ============================================================
# FINAL SUMMARY
# ============================================================

GAIN_DEFAULT=$(calc_gain "$DEFAULT_BPS" "$OPT_BPS")

echo
echo
echo "============================================================"
echo "FINAL SUMMARY"
echo "============================================================"
echo

printf "%-34s %s\n" \
    "Collection:" \
    "$FASTA"

printf "%-34s %s\n" \
    "Number of genomes:" \
    "$N_GENOMES"

printf "%-34s %s\n" \
    "Original sequence symbols:" \
    "$RAW_SYMBOLS"

printf "%-34s %s\n" \
    "ACGT symbols:" \
    "$INPUT_BYTES"

printf "%-34s %s (%s%%)\n" \
    "Removed non-ACGT symbols:" \
    "$REMOVED" \
    "$REMOVED_PERCENT"

printf "%-34s %.2f MB\n" \
    "Filtered collection size:" \
    "$(awk -v n="$INPUT_BYTES" 'BEGIN {print n/1000000}')"

printf "%-34s %.2f MB\n" \
    "Optimization sample:" \
    "$(awk -v n="$ACTUAL_SAMPLE_BYTES" 'BEGIN {print n/1000000}')"

echo

printf "%-34s %s BPS\n" \
    "JARVIS3 default (L7):" \
    "$DEFAULT_BPS"

printf "%-34s %s BPS\n" \
    "OptJV3:" \
    "$OPT_BPS"

printf "%-34s %s %%\n" \
    "Improvement vs default:" \
    "$GAIN_DEFAULT"


if [[ -n "$MANUAL_BPS" ]]; then

    GAIN_MANUAL=$(calc_gain "$MANUAL_BPS" "$OPT_BPS")

    printf "%-34s %s BPS\n" \
        "Manual/expert JARVIS3:" \
        "$MANUAL_BPS"

    printf "%-34s %s %%\n" \
        "Improvement vs manual:" \
        "$GAIN_MANUAL"
fi


echo
echo "Best OptJV3 parameters:"
echo "$OPT_PARAMS"

echo
echo "Output files:"
echo "  $RESULTS"
echo "  $FILTER_REPORT"
echo "  $BEST"
echo "  $HISTORY"
echo "  $OPTLOG"

echo
echo "Done."
