#!/usr/bin/env bash
#
# Run bin/run under a wall-clock cap, restarting from the newest checkpoint
# after each segment, so a long simulation survives being cut off by a time
# limit or by an external thermal guard (SIGSTOP/SIGCONT from
# scripts/watch_temperature.ps1) without needing any change to the
# simulation code itself.
#
# Once the cap is reached the segment is not killed immediately: bin/run has
# no SIGTERM handler, so cutting it off mid-step would waste all progress
# since the last checkpoint. Instead we let it keep producing checkpoints
# normally, and only once the cap has passed do we wait for whichever
# checkpoint is actually in flight *at that moment* to land, then send
# SIGTERM -- so a segment can run past its nominal cap but never throws
# away completed work, and never wastes time re-checking a checkpoint that
# already finished long before the cap was reached.
#
# Usage:
#   scripts/run_supervised.sh input_folder/beta3_frg1
#
# To ease the load on the machine over a long unattended run, once it has
# been running continuously for CONTINUOUS_HOURS_LIMIT (default: 2 days) it
# takes one long LONG_BREAK_HOURS (default: 3h) break instead of the normal
# cooldown, then resets the continuous-runtime clock.
#
# Env overrides:
#   SEGMENT_HOURS          wall-clock cap per segment (default: 6)
#   COOLDOWN_MIN           pause between segments (default: 10)
#   CONTINUOUS_HOURS_LIMIT continuous runtime before a long break (default: 48)
#   LONG_BREAK_HOURS       long break duration once that limit is hit (default: 3)
#   OMP_NUM_THREADS thread count passed to bin/run (default: 4)
#   CPU_LIST        comma-separated logical CPU list to pin bin/run to via
#                    taskset (default: per-input table in this script, else
#                    no pinning). Set this to
#                    disjoint physical cores when running multiple inputs
#                    concurrently -- OMP_PROC_BIND=TRUE + OMP_PLACES=cores
#                    always picks the first N cores in topology order with
#                    no coordination between processes, so two concurrent
#                    runs otherwise collide on the same physical cores and
#                    fight over the same execution units via hyperthreading.

set -uo pipefail

INPUT="${1:?usage: $0 <input_file> (e.g. input_folder/beta3)}"
SEGMENT_HOURS="${SEGMENT_HOURS:-6}"
COOLDOWN_MIN="${COOLDOWN_MIN:-10}"
CONTINUOUS_HOURS_LIMIT="${CONTINUOUS_HOURS_LIMIT:-48}"
LONG_BREAK_HOURS="${LONG_BREAK_HOURS:-3}"
CPU_LIST="${CPU_LIST:-}"

REPO_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")/.." && pwd)"
cd "$REPO_DIR"

RUN_NAME="$(basename "$INPUT")"

# Per-input thread count default (an explicit OMP_NUM_THREADS overrides this,
# same convention as CPU_LIST below) -- so a restart doesn't need the caller
# to remember to pass OMP_NUM_THREADS=8 for beta3 every time.
if [[ -z "${OMP_NUM_THREADS:-}" ]]; then
    case "$RUN_NAME" in
        test1e16) RUN_THREADS=8 ;;
        *) RUN_THREADS=4 ;;
    esac
else
    RUN_THREADS="$OMP_NUM_THREADS"
fi

# Per-input CPU pinning so concurrent runs never share physical cores.
# Logical CPUs 2k and 2k+1 are SMT siblings of P-core k (0-7); we take one
# logical CPU per physical core. An explicit CPU_LIST env var overrides this.
# To add another concurrent run, give it unused physical cores here
# (e.g. E-cores are logical CPUs 16-23).
if [[ -z "$CPU_LIST" ]]; then
    case "$RUN_NAME" in
        test1e15) CPU_LIST="0,2,4,6" ;;
        # 4 P-cores (own) + 4 E-cores (17,19,21,23 -- 16,18 stay reserved for
        # beta4_St01_new) = 8 physical cores total. P/E cores run at
        # different speeds, but the OMP loops here use dynamic scheduling.
        test1e18) CPU_LIST="8,10,12,14" ;;
        test1e17) CPU_LIST="17,19,21,23" ;;
        *) echo "[run_supervised] warning: no CPU_LIST for '$RUN_NAME' (add it to the table in this script); running unpinned" >&2 ;;
    esac
fi

LOGDIR="log"
mkdir -p "$LOGDIR"

# outdir is declared in the input file as: outdir = /path/to/dir
OUTDIR="$(grep -E '^[[:space:]]*outdir[[:space:]]*=' "$INPUT" | head -1 | sed -E 's/^[^=]*=[[:space:]]*([^ /]*\/[^[:space:]]*).*/\1/')"
if [[ -z "$OUTDIR" ]]; then
    echo "[run_supervised] could not parse outdir from $INPUT" >&2
    exit 1
fi
MAGDIR="${OUTDIR}_B"

STOP_FILE="/tmp/${RUN_NAME}.stop"
PID_FILE="/tmp/${RUN_NAME}.pid"

echo "[run_supervised] $(date '+%F %T'): input=$INPUT outdir=$OUTDIR magdir=$MAGDIR"
echo "[run_supervised] segment cap=${SEGMENT_HOURS}h cooldown=${COOLDOWN_MIN}m threads=${RUN_THREADS} cpu_list=${CPU_LIST:-<none>}"
echo "[run_supervised] long break: ${LONG_BREAK_HOURS}h after ${CONTINUOUS_HOURS_LIMIT}h of continuous running"
echo "[run_supervised] stop file: $STOP_FILE (create it to end the loop after the current segment)"

CONTINUOUS_SINCE=$(date +%s)

latest_year() {
    find "$OUTDIR" -maxdepth 1 -regextype posix-extended -regex '.*/out[0-9]+' 2>/dev/null \
        | sed -E 's#.*/out##' \
        | sort -n | tail -1
}

# Runs one segment in the background, polls it, and enforces the wall-clock
# cap and stop-file by signalling the *actual* bin/run PID (recorded in
# PID_FILE) rather than relying on `timeout`, so an external SIGSTOP/SIGCONT
# aimed at that PID (from the Windows-side thermal guard) works transparently
# -- the deadline clock just keeps ticking in real time while it's paused.
run_segment() {
    if [[ -n "$CPU_LIST" ]]; then
        taskset -c "$CPU_LIST" "$@" &> "$LOG" &
    else
        "$@" &> "$LOG" &
    fi
    local pid=$!
    echo "$pid" > "$PID_FILE"

    local last_year
    last_year="$(latest_year || true)"
    local deadline=$(( $(date +%s) + SEGMENT_HOURS * 3600 ))
    local deadline_hit=0
    while kill -0 "$pid" 2>/dev/null; do
        if [[ -e "$STOP_FILE" ]]; then
            echo "[run_supervised] $(date '+%F %T'): stop file seen, sending SIGTERM"
            kill -TERM "$pid" 2>/dev/null
            wait "$pid" 2>/dev/null
            rm -f "$PID_FILE" "$STOP_FILE"
            return 200
        fi
        if (( ! deadline_hit )) && (( $(date +%s) >= deadline )); then
            deadline_hit=1
            echo "[run_supervised] $(date '+%F %T'): segment cap reached, waiting for the checkpoint currently being computed"
        fi
        # Track the newest checkpoint continuously (not just after the cap),
        # so that once the cap is hit we wait for whichever checkpoint is
        # actually in flight at that moment -- not one that already finished
        # long before the cap, which would otherwise let the run keep going
        # uselessly toward the *next* one only to be killed before it lands.
        local cur_year
        cur_year="$(latest_year || true)"
        if [[ -n "$cur_year" && "$cur_year" != "$last_year" ]]; then
            last_year="$cur_year"
            if (( deadline_hit )); then
                echo "[run_supervised] $(date '+%F %T'): checkpoint out${cur_year} written, sending SIGTERM"
                sleep 5
                kill -TERM "$pid" 2>/dev/null
                sleep 30
                kill -0 "$pid" 2>/dev/null && kill -KILL "$pid" 2>/dev/null
                wait "$pid" 2>/dev/null
                rm -f "$PID_FILE"
                return 124
            fi
        fi
        sleep 30
    done
    wait "$pid"
    local status=$?
    rm -f "$PID_FILE"
    return "$status"
}

while true; do
    if [[ -e "$STOP_FILE" ]]; then
        echo "[run_supervised] $(date '+%F %T'): stop file found before segment start, exiting"
        rm -f "$STOP_FILE"
        exit 0
    fi

    YEAR="$(latest_year || true)"
    TS="$(date '+%Y%m%d_%H%M%S')"
    LOG="$LOGDIR/${RUN_NAME}_${TS}.log"

    export OMP_PROC_BIND=TRUE
    export OMP_PLACES=cores
    export OMP_NUM_THREADS="$RUN_THREADS"

    if [[ -z "$YEAR" ]]; then
        echo "[run_supervised] $(date '+%F %T'): no checkpoint found, starting fresh (log: $LOG)"
        run_segment ./bin/run "$INPUT"
    else
        DUST="$OUTDIR/out${YEAR}"
        GAS="$OUTDIR/out${YEAR}_gas"
        MAG="$MAGDIR/out${YEAR}_B"
        echo "[run_supervised] $(date '+%F %T'): restarting from year $YEAR (log: $LOG)"
        run_segment ./bin/run -r "$INPUT" "$DUST" "$GAS" "$MAG" "$YEAR"
    fi
    STATUS=$?

    case "$STATUS" in
        124)
            elapsed_continuous=$(( $(date +%s) - CONTINUOUS_SINCE ))
            if (( elapsed_continuous >= CONTINUOUS_HOURS_LIMIT * 3600 )); then
                echo "[run_supervised] $(date '+%F %T'): ${CONTINUOUS_HOURS_LIMIT}h of continuous running reached, taking a ${LONG_BREAK_HOURS}h break"
                sleep "$(( LONG_BREAK_HOURS * 3600 ))"
                CONTINUOUS_SINCE=$(date +%s)
            else
                echo "[run_supervised] $(date '+%F %T'): segment timed out, cooling down ${COOLDOWN_MIN}m"
                sleep "$(( COOLDOWN_MIN * 60 ))"
            fi
            ;;
        200)
            echo "[run_supervised] $(date '+%F %T'): stopped by request"
            exit 0
            ;;
        0)
            echo "[run_supervised] $(date '+%F %T'): run finished normally"
            exit 0
            ;;
        *)
            echo "[run_supervised] $(date '+%F %T'): segment exited with status $STATUS, see $LOG -- stopping (fix and re-run to resume from the last checkpoint)"
            exit "$STATUS"
            ;;
    esac
done
