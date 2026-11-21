#!/usr/bin/env bash
#
# Run scripts/run_supervised.sh for a given input, but only during a daily
# wall-clock window (default 09:00-19:00). Reuses run_supervised.sh
# unmodified for the checkpoint-detection / segment-cap / cooldown /
# restart logic; this wrapper only adds the daily start/stop schedule on
# top of it.
#
# At the window's end (WINDOW_END), if a checkpoint happens to be actively
# being written at that moment, this waits up to GRACE_HOURS for it to land
# (mirroring run_supervised.sh's own SEGMENT_HOURS-deadline behaviour) and
# cuts the run as soon as it does -- but if nothing new lands within
# GRACE_HOURS, it is cut anyway at that point regardless of what's in
# flight. Cutting is done via run_supervised.sh's own STOP_FILE mechanism,
# which sends SIGTERM immediately (does NOT wait further) -- so once the
# grace period (if any) is over, whatever progress has been made since the
# last completed checkpoint that day is lost. In the rare case SIGTERM lands
# exactly while a checkpoint file is mid-write, that file can end up
# truncated/corrupt; delete it and the next-older out<N> is still a valid
# restart point.
#
# Usage:
#   scripts/run_daytime.sh input_folder/beta4_St01_new
#
# Env overrides:
#   WINDOW_START    daily start time, HH:MM (default 09:00)
#   WINDOW_END      daily end time, HH:MM (default 19:00)
#   GRACE_HOURS     max wait past WINDOW_END for an in-flight checkpoint to
#                   land before cutting anyway (default 1)
#   OMP_NUM_THREADS thread count passed through to run_supervised.sh/bin/run
#                   (default 2 here -- run_supervised.sh's own default is 4)
#   CPU_LIST        comma-separated logical CPU list, passed through to
#                   run_supervised.sh (see its own header comment; pick
#                   physical cores disjoint from any other concurrent run)
#   (any other scripts/run_supervised.sh env override -- SEGMENT_HOURS,
#    COOLDOWN_MIN, CONTINUOUS_HOURS_LIMIT, LONG_BREAK_HOURS -- is also
#    honored, exported through to it)

set -uo pipefail

INPUT="${1:?usage: $0 <input_file> (e.g. input_folder/beta4_St01_new)}"
WINDOW_START="${WINDOW_START:-09:00}"
WINDOW_END="${WINDOW_END:-19:00}"
GRACE_HOURS="${GRACE_HOURS:-1}"
export OMP_NUM_THREADS="${OMP_NUM_THREADS:-2}"

REPO_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")/.." && pwd)"
cd "$REPO_DIR"

RUN_NAME="$(basename "$INPUT")"
STOP_FILE="/tmp/${RUN_NAME}.stop"

# outdir is declared in the input file as: outdir = /path/to/dir
OUTDIR="$(grep -E '^[[:space:]]*outdir[[:space:]]*=' "$INPUT" | head -1 | sed -E 's/^[^=]*=[[:space:]]*([^ /]*\/[^[:space:]]*).*/\1/')"
if [[ -z "$OUTDIR" ]]; then
    echo "[run_daytime] could not parse outdir from $INPUT" >&2
    exit 1
fi

latest_year() {
    find "$OUTDIR" -maxdepth 1 -regextype posix-extended -regex '.*/out[0-9]+' 2>/dev/null \
        | sed -E 's#.*/out##' \
        | sort -n | tail -1
}

SUPERVISOR_PID=""

# Make sure a Ctrl+C / kill on *this* wrapper also stops the underlying run
# cleanly instead of leaving it running detached.
cleanup() {
    if [[ -n "$SUPERVISOR_PID" ]] && kill -0 "$SUPERVISOR_PID" 2>/dev/null; then
        echo "[run_daytime] $(date '+%F %T'): wrapper stopping, cutting the current run"
        touch "$STOP_FILE"
        wait "$SUPERVISOR_PID" 2>/dev/null
    fi
}
trap cleanup EXIT INT TERM

# Next occurrence (today if still in the future, else tomorrow) of HH:MM, as
# a unix timestamp.
next_time() {
    local hhmm="$1"
    local today
    today="$(date -d "today ${hhmm}" +%s)"
    if (( today > $(date +%s) )); then
        echo "$today"
    else
        date -d "tomorrow ${hhmm}" +%s
    fi
}

echo "[run_daytime] $(date '+%F %T'): input=$INPUT outdir=$OUTDIR window=${WINDOW_START}-${WINDOW_END} grace=${GRACE_HOURS}h threads=${OMP_NUM_THREADS}"

while true; do
    START_TS="$(date -d "today ${WINDOW_START}" +%s)"
    END_TS="$(date -d "today ${WINDOW_END}" +%s)"
    NOW_TS="$(date +%s)"

    if (( NOW_TS < START_TS )); then
        echo "[run_daytime] $(date '+%F %T'): before today's window, sleeping until ${WINDOW_START}"
        sleep "$(( START_TS - NOW_TS ))"
        continue
    fi
    if (( NOW_TS >= END_TS )); then
        NEXT_START="$(next_time "$WINDOW_START")"
        echo "[run_daytime] $(date '+%F %T'): after today's window, sleeping until $(date -d @"$NEXT_START" '+%F %T')"
        sleep "$(( NEXT_START - NOW_TS ))"
        continue
    fi

    rm -f "$STOP_FILE"
    echo "[run_daytime] $(date '+%F %T'): starting run_supervised.sh (window ends $(date -d "today ${WINDOW_END}" '+%F %T'))"
    bash scripts/run_supervised.sh "$INPUT" &
    SUPERVISOR_PID=$!

    REMAINING=$(( END_TS - $(date +%s) ))
    (( REMAINING > 0 )) && sleep "$REMAINING"

    if kill -0 "$SUPERVISOR_PID" 2>/dev/null; then
        # Window end reached and the run is still going. If a checkpoint is
        # actively being written right now, give it up to GRACE_HOURS to
        # land before cutting anyway (mirrors run_supervised.sh's own
        # SEGMENT_HOURS-deadline behaviour, just bounded); cut as soon as it
        # lands rather than always waiting the full grace period.
        LAST_YEAR="$(latest_year || true)"
        GRACE_DEADLINE=$(( $(date +%s) + GRACE_HOURS * 3600 ))
        echo "[run_daytime] $(date '+%F %T'): window end, waiting up to ${GRACE_HOURS}h for the in-flight checkpoint (currently out${LAST_YEAR:-none})"
        while kill -0 "$SUPERVISOR_PID" 2>/dev/null; do
            CUR_YEAR="$(latest_year || true)"
            if [[ -n "$CUR_YEAR" && "$CUR_YEAR" != "$LAST_YEAR" ]]; then
                echo "[run_daytime] $(date '+%F %T'): checkpoint out${CUR_YEAR} landed, cutting now"
                break
            fi
            if (( $(date +%s) >= GRACE_DEADLINE )); then
                echo "[run_daytime] $(date '+%F %T'): grace period exhausted, cutting now (mid-step progress lost)"
                break
            fi
            sleep 15
        done
        if kill -0 "$SUPERVISOR_PID" 2>/dev/null; then
            touch "$STOP_FILE"
            wait "$SUPERVISOR_PID" 2>/dev/null
        fi
    else
        echo "[run_daytime] $(date '+%F %T'): run_supervised.sh already exited before window end"
    fi
    SUPERVISOR_PID=""

    NEXT_START="$(next_time "$WINDOW_START")"
    echo "[run_daytime] $(date '+%F %T'): sleeping until $(date -d @"$NEXT_START" '+%F %T')"
    sleep "$(( NEXT_START - $(date +%s) ))"
done
