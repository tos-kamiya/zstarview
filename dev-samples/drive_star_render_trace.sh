#!/usr/bin/env bash
# Drive a real zstarview window through idle / hover / rotate phases while
# trace_star_render_cache.py prints per-interval star cache counters.
#
# Requires an X11 or XWayland display and xdotool.  The config directory is
# isolated so the user's ~/.config/zstarview is not modified.
#
# Usage:
#   dev-samples/drive_star_render_trace.sh OUT_DIR [extra zstarview args...]
# Phase start times are written to OUT_DIR/phases.log and the trace to
# OUT_DIR/trace.log.
set -euo pipefail

out_dir=${1:?usage: $0 OUT_DIR [zstarview args...]}
shift
mkdir -p "$out_dir/xdg-config"
startup_wait=${STARTUP_WAIT:-75}
phase() { echo "$(date +%H:%M:%S) $*" | tee -a "$out_dir/phases.log"; }
: > "$out_dir/phases.log"

XDG_CONFIG_HOME="$out_dir/xdg-config" QT_QPA_PLATFORM=xcb \
  uv run -p .venv/bin/python dev-samples/trace_star_render_cache.py \
  --trace-interval 5 --window-frame window --window-geometry 100,100,1600,900 \
  "$@" > "$out_dir/trace.log" 2>&1 &
app_pid=$!
trap 'kill $app_pid 2>/dev/null || true' EXIT

phase "startup (waiting ${startup_wait}s)"
sleep "$startup_wait"
win=$(xdotool search --sync --onlyvisible --pid "$(pgrep -f trace_star_render_cache.py | tail -1)" | tail -1)
geom=$(xdotool getwindowgeometry --shell "$win")
eval "$geom"
cx=$((X + WIDTH / 2)); cy=$((Y + HEIGHT / 2))
xdotool mousemove "$cx" "$cy"
sleep 5

phase "idle"
sleep 20

phase "hover"
for i in $(seq 0 79); do
  xdotool mousemove $((X + 200 + i * 15)) $((Y + 150 + (i % 20) * 25))
  sleep 0.15
done
sleep 3

phase "rotate"
xdotool windowfocus "$win" 2>/dev/null || true
for _ in $(seq 1 30); do xdotool key --window "$win" Right; sleep 0.1; done
for _ in $(seq 1 15); do xdotool key --window "$win" Up; sleep 0.1; done

phase "settle"
sleep 20

phase "end"
kill -INT "$app_pid" 2>/dev/null || true
sleep 5
