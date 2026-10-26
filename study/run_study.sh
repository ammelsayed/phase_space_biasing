#!/bin/bash
# Run the phase-space biasing study: reference vs bias modules for V->leptons.
# Generates one MG5 script variant per (sample, module) in runs/ (gitignored),
# runs it with /data/ammelsayed/runmg5.sh and copies the unweighted LHE aside
# for the analyzer (the next run of the same sample overwrites the process dir).
#
# Module matrix (ptj_bias/ht_bias/ptl1_met_ht_bias are not run: the study uses
# V->leptons only, no final-state jets; Z has no MET so no MET-containing
# module is run on Z):
#   Z: reference, ptl1_bias, xsec_flat_bias
#   W: reference, ptl1_bias, met_bias, ptl1_met_bias, xsec_flat_bias
#
# Usage: study/run_study.sh     (run from anywhere; paths are resolved from it)

set -u
ROOT="$(cd "$(dirname "$0")/.." && pwd)"
RUNMG5=/data/ammelsayed/runmg5.sh
RUNS="$ROOT/runs"
mkdir -p "$RUNS"

run_one () {  # $1 = base script, $2 = module (or None), $3 = tag
  local base=$1 mod=$2 tag=$3
  local proc=${base%.mg5}
  local script="$RUNS/${tag}.mg5"
  cp "$ROOT/$base" "$script"
  if [ "$mod" != "None" ]; then
    sed -i "s|^set ebeam 6800|set bias_module $ROOT/mg5_bias_modules/$mod\nset ebeam 6800|" "$script"
  fi
  echo "=== [$(date +%H:%M:%S)] launching $tag (module: $mod) ==="
  bash "$RUNMG5" "$script"
  local pid t=0
  pid=$(cat "$RUNS/LogFiles/${tag}.pid")
  while kill -0 "$pid" 2>/dev/null; do
    sleep 30; t=$((t+30))
    if [ "$t" -gt 10800 ]; then echo "TIMEOUT waiting for $tag"; return 1; fi
  done
  echo "=== [$(date +%H:%M:%S)] $tag finished (~$((t/60)) min) ==="
  # keep the cross section and the module defaults actually used
  grep -A2 "Cross-section" "$RUNS/LogFiles/${tag}.log" | tail -3 \
      > "$RUNS/xsec_${tag}.txt" || true
  cat "$RUNS/$proc/Source/BIAS/bias.inc" > "$RUNS/biasinc_${tag}.txt" 2>/dev/null || true
  # save this run's events and MA5 report before the next run of the
  # same sample wipes the process directory
  local ev
  for ev in "$RUNS/$proc/Events"/run_*/unweighted_events.lhe*; do
    [ -e "$ev" ] && cp "$ev" "$RUNS/events_${tag}.lhe${ev##*.lhe}"
  done
  local ma5
  for ma5 in "$RUNS/$proc/Events"/run_*/tag_1_MA5_parton_analysis_*.pdf; do
    [ -e "$ma5" ] && cp "$ma5" "$RUNS/ma5_${tag}.pdf"
  done
  tail -3 "$RUNS/LogFiles/${tag}.log"
}

run_one z012j_MLM.mg5 None           z012j_ref
run_one z012j_MLM.mg5 ptl1_bias      z012j_ptl1
run_one z012j_MLM.mg5 xsec_flat_bias z012j_xsecflat
run_one w012j_MLM.mg5 None           w012j_ref
run_one w012j_MLM.mg5 ptl1_bias      w012j_ptl1
run_one w012j_MLM.mg5 met_bias       w012j_met
run_one w012j_MLM.mg5 ptl1_met_bias  w012j_ptl1_met
run_one w012j_MLM.mg5 xsec_flat_bias w012j_xsecflat

echo "ALL_RUNS_DONE"
