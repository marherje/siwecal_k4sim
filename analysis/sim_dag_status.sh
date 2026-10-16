#!/bin/bash
# Status of a make_sim_dag.py campaign: DAG node counts, condor jobs, groups finished on EOS, quota.
#   bash analysis/sim_dag_status.sh [campaign]      (default final_v7)
C=${1:-final_v7}
R="$(cd "$(dirname "${BASH_SOURCE[0]}")/.." && pwd)"
D=$R/analysis/condor_$C; O=/eos/experiment/drdcalo/siw-ecal/TB2026-06/Simulation/Processed/$C
out=$D/$C.dag.dagman.out
echo "== $C  $(date '+%F %T')"
awk '/Done +Pre +Queued +Post +Ready +Un-Ready +Failed/{getline; getline; l=$0} END{print "nodes  done pre queued post ready unready failed:"; print "      " l}' "$out"
grep -q "EXITING WITH STATUS" "$out" && grep "EXITING WITH STATUS" "$out" | tail -1
condor_q -constraint "JobBatchName==\"$C\" && JobUniverse==5" -af JobStatus 2>/dev/null \
  | awk '{st[$1]++} END{printf "jobs   running %d  idle %d  held %d\n", st[2], st[1], st[5]}'
ngrp=$(grep -c "^JOB g_" "$D/$C.dag"); nok=0
for d in "$O"/pid/*_nodigi; do [ -d "$d" ] && ls "$d"/*.valtree.root >/dev/null 2>&1 && nok=$((nok+1)); done
echo "groups finished on EOS: $nok / $ngrp"
grep -E "Node .* (failed|returned)" "$out" | grep -viE "returned 0|succeeded" | tail -5
eos root://eosproject.cern.ch quota /eos/experiment/drdcalo/siw-ecal/ 2>/dev/null | awk '/project/{print "quota  logical used " $4 " " $5 " of " $10 " " $11}'
