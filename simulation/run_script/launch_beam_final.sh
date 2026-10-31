#!/bin/bash
# Prepare (and, with SUBMIT=1, submit) the final_v6 ddsim campaign: one GPS beam sample per group of
# analysis/final_v6_groups.txt, at the group's measured beam centre, with the beam size measured in data at that
# energy (analysis/beam_spread.txt, analysis/beam_spread.py), split into chunks with their own seeds.
#
#   bash simulation/run_script/launch_beam_final.sh              # writes the .sub + job list, submits nothing
#   SUBMIT=1 bash simulation/run_script/launch_beam_final.sh     # ... and condor_submit
#
# Re-running skips chunks already on EOS, so it also resubmits the failed ones.
#
# Output: Generated/final_v6/output_beam_<p>_<E>GeV_<pos>_xy_<x>_<y>_sigx<sx>_sigy<sy>_sigE<sE><TAG>_c<NNN>.edm4hep.root
#   <pos> = table position (E1..E5, U1, S1, T<H>_<V>, see the groups file), so a file is found by what the
#   run summary calls it; x_y and the sigmas stay in the name to keep it self-describing.
#
# Environment:
#   NEVENTS        events per sample (default 50000)
#   CHUNK_EVENTS   events per condor job (default 2000: ~4 h and ~5 GB at 99 GeV, ~25 min at 7.5 GeV)
#   GROUPS_FILE         groups file (default analysis/final_v6_groups.txt)
#   SIGMAS         per-energy beam sigma table (default analysis/beam_spread.txt)
#   ONLY_SET / ONLY_POS / ONLY_E   restrict to one threshold set / table position / energy
#   PARTICLE (e-), SIGMA_E (0.02), PHYSICS_LIST (QGSP_BERT_EMZ), COMPACT_FILE (SND_compact_beamline.xml),
#   BEAM_Z_MM (-65500: the H2-like line of the default compact since 2026-10-26; -2650 with SND_compact_beamline_v7.xml), TAG (_beamline_real_w1924_mat_emz): the final_v5 settings
#   OUT_DIR        EOS output dir (default .../Simulation/Generated/final_v6)
#   SIGMA_X SIGMA_Y  beam sigma [mm] for every sample, instead of the per-energy table (e.g. pi-/mu- points)
set -uo pipefail
R="$(cd "$(dirname "${BASH_SOURCE[0]}")/../.." && pwd)"
NEVENTS=${NEVENTS:-50000}
CHUNK_EVENTS=${CHUNK_EVENTS:-2000}
GROUPS_FILE=${GROUPS_FILE:-$R/analysis/final_v6_groups.txt}
SIGMAS=${SIGMAS:-$R/analysis/beam_spread.txt}
PARTICLE=${PARTICLE:-e-}
SIGMA_E=${SIGMA_E:-0.02}
export PHYSICS_LIST=${PHYSICS_LIST:-QGSP_BERT_EMZ}
export COMPACT_FILE=${COMPACT_FILE:-SND_compact_beamline.xml}
export BEAM_Z_MM=${BEAM_Z_MM:--65500}
TAG=${TAG-_beamline_real_w1924_mat_emz}
OUT_DIR=${OUT_DIR:-/eos/experiment/drdcalo/siw-ecal/TB2026-06/Simulation/Generated/final_v6}
STEER=$R/simulation/run_script/steer/final_v6
NCHUNKS=$(( (NEVENTS + CHUNK_EVENTS - 1) / CHUNK_EVENTS ))

[ -s "$GROUPS_FILE" ] || { echo "no groups file $GROUPS_FILE" >&2; exit 1; }
[ -s "$SIGMAS" ] || { echo "no sigma table $SIGMAS (run analysis/beam_spread.py)" >&2; exit 1; }
[ -s "$R/simulation/geometry/$COMPACT_FILE" ] || { echo "no geometry $COMPACT_FILE" >&2; exit 1; }
mkdir -p "$STEER/logs"

# per-energy sigma, rounded to 0.1 mm: "E sx sy" lines of the 'per energy' block
sigma_for() {
  awk -v E="$1" '/^# per energy/{on=1; next} on && !/^#/ && NF>=4 && $1+0==E+0 {printf "%.1f %.1f\n", $3, $4; f=1}
                 END{if(!f) exit 1}' "$SIGMAS"
}

# seed offset per particle, so e-, pi- and mu- samples of one groups line never share a seed
case "$PARTICLE" in e-) poff=6000000 ;; pi-) poff=7000000 ;; mu-) poff=8000000 ;; *) poff=9000000 ;; esac
list=$STEER/jobs_$(date +%Y%m%d_%H%M%S).txt; : > "$list"
ngroups=0; nskip=0; gidx=0
while read -r th pos E xy runs; do
  gidx=$((gidx + 1))
  [ -n "${ONLY_SET:-}" ] && [ "$th" != "$ONLY_SET" ] && continue
  [ -n "${ONLY_POS:-}" ] && [ "$pos" != "$ONLY_POS" ] && continue
  [ -n "${ONLY_E:-}" ] && [ "$E" != "$ONLY_E" ] && continue
  if [ -n "${SIGMA_X:-}" ]; then sx=$SIGMA_X; sy=${SIGMA_Y:?SIGMA_Y with SIGMA_X}
  else read -r sx sy < <(sigma_for "$E") || { echo "no beam sigma for $E GeV in $SIGMAS" >&2; exit 1; }; fi
  x=${xy%_*}; y=${xy#*_}
  point=output_beam_${PARTICLE}_${E}GeV_${pos}_xy_${x}_${y}_sigx${sx}_sigy${sy}_sigE${SIGMA_E}${TAG}
  ngroups=$((ngroups + 1))
  printf "%-6s %-11s %4s GeV  (%4s,%4s) mm  sigma %4s x %4s mm  %d x %d events\n" \
         "$th" "$pos" "$E" "$x" "$y" "$sx" "$sy" "$NCHUNKS" "$CHUNK_EVENTS"
  for ((c = 0; c < NCHUNKS; c++)); do
    n=$CHUNK_EVENTS; (( (c + 1) * CHUNK_EVENTS > NEVENTS )) && n=$(( NEVENTS - c * CHUNK_EVENTS ))
    out=$OUT_DIR/$(printf "%s_c%03d.edm4hep.root" "$point" "$c")
    if [ -s "$out" ]; then nskip=$((nskip + 1)); continue; fi
    # deterministic, collision-free seed: campaign offset + group line + chunk
    seed=$(( poff + gidx * 1000 + c ))
    echo "$n $PARTICLE $E $x $y $sx $sy $SIGMA_E $seed $out $(basename "$out" .edm4hep.root)" >> "$list"
  done
done < <(grep -v '^#' "$GROUPS_FILE" | awk 'NF>=5')

njobs=$(wc -l < "$list")
echo "--- $ngroups samples, $njobs chunk jobs to run, $nskip already on EOS -> $list"
[ "$njobs" -gt 0 ] || { rm -f "$list"; exit 0; }

sub=$STEER/final_v6.sub
cat > "$sub" <<EOF
universe                = vanilla
executable              = $R/simulation/run_script/ddsim_beam_chunk.sh
arguments               = $R \$(N) \$(P) \$(E) \$(X) \$(Y) \$(SX) \$(SY) \$(SE) \$(SEED) \$(OUT)
environment             = "COMPACT_FILE=$COMPACT_FILE BEAM_Z_MM=$BEAM_Z_MM PHYSICS_LIST=$PHYSICS_LIST"
log                     = $STEER/logs/final_v6.log
output                  = $STEER/logs/\$(NAME).out
error                   = $STEER/logs/\$(NAME).err
request_cpus            = 1
request_memory          = 3 GB
request_disk            = 10 GB
should_transfer_files   = YES
when_to_transfer_output = ON_EXIT
transfer_output_files   = ""
+JobFlavour             = "workday"
queue N, P, E, X, Y, SX, SY, SE, SEED, OUT, NAME from $list
EOF
echo "--- submit file: $sub"

if [ "${SUBMIT:-0}" = 1 ]; then
  # plain mkdir per missing level: mkdir -p walks into EOS parents it cannot touch
  d=$OUT_DIR; miss=(); while [ ! -d "$d" ]; do miss=("$d" "${miss[@]}"); d=$(dirname "$d"); done
  for d in "${miss[@]}"; do mkdir "$d" || exit 1; done
  (cd "$STEER" && condor_submit "$sub")
else
  echo "--- not submitted (SUBMIT=1 to submit, or: condor_submit $sub)"
fi
