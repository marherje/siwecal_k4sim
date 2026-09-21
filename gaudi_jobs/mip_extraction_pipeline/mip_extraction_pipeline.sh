#!/bin/bash
# MIP extraction: step 0 of the digitisation.
#
#   mip_extraction_pipeline.sh <muon ddsim file> [output yaml] [geometry label]
#
# Fits the per-layer MIP peak of the simulated stack (MIPExtractor) and writes it
# to mappings/mip_values_sim.yml, which gaudi_jobs/pid2026_common/job3_digitize.py
# reads -- job3 refuses to run without it, so the values are never pasted by hand.
#
# The old three-step form (job1_shuffler -> job2_splitter -> job3_mipextract on
# timewindows.edm4hep.root) is still there for pile-up studies; for the MIP scale a
# plain muon ddsim file is the right input, with no shuffling in between.
set -eo pipefail
REPO=$(cd "$(dirname "${BASH_SOURCE[0]}")/../.." && pwd)
INPUT=${1:?muon ddsim file}
OUT=${2:-$REPO/mappings/mip_values_sim.yml}
GEOM=${3:-simulation/geometry/SND_compact.xml}
WORK=$(mktemp -d "${TMPDIR:-/tmp}/mipextract_XXXX")
cd "$WORK"
INPUT_FILE="$INPUT" INPUT_COLLECTION=${INPUT_COLLECTION:-SiPadHits} MIP_VALUES_OUT=mip_values.txt \
  MIP_HISTOS_OUT=mip_extraction.root k4run "$REPO/gaudi_jobs/mip_extraction_pipeline/job3_mipextract.py" > mipextract.log 2>&1 \
  || { tail -20 mipextract.log; exit 1; }
python3 - "$INPUT" "$OUT" "$GEOM" <<'PY'
import sys, os, datetime
src, out, geom = sys.argv[1:4]
values = []
for line in open("mip_values.txt"):
    if line.startswith("#") or not line.strip():
        continue
    parts = line.split()
    if len(parts) == 2:
        values.append(float(parts[1]))
assert len(values) == 15, f"expected 15 layers, got {len(values)}"
assert all(v > 0 for v in values), f"a layer has no MIP value: {values}"
with open(out, "w") as fh:
    fh.write("# Per-layer MIP scale of the SIMULATED stack [GeV per MIP], the peak of the raw\n"
             "# Geant4 energy per hit of a muon sample (MIPExtractor, Landau x Gaussian).\n"
             "# Written by gaudi_jobs/mip_extraction_pipeline/mip_extraction_pipeline.sh;\n"
             "# read by gaudi_jobs/pid2026_common/job3_digitize.py.  Do not edit by hand:\n"
             "# rerun the pipeline when the geometry (sensor thickness) changes.\n")
    fh.write(f"source_sample: {os.path.basename(src)}\n")
    fh.write(f"geometry: {geom}\n")
    fh.write(f"fit: langau peak (MIPExtractor FitMode 3)\n")
    fh.write(f"date: {datetime.date.today().isoformat()}\n")
    fh.write("mip_gev:\n")
    for i, v in enumerate(values):
        fh.write(f"  - {v:.6e}   # layer {i:2d}\n")
print("layer  MIP [GeV]   / layer 1")
for i, v in enumerate(values):
    print(f"{i:5d}  {v:.4e}   {v / values[1]:.3f}")
print(f"written {out}")
PY
cp mip_extraction.root "${MIP_HISTOS_KEEP:-/dev/null}" 2>/dev/null || true
cd / && rm -rf "$WORK"
