#!/usr/bin/env python3
"""Write (and with --submit, submit) the HTCondor DAG that regenerates the data valtrees with barxy_min_border_distance.

One node per data run of analysis/beam_points.txt (analysis/border_data_node.sh): the PID stage of the current
siwecal-tb2026 build on Reconstruction/<run>/ecal_<run>.root (read only), with the Reconstruction settings (hit MIP
cut 0.5, mip05_/mip1_ blocks), EDM4hep + valtree to the staging <out>/<run>/ecal_<run>.{edm4hep,valtree}.root, which
then replace those of Reconstruction/<run>/ (the new data reconstruction; final_v7 stays for the simulation only).

    python3 analysis/make_border_data_dag.py --name reco_border   # writes analysis/condor_reco_border/
    python3 analysis/make_border_data_dag.py --only TB2026CERN_run_000020 --name border_data_test

Runs whose EDM4hep and valtree are both already in <out> are written as DONE.
"""
import argparse
import math
import os
import subprocess
import sys

REPO = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
B = "/eos/experiment/drdcalo/siw-ecal/TB2026-06"
# Memory of the PID stage scales with the run's size: run 20 (1.24 GB ecal tree) peaked at 3181 MB. Requested with a
# 30 % margin, and with one core per 2 GB, the CERN slot shape (a single-core job asking for more waits for a
# drained node).
REF_BYTES, REF_MB, MARGIN, MB_PER_CPU = 1.24e9, 3181, 1.3, 2000


def resources(path):
    mb = max(2000, 500 * math.ceil(REF_MB * os.path.getsize(path) / REF_BYTES * MARGIN / 500))
    return mb, math.ceil(mb / MB_PER_CPU)


def run_name(r):
    return f"TB2026CERN_eudaq_run_{int(r[1:]):06d}" if r.startswith("e") else f"TB2026CERN_run_{int(r):06d}"


def read_runs(path):
    runs = []
    for line in open(path):
        f = line.split("#")[0].split()
        if len(f) < 9:
            continue
        for r in f[8].split(","):
            if run_name(r) not in runs:
                runs.append(run_name(r))
    return sorted(runs)


def main():
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--points", default=os.path.join(REPO, "analysis", "beam_points.txt"))
    ap.add_argument("--reco", default=f"{B}/Reconstruction")
    ap.add_argument("--out", default=f"{B}/Simulation/tests/border_cut/data")
    ap.add_argument("--name", default="border_data", help="DAG name; analysis/condor_<name>/")
    ap.add_argument("--only", default=None, help="comma-separated run names")
    ap.add_argument("--submit", action="store_true")
    args = ap.parse_args()

    runs = read_runs(args.points)
    if args.only:
        runs = [r for r in runs if r in set(args.only.split(","))]
    missing = [r for r in runs if not os.path.exists(f"{args.reco}/{r}/ecal_{r}.root")]
    if missing:
        sys.exit(f"missing ecal trees: {missing}")
    dag_dir = os.path.join(REPO, "analysis", f"condor_{args.name}")
    os.makedirs(os.path.join(dag_dir, "logs"), exist_ok=True)

    sub = os.path.join(dag_dir, "node.sub")
    with open(sub, "w") as f:
        f.write(f"""universe                = vanilla
executable              = /bin/bash
arguments               = {REPO}/analysis/border_data_node.sh $(RUN) {args.reco} {args.out}
log                     = {dag_dir}/logs/dag.log
output                  = {dag_dir}/logs/$(RUN).out
error                   = {dag_dir}/logs/$(RUN).err
request_cpus            = $(CPUS)
request_memory          = $(MEM)
request_disk            = 10000M
should_transfer_files   = NO
getenv                  = False
+JobFlavour             = "workday"
queue
""")
    lines, ndone = [], 0
    for r in runs:
        done = all(os.path.exists(f"{args.out}/{r}/ecal_{r}.{k}.root") for k in ("edm4hep", "valtree"))
        ndone += done
        mb, cpus = resources(f"{args.reco}/{r}/ecal_{r}.root")
        lines += [f"JOB {r} {sub}" + (" DONE" if done else ""), f'VARS {r} RUN="{r}" MEM="{mb}" CPUS="{cpus}"',
                  f"RETRY {r} 1"]
    dag = os.path.join(dag_dir, f"{args.name}.dag")
    with open(dag, "w") as f:
        f.write("\n".join(lines) + "\n")
    print(f"--- {len(runs)} runs ({ndone} already done) -> {dag}\n--- output: {args.out}")
    if args.submit:
        if not os.path.isdir(args.out):
            sys.exit(f"create {args.out} first (plain mkdir per level: mkdir -p fails on EOS)")
        sys.exit(subprocess.call(["condor_submit_dag", "-batch-name", args.name, os.path.basename(dag)], cwd=dag_dir))


if __name__ == "__main__":
    main()
