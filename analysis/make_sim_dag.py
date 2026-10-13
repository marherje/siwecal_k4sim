#!/usr/bin/env python3
"""Write (and with --submit, submit) the HTCondor DAG of a simulation campaign over analysis/beam_points.txt.

Per beam point: NCHUNKS chunk nodes (simulation/run_script/sim_digi_chunk.sh: ddsim + digitised and undigitised
ecal trees, staged to <out>/chunks/<name>/), then one group node (analysis/sim_reco_group.sh: hadd, ECalPid with the
0.5 MIP hit cut, valtree). Geometry/physics as final_v5/v6: SND_compact_beamline.xml, QGSP_BERT_EMZ, beam born at
z = -2650 mm, e- with sigma_E 2 %. Digitiser: the job3_digitize.py defaults (EOS Calib/MuonCalib, th210 tables,
PED_NOISE=1), reconstructed with the th210 MIP table like the data in Reconstruction.

    python3 analysis/make_sim_dag.py                       # writes analysis/condor_<campaign>/
    python3 analysis/make_sim_dag.py --submit              # ... and condor_submit_dag
    python3 analysis/make_sim_dag.py --only th230_E4_e20 --nevents 20 --chunk-events 10 --campaign test  # small test

Chunks whose two trees are already on EOS (or groups whose valtrees are) are written as DONE, so re-running the
script after a failure resubmits only what is missing. Seeds: SEED_OFFSET + point line * 1000 + chunk.
"""
import argparse
import os
import subprocess
import sys

REPO = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
SIM = "/eos/experiment/drdcalo/siw-ecal/TB2026-06/Simulation"


def read_points(path):
    points = []
    for line in open(path):
        if line.startswith("#") or len(line.split()) < 9:
            continue
        name, set_, label, E, x, y, sx, sy, runs = line.split()[:9]
        points.append(dict(name=name, set=set_, label=label, E=E, x=x, y=y, sx=sx, sy=sy, runs=runs))
    return points


def eos_mkdir(path):
    """Plain mkdir per missing level (mkdir -p walks into EOS parents it cannot touch)."""
    missing = []
    while not os.path.isdir(path):
        missing.insert(0, path)
        path = os.path.dirname(path)
    for d in missing:
        os.mkdir(d)


def main():
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--points", default=os.path.join(REPO, "analysis", "beam_points.txt"))
    ap.add_argument("--campaign", default="final_v7")
    ap.add_argument("--out", default=None, help="EOS output base (default Simulation/Processed/<campaign>)")
    ap.add_argument("--nevents", type=int, default=50000)
    ap.add_argument("--chunk-events", type=int, default=2000)
    ap.add_argument("--particle", default="e-")
    ap.add_argument("--sigma-e", default="0.02")
    ap.add_argument("--seed-offset", type=int, default=10000000)
    ap.add_argument("--only", default=None, help="comma-separated point names")
    ap.add_argument("--submit", action="store_true")
    args = ap.parse_args()

    out = args.out or os.path.join(SIM, "Processed", args.campaign)
    dag_dir = os.path.join(REPO, "analysis", f"condor_{args.campaign}")
    os.makedirs(os.path.join(dag_dir, "logs"), exist_ok=True)
    nch = (args.nevents + args.chunk_events - 1) // args.chunk_events
    points = read_points(args.points)
    only = set(args.only.split(",")) if args.only else None

    chunk_sub = os.path.join(dag_dir, "chunk.sub")
    with open(chunk_sub, "w") as f:
        f.write(f"""universe                = vanilla
executable              = {REPO}/simulation/run_script/sim_digi_chunk.sh
arguments               = {REPO} $(N) {args.particle} $(E) $(X) $(Y) $(SX) $(SY) {args.sigma_e} $(SEED) $(SET) $(REMOTE)
environment             = "COMPACT_FILE=SND_compact_beamline.xml BEAM_Z_MM=-2650 PHYSICS_LIST=QGSP_BERT_EMZ RECO_TABLE=th210"
log                     = {dag_dir}/logs/dag.log
output                  = {dag_dir}/logs/$(NAME).out
error                   = {dag_dir}/logs/$(NAME).err
request_cpus            = 1
request_memory          = 4 GB
request_disk            = 15 GB
should_transfer_files   = YES
when_to_transfer_output = ON_EXIT
transfer_output_files   = ""
+JobFlavour             = "$(FLAVOUR)"
queue
""")
    group_sub = os.path.join(dag_dir, "group.sub")
    with open(group_sub, "w") as f:
        f.write(f"""universe                = vanilla
executable              = /bin/bash
arguments               = {REPO}/analysis/sim_reco_group.sh $(NAME) {out} {nch}
log                     = {dag_dir}/logs/dag.log
output                  = {dag_dir}/logs/group_$(NAME).out
error                   = {dag_dir}/logs/group_$(NAME).err
request_cpus            = 1
request_memory          = 8000
request_disk            = 20000M
should_transfer_files   = NO
getenv                  = False
+JobFlavour             = "workday"
queue
""")

    if args.submit:
        for d in (out, os.path.join(out, "chunks")):
            eos_mkdir(d)
    lines = []
    nchunks = ndone = 0
    for i, p in enumerate(points, 1):
        if only and p["name"] not in only:
            continue
        name = p["name"]
        cdir = os.path.join(out, "chunks", name)
        if args.submit:
            eos_mkdir(cdir)
        group_done = all(any(f.endswith(".valtree.root") for f in os.listdir(os.path.join(out, "pid", f"{name}_{k}")))
                         if os.path.isdir(os.path.join(out, "pid", f"{name}_{k}")) else False
                         for k in ("digi", "nodigi"))
        flavour = "tomorrow" if float(p["E"]) >= 52 else "workday"
        parents = []
        for c in range(nch):
            n = min(args.chunk_events, args.nevents - c * args.chunk_events)
            node = f"c_{name}_{c:03d}"
            remote = os.path.join(cdir, f"ecal_{name}_c{c:03d}")
            seed = args.seed_offset + i * 1000 + c
            lines.append(f"JOB {node} {chunk_sub}" + (" DONE" if group_done or (
                os.path.exists(remote + ".root") and os.path.exists(remote + "_simple.root")) else ""))
            lines.append(f'VARS {node} N="{n}" E="{p["E"]}" X="{p["x"]}" Y="{p["y"]}" SX="{p["sx"]}" SY="{p["sy"]}" '
                         f'SEED="{seed}" SET="{p["set"]}" REMOTE="{remote}" NAME="{node}" FLAVOUR="{flavour}"')
            lines.append(f"RETRY {node} 2")
            parents.append(node)
            nchunks += 1
        g = f"g_{name}"
        lines.append(f"JOB {g} {group_sub}" + (" DONE" if group_done else ""))
        lines.append(f'VARS {g} NAME="{name}"')
        lines.append(f"RETRY {g} 1")
        lines.append(f"PARENT {' '.join(parents)} CHILD {g}")
        ndone += group_done
        print(f"{name:<30} {p['E']:>4} GeV ({p['x']:>4},{p['y']:>4}) mm  sigma {p['sx']} x {p['sy']} mm  "
              f"{nch} x {args.chunk_events}{'  [done]' if group_done else ''}")
    dag = os.path.join(dag_dir, f"{args.campaign}.dag")
    with open(dag, "w") as f:
        f.write("\n".join(lines) + "\n")
    print(f"--- {nchunks} chunk nodes, {len(lines) and nchunks // nch} groups ({ndone} already done) -> {dag}\n"
          f"--- output: {out}")
    if args.submit:
        sys.exit(subprocess.call(["condor_submit_dag", "-batch-name", args.campaign, os.path.basename(dag)],
                                 cwd=dag_dir))


if __name__ == "__main__":
    main()
