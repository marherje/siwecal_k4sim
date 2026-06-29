# Simulation run scripts

## Output location

Simulation output (`.edm4hep.root` files produced by `ddsim`) is written directly to EOS:

```
/eos/experiment/drdcalo/siw-ecal/TB2026-06/Simulation/
```

The local `data/` subdirectory is no longer used and should not be created.

## Scripts

| Script | Description |
|--------|-------------|
| `generic_condor_beam.sh` | HTCondor job for beam-test simulation (GPS-driven) |
| `generic_condor_PG.sh`   | HTCondor job for particle-gun simulation |
| `launch_beam.sh`         | Wrapper to submit beam jobs |
| `launch_PG.sh`           | Wrapper to submit PG jobs |
| `clean.sh`               | Remove generated steer/log files |
