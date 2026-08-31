# HT Quick Start

This page is the shortest path from an installed RS-utils release candidate to
a controlled high-throughput (HT) pilot. Read the support matrix before
scaling a campaign.

> **Pilot qualification**
>
> **RS-utils 0.11.0rc1 is ready for controlled pilot high-throughput RS-LMTO
> workflows within the explicitly documented support matrix.** It is a
> release candidate: some workflows remain `RESTRICTED`, `EXPERIMENTAL`, or
> `UNSUPPORTED` according to the [HT readiness report](HT_READINESS_REPORT).

## 1. Install a release wheel

For a user or cluster deployment, retrieve or build the release wheel and
install it into a dedicated environment. Do not use an editable source
checkout as the normal deployment artifact.

```bash
python -m venv rsutil-0.11rc1
source rsutil-0.11rc1/bin/activate

python -m pip install /path/to/rsutil-0.11.0rc1-<python>-<abi>-<platform>.whl
python -c "import rsutil; print(rsutil.__version__)"
```

The version check must print `0.11.0rc1`. Developers working on the source
tree can use `python -m pip install -e '.[test]'`; the reproducible build and
test procedure is documented in [HT testing](HT_TESTING).

## 2. Requirements

- Python 3.9 or newer is supported.
- The wheel includes the compiled Fortran component. A build from source
  needs a working CMake and Fortran compiler/runtime; the cluster runtime
  must be compatible with the wheel’s platform and compiler runtime.
- Core runtime dependencies are NumPy, pymatgen, and f90nml. `ase`, band/DOS
  analysis, visualization, and test tooling are optional extras.
- Slurm is required only for the rendered Slurm helper path. Campaign
  generation, status, and collection use local filesystem evidence and do not
  require a scheduler.

See [HT testing](HT_TESTING) for detailed build, wheel, and test requirements.

## 3. Legacy single-calculation workflows

The established entry points remain available for individual runs:

```bash
rsutil -S                         # setup in the current working directory
rsutil -J --cutoff 3.0 --unique  # generate exchange jobs
```

Use the existing {doc}`cli_section` for the full option set. HT
campaigns provide deterministic preparation around these established
scientific operations; they do not replace the legacy entry points.

## 4. Campaign lifecycle

The intended lifecycle is:

```text
campaign specification -> validate -> create -> render scheduler script
    -> submit -> status -> collect
```

Validate and create a representative campaign from the checked-in example:

```bash
rsutil campaign validate examples/campaigns/ht18_pilot.json
rsutil campaign create examples/campaigns/ht18_pilot.json \
    --destination pilot-runs
```

Render a site-specific Slurm array for review. The executable, launcher,
modules, partition, account, and resource choices belong to the cluster
administrator:

```bash
rsutil campaign render-slurm pilot-runs/campaign_index.json \
    --executable /path/to/rs-lmto \
    --launcher 'srun --exclusive' \
    --partition compute --account project --time 02:00:00 \
    --resource cpus_per_task=8 \
    --env-preamble 'module load gcc'

cd pilot-runs
sbatch submit_array.sbatch
cd ..

rsutil campaign status pilot-runs
rsutil campaign collect pilot-runs --csv pilot-runs/results.csv
```

`create` prepares setup-owned run directories but does not execute RS-LMTO.
`render-slurm` writes `slurm_run_list.tsv` and `submit_array.sbatch` but does
not submit them. Inspect the rendered files before submitting. `status` and
`collect` do not query Slurm, retry jobs, restart runs, or manage dependencies.
See [HT campaigns](HT_CAMPAIGNS) for the full schema and state table.

## 5. Interpret statuses scientifically

Execution success, parser validity, and SCF convergence are independent:

1. A zero process exit proves only that the launched command exited normally.
2. Parser validity proves that an output matches a supported result format.
3. Scientific convergence requires explicit SCF convergence evidence.

A parseable result is not necessarily converged. The convenience
`COMPLETED_VALID` state is reserved for an entry with zero execution exit,
valid parsed output, and explicit convergence. Non-converged or
convergence-unknown results remain collectable but are not successful
scientific results. See the [campaign state table](HT_CAMPAIGNS).

## 6. Exchange workflow

The supported exchange sequence is:

```text
converged SCF/reference -> -J exchange-job generation
    -> Atom_* calculations -> exchange-result parsing
```

RS-utils preserves the canonical periodic identity `PairKey(i,j,R)`, source
site identity, and RS type identity. Chemically identical atoms can remain
distinct RS types when the RS input assigns them separately through `izp`.
The exchange service and parsers also retain explicit source units and
conservative magnetic-symmetry decisions.

Read [exchange jobs](HT_EXCHANGE_JOBS), [exchange results](HT_EXCHANGE_RESULTS),
and [exchange symmetry](HT_EXCHANGE_SYMMETRY) before interpreting or reducing
large pair sets.

## 7. Pilot recommendation

> **Start small.** `0.11.0rc1` is a release candidate. Begin with a
> representative pilot campaign and inspect convergence, manifests, and
> parsed results before scaling to a large production database.

## 8. Deployment and provenance

The recommended cluster shape is:

```text
release wheel -> dedicated venv/environment module -> fixed RS-LMTO executable
    -> campaign manifest -> Slurm array
```

Preserve the following with every real campaign:

- RS-utils version and, when known, the Git commit;
- wheel filename and SHA-256 hash;
- RS-LMTO executable version and provenance;
- `campaign_index.json`, run manifests, and execution records;
- raw outputs, parsed result sidecars, and collection summaries.

This is deployment guidance only. The campaign layer does not add a
provenance server, restart manager, scheduler daemon, or automatic
SCF-to-exchange orchestration.

## Support boundaries

The [HT readiness report](HT_READINESS_REPORT) is authoritative for the
support matrix. Its labels mean:

- `READY`: validated within the stated boundary;
- `RESTRICTED`: usable only with the stated qualification or limited domain;
- `EXPERIMENTAL`: available for evaluation, not a general production claim;
- `UNSUPPORTED`: not a validated workflow or deliberately fail-closed.

Existing UppASD export functionality remains available, but automated UppASD
integration is not currently part of the validated `0.11.0rc1` HT campaign
workflow. In particular, this release does not claim an automatic
`SCF -> exchange -> UppASD` pipeline.
