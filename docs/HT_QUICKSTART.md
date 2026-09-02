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
prepared RS source -> SCF campaign -> converged canonical result
    -> explicit exchange stage -> canonical exchange result
```

Scientific identity, operational resources, source provenance, and exchange
stage provenance are separate: setup/source content determines the scientific
entry identity, scheduler allocation is operational metadata, CIF hashes are
campaign-level source evidence, and an exchange stage records the parent
result and policy used to prepare its jobs.

### Scaffold a prepared source collection

For a collection of prepared RS-LMTO source bundles, use an external
converter workflow:

1. Generate `input.nml` and `lattice.nml` using `cif2cell` or another trusted
   tool.
2. Keep `source.cif` alongside the source bundle when possible.
3. Run `rsutil campaign scaffold`.
4. Inspect the discovery and validation summary.
5. Run `rsutil campaign validate`.
6. Run `rsutil campaign create`.

```bash
rsutil campaign scaffold examples/campaigns/scaffold_sources \
    --name scaffolded-pilot \
    --output examples/campaigns/scaffolded_pilot.json
rsutil campaign validate examples/campaigns/scaffolded_pilot.json
rsutil campaign create examples/campaigns/scaffolded_pilot.json \
    --destination scaffolded-pilot-runs
```

An existing scaffold output is rejected by default. Use `--overwrite` only
for a deliberate replacement:

```bash
rsutil campaign scaffold examples/campaigns/scaffold_sources \
    --name scaffolded-pilot \
    --output examples/campaigns/scaffolded_pilot.json --overwrite
```

RS-utils does not require `cif2cell` and does not call it during scaffolding.
Scaffolding discovers existing prepared bundles only. In the campaign model,
`structures` is the source/material catalog, while `entries` is the
calculation matrix. Later, one structure can therefore be expanded into
multiple lattice scans or calculation recipes.

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

Entries may optionally carry operational resource requirements, for example
`"resources": {"ntasks": 8, "memory": "4G", "time": "02:00:00"}`.
With entry overrides, rendering creates one Slurm array per homogeneous
resource class and writes `resource_plan.json` plus `submit_all.sh`. Resource
allocation is operational policy, not scientific calculation identity. RS-utils
does not yet infer an optimal MPI size; use early HT campaigns to collect
scaling data before introducing an automatic policy.

`create` prepares setup-owned run directories but does not execute RS-LMTO.
For a homogeneous campaign, `render-slurm` writes the familiar
`slurm_run_list.tsv` and `submit_array.sbatch`; resource-aware campaigns put
these files under their class directories. It does not submit them. Inspect
the rendered files before submitting. `status` and `collect` do not query
Slurm, retry jobs, restart runs, or manage dependencies.
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

Exchange is an explicit second campaign stage. The supported sequence is:

```text
SCF campaign
    ↓
inspect convergence
    ↓
explicit exchange prepare
    ↓
render Slurm array
    ↓
execute Atom_* jobs
    ↓
exchange status
    ↓
exchange collect
```

By default, only successfully executed, parse-valid, explicitly converged
parent calculations are eligible for exchange preparation. Full pairs are
the conservative starting point; symmetry reduction requires an explicitly
supported magnetic context.

```bash
rsutil campaign exchange prepare pilot-runs/campaign_index.json \
    --destination pilot-runs/exchange-stage \
    --cutoff 3.0 --full-pairs

# Use the parent directory of the manifest path printed by prepare.
rsutil campaign exchange render-slurm \
    pilot-runs/exchange-stage/<stage-id> \
    --executable /path/to/rs-lmto-exchange \
    --launcher 'srun --exclusive' \
    --partition compute --account project --time 02:00:00 \
    --resource cpus_per_task=8

cd pilot-runs/exchange-stage/<stage-id>
sbatch exchange_submit_array.sbatch
cd ../../..

rsutil campaign exchange status pilot-runs/exchange-stage/<stage-id>
rsutil campaign exchange collect pilot-runs/exchange-stage/<stage-id>
```

Each stage is stored under its own `<stage-id>/exchange_stage_manifest.json`,
with its own jobs and artifacts. A changed parent result or exchange policy
therefore creates a new persistent stage; identical preparation reuses valid
exchange-job directories using the canonical exchange-job checks.

Use `--unique` only with an explicitly supported magnetic context, for
example `--magnetic-reference collinear_fm` and any required axis. RS-utils
does not yet implement a generic automatic SCF→exchange DAG or retry engine.

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
