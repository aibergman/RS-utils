# Deterministic HT campaigns

HT-16 adds a small JSON campaign layer. It generates setup inputs and
provenance directories; it does not execute RS-LMTO or act as a scheduler.

## Specification

The primary format is JSON with schema `rsutil-campaign-v1`:

```json
{
  "schema_version": "rsutil-campaign-v1",
  "name": "pilot",
  "destination_naming": "runs/{name}--{entry_id}",
  "structures": {
    "bcc-fe": "../../examples/bccFe/Dist",
    "fen": "../../examples/FeN"
  },
  "setup_template": {"rec_method": "block", "rec_steps": 20},
  "entries": [
    {
      "name": "bcc-fe-reference",
      "structure": "bcc-fe",
      "setup": {"seed": 11},
      "calculation": {"executable": "rs-lmto", "mode": "scf"}
    },
    {
      "name": "fen-reference",
      "structure": "fen",
      "setup": {"seed": 12},
      "calculation": {"executable": "rs-lmto", "mode": "scf"}
    }
  ]
}
```

An entry can use `source` instead of `structure`. A directory source must
contain `lattice.nml` and `input.nml`; an explicit source mapping can provide
`lattice` and `input` file paths. Per-entry setup values override
`setup_template`. `exchange` is retained as optional metadata for later
post-SCF job generation and is not executed by HT-16.

For source trees that already contain prepared RS-LMTO bundles, scaffold a
conservative SCF campaign without converting CIF files or generating new
namelists:

```bash
rsutil campaign scaffold sources/ \
  --name initial-materials \
  --output campaign.json
```

The command recursively finds directories containing both `input.nml` and
`lattice.nml`, uses relative source paths, and validates the in-memory
campaign with the same validator used by `campaign validate`. `source.cif` is
optional; when present, its SHA-256 is recorded in the entry's generic
`metadata.structure_source` field (`file` is relative to the referenced
structure directory). Incomplete candidates are reported and
fail scaffolding so a production campaign cannot silently omit a source
bundle. If a bundle contains multiple `.cif` files but no canonical
`source.cif`, none is selected and the ambiguity is reported.

`source.cif` provenance establishes which crystallographic source file
accompanied the prepared RS-LMTO input. It does not prove that the RS-LMTO
input can be regenerated from that CIF unless the external conversion workflow
separately records its own provenance. RS-utils does not infer converter options,
versions, symmetry tolerances, or chemical agreement between the CIF and
namelists.

When distortion or SQS is requested, the campaign requires the corresponding
seed (`seed` for lattice/basis distortion and `sqs_seed` for SQS). This keeps
campaign generation reproducible.

## IDs, paths, and collisions

An entry ID is the SHA-256 scientific hash of its source file contents,
normalized setup, calculation metadata, and exchange settings. It is
independent of entry order and human name. The default run path is
`runs/{name}--{entry_id}`. The run manifest and campaign index store only
relative/logical paths, so the complete tree can be relocated.

`rsutil campaign create` validates every entry before creating the campaign
root. Each setup run is then committed atomically by `SetupService`. An
existing run is reused only after the current source files and setup
specification have been regenerated in an isolated candidate directory. The
candidate's immutable scientific output bytes and hashes are the authority;
the existing run manifest is historical metadata and its mutable hashes cannot
authorize a coordinated file-and-manifest edit. Each generated setup file is
listed under `campaign.materialized_inputs`; missing or changed files fail
reuse with an explicit integrity error. Unexpected `.nml` files are also
rejected as unrecorded scientific inputs. Runtime artifacts such as
`report.out`, logs, and `rsutil_execution.json` are not part of generation
identity and may be added after generation. Seeded stochastic setups may be
reused; unseeded stochastic generation is rejected conservatively. An
existing directory without a matching owned manifest is an explicit collision;
it is never silently overwritten. A modified specification under the same
human name gets a different entry ID and therefore a different default
directory.

```bash
rsutil campaign validate examples/campaigns/pilot.json
rsutil campaign create examples/campaigns/pilot.json --destination pilot-runs
```

The generated top-level `campaign_index.json` is
`rsutil-campaign-index-v1`; each `run_manifest.json` retains the existing
`rsutil-run-v1` setup manifest and adds a versioned
`rsutil-campaign-run-v2` campaign identity, materialized-input hashes, and
external-run metadata. Its integrity hash covers the campaign identity and
setup scientific identity fields, but it is only a tamper-evident record of
history—not the root of trust for reuse. Reuse uses the explicit
`fresh_regeneration_v1` contract described above. Older campaign manifests
without this explicit contract are not reused implicitly.

## Slurm arrays, status, and collection

HT-17 consumes the generated index and writes two transparent artifacts. The
run list is a stable tab-separated table with zero-based array index, entry
ID, run path, and manifest path. The script uses that table, quotes paths, and
runs one user-selected command in the selected run directory. The wrapper
writes `rsutil_execution.json` atomically, first recording a running record and
then the final return code, timestamps, command, working directory, host,
campaign identity, and available Slurm IDs. It does not submit the job or
select a site-specific module/MPI stack:

```bash
rsutil campaign render-slurm pilot-runs/campaign_index.json \
  --executable /opt/rs-lmto/bin/rs-lmto \
  --launcher 'srun --exclusive' \
  --partition compute --account project --time 02:00:00 \
  --resource cpus_per_task=8 --env-preamble 'module load gcc'
```

This writes `slurm_run_list.tsv` and `submit_array.sbatch` in the campaign
root by default. Run the rendered script through Slurm, for example with
`sbatch submit_array.sbatch`, when the site configuration is ready.

Filesystem status is intentionally conservative and does not query Slurm. It
retains three independent axes for each entry:

| Axis | Values | Meaning |
| --- | --- | --- |
| Execution | `NOT_RUN`, `RUNNING_OR_UNKNOWN`, `EXITED_ZERO`, `EXITED_NONZERO`, `SIGNALED` | What the generated execution record proves about the process. |
| Output/parser | `NO_OUTPUT`, `PARSE_VALID`, `PARSE_INVALID`, `PARSE_UNSUPPORTED` | Whether the expected result output exists and is supported/valid. |
| Scientific | `CONVERGED`, `NOT_CONVERGED`, `UNKNOWN`, `NOT_APPLICABLE` | What the parsed SCF evidence establishes about scientific convergence. |

The convenience `status` is `COMPLETED_VALID` only when execution exited zero,
the output is parser-valid, and SCF convergence is explicit. A parseable
non-converged or convergence-unknown result is collected as data but is not a
successful scientific result. A report without execution evidence is legacy
or unknown and is never promoted to success. A stale parseable report cannot
override a current non-zero execution record.

```bash
rsutil campaign status pilot-runs
```

It recognizes explicit failure markers, explicit lock markers, execution
evidence, and valid or malformed `report.out` parser evidence. A scheduler
live state is never inferred from directory age or the presence of a setup
manifest.

Collection parses each entry independently and writes the versioned
`campaign_results.json` aggregate. Full SCF and exchange models remain in
structured JSON sidecars under `collected/`; optional scalar summaries can be
exported as CSV or JSONL:

```bash
rsutil campaign collect pilot-runs --csv pilot-runs/results.csv --jsonl pilot-runs/results.jsonl
```

Parser and validation errors are recorded on their entry and do not prevent
other entries from being collected. Collection preserves the three status
axes and exposes `successful` only for explicit converged successes; the
`failures` view includes failed, invalid, non-converged, and convergence-
unknown completed entries. Raw file references and SHA-256 hashes, campaign
entry IDs, and both index and manifest-content hashes are retained.
Exchange sub-job manifests, when present, are collected by their canonical
pair lists and parent campaign entry ID. Retries, restarts, dependency graphs,
and scheduler daemons remain outside this layer and belong at a future
AiiDA/workflow-engine boundary.
