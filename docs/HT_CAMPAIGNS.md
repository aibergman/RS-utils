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
existing run is reused only when its campaign manifest contains the same
entry and manifest hashes. An existing directory without a matching owned
manifest is an explicit collision; it is never silently overwritten. A
modified specification under the same human name gets a different entry ID
and therefore a different default directory.

```bash
rsutil campaign validate examples/campaigns/pilot.json
rsutil campaign create examples/campaigns/pilot.json --destination pilot-runs
```

The generated top-level `campaign_index.json` is
`rsutil-campaign-index-v1`; each `run_manifest.json` retains the existing
`rsutil-run-v1` setup manifest and adds campaign identity and external-run
metadata.

## Slurm arrays, status, and collection

HT-17 consumes the generated index and writes two transparent artifacts. The
run list is a stable tab-separated table with zero-based array index, entry
ID, run path, and manifest path. The script uses that table, quotes paths, and
executes one user-selected command in the selected run directory. It does not
submit the job or select a site-specific module/MPI stack:

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

Filesystem status is intentionally conservative and does not query Slurm:

```bash
rsutil campaign status pilot-runs
```

It recognizes explicit failure markers, explicit lock markers, valid or
malformed `report.out` parser evidence, and otherwise reports a generated or
unknown state. A scheduler live state is never inferred from directory age or
the presence of a setup manifest.

Collection parses each entry independently and writes the versioned
`campaign_results.json` aggregate. Full SCF and exchange models remain in
structured JSON sidecars under `collected/`; optional scalar summaries can be
exported as CSV or JSONL:

```bash
rsutil campaign collect pilot-runs --csv pilot-runs/results.csv --jsonl pilot-runs/results.jsonl
```

Parser and validation errors are recorded on their entry and do not prevent
other entries from being collected. Raw file references and SHA-256 hashes,
campaign entry IDs, and both index and manifest-content hashes are retained.
Exchange sub-job manifests, when present, are collected by their canonical
pair lists and parent campaign entry ID. Retries, restarts, dependency graphs,
and scheduler daemons remain outside this layer and belong at a future
AiiDA/workflow-engine boundary.
