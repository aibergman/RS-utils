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
