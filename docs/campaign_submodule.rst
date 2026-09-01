.. _campaign_submodule:

The ``campaign`` submodule
==========================

The :mod:`rsutil.campaign` package prepares reproducible high-throughput
(HT) runs around the existing RS-LMTO setup operation.  It is deliberately a
filesystem and artifact layer: it validates inputs, materializes isolated run
directories, writes manifests, renders an optional Slurm array helper, and
collects local evidence.  It does **not** submit jobs, query a scheduler,
retry runs, or provide an SCF-to-exchange workflow engine.

The supported lifecycle is:

.. code-block:: text

    JSON/YAML specification -> validate -> create -> render-slurm
        -> user/site submission -> status -> collect

Specification format
--------------------

The primary schema is ``rsutil-campaign-v1``.  JSON is always supported;
YAML is accepted when PyYAML is installed.  A campaign contains a name,
one or more entries, optional shared setup/calculation templates, and a
destination naming policy.  Each entry identifies a source directory or an
explicit pair of ``lattice.nml`` and ``input.nml`` paths:

.. code-block:: json

    {
      "schema_version": "rsutil-campaign-v1",
      "name": "pilot",
      "destination_naming": "runs/{name}--{entry_id}",
      "structures": {"bcc-fe": "../../examples/bccFe/Dist"},
      "setup_template": {"rec_method": "block", "rec_steps": 20},
      "entries": [
        {
          "name": "bcc-fe-reference",
          "structure": "bcc-fe",
          "setup": {"seed": 11},
          "calculation": {"executable": "rs-lmto", "mode": "scf"}
        }
      ]
    }

Entry setup values override the campaign template.  ``exchange`` and
``calculation`` are metadata for downstream tooling; campaign creation does
not execute them.  Distortions require ``setup.seed`` and SQS requires
``setup.sqs_seed`` so that a campaign has a stable scientific identity.

Identity, paths, and reuse
--------------------------

An entry ID is a SHA-256 scientific hash of the source file hashes, normalized
setup, calculation metadata, exchange metadata, and entry metadata.  It is
independent of entry order and human-readable names.  Campaign-level metadata
is not included in entry IDs.  The generated index is
``campaign_index.json`` (schema ``rsutil-campaign-index-v1``); every run keeps
the setup ``run_manifest.json`` plus a campaign identity and materialized-input
hashes.

Scaffolding records an optional ``source.cif`` snapshot in campaign metadata,
under ``metadata.structure_sources`` keyed by structure alias.  This records
the crystallographic source file present when scaffolding ran; it is
historical/source provenance only and does not assert that the CIF
deterministically generated ``input.nml`` or ``lattice.nml``.  The prepared
RS-LMTO files remain the authoritative campaign inputs.  ``campaign validate``
and ``campaign create`` do not revalidate or require ``source.cif`` after
scaffolding.

Generation is preflighted and committed atomically.  Existing runs are reused
only when current source/setup output is regenerated in isolation and the
scientific output bytes match the recorded identity.  A directory without a
matching owned manifest is a collision and is never silently overwritten.
Use ``collision_policy="error"`` to reject all existing run directories.

Slurm rendering
---------------

``render_slurm_array`` writes a stable tab-separated run list and an executable
``sbatch`` script.  The script reads ``SLURM_ARRAY_TASK_ID``, changes into the
selected run directory, records atomic ``rsutil_execution.json`` evidence, and
executes the user-selected command.  Site-specific executable, launcher,
modules, partition, account, and resources are caller inputs:

.. code-block:: python

    from rsutil.campaign import load_campaign_index, render_slurm_array

    index = load_campaign_index("pilot-runs/campaign_index.json")
    rendered = render_slurm_array(
        index,
        executable="/opt/rs-lmto/bin/rs-lmto",
        launcher="srun --exclusive",
        campaign_root="pilot-runs",
        partition="compute",
        account="project",
        time_limit="02:00:00",
        resources={"cpus_per_task": 8},
        environment_preamble="module load gcc",
    )
    print(rendered.script_path)

The function renders files only; submit the resulting script separately after
reviewing it.

Status and collection
---------------------

Status uses local filesystem evidence only.  It retains three independent
axes so that process execution, parser validity, and scientific convergence
are not conflated:

.. list-table:: Campaign evidence axes
   :header-rows: 1

   * - Axis
     - Values
     - Meaning
   * - Execution
     - ``NOT_RUN``, ``RUNNING_OR_UNKNOWN``, ``EXITED_ZERO``, ``EXITED_NONZERO``, ``SIGNALED``
     - What the execution record proves about the process.
   * - Output/parser
     - ``NO_OUTPUT``, ``PARSE_VALID``, ``PARSE_INVALID``, ``PARSE_UNSUPPORTED``
     - Whether a supported and valid result is available.
   * - Scientific
     - ``CONVERGED``, ``NOT_CONVERGED``, ``UNKNOWN``, ``NOT_APPLICABLE``
     - What explicit result evidence establishes about SCF convergence.

``COMPLETED_VALID`` is reserved for zero process exit, valid parsed output,
and explicit convergence.  ``collect_campaign`` records malformed or missing
results per entry and continues with other entries.  Full SCF and exchange
models are preserved as JSON sidecars; CSV and JSONL are optional scalar
summaries.

Command-line entry points
-------------------------

.. code-block:: bash

    rsutil campaign validate examples/campaigns/ht18_pilot.json
    rsutil campaign create examples/campaigns/ht18_pilot.json --destination ht18-runs
    rsutil campaign render-slurm ht18-runs/campaign_index.json \
        --executable /opt/rs-lmto/bin/rs-lmto --launcher 'srun --exclusive' \
        --resource cpus_per_task=8 --time 02:00:00
    rsutil campaign status ht18-runs
    rsutil campaign collect ht18-runs --csv ht18-runs/results.csv

See :doc:`HT_CAMPAIGNS` for the file-level contract and :doc:`HT_READINESS_REPORT`
for qualification boundaries.

Public API
----------

.. automodule:: rsutil.campaign.models
   :members:
   :show-inheritance:

.. automodule:: rsutil.campaign.generation
   :members:
   :show-inheritance:

.. automodule:: rsutil.campaign.ht17
   :members:
   :show-inheritance:
