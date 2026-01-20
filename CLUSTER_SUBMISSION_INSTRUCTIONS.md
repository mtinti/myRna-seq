# Dundee cluster batch submission instructions (no shared filesystem)

These instructions submit **myRna-seq** as a single batch job using `qsub` and a
wrapper script. The wrapper stages the workflow into node-local scratch, runs
Snakemake, and syncs results back at the end. The wrapper is configured for the
Dundee cluster’s SGE directives shown at the top of the script.

## 1) Review the wrapper script

Open `submit_snakemake_cluster.sh` and update resource requests (cores, memory,
walltime) and any path overrides.

Before submitting, make sure your pipeline configuration is complete in
`config.yaml` and that it references a properly paired `samples.csv` file. The
required `samples.csv` format and configuration instructions are documented in
the GitHub README for this repository.

Key variables you can customize inside the script:

- `CORES`: Snakemake cores to use (defaults to 40 to match `-pe smp 40`).
- `CONFIG_OVERRIDES`: Runtime overrides for `processing_dir`, `results_dir`,
  `benchmark_dir`, and any reference/sample paths.

## 2) Submit the job

From the directory containing the pipeline repository:

```bash
qsub submit_snakemake_cluster.sh
```

Check the job queue with `qstat` and review the job’s stdout/stderr log once it
finishes.

## 3) Retrieve results

The wrapper script syncs `results/`, `benchmarks/`, and Snakemake logs back to
your submission directory by default. Adjust the `RESULTS_DEST` variable in the
wrapper if you want outputs elsewhere.

---

### Notes

- Because the cluster has **no shared filesystem**, the workflow runs on a
  **single node**. Do not use Snakemake cluster submission.
- If the compute node cannot access your submission directory, update the
  staging and sync commands in the wrapper to use your preferred transfer
  method (e.g., `scp` from a data transfer node).
