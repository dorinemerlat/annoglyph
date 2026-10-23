Column suffixes specify hours, minutes, seconds or GiB. Missing values are blank, not zero.
Only actual task executions contribute to CPU/runtime totals; CACHED records are counted separately.
storeDir-skipped processes may be absent entirely: supply earlier traces to recover their measurements.
Repeated records of one execution are deduplicated; actual retries and recomputations remain included.
CPU hours use Slurm TotalCPU when requested and available; otherwise realtime * %cpu / 100 is an estimate.
CPU allocation (allocated_core_hours) differs from consumed CPU. Metrics include monitoring/sampling uncertainty.
observed_span_hours spans earliest submission to latest completion, including pauses between resumes.
observed_active_hours is the union of task start/complete intervals, including job wrapper overhead.
cumulative_task_hours sums task runtime and includes parallel tasks; it is not pipeline elapsed time.
A fresh-run wall time/critical path cannot be reconstructed exactly from traces alone.
max_task_peak_rss_gib is the largest single-task memory peak, not the simultaneous memory of the pipeline.
Controller resources, untraced tasks and missing metrics are not included; consult coverage columns.
Plot selection tools uses annotation module names; thresholds selects each metric separately; all includes every measured module.
In thresholds mode one task exceeding the cutoff selects its whole module, including failures.
Zero thresholds include all measured modules. Missing measurements are never plotted as zero.
Each metric has its own Y scale (auto/linear/log/symlog). Auto switches at a 100-fold positive range, including requested limits.
Log scaling with zero values uses symlog to preserve those observations. The actual scale is recorded in plots_index.tsv.
Clean plot labels do not change process identities; module_labels.tsv maps labels to full process names.
Requested limits come from trace time/cpus/memory. With --show-requested each module maximum is overlaid on time/memory/cores plots.
No requested limit is drawn on consumed CPU hours: reserved CPU times walltime would only be a maximum budget.
Mean cores and efficiency use the matching CPU interval: Slurm job elapsed time for Slurm CPU, Nextflow realtime otherwise.
Do not mix different experiments or changed scientific inputs/options under one run_group for a clean benchmark.
