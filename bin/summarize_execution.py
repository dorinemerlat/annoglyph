#!/usr/bin/env python3
"""Unify Nextflow traces across resumed runs without double-counting tasks."""
import argparse
import csv
import fnmatch
import json
import re
import random
import subprocess
import statistics
from collections import defaultdict
from datetime import datetime
from pathlib import Path

MISSING = {'', '-', 'NA', 'null', 'None'}


def number(value):
    return None if str(value).strip() in MISSING else float(str(value).rstrip('%'))


def seconds(value):
    value = str(value).strip()
    if value in MISSING:
        return None
    # Nextflow raw trace durations are milliseconds.
    if re.fullmatch(r'\d+(?:\.\d+)?', value):
        return float(value) / 1000
    units = {'d': 86400, 'h': 3600, 'm': 60, 's': 1, 'ms': .001}
    matches = list(re.finditer(r'(\d+(?:\.\d+)?)\s*(ms|d|h|m|s)', value))
    if not matches or re.sub(r'(\d+(?:\.\d+)?)\s*(ms|d|h|m|s)|\s+', '', value):
        raise ValueError(f'Invalid Nextflow duration: {value}')
    return sum(float(m[1]) * units[m[2]] for m in matches)


def memory_bytes(value):
    value = str(value).strip()
    if value in MISSING:
        return None
    match = re.fullmatch(r'([\d.]+)\s*([KMGTPE]?)(?:i?B)?', value, re.I)
    if not match:
        raise ValueError(f'Invalid memory value: {value}')
    return float(match[1]) * 1024 ** (' KMGTPE'.index(match[2].upper()) if match[2] else 0)


def timestamp(value):
    value = str(value).strip()
    if value in MISSING:
        return None
    if re.fullmatch(r'\d+', value):
        return datetime.fromtimestamp(int(value) / 1000).timestamp()
    return datetime.fromisoformat(value).timestamp()


def slurm_seconds(value):
    days, _, clock = value.rpartition('-')
    total = float(days or 0) * 86400
    for index, part in enumerate(reversed(clock.split(':'))):
        total += float(part) * 60 ** index
    return total


def load_accounting(rows):
    ids = sorted({r.get('native_id', '') for r in rows if re.fullmatch(r'\d+', r.get('native_id', '')) and r['status'] != 'CACHED'})
    result = {}
    for offset in range(0, len(ids), 200):
        output = subprocess.run(['sacct', '-j', ','.join(ids[offset:offset + 200]), '--parsable2', '--noheader', '--units=K', '--format=JobIDRaw,TotalCPU,ElapsedRaw,AllocCPUS,MaxRSS'], check=True, text=True, capture_output=True).stdout
        for line in output.splitlines():
            job, cpu, elapsed, cpus, rss = line.split('|')[:5]
            if job in ids:
                result.setdefault(job, {}).update(cpu_seconds=slurm_seconds(cpu), elapsed_seconds=float(elapsed), allocated_cpus=float(cpus))
            elif job.endswith('.batch') and job[:-6] in ids:
                result.setdefault(job[:-6], {})['peak_rss_bytes'] = memory_bytes(rss)
    return result


def read_traces(paths):
    unique, cached, seen = {}, [], 0
    for path in sorted(paths):
        with path.open() as handle:
            reader = csv.DictReader(handle, delimiter='\t')
            required = {'process', 'status', 'hash', 'realtime', 'cpus', '%cpu', 'peak_rss'}
            if not required.issubset(reader.fieldnames or []):
                raise ValueError(f'{path}: missing columns {sorted(required - set(reader.fieldnames or []))}')
            for line, row in enumerate(reader, 2):
                if None in row or any(value is None for value in row.values()):
                    raise ValueError(f'{path}:{line}: incomplete trace row')
                seen += 1
                if row['status'] == 'CACHED':
                    cached.append(row)
                    continue
                # Native job IDs identify executions, while task_id resets on resume.
                native = row.get('native_id', '')
                if native not in MISSING:
                    key = (row['process'], row['hash'], native, row.get('start', ''), row.get('attempt', ''))
                else:
                    key = (row['process'], row['hash'], row.get('start', ''), row.get('submit', ''), row.get('attempt', ''))
                    if not any(v not in MISSING for v in key[2:4]):
                        key += (str(path), line)  # Insufficient identity: do not silently collapse tasks.
                if key in unique:
                    previous = unique[key]
                    # Prefer a final/completed record and fill missing metrics from other traces.
                    rank = {'COMPLETED': 4, 'FAILED': 3, 'ABORTED': 2}
                    if rank.get(row['status'], 0) > rank.get(previous['status'], 0):
                        previous, row = row, previous
                    unique[key] = {k: row.get(k, '') if v in MISSING else v for k, v in previous.items()}
                else:
                    unique[key] = row
    return list(unique.values()), cached, seen


def task_metrics(row, accounting):
    runtime = seconds(row.get('realtime', '-'))
    cpus = number(row.get('cpus', '-'))
    cpu_percent = number(row.get('%cpu', '-'))
    peak = memory_bytes(row.get('peak_rss', '-'))
    slurm = accounting.get(row.get('native_id'), {})
    if 'cpu_seconds' in slurm:
        cpu, source = slurm['cpu_seconds'], 'slurm_totalcpu'
    elif runtime is not None and cpu_percent is not None:
        cpu, source = runtime * cpu_percent / 100, 'nextflow_percent_cpu_estimate'
    else:
        cpu, source = None, 'unavailable'
    if peak is None:
        peak = slurm.get('peak_rss_bytes')
    cpu_interval = slurm.get('elapsed_seconds', runtime)
    requested_cpus = cpus
    cpus = slurm.get('allocated_cpus', cpus)
    return dict(process=row['process'], name=row.get('name', ''), tag=row.get('tag', ''), hash=row['hash'],
                native_id=row.get('native_id', ''), status=row['status'], attempt=row.get('attempt', ''),
                runtime_seconds=runtime, time_requested_seconds=seconds(row.get('time', '-')), submit_to_complete_seconds=seconds(row.get('duration', '-')),
                cpus_requested=requested_cpus, cpus_allocated=cpus, cpu_seconds=cpu, cpu_source=source,
                cpu_interval_seconds=cpu_interval, scheduler_elapsed_seconds=slurm.get('elapsed_seconds'),
                allocated_core_seconds=cpu_interval * cpus if cpu_interval is not None and cpus is not None else None,
                mean_cpu_cores=cpu / cpu_interval if cpu is not None and cpu_interval else None,
                cpu_efficiency_percent=100 * cpu / (cpu_interval * cpus) if cpu is not None and cpu_interval and cpus else None,
                peak_rss_gib=peak / 1024 ** 3 if peak is not None else None,
                memory_requested_gib=(memory_bytes(row.get('memory', '-')) or 0) / 1024 ** 3 if row.get('memory', '-') not in MISSING else None,
                start=row.get('start', ''), complete=row.get('complete', ''), submit=row.get('submit', ''))


def aggregate(tasks, cached_count=0):
    known_cpu = [t['cpu_seconds'] for t in tasks if t['cpu_seconds'] is not None]
    allocated = [t['allocated_core_seconds'] for t in tasks if t['allocated_core_seconds'] is not None]
    peaks = [t['peak_rss_gib'] for t in tasks if t['peak_rss_gib'] is not None]
    runtimes = [t['runtime_seconds'] for t in tasks if t['runtime_seconds'] is not None]
    paired = [t for t in tasks if t['cpu_seconds'] is not None and t['cpu_interval_seconds'] is not None]
    paired_allocated = [t for t in paired if t['allocated_core_seconds'] is not None]
    paired_runtime = sum(t['cpu_interval_seconds'] for t in paired)
    paired_core_seconds = sum(t['allocated_core_seconds'] for t in paired_allocated)
    return dict(executed_tasks=len(tasks), completed_tasks=sum(t['status'] == 'COMPLETED' for t in tasks),
                failed_tasks=sum(t['status'] == 'FAILED' for t in tasks), aborted_tasks=sum(t['status'] == 'ABORTED' for t in tasks),
                cached_trace_records=cached_count, runtime_measured_tasks=len(runtimes),
                cumulative_task_hours=sum(runtimes) / 3600,
                runtime_mean_minutes=statistics.mean(runtimes) / 60 if runtimes else None,
                runtime_median_minutes=statistics.median(runtimes) / 60 if runtimes else None,
                runtime_max_minutes=max(runtimes) / 60 if runtimes else None,
                mean_cpu_cores=sum(t['cpu_seconds'] for t in paired) / paired_runtime if paired_runtime else None,
                cpu_efficiency_percent=100 * sum(t['cpu_seconds'] for t in paired_allocated) / paired_core_seconds if paired_core_seconds else None,
                cpu_measured_tasks=len(known_cpu),
                cpu_hours=sum(known_cpu) / 3600 if known_cpu else None,
                failed_or_aborted_cpu_hours=sum(t['cpu_seconds'] for t in tasks if t['status'] in ('FAILED', 'ABORTED') and t['cpu_seconds'] is not None) / 3600,
                cpu_slurm_tasks=sum(t['cpu_source'] == 'slurm_totalcpu' for t in tasks),
                cpu_estimated_tasks=sum(t['cpu_source'] == 'nextflow_percent_cpu_estimate' for t in tasks),
                allocated_core_hours=sum(allocated) / 3600 if allocated else None,
                memory_measured_tasks=len(peaks), max_task_peak_rss_gib=max(peaks) if peaks else None,
                max_time_requested_hours=max((t['time_requested_seconds'] / 3600 for t in tasks if t['time_requested_seconds'] is not None), default=None),
                max_cpus_requested=max((t['cpus_requested'] for t in tasks if t['cpus_requested'] is not None), default=None),
                max_memory_requested_gib=max((t['memory_requested_gib'] for t in tasks if t['memory_requested_gib'] is not None), default=None),
                max_cpus_per_task=max((t['cpus_allocated'] for t in tasks if t['cpus_allocated'] is not None), default=None))


def observed_timing(tasks):
    intervals, bounds = [], []
    for task in tasks:
        start, end, submit = (timestamp(task[k]) for k in ('start', 'complete', 'submit'))
        if submit is not None:
            bounds.append(submit)
        if end is not None:
            bounds.append(end)
        if start is not None and end is not None and end >= start:
            intervals.append((start, end))
    merged = []
    for start, end in sorted(intervals):
        if merged and start <= merged[-1][1]:
            merged[-1][1] = max(merged[-1][1], end)
        else:
            merged.append([start, end])
    return dict(observed_span_hours=(max(bounds) - min(bounds)) / 3600 if len(bounds) >= 2 else None,
                observed_active_hours=sum(end - start for start, end in merged) / 3600,
                tasks_with_complete_intervals=len(intervals))


def write_table(path, rows):
    with path.open('w') as handle:
        writer = csv.DictWriter(handle, fieldnames=list(rows[0]), delimiter='\t', lineterminator='\n')
        writer.writeheader()
        for row in rows:
            writer.writerow({k: '' if v is None else round(v, 6) if isinstance(v, float) else v for k, v in row.items()})


PLOT_METRICS = {
    'runtime': ('runtime_seconds', 1 / 60, 'Task runtime (minutes)', 30.0),
    'cpu': ('cpu_seconds', 1 / 3600, 'Consumed CPU (CPU hours)', 10.0),
    'memory': ('peak_rss_gib', 1, 'Peak resident memory (GiB)', 8.0),
    'cores': ('mean_cpu_cores', 1, 'Mean CPU use (equivalent cores)', 4.0),
}


def plot_groups(tasks, field, factor, threshold):
    grouped = defaultdict(list)
    for task in tasks:
        if task.get(field) is not None:
            grouped[task['process']].append((task[field] * factor, task['status']))
    return {process: observations for process, observations in sorted(grouped.items())
            if threshold == 0 or any(value > threshold for value, status in observations)}


ANNOTATION_TOOLS = (
    'REPEATMASKER', 'REPEATMASKER_REPEATMASKER', 'REPEATMODELER', 'REPEATMODELER_REPEATMODELER',
    'MCHELPER', 'GALBA', 'BRAKER2', 'BRAKER3', 'BLASTP', 'BLAST_BLASTP', 'BLAST_BLASTX',
    'TOXPROT_BLASTP', 'INTERPROSCAN', 'EGGNOGMAPPER', 'EGGNOG_MAPPER', 'SIGNALP',
    'TRNASCAN*', 'BARRNAP', 'BARNNAP', 'INFERNAL', 'INFERNAL_CM*', 'CMSCAN', 'CMSEARCH',
)


def select_plot_tasks(tasks, selection, tools):
    if selection != 'tools':
        return tasks
    return [task for task in tasks if any(fnmatch.fnmatchcase(task['process'].split(':')[-1], pattern)
                                          or fnmatch.fnmatchcase(task['process'], pattern) for pattern in tools)]


MODULE_LABELS = {
    'REPEATMASKER_REPEATMASKER': 'RepeatMasker', 'REPEATMASKER': 'RepeatMasker',
    'REPEATMODELER_REPEATMODELER': 'RepeatModeler', 'REPEATMODELER': 'RepeatModeler',
    'REPEATMODELER_BUILDDATABASE': 'RepeatModeler (database)',
    'MCHELPER': 'MCHelper', 'GALBA': 'GALBA', 'BRAKER2': 'BRAKER2', 'BRAKER3': 'BRAKER3',
    'BLASTP': 'BLASTp', 'BLAST_BLASTP': 'BLASTp', 'BLAST_BLASTX': 'BLASTx',
    'TOXPROT_BLASTP': 'BLASTp (ToxProt)', 'INTERPROSCAN': 'InterProScan',
    'EGGNOGMAPPER': 'eggNOG-mapper', 'EGGNOG_MAPPER': 'eggNOG-mapper', 'SIGNALP': 'SignalP',
    'TRNASCANSE': 'tRNAscan-SE', 'BARRNAP': 'Barrnap', 'BARNNAP': 'Barrnap',
    'INFERNAL': 'Infernal', 'INFERNAL_CMSCAN': 'Infernal (cmscan)', 'INFERNAL_CMSEARCH': 'Infernal (cmsearch)',
    'BUSCO_BUSCO': 'BUSCO', 'BUSCO_PROTEINS': 'BUSCO (proteins)', 'OMARK': 'OMArk',
}
HIERARCHY_LABELS = {
    'GENES_ANNOTATION': 'Protein-coding genes', 'FUNCTIONAL_ANNOTATION': 'Functional annotation',
    'STRUCTURAL_ANNOTATION': 'Gene prediction', 'NCGENES_ANNOTATION': 'Noncoding RNA',
    'REPEATS_ANNOTATION': 'Repeats', 'PREPROCESSING': 'Preprocessing', 'POSTPROCESSING': 'Postprocessing',
    'RNASEQ_ASSEMBLY': 'RNA-seq', 'ANNOTATION_QC_AND_CONTAMINATION_FILTERING': 'Quality control',
    'ANNOTATION_COMPARISON_AND_SELECTION': 'Comparison and selection',
}


def module_label(process):
    name = process.split(':')[-1]
    return MODULE_LABELS.get(name, name.replace('_', ' '))


def y_axis_scale(values):
    positive = [value for value in values if value > 0]
    if not positive:
        return 'symlog', 1
    if min(values) <= 0:
        return 'symlog', min(positive) / 2
    return 'log', None


def make_plots(tasks, output, thresholds, selection, tools, show_requested=False):
    try:
        import matplotlib
        matplotlib.use('Agg')
        import matplotlib.pyplot as plt
        from matplotlib.lines import Line2D
        import seaborn as sns
    except ImportError as exc:
        raise RuntimeError('--plots requires seaborn and matplotlib') from exc
    output.mkdir(exist_ok=True)
    tasks = select_plot_tasks(tasks, selection, tools)
    # A union across metrics gives every panel the same module order and groups.
    selected = set()
    for metric, (field, factor, label, default) in PLOT_METRICS.items():
        selected.update(plot_groups(tasks, field, factor, thresholds[metric] if selection == 'thresholds' else 0))
        for extension in ('png', 'svg'):
            (output / f'{metric}.{extension}').unlink(missing_ok=True)
    processes = sorted(selected)
    write_table(output / 'module_labels.tsv', [dict(process=p, label=module_label(p)) for p in processes])
    if not processes:
        for extension in ('png', 'svg'):
            (output / f'execution_metrics.{extension}').unlink(missing_ok=True)
        write_table(output / 'plots_index.tsv', [dict(metric=m, selection=selection, threshold=thresholds[m] if selection == 'thresholds' else 0, scale='', modules=0, status='no measured module selected') for m in PLOT_METRICS])
        return
    # Empty categorical slots separate subworkflows without breaking the shared X axis.
    order, blocks, previous = [], [], None
    for process in processes:
        parts = process.split(':')
        hierarchy = ' / '.join(parts[1:-1]) if len(parts) > 2 else ' / '.join(parts[:-1]) or 'Processes'
        if hierarchy != previous:
            if previous is not None:
                order.append(f'__gap_{len(blocks)}')
            blocks.append([hierarchy, len(order), len(order)])
        blocks[-1][2] = len(order)
        order.append(process)
        previous = hierarchy
    sns.set_theme(style='whitegrid', context='paper')
    fig, axes = plt.subplots(4, 1, sharex=True, figsize=(max(11, len(order) * .5), 13))
    index = []
    titles = {'runtime': 'Task runtime', 'cpu': 'Consumed CPU', 'memory': 'Peak resident memory', 'cores': 'Mean CPU use'}
    units = {'runtime': 'Minutes', 'cpu': 'CPU hours', 'memory': 'GiB', 'cores': 'Equivalent CPU cores'}
    for ax, (metric, (field, factor, label, default)) in zip(axes, PLOT_METRICS.items()):
        groups = plot_groups(tasks, field, factor, 0)
        groups = {p: groups[p] for p in processes if p in groups}
        observations = [(p, value, status) for p, values in groups.items() for value, status in values]
        axis_values = [value for p, value, status in observations]
        if observations:
            sns.boxplot(x=[p for p, v, s in observations], y=[v for p, v, s in observations],
                        order=order, ax=ax, width=.5, showfliers=False, color='#d4e4ec', linewidth=1)
            for completed, color, marker in [(True, '#28759d', 'o'), (False, '#b44435', 'x')]:
                subset = [(p, v) for p, v, status in observations if (status == 'COMPLETED') == completed]
                if subset:
                    sns.stripplot(x=[p for p, v in subset], y=[v for p, v in subset], order=order,
                                  ax=ax, jitter=False, color=color, marker=marker, size=4, alpha=.7, zorder=3)
        if show_requested and metric != 'cpu':
            requested_field, requested_factor = {
                'runtime': ('time_requested_seconds', 1 / 60),
                'memory': ('memory_requested_gib', 1), 'cores': ('cpus_requested', 1),
            }[metric]
            for process in processes:
                limits = [t[requested_field] * requested_factor for t in tasks
                          if t['process'] == process and t.get(requested_field) is not None]
                if limits:
                    limit = max(limits)
                    axis_values.append(limit)
                    position = order.index(process)
                    ax.hlines(limit, position - .3, position + .3, color='#9857a0', linestyle='--', linewidth=1.5)
        for block_index, (hierarchy, start, end) in enumerate(blocks):
            ax.axvspan(start - .5, end + .5, color='#28759d', alpha=.04 if block_index % 2 == 0 else .08)
        scale, linthresh = y_axis_scale(axis_values)
        ax.set_yscale(scale, **({'linthresh': linthresh} if scale == 'symlog' else {}))
        if scale == 'symlog':
            ax.set_ylim(bottom=0)
        ax.set_ylabel(titles[metric] + '\n' + units[metric], fontweight='bold', rotation=90, labelpad=12)
        ax.set_xlabel('')
        ax.grid(axis='x', visible=False)
        if selection == 'thresholds' and thresholds[metric] > 0:
            ax.axhline(thresholds[metric], color='gray', linestyle='--', linewidth=.8)
        index.append(dict(metric=metric, selection=selection, threshold=thresholds[metric] if selection == 'thresholds' else 0,
                          scale=scale, modules=len(processes), measured_modules=len(groups), status='plotted', figure='execution_metrics'))
    for hierarchy, start, end in blocks:
        axes[0].text((start + end) / 2, 1.05,
                     '\n'.join(HIERARCHY_LABELS.get(part, part.replace('_', ' ').capitalize()) for part in hierarchy.split(' / ')),
                     ha='center', va='bottom', transform=axes[0].get_xaxis_transform(), fontsize=8, rotation=25)
    axes[-1].set_xticks(range(len(order)))
    axes[-1].set_xticklabels([module_label(p) if p in selected else '' for p in order], rotation=65, ha='right', fontsize=8)
    axes[-1].set_xlabel('Modules, grouped by workflow hierarchy')
    handles = [Line2D([], [], color='#28759d', marker='o', linestyle='', label='Completed'),
               Line2D([], [], color='#b44435', marker='x', linestyle='', label='Failed / aborted / other')]
    if show_requested:
        handles.append(Line2D([], [], color='#9857a0', linestyle='--', label='Maximum requested limit'))
    fig.legend(handles=handles, loc='upper center', ncol=len(handles), bbox_to_anchor=(.5, .995))
    fig.subplots_adjust(bottom=.21, top=.86, left=.11, right=.98, hspace=.13)
    for extension in ('png', 'svg'):
        fig.savefig(output / f'execution_metrics.{extension}', dpi=180)
    plt.close(fig)
    write_table(output / 'plots_index.tsv', index)


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('traces', nargs='*', type=Path, help='Explicit trace TSV files, including older unnamed traces.')
    parser.add_argument('--trace-dir', type=Path, help='Directory containing execution_trace_*.tsv files.')
    parser.add_argument('--run-group', help='Exact run_group used when launching Nextflow.')
    parser.add_argument('--output', type=Path, default=Path('execution_stats'))
    parser.add_argument('--slurm-accounting', action='store_true', help='Query sacct for consumed CPU; otherwise estimate from Nextflow %%cpu.')
    parser.add_argument('--plots', action='store_true', help='Generate a shared-X PNG/SVG metric figure (requires seaborn and matplotlib).')
    parser.add_argument('--show-requested', action='store_true', help='Overlay maximum requested time, CPU and memory per module on corresponding plots.')
    parser.add_argument('--plot-selection', choices=['tools', 'thresholds', 'all'], default='tools', help='Tools: annotation tool list; thresholds: resource thresholds; all: every measured module.')
    parser.add_argument('--plot-tool', action='append', help='Module name or glob; repeat to replace the default annotation tool list in tools mode.')
    for metric, (_, _, _, default) in PLOT_METRICS.items():
        parser.add_argument(f'--{metric}-threshold', type=float, default=default, help=f'Plot selection threshold in the metric unit (default {default}); 0 includes all measured modules.')
    args = parser.parse_args()
    thresholds = {metric: getattr(args, metric + '_threshold') for metric in PLOT_METRICS}
    if any(value < 0 for value in thresholds.values()):
        parser.error('Plot thresholds must be nonnegative')
    if args.run_group and not re.fullmatch(r'[A-Za-z0-9][A-Za-z0-9_.-]*', args.run_group):
        parser.error('run-group must start with a letter/digit and contain only letters, digits, _, . or -')
    paths = set(args.traces)
    if args.trace_dir:
        candidates = args.trace_dir.glob('execution_trace_*.tsv')
        # The date boundary avoids mixing e.g. group a with group a_extra.
        pattern = re.compile(r'^execution_trace_' + re.escape(args.run_group) + r'_\d{8}-\d{6}(?:-\d+)?\.tsv$') if args.run_group else None
        paths.update(p for p in candidates if pattern is None or pattern.match(p.name))
    if not paths:
        parser.error('No trace files selected')
    rows, cached, seen = read_traces(paths)
    if not rows and not cached:
        parser.error('Selected traces contain no tasks')
    accounting = load_accounting(rows) if args.slurm_accounting else {}
    tasks = [task_metrics(row, accounting) for row in rows]
    grouped = defaultdict(list)
    for task in tasks:
        grouped[task['process']].append(task)
    cached_by_process = defaultdict(int)
    for row in cached:
        cached_by_process[row['process']] += 1
    summary = dict(run_group=args.run_group or '', trace_files=len(paths), input_trace_records=seen,
                   duplicate_execution_records=seen - len(cached) - len(tasks), **aggregate(tasks, len(cached)), **observed_timing(tasks))
    args.output.mkdir(parents=True, exist_ok=True)
    write_table(args.output / 'execution_summary.tsv', [summary])
    process_rows = [dict(process=process, **aggregate(grouped[process], cached_by_process[process])) for process in sorted(set(grouped) | set(cached_by_process))]
    write_table(args.output / 'process_summary.tsv', process_rows)
    if tasks:
        write_table(args.output / 'tasks.tsv', tasks)
    else:
        (args.output / 'tasks.tsv').unlink(missing_ok=True)
    if args.plots:
        make_plots(tasks, args.output / 'plots', thresholds, args.plot_selection, args.plot_tool or ANNOTATION_TOOLS, args.show_requested)
    notes = [
        'Column suffixes specify hours, minutes, seconds or GiB. Missing values are blank, not zero.',
        'Only actual task executions contribute to CPU/runtime totals; CACHED records are counted separately.',
        'storeDir-skipped processes may be absent entirely: supply earlier traces to recover their measurements.',
        'Repeated records of one execution are deduplicated; actual retries and recomputations remain included.',
        'CPU hours use Slurm TotalCPU when requested and available; otherwise realtime * %cpu / 100 is an estimate.',
        'CPU allocation (allocated_core_hours) differs from consumed CPU. Metrics include monitoring/sampling uncertainty.',
        'observed_span_hours spans earliest submission to latest completion, including pauses between resumes.',
        'observed_active_hours is the union of task start/complete intervals, including job wrapper overhead.',
        'cumulative_task_hours sums task runtime and includes parallel tasks; it is not pipeline elapsed time.',
        'A fresh-run wall time/critical path cannot be reconstructed exactly from traces alone.',
        'max_task_peak_rss_gib is the largest single-task memory peak, not the simultaneous memory of the pipeline.',
        'Controller resources, untraced tasks and missing metrics are not included; consult coverage columns.',
        'Plot selection tools uses annotation module names; thresholds selects each metric separately; all includes every measured module.',
        'In thresholds mode one task exceeding the cutoff selects its whole module, including failures.',
        'Zero thresholds include all measured modules. Missing measurements are never plotted as zero.',
        'Y axes use logarithmic scaling, with only the metric unit in the axis label.',
        'Log scaling with zero values uses symlog to preserve those observations. The actual scale is recorded in plots_index.tsv.',
        'Clean plot labels do not change process identities; module_labels.tsv maps labels to full process names.',
        'Requested limits come from trace time/cpus/memory. With --show-requested each module maximum is overlaid on time/memory/cores plots.',
        'No requested limit is drawn on consumed CPU hours: reserved CPU times walltime would only be a maximum budget.',
        'Mean cores and efficiency use the matching CPU interval: Slurm job elapsed time for Slurm CPU, Nextflow realtime otherwise.',
        'Do not mix different experiments or changed scientific inputs/options under one run_group for a clean benchmark.',
    ]
    (args.output / 'summary_metadata.json').write_text(json.dumps(dict(traces=[str(p.resolve()) for p in sorted(paths)], plot_selection=args.plot_selection if args.plots else None, show_requested=args.show_requested,  plot_tools=args.plot_tool or list(ANNOTATION_TOOLS), plot_thresholds=thresholds, notes=notes), indent=2) + '\n')
    (args.output / 'README.txt').write_text('\n'.join(notes) + '\n')
    display = lambda value: 'unavailable' if value is None else f'{value:.2f}' if isinstance(value, float) else str(value)
    report = [f'# Execution statistics: {args.run_group or "selected traces"}', '',
              'Includes actual executions across selected traces, including failed attempts and recomputations.', '',
              '| Metric | Value |', '| --- | ---: |']
    for label, key in [('Trace files', 'trace_files'), ('Task executions', 'executed_tasks'),
                       ('Completed tasks', 'completed_tasks'), ('Failed tasks', 'failed_tasks'), ('Aborted tasks', 'aborted_tasks'),
                       ('Consumed CPU hours', 'cpu_hours'), ('Failed/aborted CPU hours', 'failed_or_aborted_cpu_hours'),
                       ('Observed span, hours (includes resume gaps)', 'observed_span_hours'),
                       ('Observed active hours (union of task intervals)', 'observed_active_hours'),
                       ('Maximum single-task RSS, GiB', 'max_task_peak_rss_gib')]:
        report.append(f'| {label} | {display(summary[key])} |')
    report += ['', f'CPU coverage: {summary["cpu_measured_tasks"]}/{len(tasks)} executions; '
               f'{summary["cpu_slurm_tasks"]} from Slurm, {summary["cpu_estimated_tasks"]} estimated. '
               f'Runtime coverage: {summary["runtime_measured_tasks"]}/{len(tasks)}.', '',
               '[Global table](execution_summary.tsv) · [Process table](process_summary.tsv) · [Individual tasks](tasks.tsv)', '',
               'These timings are not a reconstruction of an uninterrupted fresh run. '
               'Processes skipped through storeDir require their original traces. '
               'Controller resources are excluded. See [measurement notes](README.txt).', '',
               '## Processes', '', '| Process | Executions | CPU hours | Max runtime, min | Peak RSS, GiB | Max requested time, h | Max requested CPUs | Max requested memory, GiB |',
               '| --- | ---: | ---: | ---: | ---: | ---: | ---: | ---: |']
    for row in sorted(process_rows, key=lambda row: row['cpu_hours'] or 0, reverse=True):
        report.append('| ' + ' | '.join(display(row[key]) for key in ['process', 'executed_tasks', 'cpu_hours', 'runtime_max_minutes', 'max_task_peak_rss_gib', 'max_time_requested_hours', 'max_cpus_requested', 'max_memory_requested_gib']) + ' |')
    report += ['', '## Figures', '', 'When enabled: [execution metrics](plots/execution_metrics.png) '
               '([SVG](plots/execution_metrics.svg)), with four panels sharing the same X axis. '
               'See plots_index.tsv for metric coverage.', '']
    (args.output / 'README.md').write_text('\n'.join(report))

    print(f'{len(tasks)} executions, {len(cached)} cached records, {len(paths)} traces -> {args.output}')


if __name__ == '__main__':
    main()
