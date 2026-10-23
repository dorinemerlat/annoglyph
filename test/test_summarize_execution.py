import csv
import importlib.util
import tempfile
import unittest
from pathlib import Path

spec = importlib.util.spec_from_file_location('summary', Path(__file__).resolve().parents[1] / 'bin/summarize_execution.py')
m = importlib.util.module_from_spec(spec)
spec.loader.exec_module(m)


class ExecutionSummaryTests(unittest.TestCase):
    def row(self, **changes):
        row = dict(process='ROOT:SUB:MODULE', status='COMPLETED', hash='aa/bb', native_id='10',
                   attempt='1', realtime='1m 30s', cpus='2', peak_rss='1.5 GB', memory='4 GB',
                   start='2026-10-23 10:00:00', complete='2026-10-23 10:01:30', submit='2026-10-23 09:59:50', duration='1m 40s', time='2h')
        row['%cpu'] = '150%'
        row.update(changes)
        return row

    def test_units(self):
        self.assertEqual(m.seconds('1h 2m 3.5s'), 3723.5)
        self.assertEqual(m.seconds('250ms'), .25)
        self.assertEqual(m.seconds('250'), .25)
        self.assertEqual(m.memory_bytes('1.5 GB'), 1.5 * 1024 ** 3)
        self.assertEqual(m.memory_bytes('512K'), 512 * 1024)
        self.assertEqual(m.slurm_seconds('1-02:03:04.5'), 93784.5)
        self.assertIsNone(m.seconds('-'))
        with self.assertRaises(ValueError): m.seconds('unknown')

    def test_consumed_and_allocated_cpu_are_distinct(self):
        task = m.task_metrics(self.row(), {})
        self.assertEqual(task['cpu_seconds'], 135)
        self.assertEqual(task['allocated_core_seconds'], 180)
        self.assertEqual(task['mean_cpu_cores'], 1.5)
        self.assertEqual(task['cpu_efficiency_percent'], 75)
        summary = m.aggregate([task])
        self.assertEqual(summary['max_time_requested_hours'], 2)
        self.assertEqual(summary['max_cpus_requested'], 2)
        self.assertEqual(summary['max_memory_requested_gib'], 4)
        measured = m.task_metrics(self.row(), {'10': {'cpu_seconds': 120}})
        self.assertEqual(measured['cpu_seconds'], 120)
        self.assertEqual(measured['cpu_source'], 'slurm_totalcpu')
        measured = m.task_metrics(self.row(), {'10': {'cpu_seconds': 120, 'elapsed_seconds': 120, 'allocated_cpus': 2}})
        self.assertEqual(measured['mean_cpu_cores'], 1)
        self.assertEqual(measured['cpu_efficiency_percent'], 50)
        self.assertEqual(measured['allocated_core_seconds'], 240)

    def test_dedup_preserves_retries_and_ignores_cached_cpu(self):
        with tempfile.TemporaryDirectory() as directory:
            a, b = [Path(directory) / x for x in ['a.tsv', 'b.tsv']]
            first = self.row()
            retry = self.row(native_id='11', attempt='2', status='FAILED')
            cached = self.row(status='CACHED')
            for path, rows in [(a, [first, retry]), (b, [first, cached])]:
                with path.open('w') as handle:
                    writer = csv.DictWriter(handle, fieldnames=list(first), delimiter='\t')
                    writer.writeheader(); writer.writerows(rows)
            rows, cache, seen = m.read_traces([a, b])
            self.assertEqual((len(rows), len(cache), seen), (2, 1, 4))
            summary = m.aggregate([m.task_metrics(row, {}) for row in rows], len(cache))
            self.assertEqual(summary['executed_tasks'], 2)
            self.assertEqual(summary['cpu_hours'], 270 / 3600)
            self.assertEqual(summary['failed_tasks'], 1)

    def test_parallel_intervals_are_not_added_for_elapsed(self):
        a = m.task_metrics(self.row(), {})
        b = m.task_metrics(self.row(start='2026-10-23 10:01:00', complete='2026-10-23 10:02:00'), {})
        timing = m.observed_timing([a, b])
        self.assertAlmostEqual(timing['observed_active_hours'], 120 / 3600)
        self.assertAlmostEqual(timing['observed_span_hours'], 130 / 3600)

    def test_missing_metrics_stay_unknown(self):
        task = m.task_metrics(self.row(realtime='-', peak_rss='-'), {})
        self.assertIsNone(task['cpu_seconds'])
        self.assertIsNone(task['peak_rss_gib'])
        summary = m.aggregate([task])
        self.assertEqual(summary['cpu_measured_tasks'], 0)
        self.assertIsNone(summary['cpu_hours'])

    def test_labels_and_y_scales(self):
        self.assertEqual(m.module_label('ROOT:REPEATS:REPEATMASKER_REPEATMASKER'), 'RepeatMasker')
        self.assertEqual(m.module_label('ROOT:GENES:INTERPROSCAN'), 'InterProScan')
        self.assertEqual(m.y_axis_scale([1, 10]), ('log', None))
        self.assertEqual(m.y_axis_scale([.01, 10]), ('log', None))
        self.assertEqual(m.y_axis_scale([0, .01, 10]), ('symlog', .005))
        self.assertEqual(m.y_axis_scale([0, 0]), ('symlog', 1))


    def test_plot_selection_modes(self):
        tasks = [dict(process='ROOT:GENES:BRAKER2'), dict(process='ROOT:GENES:AGAT_STATISTICS'), dict(process='ROOT:REPEATS:REPEATMODELER_REPEATMODELER')]
        self.assertEqual(len(m.select_plot_tasks(tasks, 'tools', m.ANNOTATION_TOOLS)), 2)
        self.assertEqual(len(m.select_plot_tasks(tasks, 'all', m.ANNOTATION_TOOLS)), 3)
        self.assertEqual(len(m.select_plot_tasks(tasks, 'thresholds', m.ANNOTATION_TOOLS)), 3)
        self.assertEqual(m.select_plot_tasks(tasks, 'tools', ['AGAT*']), [tasks[1]])

    def test_plot_threshold_selects_whole_module(self):
        tasks = [dict(process='A', value=0, status='COMPLETED'), dict(process='A', value=10, status='FAILED'),
                 dict(process='B', value=2, status='COMPLETED'), dict(process='C', value=None, status='COMPLETED')]
        self.assertEqual(list(m.plot_groups(tasks, 'value', 1, 5)), ['A'])
        self.assertEqual(len(m.plot_groups(tasks, 'value', 1, 5)['A']), 2)
        self.assertEqual(list(m.plot_groups(tasks, 'value', 1, 0)), ['A', 'B'])
        self.assertEqual(m.plot_groups(tasks, 'value', 1, 10), {})


if __name__ == '__main__': unittest.main()
