"""Regression coverage for SRA ESummary records used by minimal integration tests."""
import importlib.util
from pathlib import Path
import unittest
import xml.etree.ElementTree as ET

spec = importlib.util.spec_from_file_location('generate_runlist', Path(__file__).resolve().parents[1] / 'bin/generate_runlist.py')
runlist = importlib.util.module_from_spec(spec)
spec.loader.exec_module(runlist)


def summary(runs):
    record = ET.Element('DocSum')
    ET.SubElement(record, 'Item', Name='ExpXml').text = (
        '<Experiment acc="SRX123"/><Study acc="SRP123" name="Study"/>'
        '<Sample acc="SRS123" name="Sample"/>'
        '<Summary><Platform instrument_model="Illumina NovaSeq">ILLUMINA</Platform></Summary>'
        '<LIBRARY_LAYOUT><PAIRED/></LIBRARY_LAYOUT>')
    ET.SubElement(record, 'Item', Name='Runs').text = runs
    return record


class SummaryRunsTests(unittest.TestCase):
    def test_single_run_root_is_retained_with_metadata(self):
        rows = runlist.parse_docsum(summary('<Run acc="SRR123" total_spots="100" total_bases="30000"/>'))
        self.assertEqual(len(rows), 1)
        self.assertEqual(rows[0]['Run_acc'], 'SRR123')
        self.assertEqual(rows[0]['experiment_accession'], 'SRX123')
        self.assertEqual(rows[0]['study_accession'], 'SRP123')
        self.assertEqual(rows[0]['sample_accession'], 'SRS123')
        self.assertEqual(rows[0]['bool:paired'], 1)
        self.assertEqual(rows[0]['instrument_model'], 'Illumina NovaSeq')
        self.assertEqual(rows[0]['avg_len'], 300)

    def test_multiple_sibling_runs_are_retained(self):
        rows = runlist.parse_docsum(summary('<Run acc="SRR1"/><Run acc="SRR2"/>'))
        self.assertEqual([r['Run_acc'] for r in rows], ['SRR1', 'SRR2'])

    def test_filters_apply_to_parsed_runs(self):
        rows = runlist.parse_docsum(summary(
            '<Run acc="SRR1" total_bases="500"/><Run acc="SRR2" total_bases="2000"/>'))
        selected = runlist.select_runs(rows, limit=3, max_bases=1000, only_paired=False)
        self.assertEqual([r['Run_acc'] for r in selected], ['SRR1'])


if __name__ == '__main__':
    unittest.main()
