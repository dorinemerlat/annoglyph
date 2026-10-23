import importlib.util
import json
import sqlite3
import tempfile
import unittest
from pathlib import Path

spec = importlib.util.spec_from_file_location('filter_evidence', Path(__file__).resolve().parents[1] / 'bin/filter_protein_evidence.py')
m = importlib.util.module_from_spec(spec)
spec.loader.exec_module(m)


class EvidenceFilteringTests(unittest.TestCase):
    def setUp(self):
        self.tmp = tempfile.TemporaryDirectory()
        self.addCleanup(self.tmp.cleanup)
        self.root = Path(self.tmp.name)
        self.db = self.root / 'taxa.sqlite'
        with sqlite3.connect(self.db) as db:
            db.execute('CREATE TABLE species(taxid INT PRIMARY KEY, spname TEXT, rank TEXT, track TEXT)')
            db.executemany('INSERT INTO species VALUES(?,?,?,?)', [
                (10, 'GenusA', 'genus', '10,1'), (11, 'GenusA species', 'species', '11,10,1'),
                (12, 'GenusA other', 'species', '12,10,1'),
                (20, 'GenusB', 'genus', '20,1'), (21, 'GenusB species', 'species', '21,20,1'),
            ])
        self.fasta = self.root / 'evidence.fa'
        self.fasta.write_text('>sp|A|A OX=11\nMAAA\n>12_0:gene\nMCCC\n>sp|B|B OX=21\nMDDD\n>21_0:gene\nMEEE\n>unknown\nMFFF\n>999_0:gene\nMGGG\n')

    def run_filter(self, genus):
        output = self.root / (genus + '.fa')
        report = self.root / (genus + '.json')
        counts = m.filter_evidence(self.fasta, m.Taxonomy(self.db), genus, output, report)
        return output.read_text(), counts, json.loads(report.read_text())

    def test_excludes_entire_genus_from_both_sources_and_unknowns(self):
        fasta, counts, report = self.run_filter('GenusA')
        self.assertEqual(counts, dict(total=6, retained=2, excluded_target_genus=2, excluded_unresolved_taxonomy=2))
        self.assertIn('OX=21', fasta)
        self.assertIn('21_0:gene', fasta)
        self.assertNotIn('OX=11', fasta)
        self.assertNotIn('12_0:gene', fasta)
        self.assertEqual(report['excluded_genus'], 'GenusA')

    def test_different_genomes_get_different_evidence_and_source_unchanged(self):
        before = self.fasta.read_bytes()
        a, _, _ = self.run_filter('GenusA')
        b, _, _ = self.run_filter('GenusB')
        self.assertNotEqual(a, b)
        self.assertIn('OX=11', b)
        self.assertNotIn('OX=21', b)
        self.assertEqual(before, self.fasta.read_bytes())

    def test_empty_result_fails_with_report(self):
        self.fasta.write_text('>11_0:gene\nMAAA\n')
        report = self.root / 'empty.json'
        with self.assertRaisesRegex(ValueError, 'No usable protein evidence'):
            m.filter_evidence(self.fasta, m.Taxonomy(self.db), 'GenusA', self.root / 'empty.fa', report)
        self.assertEqual(json.loads(report.read_text())['retained'], 0)

    def test_taxid_formats(self):
        self.assertEqual(m.protein_taxid('sp|P12345|NAME OS=GenusA species OX=11 GN=x'), 11)
        self.assertEqual(m.protein_taxid('21_0:000abc'), 21)
        self.assertIsNone(m.protein_taxid('arbitrary_user_header'))


if __name__ == '__main__':
    unittest.main()
