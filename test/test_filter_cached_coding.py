import sys
import tempfile
import unittest
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parents[1] / 'bin'))
from filter_cached_coding import filter_models
from annoglyph_gff import read_gff


class CachedCodingTests(unittest.TestCase):
    def run_filter(self, root, text):
        source = root / 'source.gff3'
        source.write_text(text)
        result = filter_models(source, root / 'out.gff3', root / 'excluded.gff3',
                               root / 'report.tsv', {}, {b'a': (100, ''), b'b': (100, '')})
        return result, (root / 'out.gff3').read_text()

    def test_cross_contig_transcript_removed_and_valid_model_unchanged(self):
        with tempfile.TemporaryDirectory() as temp:
            root = Path(temp)
            valid = ('b\tAUGUSTUS\tgene\t5\t20\t.\t+\t.\tID=good\n'
                     'b\tAUGUSTUS\tmRNA\t5\t20\t.\t+\t.\tID=good.t;Parent=good\n'
                     'b\tAUGUSTUS\tCDS\t5\t20\t.\t+\t0\tID=good.c;Parent=good.t\n')
            bad = ('a\tAUGUSTUS\tgene\t1\t90\t.\t+\t.\tID=bad\n'
                   'a\tAUGUSTUS\tmRNA\t1\t90\t.\t+\t.\tID=bad.t;Parent=bad\n'
                   'b\tAUGUSTUS\tCDS\t5\t20\t.\t+\t0\tID=bad.c;Parent=bad.t\n')
            report, output = self.run_filter(root, bad + valid)
            self.assertEqual(output, valid)
            self.assertEqual(report['excluded_transcripts'], 1)
            self.assertIn('children_on_different_contigs', report['reasons']['bad.t'])
            self.assertIn(bad, (root / 'excluded.gff3').read_text())

    def test_valid_models_preserved_verbatim(self):
        with tempfile.TemporaryDirectory() as temp:
            text = 'a\tAUGUSTUS\tgene\t1\t99\t.\t+\t.\tID=g\n'
            report, output = self.run_filter(Path(temp), text)
            self.assertEqual(output, text)
            self.assertEqual(report['excluded_feature_rows'], 0)

    def test_outside_assembly_model_excluded_without_coordinate_clipping(self):
        with tempfile.TemporaryDirectory() as temp:
            text = ('a\tAUGUSTUS\tgene\t5\t150\t.\t+\t.\tID=g\n'
                    'a\tAUGUSTUS\tmRNA\t5\t150\t.\t+\t.\tID=t;Parent=g\n'
                    'a\tAUGUSTUS\tCDS\t5\t20\t.\t+\t0\tID=c;Parent=t\n'
                    'a\tAGAT\tthree_prime_UTR\t148\t150\t.\t+\t.\tID=u;Parent=t\n')
            report, output = self.run_filter(Path(temp), text)
            self.assertEqual(output, '')
            self.assertEqual(report['excluded_transcripts'], 1)
            self.assertEqual(report['excluded_feature_rows'], 4)

    def test_shared_cds_retained_for_valid_transcript(self):
        with tempfile.TemporaryDirectory() as temp:
            text = ('a\tAUGUSTUS\tgene\t1\t90\t.\t+\t.\tID=g\n'
                    'a\tAUGUSTUS\tmRNA\t1\t90\t.\t+\t.\tID=bad;Parent=g\n'
                    'a\tAUGUSTUS\tmRNA\t1\t90\t.\t+\t.\tID=good;Parent=g\n'
                    'a\tAUGUSTUS\tCDS\t5\t20\t.\t+\t0\tID=c;Parent=bad,good\n'
                    'b\tAUGUSTUS\tCDS\t5\t20\t.\t+\t0\tID=other;Parent=bad\n')
            report, output = self.run_filter(Path(temp), text)
            self.assertIn('ID=c;Parent=good', output)
            self.assertNotIn('ID=other', output)
            self.assertEqual(report['shared_feature_parent_references_updated'], 1)

    def test_unknown_contig_stops_import_instead_of_discarding_models(self):
        with tempfile.TemporaryDirectory() as temp:
            with self.assertRaisesRegex(ValueError, 'check the assembly/ID mapping'):
                self.run_filter(Path(temp), 'unknown\tAUGUSTUS\tgene\t1\t20\t.\t+\t.\tID=g\n')

    def test_reversed_legacy_intron_excludes_model_but_normal_gff_reader_is_strict(self):
        with tempfile.TemporaryDirectory() as temp:
            root = Path(temp)
            text = ('a\tAUGUSTUS\tgene\t1\t90\t.\t+\t.\tID=g\n'
                    'a\tAUGUSTUS\tmRNA\t1\t90\t.\t+\t.\tID=t;Parent=g\n'
                    'a\tAUGUSTUS\tintron\t21\t20\t.\t+\t.\tID=i;Parent=t\n')
            report, output = self.run_filter(root, text)
            self.assertEqual(output, '')
            self.assertEqual(report['excluded_transcripts'], 1)
            with self.assertRaisesRegex(ValueError, 'invalid coordinates'):
                read_gff(root / 'source.gff3')


if __name__ == '__main__':
    unittest.main()
