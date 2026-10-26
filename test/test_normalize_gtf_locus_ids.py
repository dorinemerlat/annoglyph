import sys
from pathlib import Path
import unittest

sys.path.insert(0, str(Path(__file__).resolve().parents[1] / 'bin'))
from normalize_gtf_locus_ids import normalize_lines


class LocusIdentifierTests(unittest.TestCase):
    def row(self, sequence, strand, gene='g1', transcript='g1.t1'):
        return f'{sequence}\tAUGUSTUS\tCDS\t1\t9\t.\t{strand}\t0\ttranscript_id "{transcript}"; gene_id "{gene}";\n'

    def test_unique_identifiers_are_unchanged(self):
        rows = ['chr1\tAUGUSTUS\tgene\t1\t9\t.\t+\t.\tg1\n',
                'chr1\tAUGUSTUS\ttranscript\t1\t9\t.\t+\t.\tg1.t1\n',
                self.row('chr1', '+')]
        self.assertEqual(normalize_lines(rows)[-1], rows[-1])

    def test_contig_collisions_are_disambiguated_without_changing_coordinates(self):
        rows = [self.row('chr1', '+'), self.row('chr2', '+')]
        output = [row for row in normalize_lines(rows) if '\tCDS\t' in row]
        self.assertIn('gene_id "g1_locus1"', output[0])
        self.assertIn('transcript_id "g1_locus2.t1"', output[1])
        self.assertEqual([row.split('\t')[:8] for row in output], [row.split('\t')[:8] for row in rows])
        self.assertEqual(normalize_lines(normalize_lines(rows)), normalize_lines(rows))

    def test_strand_collisions_are_disambiguated(self):
        output = [row for row in normalize_lines([self.row('chr1', '+'), self.row('chr1', '-')]) if '\tCDS\t' in row]
        self.assertIn('g1_locus1.t1', output[0])
        self.assertIn('g1_locus2.t1', output[1])

    def test_augustus_parent_rows_follow_child_identifiers(self):
        rows = [self.row('chr1', '+'), self.row('chr2', '+'),
                'chr2\tAUGUSTUS\tgene\t1\t9\t.\t+\t.\tg1\n',
                'chr2\tAUGUSTUS\ttranscript\t1\t9\t.\t+\t.\tg1.t1\n']
        output = normalize_lines(rows)
        self.assertTrue(any('gene_id "g1_locus2";' in row for row in output))
        self.assertTrue(any('transcript_id "g1_locus2.t1";' in row for row in output))

    def test_nonoverlapping_transcripts_are_separate_loci(self):
        rows = [self.row('chr1', '+'), self.row('chr1', '+', transcript='g1.t2').replace('\t1\t9\t', '\t20\t28\t')]
        output = normalize_lines(rows)
        self.assertEqual(sum('\tgene\t' in row for row in output), 2)

    def test_overlapping_isoforms_remain_in_one_gene(self):
        rows = [self.row('chr1', '+'), self.row('chr1', '+', transcript='g1.t2').replace('\t1\t9\t', '\t5\t13\t')]
        output = normalize_lines(rows)
        self.assertEqual(sum('\tgene\t' in row for row in output), 1)
        self.assertTrue(all('locus' not in row for row in output))

    def test_declared_gene_preserves_nonoverlapping_isoforms(self):
        rows = ['chr1\tAUGUSTUS\tgene\t1\t28\t.\t+\t.\tg1\n',
                self.row('chr1', '+'), self.row('chr1', '+', transcript='g1.t2').replace('\t1\t9\t', '\t20\t28\t')]
        output = normalize_lines(rows)
        self.assertEqual(sum('\tgene\t' in row for row in output), 1)

    def test_declared_transcript_preserves_a_gap_between_exons(self):
        rows = ['chr1\tAUGUSTUS\ttranscript\t1\t28\t.\t+\t.\tg1.t1\n',
                self.row('chr1', '+'), self.row('chr1', '+').replace('\t1\t9\t', '\t20\t28\t')]
        output = normalize_lines(rows)
        self.assertEqual(sum('\ttranscript\t' in row for row in output), 1)

    def test_missing_parents_are_reconstructed_from_child_coordinates(self):
        output = normalize_lines([self.row('chr1', '+')])
        self.assertEqual(len(output), 3)
        self.assertIn('\tgene\t1\t9\t', output[0])
        self.assertIn('transcript_id "g1.t1"; gene_id "g1";', output[1])

    def test_invalid_input_is_rejected(self):
        for rows in ([], ['invalid\n'], [self.row('chr1', '+').replace('gene_id', 'other')]):
            with self.assertRaises(ValueError):
                normalize_lines(rows)


if __name__ == '__main__':
    unittest.main()
