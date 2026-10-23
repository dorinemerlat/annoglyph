import importlib.util
from pathlib import Path
import tempfile
import unittest

spec = importlib.util.spec_from_file_location('reuse', Path(__file__).resolve().parents[1] / 'bin/prepare_cache_reuse.py')
reuse = importlib.util.module_from_spec(spec)
spec.loader.exec_module(reuse)


class CacheReuseTests(unittest.TestCase):
    def setUp(self):
        self.temp = tempfile.TemporaryDirectory()
        self.root = Path(self.temp.name)
        self.reference = self.root / 'reference.fa'
        self.reference.write_bytes(b'>JULSC2_1\nACGT\nAC\n>JULSC2_2\nNNAT\n')
        self.index = reuse.reference_index(self.reference)
        self.mapping = {b'JULSC1_1': b'JULSC2_1', b'JULSC1_2': b'JULSC2_2'}

    def tearDown(self):
        self.temp.cleanup()

    def test_masking_and_line_wrapping_do_not_change_sequence_identity(self):
        source = self.root / 'masked.fa'
        source.write_bytes(b'>JULSC1_1\nacg\ntAC\n>JULSC1_2\nNnaT\n')
        target = self.root / 'normalized.fa'
        result = reuse.normalize_masked(source, target, self.mapping, self.index)
        self.assertEqual(result['bases_verified'], 10)
        self.assertEqual(reuse.reference_index(target), self.index)
        self.assertIn(b'acg\ntAC', target.read_bytes())
        self.assertIn(b'>JULSC1_1', source.read_bytes())

    def test_same_sequence_name_does_not_hide_different_assembly(self):
        source = self.root / 'wrong.fa'
        source.write_bytes(b'>JULSC1_1\nACGTAT\n>JULSC1_2\nNNAT\n')
        target = self.root / 'normalized.fa'
        with self.assertRaisesRegex(ValueError, 'differs from assembly'):
            reuse.normalize_masked(source, target, self.mapping, self.index)
        self.assertFalse(target.exists())
        self.assertFalse(target.with_name('normalized.fa.tmp').exists())

    def test_gff_preserves_coordinates_and_gene_relationships(self):
        source = self.root / 'original.gff3'
        source.write_text('##gff-version 3\n##sequence-region JULSC1_1 1 6\n'
                          'JULSC1_1\ttool\tgene\t1\t6\t.\t+\t.\tID=g1\n'
                          'JULSC1_1\ttool\tmRNA\t1\t6\t.\t+\t.\tID=g1.t1;Parent=g1\n')
        target = self.root / 'normalized.gff3'
        result = reuse.normalize_gff(source, target, self.mapping, self.index)
        self.assertEqual(result, {'rows': 2, 'renamed_rows': 2})
        self.assertEqual(target.read_text(), source.read_text().replace('JULSC1_1', 'JULSC2_1'))
        self.assertIn('Parent=g1', target.read_text())

    def test_invalid_gff_coordinates_are_rejected(self):
        source = self.root / 'outside.gff'
        source.write_text('JULSC1_1\ttool\tgene\t1\t7\t.\t+\t.\tID=g1\n')
        with self.assertRaisesRegex(ValueError, 'outside cached assembly'):
            reuse.normalize_gff(source, self.root / 'new.gff', self.mapping, self.index)

    def test_zero_length_repeat_entries_are_audited_not_reinterpreted(self):
        source = self.root / 'repeats.gff'
        zero = 'JULSC1_1\tRepeatMasker\trepetitive_elements\t3\t2\t.\t+\t.\tID=zero\n'
        valid = 'JULSC1_1\tRepeatMasker\trepetitive_elements\t1\t2\t.\t+\t.\tID=valid\n'
        source.write_text(zero + valid)
        audit = self.root / 'excluded.gff'
        target = self.root / 'valid.gff'
        result = reuse.normalize_gff(source, target, self.mapping, self.index, audit)
        self.assertEqual(result['excluded_zero_length_rows'], 1)
        self.assertEqual(result['rows'], 1)
        self.assertEqual(audit.read_text(), zero)
        self.assertEqual(source.read_text(), zero + valid)
        self.assertIn('ID=valid', target.read_text())
        self.assertNotIn('ID=zero', target.read_text())

    def test_sidecar_translation_is_simultaneous(self):
        source = self.root / 'trna.struct'
        source.write_bytes(b'JULSC1_1.trna1 JULSC2_1.trna2 unrelated_JULSC1_1\n')
        target = self.root / 'new.struct'
        mapping = {b'JULSC1_1': b'JULSC2_1', b'JULSC2_1': b'JULSC1_1'}
        reuse.normalize_text(source, target, mapping)
        self.assertEqual(target.read_bytes(), b'JULSC2_1.trna1 JULSC1_1.trna2 unrelated_JULSC1_1\n')

    def test_unchanged_gff_uses_a_link(self):
        source = self.root / 'unchanged.gff'
        source.write_text('JULSC2_1\ttool\tgene\t1\t6\t.\t+\t.\tID=g1\n')
        target = self.root / 'linked.gff'
        reuse.normalize_gff(source, target, {}, self.index)
        self.assertTrue(target.is_symlink())
        self.assertTrue(target.samefile(source))


if __name__ == '__main__':
    unittest.main()
