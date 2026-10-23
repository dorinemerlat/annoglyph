import importlib.util
from pathlib import Path
import unittest

spec = importlib.util.spec_from_file_location('prefixer', Path(__file__).resolve().parents[1] / 'bin/prefix_gff_ids.py')
prefixer = importlib.util.module_from_spec(spec)
spec.loader.exec_module(prefixer)


class PrefixTests(unittest.TestCase):
    def test_first_attribute_and_parent_remain_linked(self):
        rows = ['chr1\ttool\tgene\t1\t9\t.\t+\t.\tID=g1\n',
                'chr1\ttool\tmRNA\t1\t9\t.\t+\t.\tID=g1.t1;Parent=g1\n',
                'chr1\ttool\tCDS\t1\t9\t.\t+\t0\tID=g1.t1.cds;Parent=g1.t1\n']
        converted = [prefixer.prefix_line(row, 'TEST') for row in rows]
        self.assertIn('ID=TEST.g1\n', converted[0])
        self.assertIn('ID=TEST.g1.t1;Parent=TEST.g1\n', converted[1])
        self.assertIn('Parent=TEST.g1.t1\n', converted[2])
        self.assertEqual([row.split('\t')[:8] for row in converted], [row.split('\t')[:8] for row in rows])

    def test_multiple_parents_and_unrelated_attributes(self):
        row = 'chr1\ttool\texon\t1\t9\t.\t+\t.\tID=g1.exon;Parent=g1.t1,g1.t2;Note=Parent=g1\n'
        result = prefixer.prefix_line(row, 'TEST')
        self.assertIn('Parent=TEST.g1.t1,TEST.g1.t2;', result)
        self.assertTrue(result.endswith('Note=Parent=g1\n'))
        self.assertEqual(prefixer.prefix_line(result, 'TEST'), result)

    def test_comment(self):
        self.assertEqual(prefixer.prefix_line('##gff-version 3\n', 'TEST'), '##gff-version 3\n')


if __name__ == '__main__':
    unittest.main()
