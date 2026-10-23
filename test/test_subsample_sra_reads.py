import gzip
import importlib.util
from pathlib import Path
import tempfile
import unittest

spec = importlib.util.spec_from_file_location('sampler', Path(__file__).resolve().parents[1] / 'bin/subsample_sra_reads.py')
sampler = importlib.util.module_from_spec(spec)
spec.loader.exec_module(sampler)


class SubsamplingTests(unittest.TestCase):
    def run_sample(self, ids, limit):
        with tempfile.TemporaryDirectory() as tmp:
            root = Path(tmp)
            inputs = []
            for mate, names in enumerate(ids, 1):
                path = root / f'reads_{mate}.fastq'
                path.write_text(''.join(f'@{name}/{mate}\nACGT\n+\nIIII\n' for name in names))
                inputs.append(path)
            count = sampler.subsample(inputs, root / 'out', limit)
            contents = [gzip.open(path, 'rt').read() for path in sorted((root / 'out').glob('*.gz'))]
            return count, contents

    def test_single_limit(self):
        count, reads = self.run_sample([['a', 'b', 'c']], 2)
        self.assertEqual(count, 2)
        self.assertNotIn('@c/', reads[0])

    def test_pairs_count_as_spots(self):
        count, reads = self.run_sample([['a', 'b', 'c']] * 2, 2)
        self.assertEqual(count, 2)
        self.assertTrue(all(text.count('\n') == 8 for text in reads))
        self.assertIn('@b/1', reads[0])
        self.assertIn('@b/2', reads[1])

    def test_short_run(self):
        self.assertEqual(self.run_sample([['a']], 5)[0], 1)

    def test_mismatched_pairs(self):
        with self.assertRaisesRegex(ValueError, 'same spot'):
            self.run_sample([['a'], ['b']], 2)

    def test_missing_mate(self):
        with self.assertRaisesRegex(ValueError, 'different numbers'):
            self.run_sample([['a', 'b'], ['a']], 5)


if __name__ == '__main__':
    unittest.main()
