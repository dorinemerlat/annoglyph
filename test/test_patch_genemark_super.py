import importlib.util
import subprocess
import tempfile
import unittest
from pathlib import Path


spec = importlib.util.spec_from_file_location('patch_genemark_super', Path(__file__).resolve().parents[1] / 'bin/patch_genemark_super.py')
patch_module = importlib.util.module_from_spec(spec)
spec.loader.exec_module(patch_module)


class GeneMarkSummaryTests(unittest.TestCase):
    def test_reports_exons_for_single_and_multi_exon_genes(self):
        with tempfile.TemporaryDirectory() as temporary:
            script = Path(temporary) / 'super.pl'
            script.write_text('''use strict;
use warnings;
open(my $OUT, '>&', \\*STDOUT) or die $!;
for my $gcount (1..3) {
    my @cds = (1) x $gcount;
    my @introns = (1) x ($gcount - 1);
    my $cds_length = 600;
    print $OUT ("## ". $gcount ."\\t". $cds_length ."\\t". (scalar @introns) ."\\n");
}
''')
            patch_module.patch_super(script)
            patched = script.read_text()
            patch_module.patch_super(script)
            self.assertEqual(script.read_text(), patched)
            result = subprocess.run(['perl', str(script)], check=True, capture_output=True, text=True)
            exon_counts = [int(line.split()[-1]) for line in result.stdout.splitlines()]
            self.assertEqual(exon_counts, [1, 2, 3])

    def test_rejects_unknown_summary_without_modifying_it(self):
        with tempfile.TemporaryDirectory() as temporary:
            script = Path(temporary) / 'super.pl'
            script.write_text('unknown upstream format\n')
            with self.assertRaisesRegex(ValueError, 'Unrecognized'):
                patch_module.patch_super(script)
            self.assertEqual(script.read_text(), 'unknown upstream format\n')


if __name__ == '__main__':
    unittest.main()
