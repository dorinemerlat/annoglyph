import sys
from pathlib import Path
import tempfile
import unittest

sys.path.insert(0, str(Path(__file__).resolve().parents[1] / 'bin'))
from annoglyph_gff import read_gff


class LegacyRepeatGffTests(unittest.TestCase):
    def test_mchelper_semicolons_are_preserved_and_escaped_on_export(self):
        original = ('SEQ_1\tRepeatMasker\trepetitive_elements\t1\t8\t4\t+\t.\t'
                    'ID=repeat1;MCHelper_struct=struct=(TElength: 8bps; TermRepeats: termLTR: 4; GAG );status=NA;\n')
        with tempfile.TemporaryDirectory() as temp:
            source = Path(temp) / 'legacy.gff'
            source.write_text(original)
            features = read_gff(source)
            self.assertEqual(features[0].attrs['MCHelper_struct'], ['struct=(TElength: 8bps; TermRepeats: termLTR: 4; GAG )'])
            self.assertEqual(features[0].attrs['status'], ['NA'])
            exported = features[0].line()
            self.assertIn('%3B', exported)
            source.write_text(exported + '\n')
            self.assertEqual(read_gff(source), features)

    def test_malformed_structural_annotations_still_fail(self):
        with tempfile.TemporaryDirectory() as temp:
            source = Path(temp) / 'bad.gff'
            source.write_text('SEQ_1\tBRAKER\tgene\t1\t8\t.\t+\t.\tID=g1;broken\n')
            with self.assertRaisesRegex(ValueError, 'malformed GFF3 attribute'):
                read_gff(source)


if __name__ == '__main__':
    unittest.main()
