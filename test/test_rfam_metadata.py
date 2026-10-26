import gzip
import json
import sys
import tempfile
import unittest
from pathlib import Path
from unittest.mock import MagicMock, patch
from urllib.error import URLError

sys.path.insert(0, str(Path(__file__).resolve().parents[1] / 'bin'))
from rfam_metadata import family_types, classify_family
import resolve_rfam_metadata as resolver
import merge_ncrna as merger
from annoglyph_gff import read_gff
from ncrna_stats import rna_kind


def family_row(accession, kind):
    return '\t'.join([accession] + [''] * 17 + [kind]) + '\n'


def infernal_row(accession):
    return ' '.join(['1', 'example_model', accession, 'contig1', '-', 'CL001', 'cm',
                     '1', '72', '100', '171', '+', 'no', '1', '0.5', '0', '50',
                     '1e-10', '!', '*', '-', '-', '-', '-', '-', '-', '-', 'description']) + '\n'


class RfamMetadataTests(unittest.TestCase):
    def test_official_types_and_future_subtype(self):
        for source, expected in [
            ('Gene; tRNA;', ('tRNA', None, 'tRNA')),
            ('Gene; rRNA;', ('rRNA', None, 'rRNA')),
            ('Cis-reg;', ('ncRNA', 'other', 'Regulatory RNA')),
            ('Gene; snoRNA; CD-box;', ('ncRNA', 'snoRNA', 'snoRNA')),
            ('Gene; miRNA;', ('ncRNA', 'miRNA', 'miRNA')),
            ('Gene; futureRNA;', ('ncRNA', 'other', 'futureRNA')),
            ('Intron;', ('ncRNA', 'autocatalytically_spliced_intron', 'Intron')),
            ('Gene;', ('ncRNA', 'other', 'Unclassified ncRNA')),
            ('', ('ncRNA', 'other', 'Unclassified ncRNA')),
        ]:
            with self.subTest(source=source):
                self.assertEqual(classify_family(source), expected)

    def test_sql_prose_encoding_does_not_break_type_column(self):
        with tempfile.TemporaryDirectory() as temp:
            path = Path(temp) / 'family.txt'
            row = family_row('RF04340', 'Cis-reg;').encode().split(b'\t')
            row[3] = b'legacy\x92description'
            path.write_bytes(b'\t'.join(row))
            self.assertEqual(family_types(path), {'RF04340': 'Cis-reg;'})

    def test_missing_family_resolves_from_official_metadata_and_caches(self):
        with tempfile.TemporaryDirectory() as temp:
            root = Path(temp)
            local = root / 'local.txt'
            local.write_text(family_row('RF00005', 'Gene; tRNA;'))
            table = root / 'hits.tblout'
            table.write_text(infernal_row('RF04340'))
            response = MagicMock()
            response.__enter__.return_value.read.return_value = gzip.compress(family_row('RF04340', 'Cis-reg;').encode())
            with patch.object(resolver, 'urlopen', return_value=response) as request:
                result = resolver.resolve(local, [table], root / 'cache', root / 'out')
            self.assertEqual(request.call_count, 1)
            self.assertEqual(result['automatically_resolved'], ['RF04340'])
            self.assertEqual(family_types(root / 'out/family.txt')['RF04340'], 'Cis-reg;')
            with patch.object(resolver, 'urlopen', side_effect=AssertionError('Must use cache')):
                second = resolver.resolve(local, [table], root / 'cache', root / 'out2')
            self.assertEqual(second['automatically_resolved'], ['RF04340'])

    def test_offline_unknown_is_retained_and_audited_in_english(self):
        with tempfile.TemporaryDirectory() as temp:
            root = Path(temp)
            families = root / 'family.txt'
            families.write_text(family_row('RF00005', 'Gene; tRNA;'))
            table = root / 'hits.tblout'
            table.write_text(infernal_row('RF99999'))
            with patch.object(resolver, 'urlopen', side_effect=URLError('offline')):
                report = resolver.resolve(families, [table], root / 'cache', root / 'resolved')
            self.assertEqual(report['unclassified'], ['RF99999'])
            empty = root / 'empty.gff3'
            empty.write_text('##gff-version 3\n')
            output, audit = root / 'output.gff3', root / 'audit.tsv'
            args = ['merge_ncrna', '--trnascan', str(empty), '--barrnap', str(empty),
                    '--infernal', str(table), '--families', str(root / 'resolved/family.txt'),
                    '--prefix', 'sample', '--output', str(output), '--audit', str(audit)]
            with patch('sys.argv', args):
                merger.main()
            features = read_gff(output)
            self.assertEqual(len(features), 1)
            self.assertEqual(features[0].kind, 'ncRNA')
            self.assertEqual(features[0].attrs['rna_type'], ['Unclassified ncRNA'])
            self.assertEqual(features[0].attrs['Dbxref'], ['Rfam:RF99999'])
            self.assertEqual(rna_kind(features[0]), 'Unclassified ncRNA')
            self.assertIn('Unclassified ncRNA', audit.read_text())

    def test_known_metadata_needs_no_network_and_regulatory_type_reaches_figures(self):
        with tempfile.TemporaryDirectory() as temp:
            root = Path(temp)
            local = root / 'local.txt'
            local.write_text(family_row('RF04340', 'Cis-reg;'))
            table = root / 'hits.tblout'
            table.write_text(infernal_row('RF04340'))
            with patch.object(resolver, 'urlopen', side_effect=AssertionError('No network needed')):
                report = resolver.resolve(local, [table], root / 'cache', root / 'out')
            self.assertFalse(report['unclassified'])
            feature = list(merger.infernal_features(table, family_types(root / 'out/family.txt')))[0]
            self.assertEqual(rna_kind(feature), 'Regulatory RNA')
            self.assertEqual(feature.attrs['ncRNA_class'], ['other'])


if __name__ == '__main__':
    unittest.main()
