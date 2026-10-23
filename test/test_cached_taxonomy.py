import importlib.util
from pathlib import Path
import tempfile
import unittest
import io
import json
from urllib.error import HTTPError, URLError
from unittest.mock import patch, MagicMock

spec = importlib.util.spec_from_file_location('taxonomy', Path(__file__).resolve().parents[1] / 'bin/generate_taxonomy_metadata.py')
taxonomy = importlib.util.module_from_spec(spec)
spec.loader.exec_module(taxonomy)


class CachedTaxonomyTests(unittest.TestCase):
    def setUp(self):
        taxonomy._last_ncbi_request = None

    def test_reads_exact_terminal_taxon_without_network(self):
        with tempfile.TemporaryDirectory() as temp:
            root = Path(temp)
            sample = root / 'sample-1'
            sample.mkdir()
            (sample / '123_lineage.csv').write_text('rank\tid\tname\ndomain\t2759\tEukaryota\ngenus\t12\tExample\nspecies\t123\tExample species\n')
            (sample / '123.mnemonic').write_text('EXASP\n')
            with patch.object(taxonomy, 'fetch', side_effect=AssertionError('Network must not be used')):
                result, mnemonic = taxonomy.cached_taxonomy(root, 'sample-1', '123')
            self.assertEqual(result['scientific_name'], 'Example species')
            self.assertEqual(result['lineage'][-1]['taxid'], '123')
            self.assertEqual(mnemonic, 'EXASP')

    def test_rejects_wrong_terminal_taxon(self):
        with tempfile.TemporaryDirectory() as temp:
            root = Path(temp)
            (root / 'sample').mkdir()
            (root / 'sample/123_lineage.csv').write_text('rank\tid\tname\nspecies\t456\tWrong species\n')
            with self.assertRaisesRegex(ValueError, 'does not end at requested taxon'):
                taxonomy.cached_taxonomy(root, 'sample', '123')

    def test_429_waits_retry_after_then_succeeds(self):
        error = HTTPError('https://example.test', 429, 'Too many requests', {'Retry-After': '12'}, io.BytesIO())
        response = MagicMock()
        response.__enter__.return_value.read.return_value = b'ok'
        with patch.object(taxonomy.urllib.request, 'urlopen', side_effect=[error, response]) as request, \
                patch.object(taxonomy.time, 'sleep') as sleep:
            self.assertEqual(taxonomy.fetch('https://example.test', 'application/xml'), b'ok')
        sleep.assert_called_once_with(12.0)
        self.assertEqual(request.call_count, 2)
        self.assertEqual(request.call_args.kwargs['timeout'], 60)

    def test_permanent_error_is_not_retried(self):
        error = HTTPError('https://example.test', 400, 'Bad taxid', {}, io.BytesIO())
        with patch.object(taxonomy.urllib.request, 'urlopen', side_effect=error) as request, \
                patch.object(taxonomy.time, 'sleep') as sleep:
            with self.assertRaises(HTTPError):
                taxonomy.fetch('https://example.test', 'application/xml')
        self.assertEqual(request.call_count, 1)
        sleep.assert_not_called()

    def test_transient_failure_has_bounded_attempts(self):
        with patch.object(taxonomy.urllib.request, 'urlopen', side_effect=URLError('temporary')) as request, \
                patch.object(taxonomy.time, 'sleep') as sleep:
            with self.assertRaises(URLError):
                taxonomy.fetch('https://example.test', 'application/xml', attempts=3)
        self.assertEqual(request.call_count, 3)
        self.assertEqual([call.args[0] for call in sleep.call_args_list], [2, 4])

    def test_persistent_cache_survives_a_later_failure(self):
        value = {'taxid': '123', 'scientific_name': 'Example species',
                 'lineage': [{'rank': 'species', 'taxid': '123', 'name': 'Example species'}]}
        with tempfile.TemporaryDirectory() as temp:
            root = Path(temp)
            self.assertEqual(taxonomy.cached_remote(root, 'taxonomy', '123', lambda: value), value)
            with self.assertRaises(RuntimeError):
                taxonomy.cached_remote(root, 'taxonomy', '456', lambda: (_ for _ in ()).throw(RuntimeError('offline')))
            with patch.object(taxonomy, 'fetch', side_effect=AssertionError('Network must not be used')):
                self.assertEqual(taxonomy.cached_remote(root, 'taxonomy', '123', lambda: taxonomy.fetch_taxonomy('123')), value)
            self.assertFalse((root / 'taxonomy_456.json').exists())

    def test_retry_after_http_date(self):
        self.assertGreater(taxonomy.retry_after_seconds('Thu, 01 Jan 2099 00:00:00 GMT'), 0)
        self.assertEqual(taxonomy.retry_after_seconds('bad value'), 0)

    def test_ncbi_requests_are_spaced_across_instances(self):
        with tempfile.TemporaryDirectory() as temp:
            root = Path(temp)
            with patch.object(taxonomy.time, 'time', return_value=100), patch.object(taxonomy.time, 'sleep'):
                taxonomy.throttle_ncbi(root)
            with patch.object(taxonomy.time, 'time', return_value=100.25), patch.object(taxonomy.time, 'sleep') as sleep:
                taxonomy.throttle_ncbi(root)
            sleep.assert_called_once_with(0.75)

    def test_full_command_uses_historical_cache_offline(self):
        with tempfile.TemporaryDirectory() as temp:
            root = Path(temp)
            sample = root / 'lineages/sample-1'
            sample.mkdir(parents=True)
            (sample / '123_lineage.csv').write_text('rank\tid\tname\nspecies\t123\tExample species\n')
            (sample / '123.mnemonic').write_text('EXASP\n')
            (root / 'input.tsv').write_text('id\tname\ttaxid\tmnemonic\nsample-1\tExample species\t123\tEXASP\n')
            argv = ['taxonomy', '--input', str(root / 'input.tsv'), '--output', str(root / 'out'),
                    '--taxonomy-cache', str(root / 'lineages'), '--response-cache', str(root / 'remote')]
            with patch('sys.argv', argv), patch.object(taxonomy, 'fetch', side_effect=AssertionError('Network must not be used')):
                taxonomy.main()
            self.assertEqual(json.loads((root / 'out/taxonomy.json').read_text())['sample-1']['taxid'], '123')


if __name__ == '__main__':
    unittest.main()
