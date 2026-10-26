import importlib.util
from pathlib import Path
import unittest
from unittest.mock import patch
from urllib.error import HTTPError

spec = importlib.util.spec_from_file_location('prepare_evidence', Path(__file__).resolve().parents[1] / 'bin/prepare_protein_evidence.py')
m = importlib.util.module_from_spec(spec)
spec.loader.exec_module(m)


class HttpRetryTests(unittest.TestCase):
    @patch.object(m.time, 'sleep')
    @patch.object(m.urllib.request, 'urlopen')
    def test_rate_limit_respects_retry_after(self, request, sleep):
        response = object()
        request.side_effect = [HTTPError('https://example.org', 429, 'limited', {'Retry-After': '30'}, None), response]
        self.assertIs(m.open_url('https://example.org'), response)
        sleep.assert_called_once_with(30)

    @patch.object(m.time, 'sleep')
    @patch.object(m.urllib.request, 'urlopen')
    def test_permanent_error_is_not_retried(self, request, sleep):
        request.side_effect = HTTPError('https://example.org', 404, 'missing', {}, None)
        with self.assertRaises(HTTPError):
            m.open_url('https://example.org')
        sleep.assert_not_called()
        self.assertEqual(request.call_count, 1)

    @patch.object(m.time, 'sleep')
    @patch.object(m.urllib.request, 'urlopen')
    def test_retries_are_bounded(self, request, sleep):
        request.side_effect = HTTPError('https://example.org', 503, 'unavailable', {}, None)
        with self.assertRaises(HTTPError):
            m.open_url('https://example.org')
        self.assertEqual(request.call_count, 6)
        self.assertEqual([call.args[0] for call in sleep.call_args_list], [10, 20, 40, 80, 160])
