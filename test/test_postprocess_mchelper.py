"""Regression cases for MCHelper outputs from repeat-poor integration fixtures."""
import importlib.util
from pathlib import Path
import tempfile
import unittest

spec = importlib.util.spec_from_file_location('postprocess_mchelper', Path(__file__).resolve().parents[1] / 'bin/postprocess_mchelper.py')
module = importlib.util.module_from_spec(spec)
spec.loader.exec_module(module)


class ClassificationTests(unittest.TestCase):
    def test_no_classified_te_is_a_valid_empty_result(self):
        with tempfile.TemporaryDirectory() as directory:
            path = Path(directory) / 'empty.classif'
            path.touch()
            frame = module.parse_classif_file(path)
            self.assertTrue(frame.empty)
            self.assertEqual(list(frame.columns), ['id', 'comment'])

    def test_nonempty_classification_is_preserved(self):
        with tempfile.TemporaryDirectory() as directory:
            path = Path(directory) / 'te.classif'
            path.write_text('Seq_name\tcoding\nHANUV_TE1_extended\tLTR\n')
            frame = module.parse_classif_file(path)
            self.assertEqual(frame.iloc[0]['id'], 'HANUV_TE1')
            self.assertEqual(frame.iloc[0]['comment'], 'extended')
            self.assertEqual(frame.iloc[0]['MCHelper_coding'], 'LTR')


if __name__ == '__main__':
    unittest.main()
