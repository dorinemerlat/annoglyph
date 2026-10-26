import sys
from pathlib import Path
import unittest

import pandas as pd

sys.path.insert(0, str(Path(__file__).resolve().parents[1] / 'bin'))
from plot_annotation_analysis import select_annotation_rows


class AnnotationSelectionTests(unittest.TestCase):
    def setUp(self):
        self.summary = pd.DataFrame([
            ['sample', 'galba', 10], ['sample', 'braker2', 20]],
            columns=['genome_id', 'method_annotation', 'gene_number'])
        self.introns = pd.DataFrame([
            ['sample', 'galba', 30], ['sample', 'braker2', 60]],
            columns=['genome_id', 'method_annotation', 'intron_length'])

    def test_selected_method_keeps_matching_introns(self):
        selection = self.summary.iloc[[1]][['genome_id', 'method_annotation']]
        summary, introns = select_annotation_rows(self.summary, self.introns, selection)
        self.assertEqual(summary.gene_number.tolist(), [20])
        self.assertEqual(introns.intron_length.tolist(), [60])

    def test_duplicate_selection_is_rejected(self):
        with self.assertRaises(ValueError):
            select_annotation_rows(self.summary, self.introns, self.summary)

    def test_missing_selection_is_rejected(self):
        selection = pd.DataFrame([['absent', 'galba']], columns=['genome_id', 'method_annotation'])
        with self.assertRaises(ValueError):
            select_annotation_rows(self.summary, self.introns, selection)


if __name__ == '__main__':
    unittest.main()
