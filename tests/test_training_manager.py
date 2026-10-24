import os
import tempfile
import unittest
from unittest.mock import MagicMock

from gtranslate.training_manager import TrainingManager


class TestGroundTruth(unittest.TestCase):

    def setUp(self):
        self.tmp_dir = tempfile.TemporaryDirectory()
        self.training_manager = TrainingManager()
        self.training_manager.logger = MagicMock()

    def tearDown(self):
        self.tmp_dir.cleanup()

    def _ground_truth(self, header, rows):
        taxonomy_file = os.path.join(self.tmp_dir.name, 'taxonomy.tsv')
        out_file = os.path.join(self.tmp_dir.name, 'ground_truth.tsv')
        with open(taxonomy_file, 'w') as fh:
            fh.write('\t'.join(header) + '\n')
            for row in rows:
                fh.write('\t'.join(row) + '\n')

        self.training_manager.select_ground_truth(taxonomy_file, out_file)

        with open(out_file) as fh:
            fh.readline()
            rows = [line.rstrip('\n').split('\t') for line in fh]
        return {row[0]: row[1] for row in rows}

    def test_taxonomy_column_first(self):
        # A taxonomy column at index 0 must not be ignored
        mycoplasma = 'd__Bacteria;p__Bacillota;c__Bacilli;o__Mycoplasmatales;f__Mycoplasmataceae'
        gracilibacteria = 'd__Bacteria;p__Patescibacteriota;c__JAEDAM01;o__Absconditabacterales'
        ecoli = 'd__Bacteria;p__Pseudomonadota;c__Gammaproteobacteria;o__Enterobacterales'
        ground_truth = self._ground_truth(['Taxonomy', 'Genome ID'],
                                          [[mycoplasma, 'G1'], [gracilibacteria, 'G2'], [ecoli, 'G3']])
        self.assertEqual(ground_truth, {'G1': '4', 'G2': '25', 'G3': '11'})

    def test_taxonomy_column_not_first(self):
        mycoplasma = 'd__Bacteria;p__Bacillota;c__Bacilli;o__Mycoplasmatales;f__Mycoplasmataceae'
        ground_truth = self._ground_truth(['Genome ID', 'Taxonomy'], [['G1', mycoplasma]])
        self.assertEqual(ground_truth, {'G1': '4'})


if __name__ == '__main__':
    unittest.main()
