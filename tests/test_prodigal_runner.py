import gzip
import os
import shutil
import tempfile
import unittest
from unittest.mock import patch, MagicMock

from gtranslate.biolib_lite.prodigal_runner import Prodigal

TEST_GENOME_DIR = os.path.join(os.path.dirname(os.path.dirname(os.path.abspath(__file__))),
                               'gtranslate', 'tests', 'data', 'genomes')
TT25_GENOME = os.path.join(TEST_GENOME_DIR, 'genome2_TT25.fna')


@unittest.skipUnless(shutil.which('prodigal'), 'prodigal is not on the system path')
class TestProdigalRunner(unittest.TestCase):

    def setUp(self):
        self.tmp_dir = tempfile.mkdtemp(prefix='gtranslate_test_')

        self.prodigal = Prodigal(cpus=1, verbose=False)
        self.prodigal.logger = MagicMock()
        self.prodigal.warnings = MagicMock()
        self.prodigal.called_genes = False
        self.prodigal.translation_table = None
        self.prodigal.meta = False
        self.prodigal.closed_ends = False
        self.prodigal.output_dir = os.path.join(self.tmp_dir, 'genes')
        os.makedirs(self.prodigal.output_dir)

    def tearDown(self):
        shutil.rmtree(self.tmp_dir)

    @staticmethod
    def _count_genes(faa_file):
        # Read as bytes as genes called on gzipped input contain non-UTF-8 headers
        with open(faa_file, 'rb') as fh:
            return sum(1 for line in fh if line.startswith(b'>'))

    def _call_genes_as_tt25(self, genome_file):
        # The models are not available during testing, so force the ensemble to predict table 25
        with patch('gtranslate.biolib_lite.prodigal_runner.TTPredictor') as mock_predictor:
            mock_predictor.return_value.predict_translation_table.return_value = ('25', 1.0, [], {}, {})
            result = self.prodigal._producer(('genome2_TT25', genome_file, None))

        aa_gene_file, is_empty = result[1], result[9]
        self.assertFalse(is_empty)
        return self._count_genes(aa_gene_file)

    def test_tt25_genes_called_on_gzipped_genome(self):
        # Genes re-called with table 25 must come from the decompressed genome, not the
        # gzip bytes (which Prodigal reads without error, calling only a handful of genes).
        gz_genome = os.path.join(self.tmp_dir, 'genome2_TT25.fna.gz')
        with open(TT25_GENOME, 'rb') as f_in, gzip.open(gz_genome, 'wb') as f_out:
            shutil.copyfileobj(f_in, f_out)

        num_genes_plain = self._call_genes_as_tt25(TT25_GENOME)
        num_genes_gz = self._call_genes_as_tt25(gz_genome)

        self.assertGreater(num_genes_plain, 100)
        self.assertEqual(num_genes_gz, num_genes_plain)


if __name__ == '__main__':
    unittest.main()
