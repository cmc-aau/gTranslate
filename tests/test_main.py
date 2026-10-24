import os
import tempfile
import unittest
from unittest.mock import patch, MagicMock

from gtranslate.config.common import CONFIG
from gtranslate.exceptions import GTranslateExit
from gtranslate.main import OptionsParser


class TestVerifyModelFiles(unittest.TestCase):

    def _make_model_dir(self, tmp_dir, missing=None):
        # Create placeholder model files, optionally leaving one out
        for model_file in CONFIG.MODEL_FILES:
            if model_file != missing:
                with open(os.path.join(tmp_dir, model_file), 'w') as fh:
                    fh.write('model')
        return tmp_dir

    def test_all_model_files_present(self):
        with tempfile.TemporaryDirectory() as tmp_dir:
            self.assertTrue(OptionsParser._verify_model_files(self._make_model_dir(tmp_dir)))

    def test_missing_model_file(self):
        with tempfile.TemporaryDirectory() as tmp_dir:
            model_dir = self._make_model_dir(tmp_dir, missing=CONFIG.MLP_MULTI_CLASS)
            with self.assertRaises(GTranslateExit) as ctx:
                OptionsParser._verify_model_files(model_dir)
            self.assertIn(CONFIG.MLP_MULTI_CLASS, str(ctx.exception))

    def test_missing_model_directory(self):
        with self.assertRaises(GTranslateExit):
            OptionsParser._verify_model_files('/nonexistent/gtranslate/models')

    @patch('gtranslate.main.TablePredictor')
    @patch('gtranslate.main.check_dependencies')
    def test_detect_table_stops_before_prodigal(self, mock_check_dependencies, mock_table_predictor):
        # A missing model must be reported before any genes are called
        with tempfile.TemporaryDirectory() as tmp_dir:
            genome_dir = os.path.join(tmp_dir, 'genomes')
            os.makedirs(genome_dir)
            with open(os.path.join(genome_dir, 'genome1.fna'), 'w') as fh:
                fh.write('>contig1\nATGC\n')
            model_dir = self._make_model_dir(os.path.join(tmp_dir), missing=CONFIG.LABEL_ENCODER)

            options = MagicMock(genome_dir=genome_dir, batchfile=None, extension='fna',
                                out_dir=os.path.join(tmp_dir, 'out'), cpus=1,
                                custom_model_path=model_dir)
            parser = OptionsParser('test')
            parser.logger = MagicMock()
            with self.assertRaises(GTranslateExit):
                parser.detect_table(options)
            mock_table_predictor.assert_not_called()


if __name__ == '__main__':
    unittest.main()
