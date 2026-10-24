import os
import tempfile
import unittest
from unittest.mock import patch, MagicMock, PropertyMock
from gtranslate.config.common import CONFIG
from gtranslate.exceptions import GTranslateExit
from gtranslate.misc import Misc

MODEL_FILES = [CONFIG.ADA_MULTI_CLASS, CONFIG.DT_MULTI_CLASS, CONFIG.KNN_MULTI_CLASS,
               CONFIG.XGB_MULTI_CLASS, CONFIG.MLP_MULTI_CLASS, CONFIG.LABEL_ENCODER]


class TestMisc(unittest.TestCase):
    def setUp(self):
        self.misc = Misc()
        # Mock the logger to avoid polluting the test output
        self.misc.logger = MagicMock()

    @patch('os.path.getsize')
    @patch('os.path.exists')
    def test_checkfile_exists_and_not_empty(self, mock_exists, mock_getsize):
        mock_exists.return_value = True
        mock_getsize.return_value = 1024  # File is > 0 bytes
        result = self.misc.checkfile('/fake/path/file.txt', 'file.txt')
        self.assertTrue(result)
        self.assertTrue(self.misc.logger.info.called)

    @patch('os.path.exists')
    def test_checkfile_missing(self, mock_exists):
        mock_exists.return_value = False
        result = self.misc.checkfile('/fake/path/missing.txt', 'missing.txt')
        self.assertFalse(result)
        self.assertTrue(self.misc.logger.warning.called)

    @patch('os.listdir')
    @patch('os.path.isdir')
    def test_checkfolder_exists_and_not_empty(self, mock_isdir, mock_listdir):
        mock_isdir.return_value = True
        mock_listdir.return_value = ['dummy_file.txt']  # Folder has contents
        result = self.misc.checkfolder('/fake/dir', 'dir')
        self.assertTrue(result)

    @patch('os.path.isdir')
    def test_checkfolder_missing(self, mock_isdir):
        mock_isdir.return_value = False
        result = self.misc.checkfolder('/fake/dir', 'dir')
        self.assertFalse(result)

    def _make_model_dir(self, tmp_dir, missing=None):
        # Create non-empty placeholder model files, optionally leaving one out
        for model_file in MODEL_FILES:
            if model_file != missing:
                with open(os.path.join(tmp_dir, model_file), 'w') as fh:
                    fh.write('model')
        return tmp_dir

    @patch.object(type(CONFIG), 'ENV_MODEL_PATH', new_callable=PropertyMock)
    @patch('gtranslate.misc.check_dependencies')
    def test_check_install_success(self, mock_check_dependencies, mock_model_path):
        # Mock that prodigal is on the PATH and all model files are present
        mock_check_dependencies.return_value = True
        with tempfile.TemporaryDirectory() as tmp_dir:
            mock_model_path.return_value = self._make_model_dir(tmp_dir)
            self.misc.check_install()

        # Verify it logged the presence of the tools
        self.assertTrue(self.misc.logger.info.called)

    @patch.object(type(CONFIG), 'ENV_MODEL_PATH', new_callable=PropertyMock)
    @patch('gtranslate.misc.check_dependencies')
    def test_check_install_missing_model_file(self, mock_check_dependencies, mock_model_path):
        mock_check_dependencies.return_value = True
        with tempfile.TemporaryDirectory() as tmp_dir:
            mock_model_path.return_value = self._make_model_dir(tmp_dir, missing=CONFIG.XGB_MULTI_CLASS)
            with self.assertRaises(GTranslateExit):
                self.misc.check_install()

    @patch.object(type(CONFIG), 'ENV_MODEL_PATH', new_callable=PropertyMock)
    @patch('gtranslate.misc.check_dependencies')
    def test_check_install_model_path_not_set(self, mock_check_dependencies, mock_model_path):
        mock_check_dependencies.return_value = True
        mock_model_path.return_value = None
        with self.assertRaises(GTranslateExit):
            self.misc.check_install()

    @patch.object(type(CONFIG), 'ENV_MODEL_PATH', new_callable=PropertyMock)
    @patch('gtranslate.misc.check_dependencies')
    def test_check_install_prodigal_missing(self, mock_check_dependencies, mock_model_path):
        # Models are present but Prodigal is not on the system path
        mock_check_dependencies.return_value = False
        with tempfile.TemporaryDirectory() as tmp_dir:
            mock_model_path.return_value = self._make_model_dir(tmp_dir)
            with self.assertRaises(GTranslateExit):
                self.misc.check_install()

if __name__ == '__main__':
    unittest.main()