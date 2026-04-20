import unittest
from unittest.mock import patch, MagicMock
from gtranslate.misc import Misc


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

    @patch('gtranslate.misc.check_dependencies')
    def test_check_install_success(self, mock_check_dependencies):
        # Mock that prodigal are on the PATH
        mock_check_dependencies.return_value = True

        # Test the function (note: your current check_install hardcodes `ok = True`)
        self.misc.check_install()

        # Verify it logged the presence of the tools
        self.assertTrue(self.misc.logger.info.called)


if __name__ == '__main__':
    unittest.main()