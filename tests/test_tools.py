import os
import tempfile
import unittest
from unittest.mock import patch, mock_open
from gtranslate.tools import (
    splitchunks, splitchunks_list, generateTempTableName,
    merge_two_dicts, confirm, sha256, file_has_checksum,
    symlink_f, remove_intermediate_files, get_genomes_size
)


class TestTools(unittest.TestCase):

    @patch('gtranslate.tools.read_fasta')
    def test_get_genomes_size(self, mock_read_fasta):
        # Mocking read_fasta to return a dictionary of sequences
        mock_read_fasta.return_value = {
            'seq1': 'ATGC',  # length 4
            'seq2': 'ATGCATGC'  # length 8
        }
        size = get_genomes_size('/fake/genome.fna')
        self.assertEqual(size, 12)

    def test_splitchunks_list(self):
        data = [1, 2, 3, 4, 5]
        # Split a list of 5 into 2 chunks -> sizes should be 3 and 2
        chunks = list(splitchunks_list(data, 2))
        self.assertEqual(len(chunks), 2)
        self.assertEqual(chunks[0], [1, 2, 3])
        self.assertEqual(chunks[1], [4, 5])

    def test_splitchunks_dict(self):
        data = {'a': 1, 'b': 2, 'c': 3, 'd': 4, 'e': 5}
        # Split a dict of 5 into 2 chunks -> sizes should be 3 and 2
        chunks = list(splitchunks(data, 2))
        self.assertEqual(len(chunks), 2)
        self.assertIn('a', chunks[0])
        self.assertIn('e', chunks[1])

    def test_generateTempTableName(self):
        name1 = generateTempTableName()
        name2 = generateTempTableName()
        self.assertTrue(name1.startswith('TEMP'))
        self.assertNotEqual(name1, name2)  # Ensure randomness
        self.assertGreater(len(name1), 10)

    def test_merge_two_dicts(self):
        dict1 = {'a': 1, 'b': 2}
        dict2 = {'b': 3, 'c': 4}
        merged = merge_two_dicts(dict1, dict2)

        # Check dict2 overwrote dict1 for key 'b'
        self.assertEqual(merged, {'a': 1, 'b': 3, 'c': 4})
        # Check original dict wasn't mutated
        self.assertEqual(dict1, {'a': 1, 'b': 2})

    @patch('builtins.input', return_value='y')
    def test_confirm_yes(self, mock_input):
        self.assertTrue(confirm("Do you agree?"))

    @patch('builtins.input', return_value='n')
    def test_confirm_no(self, mock_input):
        self.assertFalse(confirm("Do you agree?"))

    def test_sha256(self):
        # 1. Create a real temporary file
        with tempfile.NamedTemporaryFile(delete=False) as tmp:
            tmp.write(b"dummy content")
            tmp_path = tmp.name

        try:
            result_hash = sha256(tmp_path)

            expected_hash_256 = "bf0ecbdb9b814248d086c9b69cf26182d9d4138f2ad3d0637c4555fc8cbf68e5"
            self.assertEqual(result_hash, expected_hash_256)

        finally:
            if os.path.exists(tmp_path):
                os.remove(tmp_path)

    @patch('os.path.isfile')
    @patch('gtranslate.tools.sha256')
    @patch('builtins.open', new_callable=mock_open, read_data="fake_hash")
    def test_file_has_checksum(self, mock_file, mock_sha256, mock_isfile):
        mock_isfile.return_value = True
        mock_sha256.return_value = "fake_hash"

        result = file_has_checksum("fake_file.txt")
        self.assertTrue(result)

    @patch('os.symlink')
    @patch('os.remove')
    @patch('os.path.isfile')
    def test_symlink_f(self, mock_isfile, mock_remove, mock_symlink):
        mock_isfile.return_value = True  # Pretend the destination file already exists
        symlink_f('src.txt', 'dst.txt', force=True)

        mock_remove.assert_called_once_with('dst.txt')
        mock_symlink.assert_called_once_with('src.txt', 'dst.txt')

    @patch('shutil.rmtree')
    @patch('os.path.isdir')
    @patch('os.path.exists')
    def test_remove_intermediate_files(self, mock_exists, mock_isdir, mock_rmtree):
        mock_exists.return_value = True
        mock_isdir.return_value = True

        remove_intermediate_files('/output/dir', 'workflow_name')
        # Check if rmtree was triggered to clean the directory
        self.assertTrue(mock_rmtree.called)


if __name__ == '__main__':
    unittest.main()