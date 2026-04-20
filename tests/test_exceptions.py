import unittest
from gtranslate.exceptions import (
    GTranslateException,
    GTranslateExit,
    GenomeNameInvalid,
    ProdigalException
)


class TestExceptions(unittest.TestCase):
    def test_base_exception(self):
        exc = GTranslateException("Base error message")
        self.assertEqual(str(exc), "Base error message")
        self.assertTrue(isinstance(exc, Exception))

    def test_gtranslate_exit(self):
        exc = GTranslateExit("Quiet exit")
        self.assertEqual(str(exc), "Quiet exit")
        self.assertTrue(isinstance(exc, Exception))

    def test_specific_exceptions(self):
        exc1 = GenomeNameInvalid("Invalid name")
        self.assertTrue(isinstance(exc1, GTranslateException))

        exc2 = ProdigalException("Prodigal failed")
        self.assertTrue(isinstance(exc2, GTranslateException))


if __name__ == '__main__':
    unittest.main()