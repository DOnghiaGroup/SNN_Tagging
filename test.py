import unittest
import ISNN_data_analysis
import os

class BaseCase(unittest.TestCase):
    def setUp(self):
        self.chem_data = 'results-unregularized-matched.fits'
        self.known_data = 'table4.dat'


class FileTest(BaseCase):
    def runTest(self):
        assert os.path.isfile(self.chem_data)
        assert os.path.isfile(self.known_data)


