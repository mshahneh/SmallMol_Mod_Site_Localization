import sys
sys.path.append('..')
import unittest
from Compound_n import Compound
import json
import os
import pickle

class TestUtils(unittest.TestCase):
    def __init__(self, *args, **kwargs):
        super(TestUtils, self).__init__(*args, **kwargs)

    # TODO: test calculate_scores_n

if __name__ == '__main__':
    unittest.main()