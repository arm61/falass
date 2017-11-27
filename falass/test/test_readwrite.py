import numpy as np
from numpy.testing import assert_equal, assert_almost_equal
from falass import readwrite

class TestFiles(object):
    def test_readPDB(self):
        pdb = readwrite.readRDB('test.pdb')
        assert_equal(pdb.cell, [67.823, 67.823, 72.389])
        assert_equal(pdb.times, [0, 10000, 20000, 30000, 40000, 50000])
        assert_equal(pdb.number_of_timesteps, 6)
        assert_equal(pdb.atoms, [['C1', 'C2', 'C3'], ['C1', 'C2', 'C3'], ['C1', 'C2', 'C3'], ['C1', 'C2', 'C3'],
                                 ['C1', 'C2', 'C3'], ['C1', 'C2', 'C3']])
        assert_equal(pdb.zpos, [[45.960, 43.960, 45.630], [49.735, 50.838, 48.750], [45.721, 45.026, 46.026],
                                [45.541, 44.266, 45.383], [46.017, 43.761, 44.251], [43.957, 43.925, 43.946]])