from numpy.testing import assert_equal
from falass import dataformat
import unittest

class TestQData(unittest.TestCase):
    def test_qdata(self):
        a = dataformat.QData(1., 2., 3., 4.)
        assert_equal(a.q, 1.)
        assert_equal(a.i, 2.)
        assert_equal(a.di, 3.)
        assert_equal(a.dq, 4.)

class TestScatLens(unittest.TestCase):
    def test_scatlens(self):
        a = dataformat.ScatLens('C1', 1., 0.)
        assert_equal(a.atom, 'C1')
        assert_equal(a.real, 1.0e-5)
        assert_equal(a.imag, 0.0)

class TestAtomPositions(unittest.TestCase):
    def test_atompositions(self):
        a = dataformat.AtomPositions('C1', 1.000)
        assert_equal(a.atom, 'C1')
        assert_equal(a.zpos, 1.000)

class TestSLDPro(unittest.TestCase):
    def test_sldpro(self):
        a = dataformat.SLDPro(1., 2., 3.)
        assert_equal(a.thick, 1.)
        assert_equal(a.real, 2.)
        assert_equal(a.imag, 3.)