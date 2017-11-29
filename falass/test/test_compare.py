from numpy.testing import assert_almost_equal
from falass import compare, dataformat
import numpy as np
import unittest

class TestCompare(unittest.TestCase):
    def test_compare(self):
        data1 = dataformat.QData(0.05, 3., 0.3, 0.05 * 0.05)
        data2 = dataformat.QData(0.25, 2., 0.2, 0.05 * 0.25)
        data3 = dataformat.QData(0.50, 1., 0.1, 0.05 * 0.50)
        data = [data1, data2, data3]
        sdata1 = dataformat.QData(0.05, 1., 0.1, 0.05 * 0.05)
        sdata2 = dataformat.QData(0.25, 2., 0.2, 0.05 * 0.25)
        sdata3 = dataformat.QData(0.50, 3., 0.3, 0.05 * 0.50)
        sdata = [sdata1, sdata2, sdata3]
        a = compare.Compare(data, sdata, 1., 0.)
        assert_almost_equal(a.exp_data[0].q, 0.05)
        assert_almost_equal(a.exp_data[0].i, 3.)
        assert_almost_equal(a.exp_data[0].di, 0.3)
        assert_almost_equal(a.exp_data[0].dq, 0.05 * 0.05)
        assert_almost_equal(a.exp_data[1].q, 0.25)
        assert_almost_equal(a.exp_data[1].i, 2.)
        assert_almost_equal(a.exp_data[1].di, 0.2)
        assert_almost_equal(a.exp_data[1].dq, 0.05 * 0.25)
        assert_almost_equal(a.exp_data[2].q, 0.50)
        assert_almost_equal(a.exp_data[2].i, 1.)
        assert_almost_equal(a.exp_data[2].di, 0.1)
        assert_almost_equal(a.exp_data[2].dq, 0.05 * 0.50)
        assert_almost_equal(a.sim_data[0].q, 0.05)
        assert_almost_equal(a.sim_data[0].i, 1.)
        assert_almost_equal(a.sim_data[0].di, 0.1)
        assert_almost_equal(a.sim_data[0].dq, 0.05 * 0.05)
        assert_almost_equal(a.sim_data[1].q, 0.25)
        assert_almost_equal(a.sim_data[1].i, 2.)
        assert_almost_equal(a.sim_data[1].di, 0.2)
        assert_almost_equal(a.sim_data[1].dq, 0.05 * 0.25)
        assert_almost_equal(a.sim_data[2].q, 0.50)
        assert_almost_equal(a.sim_data[2].i, 3.)
        assert_almost_equal(a.sim_data[2].di, 0.3)
        assert_almost_equal(a.sim_data[2].dq, 0.05 * 0.50)
        assert_almost_equal(a.scale, 1.)
        assert_almost_equal(a.background, 0.)

    def test_change_scale(self):
        data1 = dataformat.QData(0.05, 3., 0.3, 0.05 * 0.05)
        data2 = dataformat.QData(0.25, 2., 0.2, 0.05 * 0.25)
        data3 = dataformat.QData(0.50, 1., 0.1, 0.05 * 0.50)
        data = [data1, data2, data3]
        sdata1 = dataformat.QData(0.05, 1., 0.1, 0.05 * 0.05)
        sdata2 = dataformat.QData(0.25, 2., 0.2, 0.05 * 0.25)
        sdata3 = dataformat.QData(0.50, 3., 0.3, 0.05 * 0.50)
        sdata = [sdata1, sdata2, sdata3]
        a = compare.Compare(data, sdata, 1., 0.)
        assert_almost_equal(a.scale, 1.)
        a.change_scale(2.)
        assert_almost_equal(a.scale, 2.)

    def test_change_background(self):
        data1 = dataformat.QData(0.05, 3., 0.3, 0.05 * 0.05)
        data2 = dataformat.QData(0.25, 2., 0.2, 0.05 * 0.25)
        data3 = dataformat.QData(0.50, 1., 0.1, 0.05 * 0.50)
        data = [data1, data2, data3]
        sdata1 = dataformat.QData(0.05, 1., 0.1, 0.05 * 0.05)
        sdata2 = dataformat.QData(0.25, 2., 0.2, 0.05 * 0.25)
        sdata3 = dataformat.QData(0.50, 3., 0.3, 0.05 * 0.50)
        sdata = [sdata1, sdata2, sdata3]
        a = compare.Compare(data, sdata, 1., 0.)
        assert_almost_equal(a.background, 0.)
        a.change_background(2.)
        assert_almost_equal(a.background, 2.)

def test_scale_and_background():
    a = np.array([1., 2., 3.])
    scale = 2.
    background = 1.
    b = compare.scale_and_background(a, scale, background)
    assert_almost_equal(b[0], 1.098612289)
    assert_almost_equal(b[1], 1.609437912)
    assert_almost_equal(b[2], 1.945910149)
