from numpy.testing import assert_almost_equal
from falass import compare, dataformat, reflect
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

    def test_fit_noq(self):
        data11 = dataformat.QData(0.05, 0.1, 0., 0.05 * 0.05)
        data12 = dataformat.QData(0.25, 0.05, 0., 0.05 * 0.25)
        data13 = dataformat.QData(0.50, 0.01, 0., 0.05 * 0.50)
        data1 = [data11, data12, data13]
        data21 = dataformat.QData(0.05, 0.12, 0., 0.05 * 0.05)
        data22 = dataformat.QData(0.25, 0.03, 0., 0.05 * 0.25)
        data23 = dataformat.QData(0.50, 0.015, 0., 0.05 * 0.50)
        data2 = [data21, data22, data23]
        data31 = dataformat.QData(0.05, 0.14, 0., 0.05 * 0.05)
        data32 = dataformat.QData(0.25, 0.07, 0., 0.05 * 0.25)
        data33 = dataformat.QData(0.50, 0.005, 0., 0.05 * 0.50)
        data3 = [data31, data32, data33]
        layer1 = dataformat.SLDPro(1., 0., 0.)
        layer2 = dataformat.SLDPro(1., 5., 0.)
        sld = [[layer1, layer2]]
        ddata1 = dataformat.QData(0.05, 0., 0., 0.05 * 0.05)
        ddata2 = dataformat.QData(0.25, 0., 0., 0.05 * 0.25)
        ddata3 = dataformat.QData(0.50, 0., 0., 0.05 * 0.50)
        ddata = [ddata1, ddata2, ddata3]
        a = reflect.Reflect(sld, ddata)
        a.reflect = [data1, data2, data3]
        a.average_ref()
        ddata = []
        k = compare.Compare(ddata, a.averagereflect, 1, 0)
        with self.assertRaises(ValueError) as context:
            k.fit()
        self.assertTrue('No q vectors have been defined -- either read a .dat file or get q vectors.' in str(context.exception))

    def test_fit_noi(self):
        data11 = dataformat.QData(0.05, 0.1, 0., 0.05 * 0.05)
        data12 = dataformat.QData(0.25, 0.05, 0., 0.05 * 0.25)
        data13 = dataformat.QData(0.50, 0.01, 0., 0.05 * 0.50)
        data1 = [data11, data12, data13]
        data21 = dataformat.QData(0.05, 0.12, 0., 0.05 * 0.05)
        data22 = dataformat.QData(0.25, 0.03, 0., 0.05 * 0.25)
        data23 = dataformat.QData(0.50, 0.015, 0., 0.05 * 0.50)
        data2 = [data21, data22, data23]
        data31 = dataformat.QData(0.05, 0.14, 0., 0.05 * 0.05)
        data32 = dataformat.QData(0.25, 0.07, 0., 0.05 * 0.25)
        data33 = dataformat.QData(0.50, 0.005, 0., 0.05 * 0.50)
        data3 = [data31, data32, data33]
        layer1 = dataformat.SLDPro(1., 0., 0.)
        layer2 = dataformat.SLDPro(1., 5., 0.)
        sld = [[layer1, layer2]]
        ddata1 = dataformat.QData(0.05, 0., 0., 0.05 * 0.05)
        ddata2 = dataformat.QData(0.25, 0., 0., 0.05 * 0.25)
        ddata3 = dataformat.QData(0.50, 0., 0., 0.05 * 0.50)
        ddata = [ddata1, ddata2, ddata3]
        a = reflect.Reflect(sld, ddata)
        a.reflect = [data1, data2, data3]
        a.average_ref()
        ddata1 = dataformat.QData(0.05, None, 0., 0.05 * 0.05)
        ddata2 = dataformat.QData(0.25, 0., 0., 0.05 * 0.25)
        ddata3 = dataformat.QData(0.50, 0., 0., 0.05 * 0.50)
        ddata = [ddata1, ddata2, ddata3]
        k = compare.Compare(ddata, a.averagereflect, 1, 0)
        with self.assertRaises(ValueError) as context:
            k.fit()
        self.assertTrue('No experimental data has been set for comparison, please read in a a .dat file.' in str(context.exception))

    def test_return_fitted(self):
        data1 = dataformat.QData(0.05, 3., 0.3, 0.05 * 0.05)
        data2 = dataformat.QData(0.25, 2., 0.2, 0.05 * 0.25)
        data3 = dataformat.QData(0.50, 1., 0.1, 0.05 * 0.50)
        data = [data1, data2, data3]
        sdata1 = dataformat.QData(0.05, 1., 0.1, 0.05 * 0.05)
        sdata2 = dataformat.QData(0.25, 2., 0.2, 0.05 * 0.25)
        sdata3 = dataformat.QData(0.50, 3., 0.3, 0.05 * 0.50)
        sdata = [sdata1, sdata2, sdata3]
        a = compare.Compare(data, sdata, 2., 0.)
        a.return_fitted()
        assert_almost_equal(a.sim_data_fitted[0].i, 2.)
        assert_almost_equal(a.sim_data_fitted[1].i, 4.)
        assert_almost_equal(a.sim_data_fitted[2].i, 6.)
        assert_almost_equal(a.sim_data_fitted[0].di, 0.2)
        assert_almost_equal(a.sim_data_fitted[1].di, 0.4)
        assert_almost_equal(a.sim_data_fitted[2].di, 0.6)

    def test_scale_and_background(self):
        a = np.array([1., 2., 3.])
        scale = 2.
        background = 1.
        b = compare.scale_and_background(a, scale, background)
        assert_almost_equal(b[0], 1.098612289)
        assert_almost_equal(b[1], 1.609437912)
        assert_almost_equal(b[2], 1.945910149)
