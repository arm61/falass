from numpy.testing import assert_almost_equal
from falass import dataformat, reflect
import unittest


class TestReflect(unittest.TestCase):
    def test_reflect(self):
        layer1 = dataformat.SLDPro(1., 0., 0.)
        layer2 = dataformat.SLDPro(1., 5., 0.)
        sld = [[layer1, layer2]]
        data1 = dataformat.QData(0.05, 0., 0., 0.05 * 0.05)
        data2 = dataformat.QData(0.25, 0., 0., 0.05 * 0.25)
        data3 = dataformat.QData(0.50, 0., 0., 0.05 * 0.50)
        data = [data1, data2, data3]
        a = reflect.Reflect(sld, data)
        assert_almost_equal(a.sld_profile[0][0].thick, 1.)
        assert_almost_equal(a.sld_profile[0][0].real, 0.)
        assert_almost_equal(a.sld_profile[0][0].imag, 0.)
        assert_almost_equal(a.sld_profile[0][1].thick, 1.)
        assert_almost_equal(a.sld_profile[0][1].real, 5.)
        assert_almost_equal(a.sld_profile[0][1].imag, 0.)
        assert_almost_equal(a.exp_data[0].q, 0.05)
        assert_almost_equal(a.exp_data[0].dq, 0.05 * 0.05)
        assert_almost_equal(a.exp_data[1].q, 0.25)
        assert_almost_equal(a.exp_data[1].dq, 0.25 * 0.05)
        assert_almost_equal(a.exp_data[2].q, 0.50)
        assert_almost_equal(a.exp_data[2].dq, 0.50 * 0.05)

    def test_average_reflect(self):
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
        assert_almost_equal(a.averagereflect[0].i, 0.12)
        assert_almost_equal(a.averagereflect[1].i, 0.05)
        assert_almost_equal(a.averagereflect[2].i, 0.01)

