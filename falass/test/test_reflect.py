from numpy.testing import assert_almost_equal, assert_equal
from falass import dataformat, reflect
import unittest
import numpy as np


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

    def test_calc_ref_basic(self):
        layer1 = dataformat.SLDPro(1., 0., 0.)
        layer2 = dataformat.SLDPro(1., 5., 0.)
        sld = [[layer1, layer2]]
        data1 = dataformat.QData(0.05, 0., 0., 0.05 * 0.05)
        data2 = dataformat.QData(0.25, 0., 0., 0.05 * 0.25)
        data3 = dataformat.QData(0.50, 0., 0., 0.05 * 0.50)
        data = [data1, data2, data3]
        a = reflect.Reflect(sld, data)
        a.calc_ref()
        assert_equal(len(a.reflect), 1)
        assert_equal(len(a.reflect[0]), 3)

    def test_calc_reflect_noq(self):
        layer1 = dataformat.SLDPro(1., 0., 0.)
        layer2 = dataformat.SLDPro(1., 5., 0.)
        sld = [[layer1, layer2]]
        ddata = []
        a = reflect.Reflect(sld, ddata)
        with self.assertRaises(ValueError) as context:
            a.calc_ref()
        self.assertTrue('No q vectors have been defined -- either read a .dat file or get q vectors.' in str(context.exception))

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

    def test_average_reflect_noq(self):
        layer1 = dataformat.SLDPro(1., 0., 0.)
        layer2 = dataformat.SLDPro(1., 5., 0.)
        sld = [[layer1, layer2]]
        ddata = []
        a = reflect.Reflect(sld, ddata)
        with self.assertRaises(ValueError) as context:
            a.average_ref()
        self.assertTrue('No q vectors have been defined -- either read a .dat file or get q vectors.' in str(context.exception))


    def test_make_kn(self):
        layer1 = dataformat.SLDPro(5., 5., 0.)
        layer2 = dataformat.SLDPro(5., 0., 0.)
        sld_profile = [layer1, layer2]
        exp_data = [0.05, 0.25, 0.5]
        layers = np.zeros((len(sld_profile), 4))
        for i in range(0, len(sld_profile)):
            layers[i][0] = sld_profile[i].thick
            layers[i][1] = sld_profile[i].real
            layers[i][2] = sld_profile[i].imag
            layers[i][3] = 0
        qvals = np.asfarray(exp_data).ravel()
        nlayers = len(sld_profile) - 2
        npnts = qvals.size
        kn = reflect.make_kn(npnts, nlayers, layers, qvals)
        assert_almost_equal(kn, [[0.025 + 0j, 7.9266940191 + 0j], [0.125 + 0j, 7.927640133 + 0j],
                                 [0.25 + 0j, 7.93059601 + 0j]])


    def test_knext_and_rj(self):
        kn = np.array([[0.025 + 0j, 7.9266940191 + 0j], [0.125 + 0j, 7.927640133 + 0j],
                       [0.25 + 0j, 7.93059601 + 0j]])
        k = np.array([0.025 + 0j, 0.125 + 0j, 0.25 + 0j])
        idx = 1
        k_next, rj = reflect.knext_and_rj(kn, idx, k)
        assert_almost_equal(k_next, np.array([7.9266940191 + 0j,  7.927640133 + 0j, 7.93059601 + 0j]))
        assert_almost_equal(rj, np.array([-1.211500325 + 0j, -2.610174734 + 0j, -6.8181020565 + 0j]))
