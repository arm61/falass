from numpy.testing import assert_equal, assert_almost_equal
from falass import readwrite, job, sld, dataformat
import os
import unittest


class TestSLD(unittest.TestCase):
    def test_sld(self):
        self.path = os.path.dirname(os.path.abspath(__file__))
        a = readwrite.Files(os.path.join(self.path, 'test.pdb'), lgtfile=os.path.join(self.path, 'test.lgt'),
                            datfile=os.path.join(self.path, 'test3.dat'))
        a.read_pdb()
        a.read_lgt()
        a.read_dat()
        b = job.Job(a, 1., 5.)
        b.set_times(times = [0., 20000., 10000.])
        b.set_lgts()
        c = sld.SLD(b)
        assert_equal(len(c.assigned_job.files.times), 6)
        assert_equal(c.assigned_job.files.times, [0, 10000., 20000., 30000., 40000., 50000.])
        assert_equal(len(c.assigned_job.times), 3)
        assert_equal(c.assigned_job.times, [0, 10000., 20000.])
        assert_equal(c.assigned_job.layer_thickness, 1.)
        assert_equal(c.assigned_job.cut_off_size, 5.)
        return

    def test_set_sld_profile(self):
        c = None
        a = dataformat.SLDPro(5., 5., 0.)
        b = dataformat.SLDPro(5., 0., 0.)
        slda = [a, b]
        a = sld.SLD(c)
        a.set_sld_profile(slda)
        assert_almost_equal(a.sld_profile[0].thick, 5.)
        assert_almost_equal(a.sld_profile[0].real, 5.)
        assert_almost_equal(a.sld_profile[0].imag, 0.)
        assert_almost_equal(a.sld_profile[1].thick, 5.)
        assert_almost_equal(a.sld_profile[1].real, 0.)
        assert_almost_equal(a.sld_profile[1].imag, 0.)
        return

    def test_set_av_sld_profile(self):
        c = None
        a = dataformat.SLDPro(5., 5., 0.)
        b = dataformat.SLDPro(5., 0., 0.)
        slda = [a, b]
        a = dataformat.SLDPro(5., 2., 0.)
        b = dataformat.SLDPro(5., 0., 0.)
        slda_err = [a, b]
        a = sld.SLD(c)
        a.set_av_sld_profile(slda, slda_err)
        assert_almost_equal(a.av_sld_profile[0].thick, 5.)
        assert_almost_equal(a.av_sld_profile[0].real, 5.)
        assert_almost_equal(a.av_sld_profile[0].imag, 0.)
        assert_almost_equal(a.av_sld_profile[1].thick, 5.)
        assert_almost_equal(a.av_sld_profile[1].real, 0.)
        assert_almost_equal(a.av_sld_profile[1].imag, 0.)
        assert_almost_equal(a.av_sld_profile_err[0].thick, 5.)
        assert_almost_equal(a.av_sld_profile_err[0].real, 2.)
        assert_almost_equal(a.av_sld_profile_err[0].imag, 0.)
        assert_almost_equal(a.av_sld_profile_err[1].thick, 5.)
        assert_almost_equal(a.av_sld_profile_err[1].real, 0.)
        assert_almost_equal(a.av_sld_profile_err[1].imag, 0.)
        return

    def test_get_sld_profile(self):
        self.path = os.path.dirname(os.path.abspath(__file__))
        a = readwrite.Files(os.path.join(self.path, 'test.pdb'), lgtfile=os.path.join(self.path, 'test.lgt'),
                            datfile=os.path.join(self.path, 'test3.dat'))
        a.read_pdb()
        a.read_lgt()
        a.read_dat()
        b = job.Job(a, 1., 0.)
        b.set_times(times=[0., 20000., 10000.])
        b.set_lgts()
        c = sld.SLD(b)
        c.get_sld_profile()
        assert_equal(len(c.sld_profile), 3)
        for i in range(0, len(c.sld_profile)):
            assert_equal(len(c.sld_profile[i]), 4)
        assert_almost_equal(c.sld_profile[0][0].thick, 1.)
        assert_almost_equal(c.sld_profile[0][0].real, 0.)
        assert_almost_equal(c.sld_profile[0][0].imag, 0.)
        assert_almost_equal(c.sld_profile[0][1].thick, 1.)
        assert_almost_equal(c.sld_profile[0][1].real, 1e-5)
        assert_almost_equal(c.sld_profile[0][1].imag, 0.)
        assert_almost_equal(c.sld_profile[0][2].thick, 1.)
        assert_almost_equal(c.sld_profile[0][2].real, 2e-5)
        assert_almost_equal(c.sld_profile[0][2].imag, 1e-5)
        assert_almost_equal(c.sld_profile[0][3].thick, 1.)
        assert_almost_equal(c.sld_profile[0][3].real, 3e-5)
        assert_almost_equal(c.sld_profile[0][3].imag, 2e-5)
        assert_almost_equal(c.sld_profile[1][0].thick, 1.)
        assert_almost_equal(c.sld_profile[1][0].real, 0. / 2.)
        assert_almost_equal(c.sld_profile[1][0].imag, 0. / 2.)
        assert_almost_equal(c.sld_profile[1][1].thick, 1.)
        assert_almost_equal(c.sld_profile[1][1].real, 2e-5 / 2.)
        assert_almost_equal(c.sld_profile[1][1].imag, 1e-5 / 2.)
        assert_almost_equal(c.sld_profile[1][2].thick, 1.)
        assert_almost_equal(c.sld_profile[1][2].real, 3e-5 / 2.)
        assert_almost_equal(c.sld_profile[1][2].imag, 2e-5 / 2.)
        assert_almost_equal(c.sld_profile[1][3].thick, 1.)
        assert_almost_equal(c.sld_profile[1][3].real, 1e-5 / 2.)
        assert_almost_equal(c.sld_profile[1][3].imag, 0 / 2.)
        assert_almost_equal(c.sld_profile[2][0].thick, 1.)
        assert_almost_equal(c.sld_profile[2][0].real, 0. / 3.)
        assert_almost_equal(c.sld_profile[2][0].imag, 0. / 3.)
        assert_almost_equal(c.sld_profile[2][1].thick, 1.)
        assert_almost_equal(c.sld_profile[2][1].real, 3e-5 / 3.)
        assert_almost_equal(c.sld_profile[2][1].imag, 2e-5 / 3.)
        assert_almost_equal(c.sld_profile[2][2].thick, 1.)
        assert_almost_equal(c.sld_profile[2][2].real, 1e-5 / 3.)
        assert_almost_equal(c.sld_profile[2][2].imag, 0 / 3.)
        assert_almost_equal(c.sld_profile[2][3].thick, 1.)
        assert_almost_equal(c.sld_profile[2][3].real, 2e-5 / 3.)
        assert_almost_equal(c.sld_profile[2][3].imag, 1e-5 / 3.)
        return

    def test_average_sld_profile(self):
        self.path = os.path.dirname(os.path.abspath(__file__))
        a = readwrite.Files(os.path.join(self.path, 'test.pdb'), lgtfile=os.path.join(self.path, 'test.lgt'),
                            datfile=os.path.join(self.path, 'test3.dat'))
        a.read_pdb()
        a.read_lgt()
        a.read_dat()
        b = job.Job(a, 1., 0.)
        b.set_times(times=[0., 20000., 10000.])
        b.set_lgts()
        c = sld.SLD(b)
        c.get_sld_profile()
        c.average_sld_profile()
        assert_equal(len(c.av_sld_profile), 4)
        assert_almost_equal(c.av_sld_profile[0].real, 0.)
        assert_almost_equal(c.av_sld_profile[0].imag, 0.)
        assert_almost_equal(c.av_sld_profile[1].real, (1e-5 + (2e-5 / 2.) + (3e-5 / 3.)) / 3.)
        assert_almost_equal(c.av_sld_profile[1].imag, (0 + (1e-5 / 2.) + (2e-5 / 3.)) / 3.)
        assert_almost_equal(c.av_sld_profile[2].real, (2e-5 + (3e-5 / 2.) + (1e-5 / 3.)) / 3.)
        assert_almost_equal(c.av_sld_profile[2].imag, ((1e-5) + (2e-5 / 2.) + (0 / 3.)) / 3.)
        assert_almost_equal(c.av_sld_profile[3].real, (3e-5 + (1e-5 / 2.) + (2e-5 / 3.)) / 3.)
        assert_almost_equal(c.av_sld_profile[3].imag, ((2e-5) + (0 / 2.) + (1e-5 / 3.)) / 3.)
        return


    def test_get_scatlen(self):
        atom1 = dataformat.ScatLens('C1', 1.0, 0.0)
        atom2 = dataformat.ScatLens('C2', 2.0, 1.0)
        atom3 = dataformat.ScatLens('C3', 3.0, 2.0)
        array = [atom1, atom2, atom3]
        real, imag = sld.get_scatlen('C3', array)
        assert_almost_equal(real, 3.0e-5)
        assert_almost_equal(imag, 2.0e-5)
        return


    def test_get_scatlen_fail(self):
        atom1 = dataformat.ScatLens('C1', 1.0, 0.0)
        atom2 = dataformat.ScatLens('C2', 2.0, 1.0)
        atom3 = dataformat.ScatLens('C3', 3.0, 2.0)
        array = [atom1, atom2, atom3]
        with self.assertRaises(ValueError) as context:
            sld.get_scatlen('C4', array)
        self.assertTrue("Attempt to get the scattering length of the atom type {} failed. This should never "
                                     "happen. Please contact the developers".format('C4') in str(context.exception))
        return
