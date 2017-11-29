from numpy.testing import assert_equal, assert_almost_equal
from falass import readwrite, dataformat
import os
import unittest


class TestFiles(unittest.TestCase):
    def test_files_standard(self):
        pdb = readwrite.Files('test.pdb')
        assert_equal(pdb.pdbfile, 'test.pdb')
        assert_equal(pdb.lgtfile, None)
        assert_equal(pdb.datfile, None)
        assert_equal(pdb.resolution, 5.)
        assert_equal(pdb.ierror, 5.)
        assert_equal(pdb.flip, False)
        return

    def test_files_variation(self):
        pdb = readwrite.Files('test1.pdb', lgtfile='test1.lgt', datfile='test1.dat', resolution=1., ierror=1.,
                              flip=True)
        assert_equal(pdb.pdbfile, 'test1.pdb')
        assert_equal(pdb.lgtfile, 'test1.lgt')
        assert_equal(pdb.datfile, 'test1.dat')
        assert_equal(pdb.resolution, 1.)
        assert_equal(pdb.ierror, 1.)
        assert_equal(pdb.flip, True)
        return

    def test_set_pdbfile(self):
        pdb = readwrite.Files('test.pdb')
        assert_equal(pdb.pdbfile, 'test.pdb')
        pdb.set_file(pdbfile='test1.pdb')
        assert_equal(pdb.pdbfile, 'test1.pdb')
        return

    def test_set_lgtfile(self):
        pdb = readwrite.Files('test.pdb')
        assert_equal(pdb.lgtfile, None)
        pdb.set_file(lgtfile='test1.lgt')
        assert_equal(pdb.lgtfile, 'test1.lgt')
        return

    def test_set_datfile(self):
        pdb = readwrite.Files('test.pdb')
        assert_equal(pdb.datfile, None)
        pdb.set_file(datfile='test1.dat')
        assert_equal(pdb.datfile, 'test1.dat')
        return

    def test_read_pdb(self):
        self.path = os.path.dirname(os.path.abspath(__file__))
        pdb = readwrite.Files(os.path.join(self.path, 'test.pdb'))
        pdb.read_pdb()
        assert_equal(pdb.number_of_timesteps, 6)
        assert_equal(pdb.times, [0., 10000., 20000., 30000., 40000., 50000.])
        assert_equal(len(pdb.atoms), 6)
        for i in range(0, len(pdb.atoms)):
            assert_equal(len(pdb.atoms[i]), 3)
            assert_equal(pdb.atoms[i][0].atom, 'C1')
            assert_equal(pdb.atoms[i][1].atom, 'C2')
            assert_equal(pdb.atoms[i][2].atom, 'C3')
        for i in range(0, 2):
            assert_equal(pdb.atoms[i * 3][0].zpos, 01.500)
            assert_equal(pdb.atoms[i * 3][1].zpos, 02.500)
            assert_equal(pdb.atoms[i * 3][2].zpos, 03.500)
            assert_equal(pdb.atoms[i * 3 + 1][0].zpos, 03.500)
            assert_equal(pdb.atoms[i * 3 + 1][1].zpos, 01.500)
            assert_equal(pdb.atoms[i * 3 + 1][2].zpos, 02.500)
            assert_equal(pdb.atoms[i * 3 + 2][0].zpos, 02.500)
            assert_equal(pdb.atoms[i * 3 + 2][1].zpos, 03.500)
            assert_equal(pdb.atoms[i * 3 + 2][2].zpos, 01.500)
        assert_equal(pdb.cell, [[1.000, 1.000, 4.000], [2.000, 1.000, 4.000], [3.000, 1.000, 4.000],
                                [4.000, 1.000, 4.000], [5.000, 1.000, 4.000], [6.000, 1.000, 4.000]])
        return

    def test_read_pdb_flip(self):
        self.path = os.path.dirname(os.path.abspath(__file__))
        pdb = readwrite.Files(os.path.join(self.path, 'test.pdb'), flip=True)
        pdb.read_pdb()
        assert_equal(pdb.number_of_timesteps, 6)
        assert_equal(pdb.times, [0., 10000., 20000., 30000., 40000., 50000.])
        assert_equal(len(pdb.atoms), 6)
        for i in range(0, len(pdb.atoms)):
            assert_equal(len(pdb.atoms[i]), 3)
            assert_equal(pdb.atoms[i][0].atom, 'C1')
            assert_equal(pdb.atoms[i][1].atom, 'C2')
            assert_equal(pdb.atoms[i][2].atom, 'C3')
        for i in range(0, 2):
            assert_equal(pdb.atoms[i * 3][0].zpos, 02.500)
            assert_equal(pdb.atoms[i * 3][1].zpos, 01.500)
            assert_equal(pdb.atoms[i * 3][2].zpos, 00.500)
            assert_equal(pdb.atoms[i * 3 + 1][0].zpos, 00.500)
            assert_equal(pdb.atoms[i * 3 + 1][1].zpos, 02.500)
            assert_equal(pdb.atoms[i * 3 + 1][2].zpos, 01.500)
            assert_equal(pdb.atoms[i * 3 + 2][0].zpos, 01.500)
            assert_equal(pdb.atoms[i * 3 + 2][1].zpos, 00.500)
            assert_equal(pdb.atoms[i * 3 + 2][2].zpos, 02.500)
        assert_equal(pdb.cell, [[1.000, 1.000, 4.000], [2.000, 1.000, 4.000], [3.000, 1.000, 4.000],
                                [4.000, 1.000, 4.000], [5.000, 1.000, 4.000], [6.000, 1.000, 4.000]])
        return

    def test_read_lgt(self):
        self.path = os.path.dirname(os.path.abspath(__file__))
        pdb = readwrite.Files('test.pdb', lgtfile=os.path.join(self.path, 'test.lgt'))
        pdb.read_lgt()
        assert_equal(pdb.scat_lens[0].atom, 'C1')
        assert_almost_equal(pdb.scat_lens[0].real, 1e-5)
        assert_almost_equal(pdb.scat_lens[0].imag, 0)
        assert_equal(pdb.scat_lens[1].atom, 'C2')
        assert_almost_equal(pdb.scat_lens[1].real, 2e-5)
        assert_almost_equal(pdb.scat_lens[1].imag, 1e-5)
        assert_equal(pdb.scat_lens[2].atom, 'C3')
        assert_almost_equal(pdb.scat_lens[2].real, 3e-5)
        assert_almost_equal(pdb.scat_lens[2].imag, 2e-5)
        return

    def test_read_lgt_not_defined(self):
        pdb = readwrite.Files('test.pdb')
        with self.assertRaises(ValueError) as context:
            pdb.read_lgt()
        self.assertTrue("No lgtfile has been defined" in str(context.exception))
        return

    def test_read_dat2(self):
        self.path = os.path.dirname(os.path.abspath(__file__))
        pdb = readwrite.Files('test.pdb', datfile=os.path.join(self.path, 'test2.dat'))
        pdb.read_dat()
        assert_almost_equal(pdb.expdata[0].q, 0.051793)
        assert_almost_equal(pdb.expdata[0].i, 0.00034985)
        assert_almost_equal(pdb.expdata[0].di, 0.00034985 * 0.05)
        assert_almost_equal(pdb.expdata[0].dq, 0.051793 * 0.05)
        assert_almost_equal(pdb.expdata[1].q, 0.13742)
        assert_almost_equal(pdb.expdata[1].i, 1.4924e-5)
        assert_almost_equal(pdb.expdata[1].di, 1.4924e-5 * 0.05)
        assert_almost_equal(pdb.expdata[1].dq, 0.13742 * 0.05)
        assert_almost_equal(pdb.expdata[2].q, 0.38285)
        assert_almost_equal(pdb.expdata[2].i, 1.215e-6)
        assert_almost_equal(pdb.expdata[2].di, 1.215e-6 * 0.05)
        assert_almost_equal(pdb.expdata[2].dq, 0.38285 * 0.05)
        return

    def test_read_dat3(self):
        self.path = os.path.dirname(os.path.abspath(__file__))
        pdb = readwrite.Files('test.pdb', datfile=os.path.join(self.path, 'test3.dat'))
        pdb.read_dat()
        assert_almost_equal(pdb.expdata[0].q, 0.051793)
        assert_almost_equal(pdb.expdata[0].i, 0.00034985)
        assert_almost_equal(pdb.expdata[0].di, 5.0354e-6)
        assert_almost_equal(pdb.expdata[0].dq, 0.051793 * 0.05)
        assert_almost_equal(pdb.expdata[1].q, 0.13742)
        assert_almost_equal(pdb.expdata[1].i, 1.4924e-5)
        assert_almost_equal(pdb.expdata[1].di, 3.7729e-7)
        assert_almost_equal(pdb.expdata[1].dq, 0.13742 * 0.05)
        assert_almost_equal(pdb.expdata[2].q, 0.38285)
        assert_almost_equal(pdb.expdata[2].i, 1.215e-6)
        assert_almost_equal(pdb.expdata[2].di, 3.478e-7)
        assert_almost_equal(pdb.expdata[2].dq, 0.38285 * 0.05)
        return

    def test_read_dat4(self):
        self.path = os.path.dirname(os.path.abspath(__file__))
        pdb = readwrite.Files('test.pdb', datfile=os.path.join(self.path, 'test4.dat'))
        pdb.read_dat()
        assert_almost_equal(pdb.expdata[0].q, 0.051793)
        assert_almost_equal(pdb.expdata[0].i, 0.00034985)
        assert_almost_equal(pdb.expdata[0].di, 5.0354e-6)
        assert_almost_equal(pdb.expdata[0].dq, 0.0003)
        assert_almost_equal(pdb.expdata[1].q, 0.13742)
        assert_almost_equal(pdb.expdata[1].i, 1.4924e-5)
        assert_almost_equal(pdb.expdata[1].di, 3.7729e-7)
        assert_almost_equal(pdb.expdata[1].dq, 0.0007)
        assert_almost_equal(pdb.expdata[2].q, 0.38285)
        assert_almost_equal(pdb.expdata[2].i, 1.215e-6)
        assert_almost_equal(pdb.expdata[2].di, 3.478e-7)
        assert_almost_equal(pdb.expdata[2].dq, 0.0019)
        return

    def test_read_dat_not_defined(self):
        pdb = readwrite.Files('test.pdb')
        assert_equal(pdb.expdata, [])
        return

    def test_get_qs(self):
        pdb = readwrite.Files('test.pdb')
        pdb.get_qs(0.0005, 0.0020, 4)
        assert_equal(len(pdb.expdata), 4)
        assert_equal(pdb.expdata[0].q, 0.0005)
        assert_almost_equal(pdb.expdata[0].dq, 0.0005 * 0.05)
        assert_equal(pdb.expdata[1].q, 0.0010)
        assert_almost_equal(pdb.expdata[1].dq, 0.0010 * 0.05)
        assert_equal(pdb.expdata[2].q, 0.0015)
        assert_almost_equal(pdb.expdata[2].dq, 0.0015 * 0.05)
        assert_equal(pdb.expdata[3].q, 0.0020)
        assert_almost_equal(pdb.expdata[3].dq, 0.0020 * 0.05)
        return

    if __name__ == "__main__":
        unittest.main()


def test_check_duplicates_true():
    atom1 = dataformat.ScatLens('C1', 1.0, 0.0)
    atom2 = dataformat.ScatLens('C2', 2.0, 1.0)
    atom3 = dataformat.ScatLens('C3', 3.0, 2.0)
    array = [atom1, atom2, atom3]
    check = dataformat.ScatLens('C2', 2.0, 1.0)
    bool_ret = readwrite.check_duplicates(array, check.atom)
    assert_equal(bool_ret, True)
    return


def test_check_duplicates_false():
    atom1 = dataformat.ScatLens('C1', 1.0, 0.0)
    atom2 = dataformat.ScatLens('C2', 2.0, 1.0)
    atom3 = dataformat.ScatLens('C3', 3.0, 2.0)
    array = [atom1, atom2, atom3]
    check = dataformat.ScatLens('C4', 4.0, 3.0)
    bool_ret = readwrite.check_duplicates(array, check.atom)
    assert_equal(bool_ret, False)
    return


def test_line_count():
    path = os.path.dirname(os.path.abspath(__file__))
    lines = readwrite.line_count(os.path.join(path, 'test.pdb'))
    assert_equal(lines, 60)
    return
