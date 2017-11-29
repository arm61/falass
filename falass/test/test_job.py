from numpy.testing import assert_equal
from falass import readwrite, job
import os
import unittest

class TestJob(unittest.TestCase):
    def test_job(self):
        self.path = os.path.dirname(os.path.abspath(__file__))
        a = readwrite.Files(os.path.join(self.path, 'test.pdb'), lgtfile=os.path.join(self.path, 'test.lgt'),
                              datfile=os.path.join(self.path, 'test3.dat'))
        a.read_pdb()
        a.read_lgt()
        a.read_dat()
        b = job.Job(a, 1., 5.)
        assert_equal(b.layer_thickness, 1.)
        assert_equal(b.cut_off_size, 5.)
        assert_equal(b.times, [0., 10000., 20000., 30000., 40000., 50000.])

    def test_set_run(self):
        self.path = os.path.dirname(os.path.abspath(__file__))
        a = readwrite.Files(os.path.join(self.path, 'test.pdb'), lgtfile=os.path.join(self.path, 'test.lgt'),
                            datfile=os.path.join(self.path, 'test3.dat'))
        a.read_pdb()
        a.read_lgt()
        a.read_dat()
        b = job.Job(a, 1., 5.)
        self.path = os.path.dirname(os.path.abspath(__file__))
        a2 = readwrite.Files(os.path.join(self.path, 'test2.pdb'), lgtfile=os.path.join(self.path, 'test.lgt'),
                            datfile=os.path.join(self.path, 'test3.dat'))
        a2.read_pdb()
        a2.read_lgt()
        a2.read_dat()
        b.set_run(files=a2, layer_thickness=2., cut_off_size=3.)
        assert_equal(b.layer_thickness, 2.)
        assert_equal(b.cut_off_size, 3.)
        assert_equal(b.times, [0., 10000., 20000., 30000., 40000.])

    def test_set_lgts(self):
        self.path = os.path.dirname(os.path.abspath(__file__))
        a = readwrite.Files(os.path.join(self.path, 'test.pdb'), lgtfile=os.path.join(self.path, 'test.lgt'))
        a.read_pdb()
        a.read_lgt()
        b = job.Job(a, 1., 5.)
        b.set_lgts()
        assert_equal(b.new_file, False)

    def test_set_times(self):
        self.path = os.path.dirname(os.path.abspath(__file__))
        a = readwrite.Files(os.path.join(self.path, 'test.pdb'), lgtfile=os.path.join(self.path, 'test.lgt'))
        a.read_pdb()
        a.read_lgt()
        b = job.Job(a, 1., 5.)
        b.set_times([0., 20000., 10000.])
        assert_equal(len(b.times), 3)
        assert_equal(b.times, [0., 10000., 20000.])

def test_check_array_true():
    array = [0, 1, 2, 3, 4]
    check = 1
    bool_ret = job.check_array(array, check)
    assert_equal(bool_ret, True)

def test_check_array_false():
    array = [0, 1, 2, 3, 4]
    check = 6
    bool_ret = job.check_array(array, check)
    assert_equal(bool_ret, False)
