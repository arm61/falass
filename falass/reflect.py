import numpy as np
from falass import dataformat
import matplotlib.pyplot as plt
from scipy.interpolate import InterpolatedUnivariateSpline

class Reflect(object):
    def __init__(self, sld_profile, exp_data, job):
        self.sld_profile = sld_profile
        self.exp_data = exp_data
        self.job = job
        self.averagereflect = []
        self.reflect = []

    def calcRef(self):
        self.reflect = []
        for i in range(0, len(self.sld_profile)):
            refl = convolution(self.exp_data, self.sld_profile[i], self.job.layer_thickness)
            a = []
            for j in range(0, len(self.exp_data)):
                a.append(dataformat.datastruct(self.exp_data[j].q, refl[j], 0, self.exp_data[j].dq))
            self.reflect.append(a)


    def averageRef(self):
        self.averagereflect = []
        for i in range(0, len(self.reflect[0])):
            self.averagereflect.append(dataformat.datastruct(self.exp_data[i].q, 0, 0, self.exp_data[i].dq))
        for i in range(0, len(self.reflect[0])):
            for j in range(0, len(self.reflect)):
                self.averagereflect[i].i += self.reflect[j][i].i
            self.averagereflect[i].i /= len(self.job.times)
            for j in range(0, len(self.reflect)):
                self.averagereflect[i].di += np.square(self.reflect[j][i].i - self.averagereflect[i].i)
            self.averagereflect[i].di = np.sqrt(1. / (len(self.job.times) - 1) * self.averagereflect[i].di)

    def plotaverageRef(self):
        x = []
        y = []
        dy = []
        plt.rc('text', usetex=True)
        plt.rc('font', family='serif')
        for i in range(0, len(self.exp_data)):
            x.append(self.averagereflect[i].q)
            y.append(self.averagereflect[i].i)
            dy.append(self.averagereflect[i].di)
        x = np.asarray(x)
        y = np.asarray(y)
        dy = np.asarray(dy)
        plt.errorbar(x, y * x ** 4, yerr=dy * x ** 4)
        plt.xlabel('$q$ (\AA)')
        plt.ylabel('$Rq^4$ (\AA$^4$)')
        plt.yscale('log')
        plt.show()


def convolution(exp_data, sld_profile, lt):
    fwhm = 2 * np.sqrt(2 * np.log(2))

    res = exp_data[0].dq / exp_data[0].q

    if exp_data[0].dq / exp_data[0].q < 0.0005:
        return reflectivity(exp_data, sld_profile, lt)

    gnum = 51
    ggpoint = (gnum - 1) / 2

    def gauss(x, s):
        return 1. / s / np.sqrt(2 * np.pi) * np.exp(-0.5 * x ** 2 /s /s)

    q = []
    for i in range(0, len(exp_data)):
        q.append(exp_data[i].q)

    lowq = np.min(q)
    highq = np.max(q)

    start = np.log10(lowq) - 6 * res / fwhm
    finish = np.log10(highq * (1 + 6 * res / fwhm))
    interp = np.round(np.abs(1 * (np.abs(start-finish)) / (1.7 * res / fwhm / ggpoint)))

    xtemp = np.linspace(start, finish, int(interp))
    xlin = np.power(10., xtemp)

    gaussx = np.linspace(-1.7 * res, 1.7 * res, gnum)
    gaussy = gauss(gaussx, res / fwhm)

    rvals = reflectivity(xlin, sld_profile, lt)
    smeared_rvals = np.convolve(rvals, gaussy, mode = 'same')
    interpol = InterpolatedUnivariateSpline(xlin, smeared_rvals)

    smeared_output = interpol(q)
    smeared_output *= gaussx[1] - gaussx[0]
    return smeared_output


def reflectivity(exp_data, sld_profile, lt):
    bn = np.zeros((2, 2), dtype='complex')
    M = np.zeros((2, 2), dtype='complex')
    bm = np.zeros((2, 2), dtype='complex')
    refl = []
    for a in range(0, len(exp_data)):
        kz = exp_data[a] / 2 + 0j
        kn = kz
        bn[0][0] = 1 + 0j
        bn[0][1] = 0 + 0j
        bn[1][0] = 0 + 0j
        bn[1][1] = 1 + 0j
        for it in range(0, len(sld_profile) - 1):
            kn1_init_layer = sld_profile[it].real + sld_profile[it].imag * 1.j
            kn1_init_super = sld_profile[0].real + sld_profile[0].imag * 1.j
            kn1_init_4pi = 4. * np.pi + 0.j
            for_exp = -2. + 0.j
            rough = 0. + 0.j
            sld_diff = np.subtract(kn1_init_layer, kn1_init_super)
            four_pi_sld_diff = np.multiply(kn1_init_4pi, sld_diff)
            kz_squared = np.multiply(kz, kz)
            kn1_squared = np.subtract(kz_squared, four_pi_sld_diff)
            kn1 = np.sqrt(kn1_squared)
            rough_squared = np.multiply(rough, rough)
            kn1_by_rough = np.multiply(kn1, rough_squared)
            kn_by_kn1_by_rough = np.multiply(kn, kn1_by_rough)
            rough_for_exp = np.multiply(for_exp, kn_by_kn1_by_rough)
            whole_rough = np.exp(rough_for_exp)
            numerator = np.subtract(kn, kn1)
            denominator = np.add(kn, kn1)
            fraction = np.divide(numerator, denominator)
            f = np.multiply(fraction, whole_rough)
            bin_w = lt + 0j
            if it > 0:
                betan = np.multiply(kn, bin_w)
            else:
                betan = 0 + 0j
            im = 0 + 1j
            negim = 0 - 1j
            M[0][0] = np.exp(np.multiply(im, betan))
            M[0][1] = np.multiply(f, M[0][0])
            M[1][1] = np.exp(np.multiply(negim, betan))
            M[1][0] = np.multiply(f, M[1][1])
            bm = np.dot(bn, M)
            bn = bm
            kn = kn1
        r = np.abs(bn[1][0]) / np.abs(bn[0][0])
        refl.append(r)
    return refl