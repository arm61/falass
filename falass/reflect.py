import numpy as np
from falass import dataformat
import matplotlib.pyplot as plt
from scipy.interpolate import InterpolatedUnivariateSpline


class Reflect:
    """Reflectometry calculation.

    A class for the calculation of reflectometry from the sld profile defined in the falass.sld.SLD class.

    Parameters
    ----------
    sld_profile: array_like falass.dataformat.SLDPro
        An array describing the scattering length density of the simulation cell under study.
    exp_data: array_like falass.dataformat.QData
        An array giving the experimental data from the datfile.
    """
    def __init__(self, sld_profile, exp_data):
        self.sld_profile = sld_profile
        self.exp_data = exp_data
        self.averagereflect = []
        self.reflect = []

    def calc_ref(self):
        """Calculate reflectometry.

        The calculation of the reflectometry profiles based on the sld profiles calculated from each of the timesteps
        under study.
        """
        if len(self.exp_data) > 0:
            self.reflect = []
            prog = 0
            k = 0
            print("[ 0 % ]")
            for i in range(0, len(self.sld_profile)):
                k += 1
                prog_new = np.floor(k / (len(self.sld_profile)) * 100)
                if prog_new > prog + 9:
                    prog = prog_new
                    print("[{} {} % ]".format('#' * int(prog / 10), int(prog / 10) * 10))
                refl = convolution(self.exp_data, self.sld_profile[i])
                a = []
                for j in range(0, len(self.exp_data)):
                    a.append(dataformat.QData(self.exp_data[j].q, refl[j], 0, self.exp_data[j].dq))
                self.reflect.append(a)
        else:
            raise ValueError('No q vectors have been defined -- either read a .dat file or get q vectors.')

    def average_ref(self):
        """Average reflectometry profiles.

        The averaging of the reflectometry profiles as calculated by the calc_ref() function.
        """
        if len(self.exp_data) > 0:
            self.averagereflect = []
            for i in range(0, len(self.reflect[0])):
                self.averagereflect.append(dataformat.QData(self.exp_data[i].q, 0, 0, self.exp_data[i].dq))
            for i in range(0, len(self.reflect[0])):
                for j in range(0, len(self.reflect)):
                    self.averagereflect[i].i += self.reflect[j][i].i
                self.averagereflect[i].i /= len(self.reflect)
                for j in range(0, len(self.reflect)):
                    self.averagereflect[i].di += np.square(self.reflect[j][i].i - self.averagereflect[i].i)
                self.averagereflect[i].di = np.sqrt(1. / (len(self.reflect) - 1) * self.averagereflect[i].di)
        else:
            raise ValueError('No q vectors have been defined -- either read a .dat file or get q vectors.')

    def plot_ref(self, rq4=True):
        """Plot reflectometry profile.

        The plotting of the calculated reflectometry data without comparison to experimental data.

        Parameters
        ----------
        rq4: bool
            Should the data be transformed to rq4 space.
        """
        if len(self.exp_data) > 0:
            x = []
            y = []
            dy = []
            plt.rc('text', usetex=True)
            plt.rc('font', family='serif')
            if rq4:
                for i in range(0, len(self.exp_data)):
                    x.append(self.averagereflect[i].q)
                    y.append(np.log10(self.averagereflect[i].i * self.averagereflect[i].q ** 4))
                    dy.append((self.averagereflect[i].di * self.averagereflect[i].q ** 4) /
                              (self.averagereflect[i].i * np.log(10)))
            else:
                for i in range(0, len(self.exp_data)):
                    x.append(self.averagereflect[i].q)
                    y.append(np.log10(self.averagereflect[i].i))
                    dy.append(self.averagereflect[i].di / (self.averagereflect[i].i * np.log(10)))
            x = np.asarray(x)
            y = np.asarray(y)
            dy = np.asarray(dy)
            plt.errorbar(x, y, yerr=dy)
            plt.xlabel('$q$ (\AA)')
            plt.ylabel('log($Rq^4$) (\AA$^4$)')
            plt.show()
        else:
            raise ValueError('No q vectors have been defined -- either read a .dat file or get q vectors.')


def convolution(exp_data, sld_profile):
    """Convolution/smearing

    The convolution of the reflectometry data by a gaussian of constant width (a percentage of the q-vector)

    Parameters
    ----------
    exp_data: falass.dataformat.QData
        The experimental data from the datfile.
    sld_profile: falass.dataformat.SLDPro
        The SLD profile calculated from the simulation trajectory.
    lt: float
        The thickness of the layers.

    Returns
    -------
    array_like
        The smeared reflectometry profile.
    """
    fwhm = 2 * np.sqrt(2 * np.log(2))

    res = exp_data[0].dq / exp_data[0].q

    if exp_data[0].dq / exp_data[0].q < 0.0005:
        return reflectivity(exp_data, sld_profile)

    gnum = 51
    ggpoint = (gnum - 1) / 2

    def gauss(x, s):
        return 1. / s / np.sqrt(2 * np.pi) * np.exp(-0.5 * x ** 2 / s / s)

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

    rvals = reflectivity(xlin, sld_profile)
    smeared_rvals = np.convolve(rvals, gaussy, mode='same')
    interpol = InterpolatedUnivariateSpline(xlin, smeared_rvals)

    smeared_output = interpol(q)
    smeared_output *= gaussx[1] - gaussx[0]
    return smeared_output


def reflectivity(exp_data, sld_profile):
    """Abeles optical matrix formalism.

    The calculation of the reflectometry using the Abeles optical matrix method.

    Parameters
    ----------
    exp_data: falass.dataformat.QData
        The experimental data from the datfile.
    sld_profile: falass.dataformat.SLDPro
        The SLD profile calculated from the simulation trajectory.
    lt: float
        The thickness of the layers.

    Returns
    -------
    array_like
        The reflectometry profile.
    """
    layers = np.zeros((len(sld_profile), 4))
    for i in range(0, len(sld_profile)):
        layers[i][0] = sld_profile[i].thick
        layers[i][1] = sld_profile[i].real
        layers[i][2] = sld_profile[i].imag
        layers[i][3] = 0
    qvals = np.asfarray(exp_data).ravel()
    nlayers = len(sld_profile) - 2
    npnts = qvals.size

    kn = np.zeros((npnts, nlayers + 2), np.complex128)

    sld = np.zeros(nlayers + 2, np.complex128)
    sld[:] += ((layers[:, 1].real - layers[0, 1]) + 1j * (layers[:, 2] - layers[0, 2]))

    kn[:] = np.sqrt(qvals[:, np.newaxis] ** 2. / 4. - 4 * np.pi * sld)

    mrtot00 = 1
    mrtot11 = 1
    mrtot10 = 0
    mrtot01 = 0
    k = kn[:, 0]

    for idx in range(1, nlayers + 2):
        k_next = kn[:, idx]
        rj = (k - k_next) / (k + k_next)
        rj *= np.exp(k * k_next)

        # work out characteristic matrix of layer
        mi00 = np.exp(k * 1j * np.fabs(layers[idx - 1, 0])) if idx - 1 else 1
        mi11 = np.exp(k * -1j * np.fabs(layers[idx - 1, 0])) if idx - 1 else 1

        mi10 = rj * mi00
        mi01 = rj * mi11

        # matrix multiply mrtot by characteristic matrix
        p0 = mrtot00 * mi00 + mrtot10 * mi01
        p1 = mrtot00 * mi10 + mrtot10 * mi11
        mrtot00 = p0
        mrtot10 = p1

        p0 = mrtot01 * mi00 + mrtot11 * mi01
        p1 = mrtot01 * mi10 + mrtot11 * mi11

        mrtot01 = p0
        mrtot11 = p1

        k = k_next

    ref = (mrtot01 * np.conj(mrtot01)) / (mrtot00 * np.conj(mrtot00))
    return np.real(np.reshape(ref, exp_data.shape))
