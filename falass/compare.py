import numpy as np
import matplotlib.pyplot as plt
from scipy.optimize import curve_fit
from falass import dataformat


class Compare():
    """Comparison and fitting.

    For the comparison and fitting of calculated and experimental reflectometry.

    Parameters
    ----------
    exp_data: array_like falass.dataformat.QData
        The experimental reflectometry data read from the datfile.
    sim_data: array_like falass.dataformat.QData
        The calculated reflectometry data from the simulation.
    scale: float
        The amount by which the calculated reflectometry should be scaled.
    background: float
        The height of the uniform background to be added to the calculated reflectometry.
    """
    def __init__(self, exp_data, sim_data, scale, background):
        self.exp_data = exp_data
        self.sim_data = sim_data
        self.scale = scale
        self.background = background
        self.sim_data_fitted = []

    def change_scale(self, scale):
        """Edit scale.

        Lets the scale factor be changed.

        Parameters
        ----------
        scale: float
            The amount by which the calculated reflectometry should be scaled.
        """
        self.scale = scale

    def change_background(self, background):
        """Edit background.

        Lets the background factor be changed.

        Parameters
        ----------
        background: float
            The height of the uniform background to be added to the calculated reflectometry.
        """
        self.background = background

    def fit(self):
        """Fit scale and background.

        Perform the fitting of the scale and background for the calculated data to the experimental data.
        Currently only a logarithmically transformed fitted can be conducted.
        """
        if len(self.exp_data) > 0:
            if self.exp_data[0].i is not None:
                y = []
                dy = []
                y2 = []
                for i in range(0, len(self.exp_data)):
                    y.append(np.log(self.exp_data[i].i))
                    dy.append(self.exp_data[i].di / (self.exp_data[i].i * np.log(10)))
                    y2.append(self.sim_data[i].i)
                popt, pcov = curve_fit(scale_and_background, y2, y, bounds=((1e-6, 0), (10, 1e-3)), sigma=dy)
                self.scale = popt[0]
                self.background = popt[1]
            else:
                raise ValueError('No experimental data has been set for comparison, please read in a a .dat file.')
        else:
            raise ValueError('No q vectors have been defined -- either read a .dat file or get q vectors.')

    def plot_compare(self, rq4=True, fitted=True):
        """Plot a comparison.

        Plotting the comparision between the calculated and experimental reflectometry.

        Parameters
        ----------
        rq4: bool
            Should the data be plotted in rq4 space.
        fitted: bool
            Should the fitted reflectometry data be used.
        """
        x = []
        y = []
        dy = []
        x2 = []
        y2 = []
        dy2 = []
        plt.rc('text', usetex=True)
        plt.rc('font', family='serif')
        if fitted:
            if len(self.sim_data_fitted) > 0:
                k = self.sim_data_fitted
            else:
                raise ValueError("The reflectometry data has not been returned yet -- please run the fit() function "
                                 "and the return_fitted().")
        else:
            k = self.sim_data
        if rq4:
            for i in range(0, len(self.exp_data)):
                x.append(k[i].q)
                y.append(k[i].i * np.power(self.exp_data[i].q, 4))
                da = k[i].di * np.power(self.exp_data[i].q, 4)
                dy.append(da)
                x2.append(self.exp_data[i].q)
                y2.append(self.exp_data[i].i * np.power(self.exp_data[i].q, 4))
                da = self.exp_data[i].di * np.power(self.exp_data[i].q, 4)
                dy2.append(da)
                plt.ylabel('$Rq^4$')
        else:
            for i in range(0, len(self.exp_data)):
                x.append(k[i].q)
                y.append(k[i].i)
                da = k[i].di
                dy.append(da)
                x2.append(self.exp_data[i].q)
                y2.append(np.log10(self.exp_data[i].i))
                da = self.exp_data[i].di
                dy2.append(da)
                plt.ylabel('log$_{10}$($R$)')
        x = np.asarray(x)
        y = np.asarray(y)
        dy = np.asarray(dy)
        x2 = np.asarray(x2)
        y2 = np.asarray(y2)
        dy2 = np.asarray(dy2)
        plt.errorbar(x, y, yerr=dy)
        plt.errorbar(x2, y2, yerr=dy2, linestyle='', marker='o')
        plt.xlabel('$q$ (\AA)')
        plt.yscale('log')
        plt.show()

    def return_fitted(self):
        """Return fitted.

        Return the fitted calculated reflectometry data for use.
        """
        self.sim_data_fitted = []
        for i in range(0, len(self.sim_data)):
            a = self.sim_data[i].i * self.scale + self.background
            b = self.sim_data[i].di * self.scale
            self.sim_data_fitted.append(dataformat.QData(self.sim_data[i].q, a, b, self.sim_data[i].dq))


def scale_and_background(sim_data, scale, background):
    """Apply scale and background.

    Apply a scale factor and uniform background to the calculated reflectometry data.

    Parameters
    ----------
    sim_data: array_like float
        The data to be scaled and have a background added.
    scale: float
        The amount by which the data should be scaled.
    background: float
        The size of the uniform background to be added.

    Returns
    -------
    array_like float
        The scaled and background added reflectometry in log space.
    """
    sim_data = np.log(sim_data * scale + background)
    return sim_data
