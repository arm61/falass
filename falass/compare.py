import numpy as np
import matplotlib.pyplot as plt
from scipy.optimize import curve_fit
from falass import dataformat

class Compare(object):
    def __init__(self, exp_data, sim_data, scale, background):
        self.exp_data = exp_data
        self.sim_data = sim_data
        self.scale = scale
        self.background = background
        self.sim_data_fitted = []

    def fit(self):
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

    def plotCompare(self, rq4=True, fitted = True):
        x = []
        y = []
        dy = []
        x2 = []
        y2 = []
        dy2 = []
        plt.rc('text', usetex=True)
        plt.rc('font', family='serif')
        if fitted:
            k = self.sim_data_fitted
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
                y.append(np.log10(k[i].i))
                dy.append(k[i].di)
                x2.append(self.exp_data[i].q)
                y2.append(np.log10(self.exp_data[i].i))
                dy2.append(self.exp_data[i].di)
                plt.ylabel('log$_{10}$($R$)')
        x = np.asarray(x)
        y = np.asarray(y)
        dy = np.asarray(dy)
        x2 = np.asarray(x2)
        y2 = np.asarray(y2)
        dy2 = np.asarray(dy2)
        plt.errorbar(x, y, yerr=dy)
        plt.errorbar(x2, y2, yerr=dy2, linestyle = '', marker='o')
        plt.xlabel('$q$ (\AA)')
        plt.yscale('log')
        plt.show()

    def returnFitted(self):
        self.sim_data_fitted = []
        for i in range(0, len(self.sim_data)):
            a = self.sim_data[i].i * self.scale + self.background
            b = self.sim_data[i].di * self.scale
            self.sim_data_fitted.append(dataformat.datastruct(self.sim_data[i].q, a, b, self.sim_data[i].dq))

def scale_and_background(sim_data, scale, background):
    sim_data = np.log(sim_data * scale + background)
    return sim_data


