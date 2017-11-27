from falass import job, readwrite, dataformat
import numpy as np
from matplotlib import rc

import matplotlib.pyplot as plt

class SLD(object):
    def __init__(self, job):
        self.job = job
        self.sld_profile = []
        self.av_sld_profile = []
        self.av_sld_profile_err = []

    def getSLDProfile(self):
        prog = 0
        self.sld_profile = []
        print("[ 0 % ]")
        for i in range(0, len(self.job.files.atoms)):
            build_sld = []
            if check_array(self.job.times, self.job.files.times[i]) is True:
                z_cut = self.job.files.cell[i][2] - self.job.cut_off_size
                number_of_bins = int(z_cut / self.job.layer_thickness)
                for j in range(0, number_of_bins):
                    build_sld.append(dataformat.sldpro(self.job.layer_thickness, 0, 0))
                self.sld_profile.append(build_sld)
        k = 0
        for i in range(0, len(self.job.files.atoms)):
            if check_array(self.job.times, self.job.files.times[i]) is True:
                z_cut = self.job.files.cell[i][2] - self.job.cut_off_size
                number_of_bins = int(z_cut / self.job.layer_thickness)
                for j in range(0, len(self.job.files.atoms[i])):
                    cutBox = int(self.job.files.cell[i][2] * self.job.layer_thickness)
                    if self.job.files.atoms[i][j].zpos < number_of_bins * self.job.layer_thickness:
                        bin = int(self.job.files.atoms[i][j].zpos / self.job.layer_thickness)
                        scatlen_to_add = get_scatlen(self.job.files.atoms[i][j].atom, self.job.files.scat_lens)
                        self.sld_profile[i][bin].real += scatlen_to_add[0]
                        self.sld_profile[i][bin].imag += scatlen_to_add[1]
                    k += 1
                    prog_new = np.floor((k) / (len(self.job.files.atoms[i]) * len(self.job.files.atoms)) * 100)
                    if prog_new > prog + 9:
                        prog = prog_new
                        print("[{} {} % ]".format('#' * int(prog / 10), int(prog / 10) * 10))
        for i in range(0, len(self.job.files.atoms)):
            if check_array(self.job.times, self.job.files.times[i]) is True:
                z_cut = self.job.files.cell[i][2] - self.job.cut_off_size
                number_of_bins = int(z_cut / self.job.layer_thickness)
                for j in range(0, number_of_bins):
                    self.sld_profile[i][j].real = self.sld_profile[i][j].real / (self.job.files.cell[i][0] * self.job.files.cell[i][1] * self.job.layer_thickness)
                    self.sld_profile[i][j].imag = self.sld_profile[i][j].imag / (self.job.files.cell[i][0] * self.job.files.cell[i][1] * self.job.layer_thickness)


    def averageSLDProfile(self):
        prog = 0
        self.av_sld_profile_err = []
        self.av_sld_profile = []
        print("[ 0 % ]")
        z_cut = self.job.files.cell[0][2] - self.job.cut_off_size
        number_of_bins = int(z_cut / self.job.layer_thickness)
        k = 0
        for j in range(0, number_of_bins):
            self.av_sld_profile.append(dataformat.sldpro(self.job.layer_thickness, 0, 0))
            for i in range(0, len(self.job.files.atoms)):
                self.av_sld_profile[j].real += self.sld_profile[i][j].real
                self.av_sld_profile[j].imag += self.sld_profile[i][j].imag
            self.av_sld_profile[j].real /= len(self.job.times)
            self.av_sld_profile[j].imag /= len(self.job.times)
            self.av_sld_profile_err.append(dataformat.sldpro(self.job.layer_thickness, 0, 0))
            for i in range(0, len(self.job.files.atoms)):
                self.av_sld_profile_err[j].real += np.square(self.sld_profile[i][j].real - self.av_sld_profile[j].real)
                self.av_sld_profile_err[j].imag += np.square(self.sld_profile[i][j].imag - self.av_sld_profile[j].imag)
                k += 1
                prog_new = np.floor((k) / (number_of_bins * len(self.job.files.atoms)) * 100)
                if prog_new > prog + 9:
                    prog = prog_new
                    print("[{} {} % ]".format('#' * int(prog / 10), int(prog / 10) * 10))
            self.av_sld_profile_err[j].real = np.sqrt(1. / (len(self.job.times) - 1)) * self.av_sld_profile_err[j].real
            self.av_sld_profile_err[j].imag = np.sqrt(1. / (len(self.job.times) - 1)) * self.av_sld_profile_err[j].imag

    def plotaverageSLDProfile(self, real=False, imag=False):
        x = []
        y = []
        dy = []
        buildx = 0
        plt.rc('text', usetex=True)
        plt.rc('font', family='serif')
        for i in range(0, len(self.av_sld_profile)):
            if real:
                y.append(self.av_sld_profile[i].real)
                dy.append(self.av_sld_profile_err[i].real)
            if imag:
                y.append(self.av_sld_profile[i].imag)
                dy.append(self.av_sld_profile_err[i].imag)
            buildx += self.av_sld_profile[i].thick
            x.append(buildx)
        plt.bar(np.asarray(x) - self.job.layer_thickness/2., np.asarray(y)*1e6, width = self.job.layer_thickness,
                yerr=np.asarray(dy)*1e6, color = 'w', edgecolor = 'k')
        plt.ylabel('SLD (1E10$^{-6}$ \AA$^{-2}$)')
        plt.xlabel('$z$ (\AA)')
        plt.xlim([0, np.amax(x)])
        plt.ylim([np.amin(np.asarray(y)*1e6), np.amax(np.asarray(y)*1e6)+1])
        plt.show()


def check_array(array, check):
    for i in range(0, len(array)):
        if check == array[i]:
            return True
    return False

def get_scatlen(atom, scat_lens):
    for i in range(0, len(scat_lens)):
        if atom == scat_lens[i].atom:
            return scat_lens[i].real, scat_lens[i].imag