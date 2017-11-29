from falass import dataformat, job
import numpy as np
import matplotlib.pyplot as plt


class SLD:
    """SLD profile calculation.

    This class enables the calculation of the SLD profile for each of the timesteps as defined in the
    falass.job.Job. Further it will then allow the calculation and plotting of the average SLD profile.

    Parameters
    ----------
    assigned_job: falass.job.Job
        The is the Job class for the particular falass run taking place. See the job.Job class for more information.
    """
    def __init__(self, assigned_job):
        self.assigned_job = assigned_job
        self.sld_profile = []
        self.av_sld_profile = []
        self.av_sld_profile_err = []

    def get_sld_profile(self):
        """Calculate SLD profile.

        This will calculate the SLD profile for each of the timesteps defined in the falass.job.Job. This is achieved
        by summing the scattering lengths for each of the atoms found in a given layer (of defined thickness). This
        total scattering length is converted to a density by division by the volume of the layer.
        """
        prog = 0
        self.sld_profile = []
        print("[ 0 % ]")
        for i in range(0, len(self.assigned_job.files.atoms)):
            build_sld = []
            if job.check_array(self.assigned_job.times, self.assigned_job.files.times[i]):
                z_cut = self.assigned_job.files.cell[i][2] - self.assigned_job.cut_off_size
                number_of_bins = int(z_cut / self.assigned_job.layer_thickness)
                for j in range(0, number_of_bins):
                    build_sld.append(dataformat.SLDPro(self.assigned_job.layer_thickness, 0, 0))
                self.sld_profile.append(build_sld)
        k = 0
        for i in range(0, len(self.assigned_job.files.atoms)):
            if job.check_array(self.assigned_job.times, self.assigned_job.files.times[i]):
                z_cut = self.assigned_job.files.cell[i][2] - self.assigned_job.cut_off_size
                number_of_bins = int(z_cut / self.assigned_job.layer_thickness)
                for j in range(0, len(self.assigned_job.files.atoms[i])):
                    if self.assigned_job.files.atoms[i][j].zpos < number_of_bins * self.assigned_job.layer_thickness:
                        bin_choose = int(self.assigned_job.files.atoms[i][j].zpos / self.assigned_job.layer_thickness)
                        scatlen_to_add = get_scatlen(self.assigned_job.files.atoms[i][j].atom,
                                                     self.assigned_job.files.scat_lens)
                        self.sld_profile[i][bin_choose].real += scatlen_to_add[0]
                        self.sld_profile[i][bin_choose].imag += scatlen_to_add[1]
                    k += 1
                    prog_new = np.floor(k / (len(self.assigned_job.files.atoms[i]) * 
                                             len(self.assigned_job.files.atoms)) * 100)
                    if prog_new > prog + 9:
                        prog = prog_new
                        print("[{} {} % ]".format('#' * int(prog / 10), int(prog / 10) * 10))
        for i in range(0, len(self.assigned_job.files.atoms)):
            if job.check_array(self.assigned_job.times, self.assigned_job.files.times[i]) is True:
                z_cut = self.assigned_job.files.cell[i][2] - self.assigned_job.cut_off_size
                number_of_bins = int(z_cut / self.assigned_job.layer_thickness)
                for j in range(0, number_of_bins):
                    self.sld_profile[i][j].real = self.sld_profile[i][j].real / (self.assigned_job.files.cell[i][0] *
                                                                                 self.assigned_job.files.cell[i][1] *
                                                                                 self.assigned_job.layer_thickness)
                    self.sld_profile[i][j].imag = self.sld_profile[i][j].imag / (self.assigned_job.files.cell[i][0] *
                                                                                 self.assigned_job.files.cell[i][1] *
                                                                                 self.assigned_job.layer_thickness)

    def average_sld_profile(self):
        """Average SLD profiles.

        Allows for the calculation of the average SLD profile across all of the timesteps that were studied.
        """
        prog = 0
        self.av_sld_profile_err = []
        self.av_sld_profile = []
        print("[ 0 % ]")
        z_cut = self.assigned_job.files.cell[0][2] - self.assigned_job.cut_off_size
        number_of_bins = int(z_cut / self.assigned_job.layer_thickness)
        k = 0
        for j in range(0, number_of_bins):
            self.av_sld_profile.append(dataformat.SLDPro(self.assigned_job.layer_thickness, 0, 0))
            for i in range(0, len(self.assigned_job.times)):
                self.av_sld_profile[j].real += self.sld_profile[i][j].real
                self.av_sld_profile[j].imag += self.sld_profile[i][j].imag
            self.av_sld_profile[j].real /= len(self.assigned_job.times)
            self.av_sld_profile[j].imag /= len(self.assigned_job.times)
            self.av_sld_profile_err.append(dataformat.SLDPro(self.assigned_job.layer_thickness, 0, 0))
            for i in range(0, len(self.assigned_job.times)):
                self.av_sld_profile_err[j].real += np.square(self.sld_profile[i][j].real - self.av_sld_profile[j].real)
                self.av_sld_profile_err[j].imag += np.square(self.sld_profile[i][j].imag - self.av_sld_profile[j].imag)
                k += 1
                prog_new = np.floor(k / (number_of_bins * len(self.assigned_job.files.atoms)) * 100)
                if prog_new > prog + 9:
                    prog = prog_new
                    print("[{} {} % ]".format('#' * int(prog / 10), int(prog / 10) * 10))
            self.av_sld_profile_err[j].real = np.sqrt(1. / (len(self.assigned_job.times) - 1)) * self.av_sld_profile_err[j].real
            self.av_sld_profile_err[j].imag = np.sqrt(1. / (len(self.assigned_job.times) - 1)) * self.av_sld_profile_err[j].imag

    def plot_sld_profile(self, real=True, imag=False):
        """Plot SLD.

        Plots the average sld profile using matplotlib.

        Parameters
        ----------
        real: bool
            Should the real SLD profile be plotted (if both real and imaginary are true the real will be plotted).
        imag: bool
            Should the imaginary SLD profile be plotted (if both real and imaginary are true the real will be plotted).
        """
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
        plt.bar(np.asarray(x) - self.assigned_job.layer_thickness/2., np.asarray(y)*1e6, 
                width=self.assigned_job.layer_thickness, yerr=np.asarray(dy)*1e6, color='w', edgecolor='k')
        plt.ylabel('SLD (10$^{-6}$ \AA$^{-2}$)')
        plt.xlabel('$z$ (\AA)')
        plt.xlim([0, np.amax(x)])
        plt.ylim([np.amin(np.asarray(y)*1e6), np.amax(np.asarray(y)*1e6)+1])
        plt.show()


def get_scatlen(atom, scat_lens):
    """Find scattering length.

    This gets the scattering length for a given atom type

    Parameters
    ----------
    atom: str
        The name of the atom type that the scattering length is needed for.
    scat_lens: array_like falass.dataformat.ScatLens
        The array of the scattering lengths that is defined in the falass.readwrite.Files class.

    Returns
    -------
    tuple_like
        The real and imaginary scattering lengths for the given atom type.
    """
    for i in range(0, len(scat_lens)):
        if atom == scat_lens[i].atom:
            return scat_lens[i].real, scat_lens[i].imag
    raise ValueError("Attempt to get the scattering length of the atom type {} failed. This should never happen. "
                     "Please contact the developers".format(atom))
