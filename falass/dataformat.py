class QData():
    """Reflectometry data.

    A class to hold the qdata information consisting of four floats associated with the q-vector, the intensity of the
    reflectometry, the uncertainty in the intensity and the uncertainty in the q-vector/resolution of the q-vector.
    """
    def __init__(self, q, i, di, dq):
        self.q = q
        self.i = i
        self.di = di
        self.dq = dq


class ScatLens():
    """Scattering lengths.

    A class to hold the scattering lengths of the different atom types consisting of a str atom type name, and two
    floats associated with the real and imaginary scattering lengths of that atom type.
    """
    def __init__(self, atom, real, imag):
        self.atom = atom
        self.real = real * 1e-5
        self.imag = imag * 1e-5


class AtomPositions():
    """z-Dimension positions.

    A class to hold the atom positions in the z-dimension, consisting of the str atom type name and a float giving the
    position in the z-dimension.
    """
    def __init__(self, atom, position):
        self.atom = atom
        self.zpos = position


class SLDPro():
    """Layer information for SLD profile.

    A class to hold the layer description of the sld profile consisting of three floats associated with the thickness,
    real scattering length density and imaginary scattering length density.
    """
    def __init__(self, thick, real, imag):
        self.thick = thick
        self.real = real
        self.imag = imag
