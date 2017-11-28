class datastruct(object):
    def __init__(self, q, i, di, dq):
        self.q = q
        self.i = i
        self.di = di
        self.dq = dq

class scatlens(object):
    def __init__(self, atom, real, imag):
        self.atom = atom
        self.real = real * 1e-5
        self.imag = imag * 1e-5

class atompositions(object):
    def __init__(self, atom, position):
        self.atom = atom
        self.zpos = position

class sldpro(object):
    def __init__(self, thick, real, imag):
        self.thick = thick
        self.real = real
        self.imag = imag
