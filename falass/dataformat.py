class datastruct(object):
    def __init__(self, q, i, di, dq):
        self.q = q
        self.i = i
        self.di = i
        self.dq = i

class scatlens(object):
    def __init__(self, atom, real, imag):
        self.atom = atom
        self.real = real * 1e-5
        self.imag = imag * 1e-5

class atompositions(object):
    def __init__(self, atom, position):
        self.atom = atom
        self.zpos = position
