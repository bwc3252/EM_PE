import numpy as np
from scipy.stats import loguniform

class Parameter:
    def __init__(self, name, llim, rlim):
        self.name = name
        self.llim = llim
        self.rlim = rlim

    def prior(self, x):
        """
        Default prior is uniform
        """
        return 1.0 / (self.rlim - self.llim)


class Distance(Parameter):
    def __init__(self):
        Parameter.__init__(self, "dist", 10.0, 100.0)


class EjectaMassRed(Parameter):
    def __init__(self):
        Parameter.__init__(self, "mej_red", 0.001, 0.1)


class EjectaMassPurple(Parameter):
    def __init__(self):
        Parameter.__init__(self, "mej_purple", 0.001, 0.1)


class EjectaMassBlue(Parameter):
    def __init__(self):
        Parameter.__init__(self, "mej_blue", 0.001, 0.1)


class EjectaMassRed(Parameter):
    def __init__(self):
        Parameter.__init__(self, "mej_red", 0.001, 0.1)


class EjectaVelocityRed(Parameter):
    def __init__(self):
        Parameter.__init__(self, "vej_red", 0.1, 0.4)


class EjectaVelocityPurple(Parameter):
    def __init__(self):
        Parameter.__init__(self, "vej_purple", 0.1, 0.4)


class EjectaVelocityBlue(Parameter):
    def __init__(self):
        Parameter.__init__(self, "vej_blue", 0.1, 0.4)


class TcRed(Parameter):
    def __init__(self):
        Parameter.__init__(self, "Tc_red", 3500.0, 4000.0)

    #def prior(self, x):
    #    """
    #    Use a log-uniform prior
    #    """
    #    return loguniform.pdf(x, self.llim, self.rlim) # FIXME should I be setting the loc and scale for this?


class TcPurple(Parameter):
    def __init__(self):
        Parameter.__init__(self, "Tc_purple", 1000.0, 1500.0)

    #def prior(self, x):
    #    """
    #    Use a log-uniform prior
    #    """
    #    return loguniform.pdf(x, self.llim, self.rlim) # FIXME should I be setting the loc and scale for this?


class TcBlue(Parameter):
    def __init__(self):
        Parameter.__init__(self, "Tc_blue", 400.0, 1000.0)

    #def prior(self, x):
    #    """
    #    Use a log-uniform prior
    #    """
    #    return loguniform.pdf(x, self.llim, self.rlim) # FIXME should I be setting the loc and scale for this?

class Sigma(Parameter):
    def __init__(self):
        Parameter.__init__(self, "sigma", 0.0, 0.5)
