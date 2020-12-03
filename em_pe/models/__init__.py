from .interpolated_model import *
from .kilonova_3c import *
from .kn_interp import *

from .parameters import *

model_dict = {
        "interpolated":interpolated,
        "kilonova_3c":kilonova_3c,
        "kn_interp":kn_interp
}

param_dict = {
        "dist":Distance,
        "mej_red":EjectaMassRed,
        "mej_purple":EjectaMassPurple,
        "mej_blue":EjectaMassBlue,
        "mej_dyn":DynamicalEjectaMass,
        "mej_wind":WindEjectaMass,
        "vej_red":EjectaVelocityRed,
        "vej_purple":EjectaVelocityPurple,
        "vej_blue":EjectaVelocityBlue,
        "vej_dyn":DynamicalEjectaVelocity,
        "vej_wind":WindEjectaVelocity,
        "Tc_red":TcRed,
        "Tc_purple":TcPurple,
        "Tc_blue":TcBlue,
        "sigma":Sigma
}
