from .interpolated_model import *
from .kilonova_3c import *

from .parameters import *

model_dict = {
        "interpolated":interpolated,
        "kilonova_3c":kilonova_3c
}

param_dict = {
        "dist":Distance,
        "mej_red":EjectaMassRed,
        "mej_purple":EjectaMassPurple,
        "mej_blue":EjectaMassBlue,
        "vej_red":EjectaVelocityRed,
        "vej_purple":EjectaVelocityPurple,
        "vej_blue":EjectaVelocityBlue,
        "Tc_red":TcRed,
        "Tc_purple":TcPurple,
        "Tc_blue":TcBlue,
        "sigma":Sigma
}
