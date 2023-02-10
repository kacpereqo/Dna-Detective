
import peptides
import math
import numpy as np

from .constants import *


class ProteinCharge():
    def __init__(self, sequence):
        self.sequence = sequence

# |------------------------------------------------------------------------------|#

    def isoelectric_point(self, scale="Rodwell", precision=0.001) -> float:

        max_ = 14.0
        min_ = 0.0
        pH = 7.0
        while not math.isclose(max_, min_):
            pH = (max_ + min_) / 2
            c = self.charge_at_ph(pH=pH, scale=scale)
            if c <= 0:
                max_ = pH
            if c >= 0:
                min_ = pH

        return pH

    # |------------------------------------------------------------------------------|#

    def charge_at_ph(self, scale="Rodwell", pH=7.0) -> float:
        pka_scale = PKA_SCALE[scale]

        net_charge = 0.0
        for aa in self.sequence:
            if aa in AMINE_ACIDS_CHARGE:
                charge = AMINE_ACIDS_CHARGE[aa]
                net_charge += charge / (1.0 + 10 ** (charge * (pH - pka_scale[aa])))

        net_charge += 1.0 / (1.0 + 10 ** (pH - pka_scale["NH2"]))
        net_charge += -1.0 / (1.0 + 10 ** (pka_scale["CO"] - pH))

        return net_charge

    # |------------------------------------------------------------------------------|#

    def net_charge(self, scale="Rodwell", start=0.0, end=14.0, step=0.5) -> float:
        net_charge = {}
        for pH in np.arange(start, end, step):
            pH = round(pH, 2)
            net_charge[pH] = format(self.charge_at_ph(scale, pH), '.3f')

        return net_charge

    # |------------------------------------------------------------------------------|#
