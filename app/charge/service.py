

import numpy as np

from .constants import *

class ProteinCharge():
    def __init__(self, sequence):
        self.sequence = sequence

#|------------------------------------------------------------------------------|#

    def isoelectric_point(self, scale="IPC_protein", precision = 0.01) -> float:
    
        pka_scale = PKA_SCALE[scale]

        pH = 6.51       
        pH_prev = 0.0
        pH_next = 14.0
        temp = 0.01


        while True:

            net_charge = 0.0
            net_charge += -1.0/(1.0+pow(10,(pka_scale['C']-pH)))
            net_charge +=  1.0/(1.0+pow(10,(pH-pka_scale['NH2'])))

            for aa, charge in AMINE_ACIDS_CHARGE.items():
                if charge == -1:
                    net_charge += -self.sequence.count(aa)/(1.0+pow(10,(pka_scale[aa]-pH)))

                elif charge == 1:
                    net_charge += self.sequence.count(aa)/(1.0 + pow(10,(pH-pka_scale[aa])))


            temp = pH

            if net_charge < 0.0:                   
                pH = pH-((pH-pH_prev)/2.0)
                pH_next = temp

            elif net_charge > 0.0:
                pH = pH + ((pH_next-pH)/2.0)
                pH_prev = temp
                        

            if (pH - pH_prev < precision) and (pH_next - pH < precision):
                return pH

    #|------------------------------------------------------------------------------|#

    def charge_at_ph(self, scale = "Rodwell", pH = 7.0) -> float:
                pka_scale = PKA_SCALE[scale]

                net_charge = 0.0
                for aa in self.sequence:
                    if aa in AMINE_ACIDS_CHARGE:
                        charge = AMINE_ACIDS_CHARGE[aa]
                        net_charge += charge / (1.0 + 10 ** (charge * (pH - pka_scale[aa])))

                net_charge += 1.0 / (1.0 + 10 ** (pH - pka_scale["NH2"]))
                net_charge += -1.0 / (1.0 + 10 ** (pka_scale["CO"] - pH))

                return net_charge

    #|------------------------------------------------------------------------------|#
    
    def net_charge(self, scale = "Rodwell", start = 0.0, end = 14.0, step = 0.5) -> float:
        net_charge = {} 
        for pH in np.arange(start, end, step):
            net_charge[pH] = self.charge_at_ph(scale, pH)
        
        return net_charge

    #|------------------------------------------------------------------------------|#
    
        