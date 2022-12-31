"""Calculate isoelectric points of polypeptides using methods of Bjellqvist.
pK values and the methods are taken from::
    * Bjellqvist, B.,Hughes, G.J., Pasquali, Ch., Paquet, N., Ravier, F.,
    Sanchez, J.-Ch., Frutiger, S. & Hochstrasser, D.F.
    The focusing positions of polypeptides in immobilized pH gradients can be
    predicted from their amino acid sequences. Electrophoresis 1993, 14,
    1023-1031.
    * Bjellqvist, B., Basse, B., Olsen, E. and Celis, J.E.
    Reference points for comparisons of two-dimensional maps of proteins from
    different human cell types defined in a pH scale where isoelectric points
    correlate with polypeptide compositions. Electrophoresis 1994, 15, 529-539.
I designed the algorithm according to a note by David L. Tabb, available at:
http://fields.scripps.edu/DTASelect/20010710-pI-Algorithm.pdf
"""

from .constants import HYDROPHOBICITY_SCALE

positive_pKs = {"Nterm": 8.6, "K": 10.8, "R": 12.5, "H": 6.5}
negative_pKs = {"Cterm": 3.6, "D": 3.9, "E": 4.1, "C": 8.5, "Y": 10.1}
pKcterminal = {"D": 1, "E":1}
pKnterminal = {"A": 1,"M": 1,"S": 1,"P": 1,"T": 1,"V":1,"E": 1,}
charged_aas = ("K", "R", "H", "D", "E", "C", "Y")

scales = {
"EMBOSS":     {'Cterm': 3.6, 'pKAsp': 3.9,  'pKGlu': 4.1, 'pKCys': 8.5, 'pKTyr': 10.1, 'pk_his': 6.5, 'Nterm': 8.6, 'pKLys': 10.8, 'pKArg': 12.5},
"DTASelect":  {'Cterm': 3.1, 'pKAsp': 4.4,  'pKGlu': 4.4, 'pKCys': 8.5, 'pKTyr': 10.0, 'pk_his': 6.5, 'Nterm': 8.0, 'pKLys': 10.0, 'pKArg': 12.0},
"Solomon":    {'Cterm': 2.4, 'pKAsp': 3.9,  'pKGlu': 4.3, 'pKCys': 8.3, 'pKTyr': 10.1, 'pk_his': 6.0, 'Nterm': 9.6, 'pKLys': 10.5, 'pKArg': 12.5}, 
"Sillero":    {'Cterm': 3.2, 'pKAsp': 4.0,  'pKGlu': 4.5, 'pKCys': 9.0, 'pKTyr': 10.0, 'pk_his': 6.4, 'Nterm': 8.2, 'pKLys': 10.4, 'pKArg': 12.0},
"Rodwell":    {'Cterm': 3.1, 'pKAsp': 3.68, 'pKGlu': 4.25,'pKCys': 8.33,'pKTyr': 10.07,'pk_his': 6.0, 'Nterm': 8.0, 'pKLys': 11.5, 'pKArg': 11.5},
"Patrickios": {'Cterm': 4.2, 'pKAsp': 4.2,  'pKGlu': 4.2, 'pKCys': 0.0, 'pKTyr':  0.0, 'pk_his': 0.0, 'Nterm': 11.2,'pKLys': 11.2, 'pKArg': 11.2},
"Wikipedia":  {'Cterm': 3.65,'pKAsp': 3.9,  'pKGlu': 4.07,'pKCys': 8.18,'pKTyr': 10.46,'pk_his': 6.04,'Nterm': 8.2, 'pKLys': 10.54,'pKArg': 12.48},
"Grimsley":   {'Cterm': 3.3, 'pKAsp': 3.5,  'pKGlu': 4.2, 'pKCys': 6.8, 'pKTyr': 10.3, 'pk_his': 6.6, 'Nterm': 7.7, 'pKLys': 10.5, 'pKArg': 12.04},
'Lehninger':  {'Cterm': 2.34,'pKAsp': 3.86, 'pKGlu': 4.25,'pKCys': 8.33,'pKTyr': 10.0, 'pk_his': 6.0, 'Nterm': 9.69,'pKLys': 10.5, 'pKArg': 12.4},
'Bjellqvist': {'Cterm': 3.55,'pKAsp': 4.05, 'pKGlu': 4.45,'pKCys': 9.0, 'pKTyr': 10.0, 'pk_his': 5.98,'Nterm': 7.5, 'pKLys': 10.0, 'pKArg': 12.0},   
'IPC_peptide':{'Cterm': 2.383, 'pKAsp': 3.887, 'pKGlu': 4.317, 'pKCys': 8.297, 'pKTyr': 10.071, 'pk_his': 6.018, 'Nterm': 9.564, 'pKLys': 10.517, 'pKArg': 12.503},    # IPC peptide
'IPC_protein':{'Cterm': 2.869, 'pKAsp': 3.872, 'pKGlu': 4.412, 'pKCys': 7.555, 'pKTyr': 10.85,  'pk_his': 5.637, 'Nterm': 9.094, 'pKLys': 9.052,  'pKArg': 11.84},     # IPC protein 
'Toseland':   {'Cterm': 3.19,'pKAsp': 3.6,  'pKGlu': 4.29,'pKCys': 6.87,'pKTyr': 9.61, 'pk_his': 6.33,'Nterm': 8.71, 'pKLys': 10.45, 'pKArg':  12},
'Thurlkill':  {'Cterm': 3.67,'pKAsp': 3.67, 'pKGlu': 4.25,'pKCys': 8.55,'pKTyr': 9.84, 'pk_his': 6.54,'Nterm': 8.0, 'pKLys': 10.4, 'pKArg': 12.0},
'Nozaki':     {'Cterm': 3.8, 'pKAsp': 4.0,  'pKGlu': 4.4, 'pKCys': 9.5, 'pKTyr': 9.6,  'pk_his': 6.3, 'Nterm': 7.5, 'pKLys': 10.4, 'pKArg': 12},   
'Dawson':     {'Cterm': 3.2, 'pKAsp': 3.9,  'pKGlu': 4.3, 'pKCys': 8.3, 'pKTyr': 10.1, 'pk_his': 6.0, 'Nterm': 8.2, 'pKLys': 10.5, 'pKArg':  12},   
          }

aaDict = {'Asp':'D', 'Glu':'E', 'Cys':'C', 'Tyr':'Y', 'His':'H', 
          'Lys':'K', 'Arg':'R', 'Met':'M', 'Phe':'F', 'Leu':'L', 
          'Val':'V', 'Ala':'A', 'Gly':'G', 'Gln':'Q', 'Asn':'N',
          'Ile':'I', 'Trp':'W', 'Ser':'S', 'Thr':'T', 'Sec':'U',
          'Pro':'P', 'Xaa':'X', 'Sec':'U', 'Pyl':'O', 'Asx':'B',
          'Xle':'J', }

acidic = ['D', 'E', 'C', 'Y']
basic = ['K', 'R', 'H']

pKcterminal = {'D': 4.55, 'E': 4.75} 
pKnterminal = {'A': 7.59, 'M': 7.0, 'S': 6.93, 'P': 8.36, 'T': 6.82, 'V': 7.44, 'E': 7.7} 

#N-teminus, middle, C-terminus
promost ={
'K':[10.00,  9.80, 10.30],
'R':[11.50, 12.50, 11.50],
'H':[ 4.89,  6.08,  6.89],
'D':[ 3.57,  4.07,  4.57],
'E':[ 4.15,  4.45,  4.75],
'C':[ 8.00,  8.28,  9.00],
'Y':[ 9.34,  9.84, 10.34],
'U':[ 5.20,  5.43,  5.60], # ref (http://onlinelibrary.wiley.com/doi/10.1002/bip.21581/pdf)
}
           
promost_mid = {
"G":[7.50, 3.70],
"A":[7.58, 3.75],
"S":[6.86, 3.61],
"P":[8.36, 3.40],
"V":[7.44, 3.69],
"T":[7.02, 3.57],
"C":[8.12, 3.10],
"I":[7.48, 3.72],
"L":[7.46, 3.73],
"J":[7.46, 3.73],
"N":[7.22, 3.64],
"D":[7.70, 3.50],
"Q":[6.73, 3.57],
"K":[6.67, 3.40],
"E":[7.19, 3.50],
"M":[6.98, 3.68],
"H":[7.18, 3.17],
"F":[6.96, 3.98],
"R":[6.76, 3.41],
"Y":[6.83, 3.60],
"W":[7.11, 3.78],
"X":[7.26, 3.57],   #avg
"Z":[6.96, 3.535],  #("E"+"Q")/2
'B':[7.46, 3.57],   #("N"+"D")/2
'U':[5.20, 5.60], 
'O':[7.00, 3.50],     
}


from collections import Counter


class ProteinPropeties():
    def __init__(self, sequence):
        self.sequence = sequence

    def isoelectric_point(self):
        seq = self.sequence
        scale = "IPC_protein"

        pKCterm = scales[scale]['Cterm']
        pKAsp = scales[scale]['pKAsp']
        pKGlu = scales[scale]['pKGlu']
        pKCys = scales[scale]['pKCys']
        pKTyr = scales[scale]['pKTyr']
        pKHis = scales[scale]['pk_his']
        pKNterm = scales[scale]['Nterm']
        pKLys = scales[scale]['pKLys'] 
        pKArg = scales[scale]['pKArg']

        pH = 6.51             #starting po pI = 6.5 - theoretically it should be 7, but average protein pI is 6.5 so we increase the probability of finding the solution
        pHprev = 0.0         
        pHnext = 14.0        
        E = 0.01             #epsilon means precision [pI = pH +- E]
        temp = 0.01
        nterm=seq[0]
                
        while 1:             #the infinite loop
            QN1=-1.0/(1.0+pow(10,(pKCterm-pH)))                                        
            QN2=-seq.count('D')/(1.0+pow(10,(pKAsp-pH)))           
            QN3=-seq.count('E')/(1.0+pow(10,(pKGlu-pH)))           
            QN4=-seq.count('C')/(1.0+pow(10,(pKCys-pH)))           
            QN5=-seq.count('Y')/(1.0+pow(10,(pKTyr-pH)))        
            QP1=seq.count('H')/(1.0+pow(10,(pH-pKHis)))            
            QP2=1.0/(1.0+pow(10,(pH-pKNterm)))                
            QP3=seq.count('K')/(1.0+pow(10,(pH-pKLys)))           
            QP4=seq.count('R')/(1.0+pow(10,(pH-pKArg)))            
            NQ=QN1+QN2+QN3+QN4+QN5+QP1+QP2+QP3+QP4
                #print NQ
            if NQ<0.0:              #we are out of range, thus the new pH value must be smaller                     
                temp = pH
                pH = pH-((pH-pHprev)/2.0)
                pHnext = temp
                print(pH)
                #print "pH: ", pH, ", \tpHnext: ",pHnext
            else:
                temp = pH
                pH = pH + ((pHnext-pH)/2.0)
                pHprev = temp
                print(pH)
                #print "pH: ", pH, ",\tpHprev: ", pHprev

            if (pH-pHprev<E) and (pHnext-pH<E): #terminal condition, finding pI with given precision
                return pH

    def hydrophobicity(self):
        hydrophobicity = 0.0
        for aa in self.sequence:
            hydrophobicity += HYDROPHOBICITY_SCALE[aa]

        return hydrophobicity/len(self.sequence)


