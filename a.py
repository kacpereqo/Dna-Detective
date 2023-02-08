def compute_pI(protein_sequence):
    amino_acids = 'RHKDESTNQCGPAVILFMYW'
    pKs = [4.5, 10.8, 13.0, 2.1, 9.6, 9.2, 8.8, 10.6, 5.0,
           2.2, 2.1, 14.0, 2.8, 3.9, 1.8, 7.6, 5.6, 9.3, 6.0]
    charges = [-4, -3, -3, -3, -3, -3, -3, -3,
               - 2, -2, -2, -2, -2, -2, -1, -1, -1, -1, 0]

    n_residues = len(protein_sequence)
    pIs = [0.0] * n_residues

    for i in range(n_residues):
        res = protein_sequence[i]
        if res not in amino_acids:
            raise ValueError(f"Amino acid {res} not recognized")
        idx = amino_acids.index(res)
        pIs[i] = pKs[idx] * charges[idx]

    pI = 7.0
    min_diff = 0.01

    for pH in range(1, 14):
        pos_charge = 0.0
        neg_charge = 0.0
        for i in range(n_residues):
            if pIs[i] < pH:
                neg_charge += 10 ** (pH - pIs[i])
            else:
                pos_charge += 10 ** (pIs[i] - pH)
        charge = pos_charge - neg_charge
        if abs(charge) < min_diff:
            pI = pH
            min_diff = abs(charge)

    return pI


print(compute_pI("MNENLFASFIAPTILGLP"))
