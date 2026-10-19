"""Liquid molar volume at 298.15 K from the group contributions of Bouteloup and Mathieu.

R. Bouteloup, D. Mathieu, Predicting dielectric constants of pure liquids: fragment-based
Kirkwood-Froehlich model applicable over a wide range of polarity, Phys. Chem. Chem. Phys. 21 (2019)
11043-11057. https://doi.org/10.1039/C9CP01704F
"""
from collections import Counter


# Atom-type and ring contributions to V_m (cm^3/mol) of R. Bouteloup, D. Mathieu, Phys. Chem. Chem. Phys.
# 21 (2019) 11043-11057, https://doi.org/10.1039/C9CP01704F (as implemented in their program smi2DC).
VM_CONTRIBUTIONS = {
    'B30': 7.61,
    'Br10': 27.52,
    'C10': 14.68,
    'C20': 9.89,
    'C21': 22.77,
    'C30': -0.0,
    'C30a': -0.0,
    'C31': 13.23,
    'C31a': 13.23,
    'C32': 27.17,
    'C40': -8.4,
    'C41': 4.46,
    'C42': 16.57,
    'C43': 29.58,
    'Cl10': 24.74,
    'F10': 17.63,
    'Ge40': 9.36,
    'I10': 35.64,
    'N10': 14.09,
    'N20': 7.42,
    'N20a': 7.42,
    'N21': 18.14,
    'N30': -3.08,
    'N30a': -3.08,
    'N31': 7.74,
    'N31a': 7.74,
    'N32': 17.81,
    'N43': 5.36,
    'O10': 14.89,
    'O20': 6.25,
    'O20a': 6.25,
    'O21': 11.78,
    'P30': 10.42,
    'P40': -1.94,
    'P41': 10.06,
    'R5': 9.41,
    'R6': 6.89,
    'R<5': 10.89,
    'R>6': 3.75,
    'S10': 25.92,
    'S20': 14.9,
    'S20a': 14.9,
    'S21': 26.14,
    'S30': 5.58,
    'S40': -3.74,
    'Se20': 19.0,
    'Se20a': 19.0,
    'Si40': 9.28,
    'Si41': 23.35,
    'Sn40': 14.47,
    'Ti40': 6.09,
    'aromat': 1.82,
}


def get_atom_code(atom):
    """ return an atom code consistent with the keys of the 'vi' dictionary """
    # be careful to take into account nD = number of deuterium atoms !
    nD = len([x for x in atom.GetNeighbors() if x.GetMass()>2 and x.GetSymbol()=='H'])
    code = atom.GetSymbol() + str(atom.GetTotalDegree()) + str(atom.GetTotalNumHs()+nD)
    code += 'a'*int(atom.GetIsAromatic())
    return code


def get_ring_descriptors(mol, maxi=6, mini=5):
    """ return dict of ring descriptors for molecule provided as input """
    dic = Counter()
    ri = mol.GetRingInfo()
    for ring in ri.AtomRings():
        size = len(ring)
        if size > maxi:
            label = 'R>'+str(maxi)
        elif size < mini:
            label = 'R<'+str(mini)
        else:
            label = 'R%u' %len(ring)
        dic[label] += 1
        # specific contribution of aromatic rings
        atoms = [mol.GetAtomWithIdx(i) for i in ring]
        if all(at.GetIsAromatic() for at in atoms):
            dic['aromat'] += 1
    return dic


def molar_volume(mol):
    """Liquid molar volume at 298.15 K in cm^3/mol from the atom-type and ring contributions of
    R. Bouteloup, D. Mathieu, Phys. Chem. Chem. Phys. 21 (2019) 11043-11057, https://doi.org/10.1039/C9CP01704F
    (contributions VM_CONTRIBUTIONS, atom codes get_atom_code, ring terms get_ring_descriptors)."""
    atoms = Counter(get_atom_code(atom) for atom in mol.GetAtoms())
    missing = set(atoms) - set(VM_CONTRIBUTIONS)
    if missing:
        raise ValueError(f"no molar-volume contribution for atom type(s) {', '.join(sorted(missing))}")
    rings = get_ring_descriptors(mol)
    return (sum(n * VM_CONTRIBUTIONS[code] for code, n in atoms.items())
            + sum(n * VM_CONTRIBUTIONS[code] for code, n in rings.items()))
