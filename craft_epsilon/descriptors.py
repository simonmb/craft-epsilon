"""Boltzmann weights, Boltzmann-averaged atomic polarizabilities and polar-group dipole vectors."""
from collections import defaultdict
import numpy as np
from rdkit import Chem
from .constants import ATOM_TYPES, AU_TO_DEBYE, HARTREE_TO_KCAL, KAPPA_D, R_KCAL, TEMPERATURE_K
from .groups import atom_type


def boltzmann_weights(energies_hartree):
    relative = (np.asarray(energies_hartree) - np.min(energies_hartree)) * HARTREE_TO_KCAL
    w = np.exp(-relative / (R_KCAL * TEMPERATURE_K))
    return w / w.sum()


def ensemble_weights(conformers):
    """Boltzmann weights within each diastereomer ensemble, equal weight between ensembles."""
    ensembles = [c.get("diastereomer", 0) for c in conformers]
    labels = sorted(set(ensembles))
    weights = np.zeros(len(conformers))
    for label in labels:
        members = [i for i, ensemble in enumerate(ensembles) if ensemble == label]
        weights[members] = boltzmann_weights(
            [conformers[i]["alpb"]["qrrho_free_energy_hartree"] for i in members]) / len(labels)
    return weights


def fragment_anchor(mol, atoms):
    """Atoms whose mean position anchors a fragment's charge dipole: the heavy atoms bonded to the
    fragment, else any bonded atoms, else the fragment itself. A fragment with a net charge then
    carries its local bond polarity instead of a term that grows with its distance from the
    molecular centroid; the net charge at the anchor stays in the remainder."""
    members = set(atoms)
    outside = sorted({n.GetIdx() for a in atoms for n in mol.GetAtomWithIdx(a).GetNeighbors()
                      if n.GetIdx() not in members})
    heavy = [i for i in outside if mol.GetAtomWithIdx(i).GetAtomicNum() > 1]
    return heavy or outside or atoms


def conformer_descriptors(record, frags):
    """Boltzmann-averaged atom-type polarizabilities and per-conformer fragment vectors."""
    conformers = record["conformers"]
    weights = ensemble_weights(conformers)
    alpha = np.zeros(len(ATOM_TYPES))
    type_index = {t: i for i, t in enumerate(ATOM_TYPES)}
    per_conformer = []
    for conformer, weight in zip(conformers, weights):
        mol = Chem.MolFromMolBlock(conformer["molblock"], removeHs=False)
        if mol is None:
            raise ValueError("molblock could not be parsed")
        alpb = conformer["alpb"]
        positions = np.asarray(mol.GetConformer().GetPositions())
        charges = np.asarray(alpb["charges_e"], dtype=float)
        # New records store the components in Debye. Older records used this key
        # for raw xTB e*bohr values, so convert those when the unit marker is absent.
        total = np.asarray(alpb["dipole_vector_debye"], dtype=float)
        if alpb.get("dipole_vector_unit") != "debye":
            total *= AU_TO_DEBYE
        if len(charges) != mol.GetNumAtoms():
            raise ValueError("charge count does not match molblock")
        for atom, value in zip(mol.GetAtoms(), alpb["alpha_atoms_au"]):
            name = atom_type(atom)
            if name not in type_index:
                raise ValueError(f"atom type {name} is not in the regression list")
            alpha[type_index[name]] += weight * float(value)
        atomic_dipoles = np.asarray(alpb["atomic_dipoles_debye"], dtype=float)
        groups = defaultdict(lambda: np.zeros(3))
        for fragment in frags:
            atoms = list(fragment.atoms)
            anchor = positions[fragment_anchor(mol, atoms)].mean(axis=0)
            groups[fragment.label] += (KAPPA_D * (charges[atoms, None] * (positions[atoms] - anchor)).sum(axis=0)
                                       + atomic_dipoles[atoms].sum(axis=0))
        remainder = total - sum(groups.values(), np.zeros(3))
        per_conformer.append((float(weight), dict(groups), remainder))
    return alpha, per_conformer
