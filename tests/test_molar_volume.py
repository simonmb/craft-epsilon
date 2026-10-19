"""Molar volume of Bouteloup and Mathieu (2019): identical to their program smi2DC."""
import pytest
from rdkit import Chem

from craft_epsilon.molar_volume import VM_CONTRIBUTIONS, get_atom_code, molar_volume


def test_reproduces_smi2dc(reference):
    for smiles, entry in reference.items():
        assert molar_volume(Chem.MolFromSmiles(smiles)) == pytest.approx(entry["molar_volume_smi2DC"], rel=1e-12)


def test_hydrogens_counted_on_heavy_atoms():
    # atom code = element + total degree + number of H: ethanol = CH3 (C43) + CH2 (C42) + OH (O21), no ring terms
    mol = Chem.MolFromSmiles("CCO")
    assert sorted(get_atom_code(atom) for atom in mol.GetAtoms()) == ["C42", "C43", "O21"]
    assert molar_volume(mol) == pytest.approx(VM_CONTRIBUTIONS["C43"] + VM_CONTRIBUTIONS["C42"] + VM_CONTRIBUTIONS["O21"])


def test_ring_and_aromatic_terms():
    benzene = Chem.MolFromSmiles("c1ccccc1")
    expected = 6 * VM_CONTRIBUTIONS["C31a"] + VM_CONTRIBUTIONS["R6"] + VM_CONTRIBUTIONS["aromat"]
    assert molar_volume(benzene) == pytest.approx(expected)


def test_volume_grows_with_chain_length():
    volumes = [molar_volume(Chem.MolFromSmiles("C" * n)) for n in range(2, 10)]
    increments = [b - a for a, b in zip(volumes, volumes[1:])]
    assert all(inc == pytest.approx(VM_CONTRIBUTIONS["C42"]) for inc in increments)   # one CH2 per step


def test_unknown_atom_type_raises():
    with pytest.raises(ValueError, match="no molar-volume contribution"):
        molar_volume(Chem.MolFromSmiles("C[Te]C"))
