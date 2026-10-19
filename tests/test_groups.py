"""Polar groups, classes and atom types (rules of the paper)."""
from collections import Counter

import pytest
from rdkit import Chem

from craft_epsilon.constants import ATOM_TYPES
from craft_epsilon.groups import atom_group, atom_type, base_label, fragments
from craft_epsilon.model import classify


def groups(smiles):
    return Counter(base_label(f.label) for f in fragments(smiles))


@pytest.mark.parametrize("smiles, expected", [
    ("CCOC(=O)C(Cl)(Cl)Cl", {"ester": 1, "chloro": 3}),        # ethyl trichloroacetate (paper figure)
    ("COC(=O)c1ccccc1O", {"ester": 1, "phenol": 1}),            # methyl salicylate
    ("COC(=O)c1ccccc1C(=O)OC", {"ester": 2}),                   # dimethyl phthalate
    ("CC(C)=O", {"ketone": 1}),
    ("CCO", {"alcohol": 1}),
    ("CCC(=O)O", {"carboxylic_acid": 1}),
    ("CN(C)C=O", {"formamide": 1}),
    ("COCCOC", {"ether": 2}),
    ("CCCCCC", {}),
    ("c1ccccc1", {}),
])
def test_polar_groups(smiles, expected):
    assert groups(smiles) == Counter(expected)


def test_groups_do_not_overlap_and_include_polar_hydrogens():
    for smiles in ("OCC(O)CO", "NCCO", "CC(=O)NCCO", "O=C(O)CC(=O)O"):
        frags = fragments(smiles)
        atoms = [a for f in frags for a in f.atoms]
        assert len(atoms) == len(set(atoms))
        full = Chem.AddHs(Chem.MolFromSmiles(smiles))
        polar_h = {a.GetIdx() for a in full.GetAtoms() if a.GetAtomicNum() == 1
                   and a.GetNeighbors()[0].GetAtomicNum() in (7, 8)}
        assert polar_h <= set(atoms)


def test_context_suffix():
    assert {f.label for f in fragments("Clc1ccccc1")} == {"chloro_arom"}
    assert {f.label for f in fragments("O=C1CCCCC1")} == {"ketone_alicyc"}
    assert {f.label for f in fragments("CCCCl")} == {"chloro"}


@pytest.mark.parametrize("smiles, cls", [
    ("CCCCCC", "non_polar"), ("C=Cc1ccccc1", "non_polar"), ("C[Si](C)(C)C", "non_polar"),
    ("CC(C)=O", "single_polar"), ("CCBr", "single_polar"),
    ("CCOC(=O)C(Cl)(Cl)Cl", "multiple_polar"), ("ClCCCl", "multiple_polar"),
    ("CCO", "H_bonding"), ("CCS", "H_bonding"), ("c1cc[nH]c1", "H_bonding"), ("CC(N)=O", "H_bonding"),
])
def test_classification(smiles, cls):
    assert classify(smiles, fragments(smiles)) == cls


def test_atom_types_and_tied_ring_variants():
    mol = Chem.AddHs(Chem.MolFromSmiles("O=C1CCCC1CC#N"))
    types = {atom_type(a) for a in mol.GetAtoms()}
    assert types <= set(ATOM_TYPES)
    assert {"C_SP2_R", "O_SP2_A", "C_SP3_R", "C_SP3_A", "C_SP_A", "N_SP_A", "H"} <= types
    assert atom_group("C_SP3_R") == atom_group("C_SP3_A") == "C_SP3"
    assert atom_group("Cl") == "Cl"
