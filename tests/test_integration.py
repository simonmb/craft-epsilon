"""End-to-end tests (opt-in): run with  pytest -m xtb  (needs xtb or downloads it)  and  pytest -m network.

The conformer ensembles depend slightly on the RDKit version and platform, so the comparison with the paper
values uses a tolerance."""
import pytest

from craft_epsilon import epsilon, predict
from craft_epsilon.model import resolve


@pytest.mark.xtb
@pytest.mark.parametrize("smiles", ["CCBr", "CC(C)S(C)(=O)=O", "CC(Cl)(Cl)Cl", "C#CCO", "CCNC=O"])
def test_full_workflow_matches_paper(smiles, reference):
    assert epsilon(smiles) == pytest.approx(reference[smiles]["epsilon_paper"], rel=0.03)


@pytest.mark.xtb
def test_predict_reports_protocol():
    r = predict("CCO")
    assert r["class"] == "H_bonding" and r["resolved_from"] == "SMILES"
    assert r["conformer_protocol"]["n_optimized"] >= 1 and r["n_conformers"] >= 1
    assert 15 < r["epsilon"] < 30


@pytest.mark.network
def test_resolve_name():
    assert resolve("ethanol") == ("CCO", "name (MoleculeResolver)")
