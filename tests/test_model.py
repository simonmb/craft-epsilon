"""The model on stored paper conformer records: regression against the paper fit, physics checks, domain checks."""
import math

import pytest

from craft_epsilon import PARAMETERS, epsilon_from_record
from craft_epsilon.groups import fragments
from craft_epsilon.model import check_domain, group_parameter, resolve


def test_reproduces_paper_predictions(reference, records):
    """Same conformer records -> same epsilon as the fit of the paper (16 liquids, all four classes)."""
    for smiles, entry in reference.items():
        result = epsilon_from_record(records[smiles], smiles)
        assert result["class"] == entry["class"], smiles
        assert result["epsilon"] == pytest.approx(entry["epsilon_paper"], rel=1e-8), smiles


def test_onsager_relation_is_satisfied(reference, records):
    for smiles in reference:
        r = epsilon_from_record(records[smiles], smiles)
        eps, eps_inf, a = r["epsilon"], r["eps_inf"], r["A"]
        lhs = (2 * eps + eps_inf) * (eps - eps_inf) / (eps * (eps_inf + 2) ** 2)
        assert lhs == pytest.approx(a, rel=1e-10), smiles
        assert eps > eps_inf > 1.0


def test_without_dipole_epsilon_equals_eps_inf(reference, records):
    """Setting all dipole scales to zero leaves only the electronic part (Clausius-Mossotti)."""
    smiles = "CCBr"
    record = records[smiles]
    saved = PARAMETERS["polar"]["rest_scale"], dict(PARAMETERS["polar"]["family_scales"])
    try:
        PARAMETERS["polar"]["rest_scale"] = 0.0
        PARAMETERS["polar"]["family_scales"].update({k: 0.0 for k in PARAMETERS["polar"]["family_scales"]})
        r = epsilon_from_record(record, smiles)
    finally:
        PARAMETERS["polar"]["rest_scale"] = saved[0]
        PARAMETERS["polar"]["family_scales"].update(saved[1])
    assert r["A"] == pytest.approx(0.0, abs=1e-12)
    assert r["epsilon"] == pytest.approx(r["eps_inf"], rel=1e-12)


def test_permittivity_decreases_with_temperature(reference, records):
    smiles = "CCNC=O"
    values = [epsilon_from_record(records[smiles], smiles, temperature=t)["epsilon"] for t in (280.0, 298.15, 320.0)]
    assert values[0] > values[1] > values[2]


def test_molar_volume_follows_thermal_expansion(reference, records):
    """V_m(T) = V_m(298.15 K) (1 + alpha_V (T - 298.15 K)) with the fixed alpha_V of the paper fit."""
    from craft_epsilon.constants import THERMAL_EXPANSION_PER_K
    smiles = "CCBr"
    r = epsilon_from_record(records[smiles], smiles, temperature=350.0)
    assert r["molar_volume_298K_cm3_mol"] == pytest.approx(reference[smiles]["molar_volume_smi2DC"], abs=0.01)
    assert r["molar_volume_cm3_mol"] == pytest.approx(
        r["molar_volume_298K_cm3_mol"] * (1 + THERMAL_EXPANSION_PER_K * (350.0 - 298.15)), rel=1e-12)


def test_weights_and_effective_dipole(reference, records):
    smiles = "CC(Br)C(C)Br"            # two diastereomer ensembles with mirror twins
    r = epsilon_from_record(records[smiles], smiles)
    assert r["n_conformers"] == len(records[smiles]["conformers"])
    assert len({c["diastereomer"] for c in records[smiles]["conformers"]}) == 2
    assert r["mu_eff_rms_D"] > 0 and r["mu_rms_xtb_D"] > 0
    assert math.isfinite(r["epsilon"])


@pytest.mark.parametrize("smiles, parameter", [
    ("CCO", "hbond:alcohol"), ("c1cc[nH]c1", "hbond:rare"),     # azole: pooled into the rare scale
    ("CC(C)=O", "family:carbonyl"), ("CCBr", "family:halogen"), ("CO[Si](C)(C)C", "family:silicon"),
])
def test_group_parameters(smiles, parameter):
    from craft_epsilon.model import classify
    frags = fragments(smiles)
    assert group_parameter(frags[0].label, classify(smiles, frags))[0] == parameter


def test_unseen_hydrogen_bonding_fragment_uses_rare_scale():
    assert group_parameter("nitrile", "H_bonding") == ("hbond:rare", PARAMETERS["hbond"]["rare_scale"])


def test_resolve_canonicalises_smiles():
    assert resolve("OCC") == ("CCO", "SMILES")
    assert resolve("C(=O)(C)C") == ("CC(C)=O", "SMILES")


@pytest.mark.parametrize("smiles, message", [
    ("[Na+].[Cl-]", "single molecules"),
    ("CC[O-]", "neutral"),
    ("C", "two carbon atoms"),
    ("C[Ge](C)(C)C", "not covered"),
])
def test_domain_errors(smiles, message):
    with pytest.raises(ValueError, match=message):
        check_domain(smiles)


def test_charge_separated_neutral_groups_are_accepted():
    check_domain("CC[N+](=O)[O-]")                               # nitroethane
    check_domain("[O-][n+]1ccccc1C")                             # 2-picoline N-oxide
