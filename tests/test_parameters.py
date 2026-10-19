"""The fitted parameters: 41 in total (3 + 21 + 17), within the bounds of the fit."""
from craft_epsilon import PARAMETERS
from craft_epsilon.groups import family

SCALE_BOUND, DIPOLE_BOUND = 10.0, 5.0
FAMILIES = {"carbonyl", "halogen", "nitro", "nitrogen", "other", "oxygen", "phosphorus", "silicon", "sulfur"}


def stage_counts():
    nonpolar = len(PARAMETERS["nonpolar"]["atom_group_scales"]) + 1
    polar = len(PARAMETERS["polar"]["atom_group_scales"]) + len(PARAMETERS["polar"]["family_scales"]) + 1
    hbond = len(PARAMETERS["hbond"]["fragment_scales"]) + 2
    return nonpolar, polar, hbond


def test_parameter_count_per_stage_and_total():
    assert stage_counts() == (3, 21, 17)
    assert sum(stage_counts()) == 41


def test_parameters_within_fit_bounds():
    scales = [*PARAMETERS["nonpolar"]["atom_group_scales"].values(), *PARAMETERS["polar"]["atom_group_scales"].values(),
              *PARAMETERS["polar"]["family_scales"].values(), *PARAMETERS["hbond"]["fragment_scales"].values(),
              PARAMETERS["hbond"]["rare_scale"]]
    assert all(0.0 <= value <= SCALE_BOUND for value in scales)
    for value in (PARAMETERS["nonpolar"]["dipole_scale"], PARAMETERS["polar"]["rest_scale"],
                  PARAMETERS["hbond"]["rest_scale"]):
        assert 0.0 <= value <= DIPOLE_BOUND


def test_stages_fit_disjoint_atom_groups_and_all_families():
    assert not set(PARAMETERS["nonpolar"]["atom_group_scales"]) & set(PARAMETERS["polar"]["atom_group_scales"])
    assert "H" not in PARAMETERS["polar"]["atom_group_scales"]          # s_H = 1 throughout
    assert set(PARAMETERS["polar"]["family_scales"]) == FAMILIES
    assert family("ester") == "carbonyl" and family("chloro_arom") == "halogen" and family("xyz") == "other"


def test_hydrogen_bonding_parameters_at_lower_bound():
    # paper: b_HB and the chloro and ether scales of the hydrogen-bonding stage converge to zero
    hbond = PARAMETERS["hbond"]
    assert hbond["rest_scale"] < 1e-10
    assert hbond["fragment_scales"]["chloro"] < 1e-10 and hbond["fragment_scales"]["ether"] < 1e-10
