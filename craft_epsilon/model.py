"""CRAFT-epsilon: classification, effective dipole, Clausius-Mossotti eps_inf and the
Onsager-Kirkwood-Froehlich relation; identifier resolution and the public functions."""
import math
import warnings
import numpy as np
from rdkit import Chem
from rdkit.Chem import rdMolDescriptors
from .conformers import find_xtb_binary, generate_conformers
from .constants import ATOM_TYPES, AU_TO_A3, DEBYE, EPS0, KB, N_A, TEMPERATURE_K, THERMAL_EXPANSION_PER_K, X_MAX
from .descriptors import conformer_descriptors
from .groups import atom_group, atom_type, base_label, family, fragments
from .molar_volume import molar_volume
from .parameters import PARAMETERS


def classify(smiles, frags):
    """Molecular class of the model: H-bond donor present -> H_bonding; else 0 / 1 / >= 2 polar groups."""
    if rdMolDescriptors.CalcNumHBD(Chem.MolFromSmiles(smiles)) > 0:
        return "H_bonding"
    return "non_polar" if not frags else "single_polar" if len(frags) == 1 else "multiple_polar"


def group_parameter(label, cls):
    """(parameter name, value) of the scale c_k of one polar group in the stage of the molecule's class."""
    if cls == "H_bonding":
        base = base_label(label)
        scales = PARAMETERS["hbond"]["fragment_scales"]
        if base in scales:
            return f"hbond:{base}", scales[base]
        return "hbond:rare", PARAMETERS["hbond"]["rare_scale"]
    fam = family(label)
    return f"family:{fam}", PARAMETERS["polar"]["family_scales"][fam]


def atom_type_scales(cls):
    """Scale s_t of every atom type for the stage of the class (ring and acyclic variants tied, H = 1)."""
    groups = dict(PARAMETERS["nonpolar"]["atom_group_scales"])
    if cls != "non_polar":
        groups.update(PARAMETERS["polar"]["atom_group_scales"])
    return np.array([groups.get(atom_group(t), 1.0) for t in ATOM_TYPES])


def epsilon_from_record(record, smiles, temperature=TEMPERATURE_K):
    """epsilon_r and intermediate quantities from a conformer record (output of generate_conformers).

    The model contains a temperature dependence through T in the orientation term and through the molar volume,
    V_m(T) = V_m(298.15 K) (1 + alpha_V (T - 298.15 K)) with the fixed alpha_V = THERMAL_EXPANSION_PER_K, as in the
    fit of the paper (data measured at other temperatures). It is meant to be used at 298.15 K only: the conformer
    ensembles, their Boltzmann weights and the group-contribution molar volume refer to 298.15 K, and the
    temperature dependence of association is not represented. The temperature argument exists to reproduce the
    fit; predict() and epsilon() always use 298.15 K."""
    mol = Chem.MolFromSmiles(smiles)
    frags = fragments(smiles)
    cls = classify(smiles, frags)
    vm_298 = molar_volume(mol)        # 298.15 K; corrected to T with the fixed alpha_V, as in the paper fit
    vm = vm_298 * (1 + THERMAL_EXPANSION_PER_K * (temperature - TEMPERATURE_K))
    alpha, per_conformer = conformer_descriptors(record, frags)

    # Clausius-Mossotti with atom-group scaled, Boltzmann-averaged atomic polarizabilities
    alpha_scaled = float(alpha @ atom_type_scales(cls))
    x = min(4 * math.pi / 3 * N_A * alpha_scaled * AU_TO_A3 * 1e-24 / vm, X_MAX)
    eps_inf = (1 + 2 * x) / (1 - x)

    # effective dipole: b mu_rest + sum_k c_k mu_k per conformer, Boltzmann mean square
    if cls == "non_polar":
        rest_scale = PARAMETERS["nonpolar"]["dipole_scale"]
        group_scale = {}
    else:
        rest_scale = PARAMETERS["hbond" if cls == "H_bonding" else "polar"]["rest_scale"]
        group_scale = {f.label: group_parameter(f.label, cls)[1] for f in frags}
    mu_sq_eff = mu_sq_raw = 0.0
    for weight, groups, remainder in per_conformer:
        vector = rest_scale * remainder + sum((group_scale[label] * v for label, v in groups.items()), np.zeros(3))
        raw = remainder + sum(groups.values(), np.zeros(3))
        mu_sq_eff += weight * float(vector @ vector)
        mu_sq_raw += weight * float(raw @ raw)

    # Onsager-Kirkwood-Froehlich relation, physical root
    a = N_A * mu_sq_eff * DEBYE ** 2 / (9 * EPS0 * KB * temperature * vm * 1e-6)
    b = eps_inf + a * (eps_inf + 2) ** 2
    eps = (b + math.sqrt(b * b + 8 * eps_inf ** 2)) / 4
    return {"epsilon": eps, "class": cls, "molar_volume_cm3_mol": vm, "molar_volume_298K_cm3_mol": vm_298,
            "eps_inf": eps_inf,
            "alpha_boltzmann_A3": float(alpha.sum()) * AU_TO_A3, "alpha_scaled_A3": alpha_scaled * AU_TO_A3,
            "mu_rms_xtb_D": math.sqrt(mu_sq_raw), "mu_eff_rms_D": math.sqrt(mu_sq_eff), "A": a,
            "rest_scale": rest_scale,
            "polar_groups": [{"label": f.label, "parameter": group_parameter(f.label, cls)[0],
                              "scale": group_parameter(f.label, cls)[1]} for f in frags],
            "n_conformers": len(record["conformers"])}


def resolve(identifier):
    """Canonical SMILES of a SMILES string or of a name (MoleculeResolver, cross-checked, neutral structure)."""
    mol = Chem.MolFromSmiles(identifier)
    if mol is not None:
        return Chem.MolToSmiles(mol), "SMILES"
    from moleculeresolver import MoleculeResolver   # imported here: only needed for names
    with warnings.catch_warnings():                  # services without API key announce themselves as warnings
        warnings.simplefilter("ignore")
        with MoleculeResolver() as mr:
            molecule = mr.find_single_molecule_crosschecked([identifier], ["name"], required_charge=0,
                                                            required_structure_type="neutral")
    smiles = getattr(molecule, "SMILES", None) if molecule else None
    mol = Chem.MolFromSmiles(smiles) if smiles else None
    if mol is None:
        raise ValueError(f"no cross-checked structure found for the name {identifier!r}")
    return Chem.MolToSmiles(mol), "name (MoleculeResolver)"


def check_domain(smiles):
    """Raise ValueError if the molecule is outside the domain of the model."""
    mol = Chem.MolFromSmiles(smiles)
    if len(Chem.GetMolFrags(mol)) > 1:
        raise ValueError("only single molecules are supported (the SMILES contains several fragments)")
    if Chem.GetFormalCharge(mol) != 0:   # charge-separated neutral groups (nitro, N-oxides) are allowed
        raise ValueError("only neutral molecules are supported")
    if sum(atom.GetAtomicNum() == 6 for atom in mol.GetAtoms()) < 2:
        raise ValueError("molecules with fewer than two carbon atoms are outside the domain of the model")
    for atom in Chem.AddHs(mol).GetAtoms():
        try:
            name = atom_type(atom)
        except ValueError as error:
            raise ValueError(f"atom {atom.GetSymbol()}{atom.GetIdx() + 1}: {error}") from error
        if name not in ATOM_TYPES:
            raise ValueError(f"atom type {name} is not covered by the model")
    molar_volume(mol)   # raises if an atom type has no molar-volume contribution


def predict(identifier, xtb_binary=None):
    """epsilon_r at 298.15 K and all intermediate quantities for a SMILES string or a name."""
    smiles, source = resolve(identifier)
    check_domain(smiles)
    binary = find_xtb_binary(xtb_binary)
    record = generate_conformers(smiles, binary)
    result = epsilon_from_record(record, smiles)
    result.update({"input": identifier, "smiles": smiles, "resolved_from": source, "temperature_K": TEMPERATURE_K,
                   "xtb": binary, "conformer_protocol": {k: v for k, v in record.items() if k != "conformers"}})
    return result


def epsilon(identifier, xtb_binary=None):
    """Static relative permittivity at 298.15 K for a SMILES string or a name."""
    return predict(identifier, xtb_binary)["epsilon"]
