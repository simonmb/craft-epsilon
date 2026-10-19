"""GFN2-xTB/ALPB conformer ensembles: RDKit ETKDGv3 embedding, diastereomer enumeration, polar-H
variants, xtb optimization / single point / qRRHO Hessian, symmetry-corrected clustering and mirror
twins (the protocol used to fit CRAFT-epsilon)."""
import copy
import hashlib
import json
import logging
import math
import os
import platform
import random
import shutil
import subprocess
import tempfile
import time
from concurrent.futures import ThreadPoolExecutor
from concurrent.futures import as_completed
from itertools import product
from pathlib import Path
from urllib.request import urlopen
from urllib.request import urlretrieve
import numpy as np
from rdkit import Chem
from rdkit.Chem import AllChem
from rdkit.Geometry import Point3D
from spyrmsd import graph as spyrmsd_graph
from .constants import (
    AU_DIPOLE_TO_DEBYE,
    BOHR_PER_ANGSTROM,
    DIPOLE_PATTERN,
    DIPOLE_RECONSTRUCTION_TOL_AU,
    EMBED_PRUNE_RMS_A,
    ENERGY_PATTERN,
    MAX_ENUMERATED_ELEMENTS,
    MIRROR_DISTINCT_RMS_A,
    N_SEEDS,
    N_WORKERS,
    POLARIZABILITY_HEADER,
    POLAR_Z,
    POSTOPT_MU_D,
    POSTOPT_RMS_A,
    QRRHO_FREE_ENERGY,
    RMSD_CHUNK,
    SOLVENT,
    STEREO_COUPLING_BONDS,
    STEREO_POLAR_BONDS,
    TEMPERATURE_K,
)


def _polar_h_groups(molecule):
    """Return rotatable terminal OH, SH, and NH2 groups with their torsion atoms."""
    sp3 = Chem.HybridizationType.SP3
    n_h_needed = {7: 2, 8: 1, 16: 1}
    groups = []
    for atom in molecule.GetAtoms():
        if atom.GetAtomicNum() not in n_h_needed:
            continue
        hydrogens = [neighbor for neighbor in atom.GetNeighbors() if neighbor.GetAtomicNum() == 1]
        heavy_neighbors = [neighbor for neighbor in atom.GetNeighbors() if neighbor.GetAtomicNum() > 1]
        if len(hydrogens) != n_h_needed[atom.GetAtomicNum()] or len(heavy_neighbors) != 1:
            continue
        parent = heavy_neighbors[0]
        staggered = parent.GetHybridization() == sp3
        if atom.GetAtomicNum() == 7 and not staggered:
            continue
        others = [neighbor for neighbor in parent.GetNeighbors() if neighbor.GetIdx() != atom.GetIdx()]
        if not others:
            continue
        reference = max(others, key=lambda neighbor: neighbor.GetAtomicNum() > 1)
        groups.append(("staggered" if staggered else "planar",
                       (reference.GetIdx(), parent.GetIdx(), atom.GetIdx(), hydrogens[0].GetIdx()),
                       [hydrogen.GetIdx() for hydrogen in hydrogens]))
    return groups


def randomize_polar_h(molecule, seed):
    """Set each rotatable polar-H torsion to a random physically distinct orientation."""
    from rdkit.Chem import rdMolTransforms

    rng = random.Random(seed)
    angles = {"planar": [0.0, 180.0], "staggered": [60.0, 180.0, 300.0]}
    groups = _polar_h_groups(molecule)
    for conformer in molecule.GetConformers():
        for kind, (a, b, c, d), _ in groups:
            rdMolTransforms.SetDihedralDeg(conformer, a, b, c, d, rng.choice(angles[kind]))
    return len(groups)


def _is_hbond_acceptor(atom):
    """Whether an atom can accept an intramolecular H bond for the seed variant."""
    atomic_number, formal_charge = atom.GetAtomicNum(), atom.GetFormalCharge()
    if atomic_number == 8:
        return formal_charge <= 0
    if atomic_number != 7 or formal_charge > 0:
        return False
    if atom.GetIsAromatic():
        return atom.GetTotalNumHs(includeNeighbors=True) == 0 and atom.GetDegree() == 2
    if atom.GetHybridization() == Chem.HybridizationType.SP3:
        return True
    return any(bond.GetBondType() in (Chem.BondType.DOUBLE, Chem.BondType.TRIPLE)
               for bond in atom.GetBonds())


def force_intramolecular_hbonds(molecule, max_bonds, min_bonds=3, max_dist=2.6,
                                min_angle=95.0, step=10.0):
    """Aim rotatable polar H atoms at feasible intramolecular acceptors."""
    from rdkit.Chem import rdMolTransforms

    groups = _polar_h_groups(molecule)
    topology = Chem.GetDistanceMatrix(molecule)
    acceptors = [atom.GetIdx() for atom in molecule.GetAtoms() if _is_hbond_acceptor(atom)]
    candidates = [(group, [acceptor for acceptor in acceptors
                           if acceptor != group[1][2] and min_bonds <= topology[group[1][2], acceptor] <= max_bonds])
                  for group in groups]
    candidates = [(group, acceptors) for group, acceptors in candidates if acceptors]
    counts = []
    for conformer in molecule.GetConformers():
        options = []
        for (_, (a, b, c, d), hydrogens), nearby_acceptors in candidates:
            start = rdMolTransforms.GetDihedralDeg(conformer, a, b, c, d)
            best = None
            for torsion in np.arange(0.0, 360.0, step):
                rdMolTransforms.SetDihedralDeg(conformer, a, b, c, d, float(torsion))
                donor = np.array(conformer.GetAtomPosition(c))
                for hydrogen in hydrogens:
                    hydrogen_position = np.array(conformer.GetAtomPosition(hydrogen))
                    for acceptor in nearby_acceptors:
                        acceptor_position = np.array(conformer.GetAtomPosition(acceptor))
                        distance = np.linalg.norm(acceptor_position - hydrogen_position)
                        v1, v2 = donor - hydrogen_position, acceptor_position - hydrogen_position
                        angle = np.degrees(np.arccos(np.clip(v1 @ v2 / (np.linalg.norm(v1) * np.linalg.norm(v2)), -1.0, 1.0)))
                        if distance <= max_dist and angle >= min_angle and (best is None or distance < best[0]):
                            best = (distance, float(torsion), acceptor)
            rdMolTransforms.SetDihedralDeg(conformer, a, b, c, d, start)
            if best is not None:
                options.append((best[0], best[1], (a, b, c, d), best[2]))
        formed = set()
        for _, torsion, (a, b, c, d), acceptor in sorted(options):
            if (acceptor, c) not in formed:
                rdMolTransforms.SetDihedralDeg(conformer, a, b, c, d, torsion)
                formed.add((c, acceptor))
        counts.append(len(formed))
    return counts


def amplify_polar_h_conformers(molecule, seed=0, n_random=2, hbond_variant=True,
                               max_bonds=5, min_groups=2):
    """Add polar-H variants for conformers with at least two rotatable polar-H groups."""
    if len(_polar_h_groups(molecule)) < min_groups or molecule.GetNumConformers() == 0:
        return []
    variants = []
    for offset in range(n_random):
        variant = Chem.Mol(molecule)
        randomize_polar_h(variant, seed + offset)
        variants.append(variant)
    if hbond_variant:
        variant = Chem.Mol(molecule)
        randomize_polar_h(variant, seed + n_random)
        force_intramolecular_hbonds(variant, max_bonds)
        variants.append(variant)
    added = []
    for variant in variants:
        for conformer in variant.GetConformers():
            added.append(molecule.AddConformer(Chem.Conformer(conformer), assignId=True))
    return added


def run_xtb(binary_path, arguments, working_directory):
    result = subprocess.run(
        [binary_path, "mol.xyz", "--gfn", "2", "--chrg", str(arguments["charge"]), "--norestart", *arguments["calculation"]],
        cwd=working_directory,
        capture_output=True,
        text=True,
        encoding="utf-8",
        errors="replace",
        timeout=1200,
        creationflags=getattr(subprocess, "CREATE_NO_WINDOW", 0),
        env={**os.environ, "OMP_NUM_THREADS": "1", "MKL_NUM_THREADS": "1", "OPENBLAS_NUM_THREADS": "1"},
    )
    if result.returncode:
        raise RuntimeError((result.stdout + result.stderr)[-1500:])
    return result.stdout


def write_xyz(molecule, conformer_id, path):
    conformer = molecule.GetConformer(conformer_id)
    with path.open("w", encoding="utf-8") as file:
        file.write(f"{molecule.GetNumAtoms()}\n\n")
        for atom in molecule.GetAtoms():
            position = conformer.GetAtomPosition(atom.GetIdx())
            file.write(f"{atom.GetSymbol()} {position.x:.8f} {position.y:.8f} {position.z:.8f}\n")


def molblock_from_xyz(molecule, xyz_path):
    lines = xyz_path.read_text(encoding="utf-8").splitlines()
    optimized = Chem.Mol(molecule)
    optimized.RemoveAllConformers()
    conformer = Chem.Conformer(optimized.GetNumAtoms())
    for index, line in enumerate(lines[2:2 + optimized.GetNumAtoms()]):
        _, x, y, z = line.split()[:4]
        conformer.SetAtomPosition(index, Point3D(float(x), float(y), float(z)))
    optimized.AddConformer(conformer, assignId=True)
    return Chem.MolToMolBlock(optimized)


def read_xtb_json(path):
    """xtbout.json written by --json; xTB does not escape the backslashes of Windows paths."""
    return json.loads(path.read_text(encoding="utf-8").replace("\\", "\\\\"))


def xyz_positions_angstrom(xyz_text, n_atoms):
    return np.array([line.split()[1:4] for line in xyz_text.splitlines()[2:2 + n_atoms]], dtype=float)


def parse_single_point(output, molecule, xtb_json, xyz_text):
    energy_matches = ENERGY_PATTERN.findall(output)
    dipole_matches = DIPOLE_PATTERN.findall(output)
    lines = output.splitlines()
    headers = [index for index, line in enumerate(lines) if POLARIZABILITY_HEADER.search(line)]
    if not energy_matches or not dipole_matches or not headers:
        raise ValueError("xTB output is missing energy, dipole, or atomic polarizability data")
    rows = lines[headers[-1] + 1:headers[-1] + 1 + molecule.GetNumAtoms()]
    atomic_polarizabilities = []
    atomic_charges = []
    for atom, row in zip(molecule.GetAtoms(), rows):
        fields = row.split()
        # xTB prints: index, Z, element, covCN, q, C6AA, alpha(0).
        # Keep the element check here so a changed table layout cannot silently
        # turn coordination numbers or another property into atomic charges.
        if (len(fields) < 7 or int(fields[0]) != atom.GetIdx() + 1 or
                int(fields[1]) != atom.GetAtomicNum() or fields[2] != atom.GetSymbol()):
            raise ValueError(f"xTB atomic data do not match the molecule: {row}")
        atomic_charges.append(float(fields[4]))
        atomic_polarizabilities.append(float(fields[-1]))
    if len(atomic_polarizabilities) != molecule.GetNumAtoms():
        raise ValueError("xTB returned an incomplete atomic polarizability table")
    # The JSON charges are the same as the table's but printed with 8 instead of 3 decimals.
    json_charges = np.asarray(xtb_json["partial charges"], dtype=float)
    atomic_dipoles_au = np.asarray(xtb_json["atomic dipole moments"], dtype=float)
    if json_charges.shape != (molecule.GetNumAtoms(),) or atomic_dipoles_au.shape != (molecule.GetNumAtoms(), 3):
        raise ValueError("xtbout.json charges or atomic dipoles do not match the molecule")
    if np.max(np.abs(json_charges - np.asarray(atomic_charges))) > 6e-4:
        raise ValueError("xtbout.json charges differ from the printed atomic charges")
    atomic_charges = json_charges.tolist()
    expected_charge = Chem.GetFormalCharge(molecule)
    charge_tolerance = max(0.02, 0.0006 * molecule.GetNumAtoms())
    if abs(sum(atomic_charges) - expected_charge) > charge_tolerance:
        raise ValueError("xTB atomic partial charges do not sum to the molecular charge")
    x, y, z, magnitude = map(float, dipole_matches[-1])
    # GFN2 dipole = point-charge part + sum of atomic dipoles (both e*bohr); a mismatch
    # means the atom order or the units of the stored arrays are wrong.
    positions_bohr = xyz_positions_angstrom(xyz_text, molecule.GetNumAtoms()) * BOHR_PER_ANGSTROM
    positions_bohr -= positions_bohr.mean(axis=0)
    reconstructed = (json_charges[:, None] * positions_bohr).sum(axis=0) + atomic_dipoles_au.sum(axis=0)
    if np.max(np.abs(reconstructed - np.array([x, y, z]))) > DIPOLE_RECONSTRUCTION_TOL_AU:
        raise ValueError("charges and atomic dipoles do not reproduce the xTB molecular dipole")
    return {
        "energy_hartree": float(energy_matches[-1]),
        "dipole_debye": magnitude,
        "dipole_vector_debye": [component * AU_DIPOLE_TO_DEBYE for component in (x, y, z)],
        "dipole_vector_unit": "debye",
        "charges_e": atomic_charges,
        "atomic_dipoles_debye": (atomic_dipoles_au * AU_DIPOLE_TO_DEBYE).tolist(),
        "alpha_atoms_au": atomic_polarizabilities,
    }


def parse_qrrho_free_energy(output):
    matches = QRRHO_FREE_ENERGY.findall(output)
    if not matches:
        raise ValueError("xTB Hessian output is missing its qRRHO total free energy")
    return float(matches[-1])


def calculate_conformer(binary_path, molecule, conformer_id, solvent):
    charge = Chem.GetFormalCharge(molecule)
    solvent_arguments = [] if solvent == "gas" else ["--alpb", solvent]
    command = {"charge": charge, "calculation": []}
    with tempfile.TemporaryDirectory() as temporary:
        working_directory = Path(temporary)
        geometry_path = working_directory / "mol.xyz"
        write_xyz(molecule, conformer_id, geometry_path)
        optimized_path = working_directory / "xtbopt.xyz"
        for attempt in range(1, 5):
            command["calculation"] = ["--opt", "tight", *solvent_arguments]
            optimization_output = run_xtb(binary_path, command, working_directory)
            if "GEOMETRY OPTIMIZATION CONVERGED" in optimization_output or not optimized_path.exists():
                break
            geometry_path.write_text(optimized_path.read_text(encoding="utf-8"), encoding="utf-8")
        if not optimized_path.exists() or "GEOMETRY OPTIMIZATION CONVERGED" not in optimization_output:
            raise RuntimeError("xTB geometry optimization did not converge")

        optimized_molblock = molblock_from_xyz(molecule, optimized_path)
        optimized_xyz = optimized_path.read_text(encoding="utf-8")
        geometry_path.write_text(optimized_xyz, encoding="utf-8")
        command["calculation"] = ["--sp", "--json", *solvent_arguments]
        single_point_output = run_xtb(binary_path, command, working_directory)
        solvent_properties = parse_single_point(
            single_point_output, molecule, read_xtb_json(working_directory / "xtbout.json"), optimized_xyz
        )
        record = {"seed_id": conformer_id, "optimization_attempts": attempt,
                  "molblock": optimized_molblock, "_optimized_xyz": optimized_xyz}
        if solvent == "gas":
            record["gas"] = solvent_properties
        else:
            record["alpb"] = solvent_properties
        return record


def calculate_qrrho(binary_path, molecule, conformer, solvent):
    solvent_arguments = [] if solvent == "gas" else ["--alpb", solvent]
    command = {
        "charge": Chem.GetFormalCharge(molecule),
        "calculation": ["--hess", *solvent_arguments, "--input", "thermo.inp"],
    }
    with tempfile.TemporaryDirectory() as temporary:
        working_directory = Path(temporary)
        (working_directory / "thermo.inp").write_text(
            f"$thermo\n temp={TEMPERATURE_K}\n sthr=50.0\n$end\n", encoding="utf-8"
        )
        (working_directory / "mol.xyz").write_text(conformer["_optimized_xyz"], encoding="utf-8")
        output = run_xtb(binary_path, command, working_directory)
    return parse_qrrho_free_energy(output)


def prepare_molecule(molblock):
    molecule = Chem.MolFromMolBlock(molblock, removeHs=False, sanitize=False)
    carbon_hydrogens = [atom.GetIdx() for atom in molecule.GetAtoms()
                        if atom.GetAtomicNum() == 1 and
                        any(neighbor.GetAtomicNum() == 6 for neighbor in atom.GetNeighbors())]
    editable = Chem.RWMol(molecule)
    for atom_index in sorted(carbon_hydrogens, reverse=True):
        editable.RemoveAtom(atom_index)
    return editable.GetMol()


def clustering_coordinates(molblock):
    """Coordinates without carbon-bound hydrogens, centred at the origin."""
    coordinates = np.asarray(prepare_molecule(molblock).GetConformer().GetPositions())
    return coordinates - coordinates.mean(axis=0)


def symmetry_mappings(molblock):
    """All atom permutations that map the clustering graph onto itself.

    The automorphisms come from spyrmsd's graph isomorphism matching (atoms labelled by
    atomic number), so symmetry-equivalent atoms such as ring flips or CF3 fluorines are
    matched. RDKit's GetBestRMS does not apply these mappings to the unsanitized,
    partly hydrogen-stripped molecules used here and falls back to the identity order.
    """
    molecule = prepare_molecule(molblock)
    atomic_numbers = np.array([atom.GetAtomicNum() for atom in molecule.GetAtoms()])
    graph = spyrmsd_graph.graph_from_adjacency_matrix(Chem.GetAdjacencyMatrix(molecule), atomic_numbers)
    mappings = []
    for reference_order, mapped_order in spyrmsd_graph.match_graphs(graph, graph):
        permutation = np.empty(len(atomic_numbers), dtype=np.intp)
        permutation[np.asarray(reference_order)] = np.asarray(mapped_order)
        mappings.append(permutation)
    return np.array(mappings)


def symmetric_rmsd(reference, coordinates, mappings, allow_reflection=False):
    """Lowest Kabsch RMSD over all symmetry mappings, evaluated in batches.

    With allow_reflection the superposition may also invert the structure, so a conformer
    and its mirror image have RMSD 0; otherwise only proper rotations are used.
    """
    reference_norm = np.sum(reference * reference)
    coordinates_norm = np.sum(coordinates * coordinates)
    best = np.inf
    for start in range(0, len(mappings), RMSD_CHUNK):
        permuted = coordinates[mappings[start:start + RMSD_CHUNK]]
        covariance = np.einsum("ni,knj->kij", reference, permuted)
        singular = np.linalg.svd(covariance, compute_uv=False)
        if not allow_reflection:
            singular[:, -1] *= np.sign(np.linalg.det(covariance))
        squared = (reference_norm + coordinates_norm - 2 * singular.sum(axis=1)) / len(reference)
        best = min(best, float(np.sqrt(max(squared.min(), 0.0))))
    return best


def mirror_rmsd(coordinates, mappings):
    """RMSD between a conformer and its own mirror image (proper rotations only)."""
    return symmetric_rmsd(coordinates, coordinates * np.array([1.0, 1.0, -1.0]), mappings)


def conformer_groups(records, mappings):
    """Energy-ordered leader clustering on dipole difference and symmetry-corrected RMSD.

    The RMSD allows reflections, so mirror-image conformers always fall into one group;
    the mirror image is restored afterwards as an explicit twin (see mirror_twin). All
    conformers share one topology, so the symmetry mappings are computed once per molecule.
    """
    records = sorted(records, key=lambda record: record["energy_hartree"])
    coordinates = [clustering_coordinates(record["molblock"]) for record in records]
    groups, group_coordinates = [], []
    for index, record in enumerate(records):
        duplicate_group = None
        for group, reference in zip(groups, group_coordinates):
            if abs(record["dipole_debye"] - group[0]["dipole_debye"]) >= POSTOPT_MU_D:
                continue
            if symmetric_rmsd(reference, coordinates[index], mappings, allow_reflection=True) < POSTOPT_RMS_A:
                duplicate_group = group
                break

        if duplicate_group is None:
            groups.append([record])
            group_coordinates.append(coordinates[index])
        else:
            duplicate_group.append(record)
    return groups


def mirror_twin(conformer, property_name):
    """Mirror image (z -> -z) of a finished conformer record; no new xTB calculation is needed.

    Energies, free energy, charges, polarizabilities and |mu| are identical for a mirror
    image; the z components of the dipole vector and of the atomic dipoles change sign.
    """
    twin = copy.deepcopy(conformer)
    lines = conformer["molblock"].split("\n")
    if "V2000" not in lines[3]:
        raise ValueError("mirror_twin expects a V2000 molblock")
    for index in range(4, 4 + int(lines[3][:3])):
        line = lines[index]
        lines[index] = f"{line[:20]}{-float(line[20:30]):10.4f}{line[30:]}"
    twin["molblock"] = "\n".join(lines)
    properties = twin[property_name]
    properties["dipole_vector_debye"][2] = -properties["dipole_vector_debye"][2]
    for dipole in properties["atomic_dipoles_debye"]:
        dipole[2] = -dipole[2]
    twin["mirror_of"] = conformer["seed_id"]
    return twin


def stereo_elements(molecule):
    """Tetrahedral centres and double bonds as (kind, index, atoms, specified)."""
    elements = []
    for info in Chem.FindPotentialStereo(Chem.Mol(molecule)):
        specified = info.specified == Chem.StereoSpecified.Specified
        if info.type == Chem.StereoType.Atom_Tetrahedral:
            elements.append(("atom", info.centeredOn, (info.centeredOn,), specified))
        elif info.type == Chem.StereoType.Bond_Double:
            bond = molecule.GetBondWithIdx(info.centeredOn)
            begin, end = bond.GetBeginAtom(), bond.GetEndAtom()
            # A double bond that is adjacent to another multiple bond is part of a
            # cumulene-like system (for example aryl-N=C=O).  It does not have an
            # independently assignable E/Z configuration.  RDKit may nevertheless
            # list its first double bond as "potential stereo"; assigning stereo
            # atoms to it makes canonicalization fail with "could not find atom2".
            def has_other_multiple_bond(atom):
                return any(
                    neighbour_bond.GetIdx() != bond.GetIdx()
                    and neighbour_bond.GetBondType() != Chem.BondType.SINGLE
                    for neighbour_bond in atom.GetBonds()
                )
            if has_other_multiple_bond(begin) or has_other_multiple_bond(end):
                continue
            # E/Z can only be set on the heavy-atom graph when both ends carry another heavy
            # atom; a bond like C=[NH] is left to the embedding.
            if begin.GetDegree() < 2 or end.GetDegree() < 2:
                continue
            elements.append(("bond", info.centeredOn, (begin.GetIdx(), end.GetIdx()), specified))
    return elements


def ring_systems(molecule):
    systems = []
    for ring in molecule.GetRingInfo().AtomRings():
        system = set(ring)
        for other in [other for other in systems if other & system]:
            system |= other
            systems.remove(other)
        systems.append(system)
    return systems


def stereo_groups(molecule, elements):
    """Unspecified stereo elements whose relative configuration can change the dipole.

    Elements at most STEREO_COUPLING_BONDS apart or in one ring system form a group. A group
    is kept when it has a relative configuration (two or more elements, or a double bond)
    and a polar atom lies within STEREO_POLAR_BONDS of it. A lone tetrahedral centre only
    has enantiomers, which the mirror twins already cover.
    """
    distance = Chem.GetDistanceMatrix(molecule)
    systems = ring_systems(molecule)

    def coupled(first, second):
        if min(distance[i, j] for i in first[2] for j in second[2]) <= STEREO_COUPLING_BONDS:
            return True
        return any(system.intersection(first[2]) and system.intersection(second[2]) for system in systems)

    groups = []
    for element in [element for element in elements if not element[3]]:
        linked = [group for group in groups if any(coupled(element, other) for other in group)]
        groups = [group for group in groups if group not in linked]
        groups.append([element] + [other for group in linked for other in group])
    polar = [atom.GetIdx() for atom in molecule.GetAtoms() if atom.GetAtomicNum() in POLAR_Z]

    def changes_dipole(group):
        has_relative = len(group) > 1 or group[0][0] == "bond"
        return has_relative and any(distance[i, j] <= STEREO_POLAR_BONDS
                                    for element in group for i in element[2] for j in polar)

    return [sorted(group, key=lambda element: (element[0], element[1])) for group in groups if changes_dipole(group)]


def invert_centres(molecule, atom_indices):
    swap = {Chem.ChiralType.CHI_TETRAHEDRAL_CW: Chem.ChiralType.CHI_TETRAHEDRAL_CCW,
            Chem.ChiralType.CHI_TETRAHEDRAL_CCW: Chem.ChiralType.CHI_TETRAHEDRAL_CW}
    for index in atom_indices:
        atom = molecule.GetAtomWithIdx(index)
        atom.SetChiralTag(swap.get(atom.GetChiralTag(), atom.GetChiralTag()))


def assign_stereo(molecule, elements, states):
    """Copy of the molecule (same atom order) with the given configuration of each element."""
    assigned = Chem.Mol(molecule)
    for (kind, index, _atoms, _specified), state in zip(elements, states):
        if kind == "atom":
            assigned.GetAtomWithIdx(index).SetChiralTag(
                Chem.ChiralType.CHI_TETRAHEDRAL_CW if state else Chem.ChiralType.CHI_TETRAHEDRAL_CCW)
        else:
            bond = assigned.GetBondWithIdx(index)
            begin, end = bond.GetBeginAtom(), bond.GetEndAtom()
            left = next(n.GetIdx() for n in begin.GetNeighbors() if n.GetIdx() != end.GetIdx())
            right = next(n.GetIdx() for n in end.GetNeighbors() if n.GetIdx() != begin.GetIdx())
            bond.SetStereoAtoms(left, right)
            bond.SetStereo(Chem.BondStereo.STEREOTRANS if state else Chem.BondStereo.STEREOCIS)
    return assigned


def diastereomer_key(assigned, inversion_sets):
    """Canonical label shared by all assignments that differ only by inverting whole groups."""
    labels = []
    for flips in product((False, True), repeat=len(inversion_sets)):
        variant = Chem.Mol(assigned)
        for flip, atoms in zip(flips, inversion_sets):
            if flip:
                invert_centres(variant, atoms)
        labels.append(Chem.CanonSmiles(Chem.MolToSmiles(variant)))
    return min(labels)


def enumerate_diastereomers(molecule):
    """Heavy-atom molecules (input atom order), one per diastereomer that can change the dipole.

    Configurations outside the enumerated groups stay unspecified and are sampled at random
    by the embedding, as before. Inverting all centres of a group gives the same diastereomer
    (its enantiomer is covered by the mirror twins), unless the molecule also has specified
    stereo, in which case only inversion of the whole molecule is treated as equivalent.
    """
    elements = stereo_elements(molecule)
    groups = stereo_groups(molecule, elements)
    selected = [element for group in groups for element in group]
    info = {"enumerated_stereo_elements": [f"{kind}{index}" for kind, index, _, _ in selected],
            "stereo_fallback": None}
    if not selected:
        return [molecule], info
    if len(selected) > MAX_ENUMERATED_ELEMENTS:
        info["stereo_fallback"] = (f"{len(selected)} stereo elements exceed {MAX_ENUMERATED_ELEMENTS}; "
                                   "configurations sampled at random")
        return [molecule], info
    specified_atoms = [index for kind, index, _, specified in elements if specified and kind == "atom"]
    if any(element[3] for element in elements):
        inversion_sets = [[index for kind, index, _, _ in selected if kind == "atom"] + specified_atoms]
    else:
        inversion_sets = [[index for kind, index, _, _ in group if kind == "atom"] for group in groups]
    diastereomers = {}
    for states in product((False, True), repeat=len(selected)):
        assigned = assign_stereo(molecule, selected, states)
        diastereomers.setdefault(diastereomer_key(assigned, inversion_sets), assigned)
    return list(diastereomers.values()), info


def embed_conformers(molecule, n_conformers):
    embedding = AllChem.ETKDGv3()
    embedding.randomSeed = 0xF00D
    embedding.useExpTorsionAnglePrefs = False
    embedding.useBasicKnowledge = True
    embedding.pruneRmsThresh = EMBED_PRUNE_RMS_A
    embedding.onlyHeavyAtomsForRMS = True
    embedding.useSymmetryForPruning = True
    conformer_ids = list(AllChem.EmbedMultipleConfs(molecule, numConfs=n_conformers, params=embedding))
    if not conformer_ids:
        embedding.useRandomCoords = True
        conformer_ids = list(AllChem.EmbedMultipleConfs(molecule, numConfs=n_conformers, params=embedding))
    return conformer_ids


def can_embed(molecule):
    """Quick test with one structure: geometrically impossible diastereomers (e.g. a
    trans-fused bridged bicycle) fail here instead of in the full multi-conformer run."""
    for random_coordinates in (False, True):
        embedding = AllChem.ETKDGv3()
        embedding.randomSeed = 0xF00D
        embedding.useRandomCoords = random_coordinates
        if AllChem.EmbedMolecule(Chem.Mol(molecule), embedding) == 0:
            return True
    return False


def calculate_qrrho_for_group(binary_path, molecule, conformers, solvent):
    property_name = "gas" if solvent == "gas" else "alpb"
    failures = {}
    for conformer in conformers:
        try:
            conformer[property_name]["qrrho_free_energy_hartree"] = calculate_qrrho(
                binary_path, molecule, conformer, solvent
            )
            return conformer, failures
        except Exception as error:
            failures[str(conformer["seed_id"])] = repr(error)
    return None, failures


# ---------------------------------------------------------------- xtb executable
# Search order: given path / command, PATH, earlier download, latest GitHub release.


XTB_RELEASE_API = "https://api.github.com/repos/grimme-lab/xtb/releases/latest"


# Downloaded xtb releases are kept next to the installed package (like the cache database of MoleculeResolver);
# if that folder is not writable, in the user's home. CRAFT_EPSILON_XTB_DIR overrides both.
PACKAGE_XTB_DIR = Path(__file__).resolve().parent / "xtb"


USER_XTB_DIR = Path.home() / ".craft_epsilon" / "xtb"


XTB_NAMES = ("xtb.exe", "xtb") if platform.system() == "Windows" else ("xtb",)


OS_ALIASES = {"linux": ["linux"], "windows": ["windows"], "darwin": ["macos", "darwin", "osx"]}


ARCH_ALIASES = {"x86_64": ["x86_64", "amd64"], "amd64": ["x86_64", "amd64"],
                "aarch64": ["aarch64", "arm64"], "arm64": ["aarch64", "arm64"]}


logger = logging.getLogger("craft_epsilon")


def xtb_cache_dirs():
    """Folders searched for downloaded xtb releases, in order: CRAFT_EPSILON_XTB_DIR if set, otherwise the
    xtb folder inside the installed package and the user's ~/.craft_epsilon/xtb."""
    if os.environ.get("CRAFT_EPSILON_XTB_DIR"):
        return [Path(os.environ["CRAFT_EPSILON_XTB_DIR"])]
    return [PACKAGE_XTB_DIR, USER_XTB_DIR]


def writable_xtb_cache_dir():
    """First folder of xtb_cache_dirs() that can be created and written to."""
    for directory in xtb_cache_dirs():
        try:
            directory.mkdir(parents=True, exist_ok=True)
            probe = directory / ".write_test"
            probe.write_text("", encoding="utf-8")
            probe.unlink()
            return directory
        except OSError:
            continue
    raise PermissionError(f"no writable folder for the xtb download among {[str(d) for d in xtb_cache_dirs()]}; "
                          "set CRAFT_EPSILON_XTB_DIR")


def cached_xtb_binary(directories=None):
    """xtb executable inside a previously downloaded release (newest release folder first), or None."""
    for directory in directories or xtb_cache_dirs():
        if not Path(directory).is_dir():
            continue
        for release in sorted((p for p in Path(directory).iterdir() if p.is_dir()), reverse=True):
            for name in XTB_NAMES:
                for candidate in sorted(release.rglob(name)):
                    if candidate.is_file() and candidate.parent.name == "bin":
                        return str(candidate)
    return None


def select_xtb_asset(release, os_name, arch):
    """Binary archive of a GitHub release (dict of the releases API) for an operating system and architecture
    (platform.system().lower(), platform.machine().lower()), or None if the release has none."""
    return next((asset for asset in release["assets"]
                 if any(x in asset["name"].lower() for x in OS_ALIASES.get(os_name, []))
                 and any(x in asset["name"].lower() for x in ARCH_ALIASES.get(arch, [arch]))
                 and asset["name"].endswith((".tar.xz", ".tar.gz", ".zip"))), None)


def download_xtb(directory=None):
    """Download the latest xtb release (grimme-lab/xtb on GitHub) for this operating system and architecture,
    check it against the published SHA-256 checksum, unpack it into directory/<tag> and return the executable.
    Default directory: the xtb folder of the installed package, else ~/.craft_epsilon/xtb (writable_xtb_cache_dir)."""
    directory = Path(directory) if directory else writable_xtb_cache_dir()
    os_name, arch = platform.system().lower(), platform.machine().lower()
    with urlopen(XTB_RELEASE_API, timeout=60) as response:
        release = json.load(response)
    assets = {asset["name"]: asset for asset in release["assets"]}
    asset = select_xtb_asset(release, os_name, arch)
    if asset is None:
        raise FileNotFoundError(f"xtb was not found and release {release['tag_name']} has no binary for "
                                f"{os_name}/{arch}; install xtb yourself (e.g. conda install -c conda-forge xtb) "
                                "or pass xtb_binary")
    target = Path(directory) / release["tag_name"]
    logger.info("xtb not found; downloading %s from github.com/grimme-lab/xtb to %s", asset["name"], target)
    with tempfile.TemporaryDirectory() as temporary:
        archive = Path(temporary) / asset["name"]
        urlretrieve(asset["browser_download_url"], archive)
        checksum = assets.get(asset["name"] + ".sha256")
        if checksum is not None:
            with urlopen(checksum["browser_download_url"], timeout=60) as response:
                expected = response.read().decode().split()[0].lower()
            actual = hashlib.sha256(archive.read_bytes()).hexdigest()
            if actual != expected:
                raise RuntimeError(f"checksum mismatch for {asset['name']}: expected {expected}, got {actual}")
        target.mkdir(parents=True, exist_ok=True)
        shutil.unpack_archive(archive, target)
    binary = cached_xtb_binary([directory])
    if binary is None:
        raise FileNotFoundError(f"no xtb executable found in the downloaded release {asset['name']}")
    if os.name != "nt":
        os.chmod(binary, os.stat(binary).st_mode | 0o111)
    logger.info("using xtb %s: %s", release["tag_name"], binary)
    return binary


def find_xtb_binary(xtb_binary=None, download=True):
    """Path of the xtb executable: the given path or command, otherwise xtb on the PATH (shutil.which), otherwise
    a release downloaded earlier (xtb_cache_dirs), otherwise (download=True) the latest release from GitHub."""
    if xtb_binary:
        path = shutil.which(str(xtb_binary)) or (str(Path(xtb_binary).resolve()) if Path(xtb_binary).is_file() else None)
        if path is None:
            raise FileNotFoundError(f"xtb executable not found: {xtb_binary}")
        return path
    for name in ("xtb", "xtb.exe", "otool_xtb", "otool_xtb.exe"):
        path = shutil.which(name)
        if path:
            return path
    path = cached_xtb_binary()
    if path:
        return path
    if not download:
        raise FileNotFoundError("xtb executable not found on the PATH; install xtb (6.7.1) or pass xtb_binary")
    return download_xtb()


def generate_conformers(smiles, binary_path, solvent=SOLVENT, workers=N_WORKERS):
    """GFN2-xTB/ALPB conformer ensemble of a molecule (the conformer protocol used to fit the model).

    Returns a record with the unique conformers (molblock, energies, dipoles, partial charges, atomic dipoles,
    atomic polarizabilities, qRRHO free energy), their diastereomer ensemble and the protocol counts."""
    started_at = time.perf_counter()
    heavy_molecule = Chem.MolFromSmiles(smiles)
    property_name = "gas" if solvent == "gas" else "alpb"
    # Diastereomers whose relative configuration can change the dipole get separate
    # ensembles; the seeds are shared between them so the total effort stays N_SEEDS.
    diastereomers, stereo_info = enumerate_diastereomers(heavy_molecule)
    if len(diastereomers) > 1:
        possible = [diastereomer for diastereomer in diastereomers if can_embed(Chem.AddHs(diastereomer))]
        diastereomers = possible or [heavy_molecule]
    seeds_per_diastereomer = math.ceil(N_SEEDS / len(diastereomers))
    molecules, jobs, n_polar_h_variants = [], [], 0
    for index, diastereomer in enumerate(diastereomers):
        molecule = Chem.AddHs(diastereomer)
        conformer_ids = embed_conformers(molecule, seeds_per_diastereomer)
        added_polar_h_ids = amplify_polar_h_conformers(molecule, seed=0xF00D + index)
        n_polar_h_variants += len(added_polar_h_ids)
        molecules.append(molecule)
        jobs += [(index, conformer_id) for conformer_id in conformer_ids + added_polar_h_ids]
    if not jobs:
        raise RuntimeError("RDKit could not generate conformers")

    conformers, failures = [], {}
    with ThreadPoolExecutor(max_workers=min(len(jobs), workers)) as pool:
        futures = {pool.submit(calculate_conformer, binary_path, molecules[index], conformer_id, solvent):
                   (job_index, index) for job_index, (index, conformer_id) in enumerate(jobs)}
        for future in as_completed(futures):
            seed_id, index = futures[future]
            try:
                record = future.result()
                record["seed_id"] = seed_id
                record["diastereomer"] = index
                conformers.append(record)
            except Exception as error:
                failures[str(seed_id)] = repr(error)
    if not conformers:
        raise RuntimeError(f"No conformers completed successfully: {failures}")
    n_optimized = len(conformers)

    mappings = symmetry_mappings(conformers[0]["molblock"])
    conformers_by_seed = {record["seed_id"]: record for record in conformers}
    qrrho_groups = []
    for index in range(len(diastereomers)):
        records_for_clustering = [
            {"seed_id": record["seed_id"], "molblock": record["molblock"],
             "energy_hartree": record[property_name]["energy_hartree"],
             "dipole_debye": record[property_name]["dipole_debye"]}
            for record in conformers if record["diastereomer"] == index]
        if records_for_clustering:
            qrrho_groups += [[conformers_by_seed[record["seed_id"]] for record in cluster]
                             for cluster in conformer_groups(records_for_clustering, mappings)]
    unique_conformers = []
    with ThreadPoolExecutor(max_workers=min(len(qrrho_groups), workers)) as pool:
        futures = [pool.submit(calculate_qrrho_for_group, binary_path, molecules[group[0]["diastereomer"]], group,
                               solvent) for group in qrrho_groups]
        for future in as_completed(futures):
            try:
                conformer, group_failures = future.result()
                failures.update({f"qrrho_{key}": value for key, value in group_failures.items()})
                if conformer is not None:
                    unique_conformers.append(conformer)
            except Exception as error:
                failures["qrrho_group"] = repr(error)
    if not unique_conformers:
        raise RuntimeError(f"No unique conformers completed qRRHO successfully: {failures}")

    # Every chiral conformer is stored together with its mirror image, every achiral one once.
    conformers = []
    for conformer in sorted(unique_conformers, key=lambda record: (
            record["diastereomer"], record[property_name]["qrrho_free_energy_hartree"])):
        conformer.pop("_optimized_xyz", None)
        distance = mirror_rmsd(clustering_coordinates(conformer["molblock"]), mappings)
        conformer["mirror_rmsd_a"] = round(distance, 4)
        conformers.append(conformer)
        if distance >= MIRROR_DISTINCT_RMS_A:
            conformers.append(mirror_twin(conformer, property_name))
    return {"smiles": smiles, "solvent": solvent, "conformers": conformers, "failures": failures,
            "n_diastereomers": len(diastereomers), "n_embedded": len(jobs), "n_polar_h_variants": n_polar_h_variants,
            "n_optimized": n_optimized, "n_unique": len(unique_conformers),
            "stereo_fallback": stereo_info.get("stereo_fallback"),
            "elapsed_seconds": time.perf_counter() - started_at}
