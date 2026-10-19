"""CRAFT-epsilon: A Conformer-Resolved Atom- and Fragment-based Theory for the Static Relative Permittivity.

CRAFT-epsilon (S. Müller, "CRAFT-epsilon: A Conformer-Resolved Atom- and Fragment-based Theory for the Static
Relative Permittivity of Pure Organic Liquids") predicts the static relative permittivity of a neutral pure organic
liquid at 298.15 K from its SMILES string or name. All 41 fitted parameters are defined in
craft_epsilon.parameters; they were fitted to the static relative permittivities of 1393 pure liquids with
GFN2-xTB/ALPB(diethyl ether) conformer ensembles.

Usage
    command line:  craft-epsilon "<SMILES or name>"           prints all relevant intermediate quantities
                   python -m craft_epsilon "<SMILES or name>"
    Python:        from craft_epsilon import epsilon, predict
                   epsilon("CCO")                              -> float, prints nothing
                   predict("ethanol")                          -> dict with all intermediate quantities

    A valid SMILES string is used directly. Anything else is treated as a name and resolved to a structure
    with MoleculeResolver (find_single_molecule_crosschecked, neutral structure, cross-checked across
    several chemical databases) [7]; this needs an internet connection.

Requirements
    Python >= 3.10, RDKit, NumPy, spyrmsd, molecule-resolver (all installed by pip install craft-epsilon);
    the xtb program (version 6.7.1 was used to fit the model) is taken from xtb_binary, from the PATH (found
    with shutil.which) or from an earlier download; if none is available, the latest xtb release for this
    operating system and architecture is downloaded once from github.com/grimme-lab/xtb (checked against its
    SHA-256 checksum) into the folder xtb inside the installed package (like the cache database of
    MoleculeResolver); if that folder is not writable, into ~/.craft_epsilon/xtb. The environment variable
    CRAFT_EPSILON_XTB_DIR overrides both.

Workflow (from the SMILES to epsilon_r)
    1. Conformers (craft_epsilon.conformers): up to 100 RDKit ETKDGv3 conformers [4] (heavy-atom RMSD pruning
       0.125 A); unspecified stereo elements whose relative configuration can change the dipole are enumerated
       as separate diastereomer ensembles; molecules with >= 2 rotatable OH/SH/NH2 groups get extra polar-H
       variants.
    2. GFN2-xTB [2] in the ALPB continuum [3] of diethyl ether: tight optimization, single point (dipole,
       partial charges, atomic dipoles, D4 atomic polarizabilities [6]), leader clustering (0.2 D, symmetry-
       corrected RMSD 0.5 A [5]) and a numerical Hessian per unique conformer for the qRRHO free energy [8].
       Chiral conformers are stored together with their mirror image.
    3. Polar groups (craft_epsilon.groups) split every conformer dipole exactly into group vectors and a
       remainder; mu_eff = b * mu_rest + sum_k c_k mu_k, Boltzmann mean square over the free-energy-weighted
       ensemble (craft_epsilon.descriptors).
    4. eps_inf from Clausius-Mossotti with atom-group scaled, Boltzmann-averaged atomic polarizabilities.
    5. Molar volume from the atom-type and ring group contributions of Bouteloup and Mathieu [1]
       (craft_epsilon.molar_volume).
    6. epsilon_r from the Onsager-Kirkwood-Froehlich relation (craft_epsilon.model)
           (2 eps + eps_inf)(eps - eps_inf) / (eps (eps_inf + 2)^2) = N_A <mu_eff^2> / (9 eps0 kB T V_m).

    Conformer ensembles depend slightly on the RDKit version and platform (the embedding of a given random seed
    can differ); the fitted parameters refer to the ensembles used in the paper.

References
    [1] R. Bouteloup, D. Mathieu, Predicting dielectric constants of pure liquids: fragment-based
        Kirkwood-Froehlich model applicable over a wide range of polarity, Phys. Chem. Chem. Phys. 21 (2019)
        11043-11057. https://doi.org/10.1039/C9CP01704F  (molar-volume group contributions used here)
    [2] C. Bannwarth, S. Ehlert, S. Grimme, J. Chem. Theory Comput. 15 (2019) 1652-1671.
        https://doi.org/10.1021/acs.jctc.8b01176  (GFN2-xTB)
    [3] S. Ehlert, M. Stahn, S. Spicher, S. Grimme, J. Chem. Theory Comput. 17 (2021) 4250-4261.
        https://doi.org/10.1021/acs.jctc.1c00471  (ALPB implicit solvation)
    [4] RDKit, https://www.rdkit.org; S. Wang, J. Witek, G. A. Landrum, S. Riniker, J. Chem. Inf. Model. 60
        (2020) 2044-2058. https://doi.org/10.1021/acs.jcim.0c00025  (ETKDGv3)
    [5] R. Meli, P. C. Biggin, J. Cheminform. 12 (2020) 49. https://doi.org/10.1186/s13321-020-00455-2  (spyrmsd)
    [6] E. Caldeweyher et al., J. Chem. Phys. 150 (2019) 154122. https://doi.org/10.1063/1.5090222  (D4)
    [7] S. Müller, J. Cheminform. 17 (2025). https://doi.org/10.1186/s13321-025-01064-7  (MoleculeResolver)
    [8] S. Grimme, Chem. Eur. J. 18 (2012) 9955-9964. https://doi.org/10.1002/chem.201200497  (qRRHO)
"""
from rdkit import RDLogger

RDLogger.DisableLog("rdApp.*")

from .model import epsilon, epsilon_from_record, predict  # noqa: E402
from .parameters import PARAMETERS  # noqa: E402

__version__ = "1.0.0"
__license__ = "MIT"
__all__ = ["epsilon", "predict", "epsilon_from_record", "PARAMETERS", "__version__"]
