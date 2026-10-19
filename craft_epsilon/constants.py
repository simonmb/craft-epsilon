"""Physical constants and the settings of the conformer workflow used to fit CRAFT-epsilon."""
import re


SOLVENT = "ether"                     # ALPB solvent of the conformer ensembles (diethyl ether)


TEMPERATURE_K = 298.15                # conformer ensembles, Boltzmann weights, molar volume and prediction


THERMAL_EXPANSION_PER_K = 1.0e-3      # alpha_V (K^-1), fixed (not fitted): V_m(T) = V_m(298.15 K) (1 + alpha_V (T - 298.15 K))


N_WORKERS = 8                         # parallel single-threaded xtb processes


N_A = 6.02214076e23                   # mol^-1 (CODATA 2018, exact)


EPS0 = 8.8541878128e-12               # F m^-1


KB = 1.380649e-23                     # J K^-1 (exact)


DEBYE = 1e-21 / 299792458.0           # C m per debye


R_KCAL = 0.0019872041                 # kcal mol^-1 K^-1


HARTREE_TO_KCAL = 627.5094740631


KAPPA_D = 4.80320471257               # D per (e angstrom)


AU_TO_A3 = 0.14818471                 # bohr^3 -> angstrom^3


AU_TO_DEBYE = 2.541746473             # D per (e bohr)


AU_DIPOLE_TO_DEBYE = AU_TO_DEBYE


BOHR_PER_ANGSTROM = 1 / 0.529177210903


DIPOLE_RECONSTRUCTION_TOL_AU = 0.002  # q*r + atomic dipoles vs printed xTB dipole (3 decimals)


X_MAX = 0.99                          # keeps (1 + 2x)/(1 - x) finite


N_SEEDS = 100


EMBED_PRUNE_RMS_A = 0.125


POSTOPT_RMS_A = 0.5


POSTOPT_MU_D = 0.2


RMSD_CHUNK = 20000                    # symmetry mappings evaluated per numpy batch


MIRROR_DISTINCT_RMS_A = 0.1           # a conformer closer than this to its own mirror image is achiral


STEREO_COUPLING_BONDS = 3             # stereo elements this close (or in one ring system) share a relative configuration


STEREO_POLAR_BONDS = 2                # a stereo group is enumerated only with a polar atom this close


MAX_ENUMERATED_ELEMENTS = 10          # more selected stereo elements than this: sample configurations at random


POLAR_Z = {7, 8, 9, 15, 16, 17, 35, 53}


ENERGY_PATTERN = re.compile(r"TOTAL ENERGY\s+(-?\d+\.\d+)\s*Eh")


DIPOLE_PATTERN = re.compile(r"molecular dipole:.*?\n\s*full:\s+(\S+)\s+(\S+)\s+(\S+)\s+(\S+)", re.DOTALL)


POLARIZABILITY_HEADER = re.compile(r"#\s+Z\s+covCN.*(?:α|alpha)\(0\)")


QRRHO_FREE_ENERGY = re.compile(r"::\s*total free energy\s+([-+\d.Ee]+)\s*Eh\s*::", re.IGNORECASE)


ATOM_TYPES = ["B_SP2_A", "Br", "C_SP2_A", "C_SP2_R", "C_SP3_A", "C_SP3_R", "C_SP_A", "Cl", "F", "H", "I",
              "N_SP2_A", "N_SP2_R", "N_SP3_A", "N_SP3_R", "N_SP_A", "O_SP2_A", "O_SP2_R", "O_SP3_A",
              "O_SP3_R", "P_SP3_A", "S_SP2_A", "S_SP2_R", "S_SP3_A", "S_SP3_R", "Si_SP3_A"]
