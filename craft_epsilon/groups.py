"""Polar groups (prioritized SMARTS rules), chemical families and atom types of CRAFT-epsilon."""
from dataclasses import dataclass
from rdkit import Chem
from .constants import POLAR_Z


PATTERNS = {
    "ketone": "[CX3](=[OX1])", "aldehyde": "[CX3H1](=[OX1])", "nitrile": "[CX2]#[NX1]",
    "furan": "[o;r5]", "alcohol": "[OX2H]", "phenol": "[OX2H][c]", "ether": "[OD2]([#6])[#6]",
    "ester": "[CX3](=[OX1])[OX2]", "carboxylic_acid": "[CX3](=[OX1])[OX2H]",
    "amide": "[NX3][CX3](=[OX1])", "amine": "[NX3]", "nitro": "[N+](=O)[O-]",
    "phosphoryl_P=O": "[#15]=[OX1]", "phosphate_phosphonate": "[#15](=[OX1])[OX2]",
    "sulfoxide": "[#16X3]=[OX1]", "sulfone": "[#16X4](=[OX1])=[OX1]",
    "thioether": "[#16X2]", "thiol": "[#16X2H]", "azole": "[n;r5]",
    "aromatic_N_heterocycle": "[n]", "aromatic_N_oxide": "[n+][O-]",
    "isocyanate": "[NX2]=[CX2]=[OX1]", "thiocyanate": "[$([#16X2][CX2]#[NX1]),$([NX2]=[CX2]=[SX1])]",
    "urea": "[NX3][CX3](=[OX1])[NX3]", "carbonate/carbamate": "[OX2][CX3](=[OX1])[OX2,NX3]",
    "nitrate_ester": "[OX2][N+](=O)[O-]", "nitrite_ester": "[OX2][NX2]=[OX1]",
    "acyl_halide": "[CX3](=[OX1])[F,Cl,Br,I]", "epoxide": "[O;r3]1CC1",
    "fluoro": "[F]", "chloro": "[Cl]", "bromo": "[Br]", "iodo": "[I]",
    "anhydride": "[CX3](=[OX1])[OX2][CX3]=[OX1]", "sulfate_ester": "[OX2][SX4](=[OX1])(=[OX1])[OX2]",
    "siloxane/silyl_ether": "[#14][OX2]", "borate_B-O": "[#5][OX2]",
}
PATTERNS = {name: Chem.MolFromSmarts(value) for name, value in PATTERNS.items()}
PRIORITY = ["urea", "carbonate/carbamate", "anhydride", "carboxylic_acid", "ester", "amide",
            "acyl_halide", "nitrate_ester", "nitrite_ester", "nitro", "phosphate_phosphonate",
            "sulfate_ester", "sulfone", "sulfoxide", "aromatic_N_oxide", "isocyanate",
            "thiocyanate", "nitrile", "aldehyde", "ketone", "phenol", "epoxide", "furan",
            "alcohol", "ether", "thiol", "thioether", "azole", "amine"]
FORMAMIDE = Chem.MolFromSmarts("[NX3][CX3H1]=[OX1]")
ACETAMIDE = Chem.MolFromSmarts("[NX3][CX3](=[OX1])[CH3]")
HB_PATTERNS = {name: Chem.MolFromSmarts(value) for name, value in {
    "carboxylic_acid": "[CX3](=O)[OX2H1]", "phosphate_acid": "[PX4](=O)[OX2H1]",
    "amide_primary": "[CX3](=[OX1])[NX3H2]", "amide_secondary": "[CX3](=[OX1])[NX3H1][#6]",
    "phenol": "[OX2H][c]", "alcohol_primary": "[CX4H2;!$(C(O)(O))][OX2H1]",
    "alcohol_secondary": "[CX4H1;!$(C(O)(O))][OX2H1]", "alcohol_tertiary": "[CX4H0;!$(C(O)(O))][OX2H1]",
    "thiol": "[#6][SX2H1]", "amine_primary": "[NX3H2;!$(NC=O)]",
    "amine_secondary": "[NX3H1;!$(NC=O)]([!#1])[!#1]"}.items()}


@dataclass(frozen=True)
class Fragment:
    label: str
    atoms: tuple
    donor_labels: tuple


def fragment_context(mol):
    atoms = {a.GetIdx() for a in mol.GetAtoms() if a.GetAtomicNum() in POLAR_Z}
    atoms.update(n.GetIdx() for i in tuple(atoms) for n in mol.GetAtomWithIdx(i).GetNeighbors() if n.GetAtomicNum() == 6)
    selected = [mol.GetAtomWithIdx(i) for i in atoms]
    if any(a.GetIsAromatic() or any(n.GetIsAromatic() for n in a.GetNeighbors()) for a in selected):
        return "_arom"
    if any(a.IsInRing() or any(n.IsInRing() for n in a.GetNeighbors()) for a in selected):
        return "_alicyc"
    return ""


def fragments(smiles):
    heavy = Chem.MolFromSmiles(smiles)
    full = Chem.AddHs(heavy)
    suffix, claimed, result = fragment_context(heavy), set(), []
    formamide = [set(m[:3]) for m in heavy.GetSubstructMatches(FORMAMIDE)]
    acetamide = [set(m[:3]) for m in heavy.GetSubstructMatches(ACETAMIDE)]
    donor_matches = [(name, next(i for i in match if heavy.GetAtomWithIdx(i).GetAtomicNum() not in (1, 6)))
                     for name, pattern in HB_PATTERNS.items() for match in heavy.GetSubstructMatches(pattern)]
    for name in PRIORITY + sorted(set(PATTERNS) - set(PRIORITY)):
        for match in heavy.GetSubstructMatches(PATTERNS[name]):
            atoms, label = set(match), name
            if atoms & claimed:
                continue
            if name == "amide":
                if any(atoms <= x for x in formamide):
                    label = "formamide"
                elif any(atoms <= x for x in acetamide):
                    label = "acetamide"
            atoms.update(n.GetIdx() for i in tuple(atoms) if full.GetAtomWithIdx(i).GetAtomicNum() in POLAR_Z
                         for n in full.GetAtomWithIdx(i).GetNeighbors() if n.GetAtomicNum() == 1)
            donors = tuple(n for n, anchor in donor_matches if anchor in atoms)
            result.append(Fragment(label + suffix, tuple(sorted(atoms)), donors))
            claimed.update(atoms)
    for atom in heavy.GetAtoms():
        if atom.GetAtomicNum() in POLAR_Z and atom.GetIdx() not in claimed:
            atoms = {atom.GetIdx()}
            atoms.update(n.GetIdx() for n in full.GetAtomWithIdx(atom.GetIdx()).GetNeighbors() if n.GetAtomicNum() == 1)
            result.append(Fragment("other" + suffix, tuple(sorted(atoms)), ()))
            claimed.update(atoms)
    if len({i for f in result for i in f.atoms}) != sum(len(f.atoms) for f in result):
        raise ValueError("overlapping fragments")
    return result


def base_label(label):
    for suffix in ("_arom", "_alicyc"):
        if label.endswith(suffix):
            return label[:-len(suffix)]
    return label


def family(label):
    base = base_label(label)
    if base.startswith(("acetamide", "formamide", "amide", "aldehyde", "ketone", "ester", "anhydride", "carbonate", "carboxylic", "acyl_halide")):
        return "carbonyl"
    if base.startswith(("ether", "furan", "alcohol", "phenol", "epoxide")):
        return "oxygen"
    if base.startswith(("amine", "nitrile", "aromatic_N", "azole", "isocyanate", "urea")):
        return "nitrogen"
    if base.startswith(("nitro", "nitrate", "nitrite")):
        return "nitro"
    if base.startswith(("fluoro", "chloro", "bromo", "iodo")):
        return "halogen"
    if base.startswith(("phosphate", "phosphoryl")):
        return "phosphorus"
    if base.startswith(("sulfate", "sulfone", "sulfoxide", "thioether", "thiol", "thiocyanate")):
        return "sulfur"
    if base.startswith("siloxane"):
        return "silicon"
    return "other"


def atom_group(type_name):
    """Tie the acyclic (_A) and ring (_R) versions of an atom type."""
    return type_name[:-2] if type_name.endswith(("_A", "_R")) else type_name


def atom_type(atom):
    symbol = atom.GetSymbol()
    if symbol == "H":
        return "H"
    if symbol in ("F", "Cl", "Br", "I"):
        return symbol
    if symbol == "Si":
        return "Si_SP3_A"
    hybrid = {Chem.HybridizationType.SP: "SP", Chem.HybridizationType.SP2: "SP2",
              Chem.HybridizationType.SP3: "SP3"}.get(atom.GetHybridization())
    if hybrid is None:
        raise ValueError(f"unsupported hybridization for {symbol}")
    return f"{symbol}_{hybrid}_{'R' if atom.IsInRing() else 'A'}"
