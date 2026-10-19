# CRAFT-ε

**CRAFT-ε: A Conformer-Resolved Atom- and Fragment-based Theory for the Static Relative Permittivity of Pure
Organic Liquids.** Predicts the static relative permittivity ε_r of a neutral pure organic liquid at 298.15 K
from its SMILES string or name, with the 41 fitted parameters of the paper and GFN2-xTB/ALPB(diethyl ether)
conformer ensembles.

## Installation

```bash
pip install craft-epsilon
```

or from a local copy of this repository:

```bash
pip install .                 # from this folder
```

The [xtb](https://github.com/grimme-lab/xtb) program (version 6.7.1 was used to fit the model) is taken from
`xtb_binary`, from the `PATH` or from an earlier download. If none is available, the latest xtb release for the
current operating system and architecture (Linux and Windows x86_64) is downloaded once from GitHub, checked
against its SHA-256 checksum and stored in the folder `xtb` inside the installed package (`craft_epsilon/xtb`,
next to the module files; like the cache database of MoleculeResolver). If that folder is not writable, it goes to
`~/.craft_epsilon/xtb`; the environment variable `CRAFT_EPSILON_XTB_DIR` overrides both. The download is not
removed by `pip uninstall`; delete the folder by hand if needed.

## Usage

Command line (prints all relevant intermediate quantities):

```bash
craft-epsilon "CCO"
craft-epsilon ethanol
python -m craft_epsilon "CC(C)=O"
```

Python (returns values, prints nothing):

```python
from craft_epsilon import epsilon, predict

epsilon("CCO")            # static relative permittivity at 298.15 K
predict("ethanol")        # dict with class, polar groups, V_m, eps_inf, dipoles, conformer counts, ...
```

A valid SMILES string is used directly; anything else is resolved as a name with MoleculeResolver
(cross-checked, neutral structure; needs an internet connection).

**Temperature.** The model contains a temperature dependence through the orientation term and the molar volume,
V_m(T) = V_m(298.15 K)·(1 + α_V(T − 298.15 K)) with a fixed α_V = 10⁻³ K⁻¹, which was used to fit data measured at
other temperatures. CRAFT-ε is meant to be used at 298.15 K only: the conformer ensembles, their Boltzmann weights
and the molar volume refer to 298.15 K, and the temperature dependence of association is not represented.
`epsilon()`, `predict()` and the command line therefore always use 298.15 K.

## Tests

```bash
pip install pytest
pytest                        # fast tests, no xtb and no internet needed
pytest -m xtb                 # end-to-end calculations with xtb
pytest -m network             # name resolution with MoleculeResolver
```

The fast tests check the 41 parameters, the molar volume against the original smi2DC, the polar groups and
classes, the domain checks, the xtb lookup and the command line, and reproduce the paper predictions of 16 liquids
from their stored conformer ensembles (`tests/data`).

## License

MIT, see `LICENSE`.

## Package layout

| module | content |
|---|---|
| `craft_epsilon.model` | classification, effective dipole, Clausius–Mossotti ε∞, Onsager–Kirkwood–Fröhlich relation, `epsilon`, `predict` |
| `craft_epsilon.parameters` | the 41 fitted parameters |
| `craft_epsilon.conformers` | conformer workflow (RDKit ETKDGv3, GFN2-xTB/ALPB optimization, single point, qRRHO) |
| `craft_epsilon.descriptors` | Boltzmann weights, polarizabilities, polar-group dipole vectors |
| `craft_epsilon.groups` | polar groups (SMARTS), chemical families, atom types |
| `craft_epsilon.molar_volume` | molar volume of Bouteloup and Mathieu |
| `craft_epsilon.cli` | `craft-epsilon` command |

Conformer ensembles depend slightly on the RDKit version and platform (the embedding of a given random seed can
differ), so predictions can deviate from the paper values by up to a few per cent for some small rings and diols.

## References

- Molar volume: R. Bouteloup, D. Mathieu, *Predicting dielectric constants of pure liquids: fragment-based
  Kirkwood–Fröhlich model applicable over a wide range of polarity*, Phys. Chem. Chem. Phys. 21 (2019)
  11043–11057, [doi:10.1039/C9CP01704F](https://doi.org/10.1039/C9CP01704F).
- GFN2-xTB: C. Bannwarth, S. Ehlert, S. Grimme, J. Chem. Theory Comput. 15 (2019) 1652–1671,
  [doi:10.1021/acs.jctc.8b01176](https://doi.org/10.1021/acs.jctc.8b01176).
- ALPB: S. Ehlert, M. Stahn, S. Spicher, S. Grimme, J. Chem. Theory Comput. 17 (2021) 4250–4261,
  [doi:10.1021/acs.jctc.1c00471](https://doi.org/10.1021/acs.jctc.1c00471).
- ETKDGv3: S. Wang, J. Witek, G. A. Landrum, S. Riniker, J. Chem. Inf. Model. 60 (2020) 2044–2058,
  [doi:10.1021/acs.jcim.0c00025](https://doi.org/10.1021/acs.jcim.0c00025).
- spyrmsd: R. Meli, P. C. Biggin, J. Cheminform. 12 (2020) 49,
  [doi:10.1186/s13321-020-00455-2](https://doi.org/10.1186/s13321-020-00455-2).
- D4: E. Caldeweyher et al., J. Chem. Phys. 150 (2019) 154122,
  [doi:10.1063/1.5090222](https://doi.org/10.1063/1.5090222).
- qRRHO: S. Grimme, Chem. Eur. J. 18 (2012) 9955–9964,
  [doi:10.1002/chem.201200497](https://doi.org/10.1002/chem.201200497).
- MoleculeResolver: S. Müller, J. Cheminform. 17 (2025),
  [doi:10.1186/s13321-025-01064-7](https://doi.org/10.1186/s13321-025-01064-7).
