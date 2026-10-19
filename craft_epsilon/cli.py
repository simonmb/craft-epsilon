"""Command line interface: craft-epsilon "<SMILES or name>" (or python -m craft_epsilon ...)."""
import logging
import sys
from .model import predict


def report(result):
    protocol = result["conformer_protocol"]
    lines = [
        "CRAFT-epsilon: static relative permittivity (41 parameters, GFN2-xTB/ALPB diethyl ether)",
        f"  input                      {result['input']}",
        f"  SMILES                     {result['smiles']}   (from {result['resolved_from']})",
        f"  class                      {result['class']}",
        "  polar groups               " + (", ".join(f"{g['label']} [{g['parameter']} = {g['scale']:.4f}]"
                                                   for g in result["polar_groups"]) or "none"),
        f"  remainder scale            {result['rest_scale']:.4f}" + ("  (nonpolar dipole scale c_np)"
                                                                    if result["class"] == "non_polar" else ""),
        f"  conformers                 {protocol['n_diastereomers']} diastereomer ensemble(s), "
        f"{protocol['n_embedded']} embedded ({protocol['n_polar_h_variants']} polar-H variants), "
        f"{protocol['n_optimized']} optimized, {protocol['n_unique']} unique, "
        f"{result['n_conformers']} incl. mirror images",
        f"  wall time                  {protocol['elapsed_seconds']:.1f} s   (xtb: {result['xtb']})",
        f"  temperature                {result['temperature_K']:.2f} K",
        f"  molar volume V_m           {result['molar_volume_cm3_mol']:.2f} cm3/mol   (Bouteloup and Mathieu 2019)",
        f"  polarizability alpha       {result['alpha_boltzmann_A3']:.3f} A3 (xTB/D4), {result['alpha_scaled_A3']:.3f} A3 scaled",
        f"  eps_inf                    {result['eps_inf']:.4f}",
        f"  rms dipole                 {result['mu_rms_xtb_D']:.3f} D (xTB), {result['mu_eff_rms_D']:.3f} D effective",
        f"  orientation term A         {result['A']:.4f}",
        f"  epsilon_r                  {result['epsilon']:.3f}",
    ]
    if protocol.get("stereo_fallback"):
        lines.append(f"  note                       {protocol['stereo_fallback']}")
    if protocol["failures"]:
        lines.append(f"  note                       {len(protocol['failures'])} xtb calculation(s) failed and were skipped")
    lines += ["", "Molar volume: R. Bouteloup, D. Mathieu, Phys. Chem. Chem. Phys. 21 (2019) 11043-11057,",
              "              https://doi.org/10.1039/C9CP01704F"]
    return "\n".join(lines)


def main(argv=None):
    """Entry point of the craft-epsilon command; returns the exit code."""
    argv = sys.argv[1:] if argv is None else argv
    logging.basicConfig(level=logging.INFO, format="%(message)s", stream=sys.stderr)   # e.g. xtb download
    if len(argv) != 1 or argv[0] in ("-h", "--help"):
        from . import __doc__ as package_doc
        print(package_doc.split("\n\nRequirements")[0])
        print('\nusage: craft-epsilon "<SMILES or name>"')
        return 2
    try:
        print(report(predict(argv[0])))
    except (ValueError, FileNotFoundError, RuntimeError) as error:
        print(f"error: {error}", file=sys.stderr)
        return 1
    return 0
