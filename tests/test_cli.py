"""Command line interface: exit codes and messages (no xtb needed)."""
from craft_epsilon.cli import main


def test_help_returns_2(capsys):
    assert main(["--help"]) == 2
    out = capsys.readouterr().out
    assert "CRAFT-epsilon" in out and 'craft-epsilon "<SMILES or name>"' in out


def test_wrong_number_of_arguments_returns_2():
    assert main([]) == 2
    assert main(["CCO", "CCC"]) == 2


def test_out_of_domain_molecule_returns_1(capsys):
    assert main(["[Na+].[Cl-]"]) == 1
    assert "only single molecules" in capsys.readouterr().err
