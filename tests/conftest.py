"""Shared fixtures. Reference data (tests/data): stored GFN2-xTB/ALPB(diethyl ether) conformer records of 16 liquids
measured at 298.15 K, the in-sample predictions of the paper fit and the molar volumes of the original smi2DC."""
import json
from pathlib import Path

import pytest

DATA = Path(__file__).parent / "data"


@pytest.fixture(scope="session")
def reference():
    return json.loads((DATA / "reference.json").read_text(encoding="utf-8"))


@pytest.fixture(scope="session")
def records(reference):
    return {smiles: json.loads((DATA / "records" / entry["record"]).read_text(encoding="utf-8"))
            for smiles, entry in reference.items()}
