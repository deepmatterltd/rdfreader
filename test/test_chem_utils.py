import pytest

from rdfreader.chem.mol import Molecule
from rdfreader.chem.utils import mol_list_to_smiles, reaction_smiles


@pytest.fixture
def molecule_list() -> list[Molecule]:
    return [
        Molecule.from_smiles("C"),
        Molecule.from_smiles("CO"),
        Molecule.from_smiles("OCO"),
    ]


@pytest.fixture
def molecule_list_smiles() -> str:
    return "C.CO.OCO"


def test_mol_list_to_smiles(molecule_list, molecule_list_smiles):
    assert mol_list_to_smiles(molecule_list) == molecule_list_smiles


def test_reaction_to_smiles(molecule_list, molecule_list_smiles):
    assert (
        reaction_smiles(molecule_list, molecule_list, molecule_list)
        == f"{molecule_list_smiles}>{molecule_list_smiles}>{molecule_list_smiles}"  # noqa: E501
    )


def test_reaction_to_smiles_no_reagents(molecule_list, molecule_list_smiles):
    assert reaction_smiles(molecule_list, molecule_list) == f"{molecule_list_smiles}>>{molecule_list_smiles}"
