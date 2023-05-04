import datetime
from pathlib import Path
from typing import Any

import pytest
from pytest_mock import MockerFixture
from pytest_mock import mocker as mocker_fixture  # noqa: F401

from rdfreader.chem.mol import Molecule


@pytest.fixture
def mocker(mocker_fixture) -> MockerFixture:  # noqa: F811
    """
    Wraps pytest_mock.mocker fixture to allow for
    easier importing.
    """
    return mocker_fixture


def get_sample_mol_block() -> str:
    """
    Load the sample mol block as a string from the test resources.
    """
    sample_mol_block_path = "test/resources/sample_mol_block.txt"
    with open(sample_mol_block_path, "r") as f:
        sample_mol_block = f.read()
    return sample_mol_block


def get_rdf_path() -> Path:
    """
    Load the sample rdf as a string from the test resources.
    """
    rdf_path = Path("test/resources/sample_rdf.rdf")
    return rdf_path


@pytest.fixture
def rdf_path() -> Path:
    return get_rdf_path()


@pytest.fixture
def first_sample_rxn() -> str:
    """
    Return the first rxn block from the sample rdf.
    """
    with open("test/resources/sample_rdf_first_rxn.rxn", "r") as f:
        first_sample_rxn = f.read()
    return first_sample_rxn


def get_sample_rdf_string() -> str:
    """
    Load the sample rdf as a string from the test resources.
    """
    sample_rdf_string_path = get_rdf_path()
    with open(sample_rdf_string_path, "r") as f:
        sample_rdf_string = f.read()
    return sample_rdf_string


def get_sample_rxn_block() -> str:
    """
    Load the sample rxn block as a string from the test resources.
    """
    sample_rxn_block_path = "test/resources/sample_rxn_block.txt"
    with open(sample_rxn_block_path, "r") as f:
        sample_rxn_block = f.read()
    return sample_rxn_block


@pytest.fixture
def sample_mol_block() -> str:
    return get_sample_mol_block()


@pytest.fixture
def sample_mol_block_lines() -> str:
    """
    Return the sample mol block split into a list of lines.
    """
    return get_sample_mol_block().split("\n")


@pytest.fixture
def sample_rxn_block() -> str:
    return get_sample_rxn_block()


@pytest.fixture
def sample_rxn_block_lines() -> str:
    """
    Return the sample rxn block split into a list of lines.
    """
    return get_sample_rxn_block().split("\n")


@pytest.fixture
def sample_molecule() -> Molecule:
    """
    Create a test molecule.
    """
    mol = Molecule()
    mol.mol_block = get_sample_mol_block()
    return mol


@pytest.fixture
def sample_molecule_metadata() -> dict[str, Any]:
    """
    Return the sample mol block metadata.
    """
    return dict(
        molecule_name="sample name",
        user_initials="II",
        program_name="PPPPPPPP",
        date_time=datetime.datetime(22, 5, 24, 14, 23),
        dimensional_codes="dd",
        scaling_factor_1=12,
        scaling_factor_2=1.12345678,
        energy=1.2345678912,
        registry_number="overflowing reg number",
        comment="sample comment",
    )


@pytest.fixture
def sample_rxn_block_metadata() -> str:
    """
    Return sample rxn block metadata.
    """
    return dict(
        reaction_name="sample reaction name",
        user_initials="IIIIII",
        program_name="PPPPPPPPP",
        date_time=datetime.datetime(2022, 5, 24, 14, 55),
        registry_number="RRRRRRR",
        comment="sample reaction comment",
        product_count=1,
        reactant_count=3,
    )
