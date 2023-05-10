import pytest
from rdkit.Chem import MolFromMolBlock

from rdfreader.chem.mol import Molecule
from rdfreader.parse.rxnblock import (
    DatumParser,
    get_rxn_block_metadata,
    mol_blocks_from_rxn_block,
    parse_dtype_string,
    preprocess_datum_string,
    validate_rxn_block,
)


def get_sample_datum_mol_block_processed():
    file_path = "test/resources/sample_datum_mol_block_processed.txt"
    with open(file_path, "r") as f:
        mol_block = f.read()
    return mol_block


@pytest.fixture
def sample_datum_mol_block():
    file_path = "test/resources/sample_datum_mol_block.txt"
    with open(file_path, "r") as f:
        mol_block = f.read()
    return mol_block


@pytest.fixture
def sample_datum_mol_block_parsed():
    file_path = "test/resources/sample_datum_mol_block_parsed.txt"
    with open(file_path, "r") as f:
        mol_block = f.read()
    return mol_block


@pytest.fixture
def empty_datum_parser():
    return


def test_get_rxn_block_metadata(sample_rxn_block, sample_rxn_block_metadata):
    """Test that the get_rxn_block_metadata function correctly parses the
    sample rxn block metadata."""
    assert sample_rxn_block_metadata == get_rxn_block_metadata(sample_rxn_block)


@pytest.mark.parametrize(
    "sample_dtype_string, expected_dtype_string",
    [
        ("$DTYPE RXN:VARIATION:STEPNO:SOLVENT(1):MOL:SYMBOL", "rxn_variation_stepno_solvent_1__mol_symbol"),
        ("$DTYPE RXN:VARIATION:PRODUCT:YIELD", "rxn_variation_product_yield"),
        ("$DTYPE RXN:CLASSIFICATION(1):MEDIUM", "rxn_classification_1__medium"),
        ("$DTYPE RXN:VARIATION:LITREF:JOURNAL_ISSN", "rxn_variation_litref_journal_issn"),
    ],
)
def test_parse_dtype_string(sample_dtype_string, expected_dtype_string):
    assert expected_dtype_string == parse_dtype_string(sample_dtype_string)


def test_datum_parser_rxnblock(sample_rxn_block):
    """Test that the DatumParser parses the sample rxn block."""
    datum_parser = DatumParser(sample_rxn_block)

    for parsed_dtype, datum in datum_parser:
        # difficult to test this as there's such a range of datatypes and
        # functions used in the parsing
        # but we can at least test that what we get is the right types
        assert isinstance(parsed_dtype, str)
        assert isinstance(datum, (str, Molecule))


@pytest.mark.parametrize(
    "expected_datum_string, sample_datum_string",
    [
        ("methanol", "$DATUM methanol"),
        ("87.0-87.0", "$DATUM 87.0-87.0"),
        ("384991457334703", "$DATUM 384991457334703"),
        ("0040-4039", "$DATUM 0040-4039"),
        ("a\nmultiline\nstring", "$DATUM a\nmultiline\nstring"),
    ],
)
def test_preprocess_datum_string(expected_datum_string, sample_datum_string):
    result = preprocess_datum_string(sample_datum_string)
    assert expected_datum_string == result


def test_validate_rxn_block(sample_rxn_block):
    """Test that the validate_rxn_block function correctly validates the sample
    rxn block."""
    assert validate_rxn_block(sample_rxn_block)


def test_validate_rxn_block_invalid():
    """Test that the validate_rxn_block function correctly validates the sample
    rxn block."""
    assert not validate_rxn_block("invalid rxn block")


def test_mol_blocks_from_rxn_block(sample_rxn_block):
    """Test that the mol_blocks_from_rxn_block function correctly parses the
    sample rxn block."""
    reactant_count: int = 3
    product_count: int = 1

    reactants, products = mol_blocks_from_rxn_block(sample_rxn_block, reactant_count, product_count)

    # verify that the correct number of mol_blocks are returned
    assert len(products) == product_count
    assert len(reactants) == reactant_count

    # verify that the mol_blocks can be parsed
    for molecules in [reactants, products]:
        for mol_block in molecules:
            rd_mol = MolFromMolBlock(mol_block)
            assert rd_mol is not None
