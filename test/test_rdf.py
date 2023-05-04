from unittest.mock import MagicMock

import pytest
from pytest_mock import MockerFixture

from rdfreader.chem.reaction import Reaction
from rdfreader.exceptions import InvalidReactionError
from rdfreader.rdf import RDFParser, parse_rdf_reg_num


@pytest.fixture
def sample_rdf_file() -> str:
    return "test/resources/spresi-100.rdf"


def test_parse_rdf_reg_num():
    reg_num = parse_rdf_reg_num("$RFMT $RIREG 4620744")
    assert reg_num == "4620744"


def test_parse_rdf(sample_rdf_file: str):
    with open(sample_rdf_file, "r") as f:
        rdf_parser = RDFParser(f)
        reactions = [reaction for reaction in rdf_parser]

    assert len(reactions) == 100

    # we can't test all 100 reactions in all this detail, so we'll just check the first two
    expected_reaction_ids = ["1274842", "808226"]
    expected_line_numbers = [3, 212]
    expected_product_counts = [1, 1]
    expected_reactant_counts = [3, 2]
    expected_catalyst_counts = [1, 1]
    expected_solvent_counts = [0, 0]
    expected_other_reagent_counts = [0, 0]

    for reaction, eid, lineno, product_count, reactant_count, catalyst_count, solvent_count, other_count in zip(
        reactions,
        expected_reaction_ids,
        expected_line_numbers,
        expected_product_counts,
        expected_reactant_counts,
        expected_catalyst_counts,
        expected_solvent_counts,
        expected_other_reagent_counts,
    ):
        assert reaction.id == eid
        assert reaction.lineno == lineno
        assert reaction.rdf_metadata == {
            "version": "1",
            "date_stamp": "02/12/04 11:58",
        }
        assert len(reaction.products) == product_count
        assert len(reaction.reactants) == reactant_count
        assert len(reaction.catalysts) == catalyst_count
        assert len(reaction.solvents) == solvent_count
        assert len(reaction.other_reagents) == other_count


@pytest.fixture
def reaction_raise_exception(mocker: MockerFixture) -> MagicMock:
    mocker.patch.object(Reaction, "__init__", side_effect=InvalidReactionError("Test exception"))


def test_parse_rdf_catches(reaction_raise_exception: MagicMock, sample_rdf_file: str):
    with open(sample_rdf_file, "r") as f:
        rdf_parser = RDFParser(f, except_on_invalid_reaction=False)
        reactions = [reaction for reaction in rdf_parser]

    for reaction in reactions:
        assert reaction is None


def test_parse_rdf_raises(reaction_raise_exception: MagicMock, sample_rdf_file: str):
    with pytest.raises(InvalidReactionError):
        with open(sample_rdf_file, "r") as f:
            rdf_parser = RDFParser(f)
            [reaction for reaction in rdf_parser]
