import datetime

import pytest

from rdfreader.parse.utils import (
    CTF_DEFAULT_LETTER_TO_FIELD_MAPPING,
    CTF_DEFAULT_MOLBLOCK_HEADER_FORMAT_STRING,
    CTF_RXNBLOCK_HEADER_FORMAT_STRING,
    _default_line_item,
    _parse_block_header_line,
    dict_elements_to_datetime,
    get_line_item,
    get_whole_line_item,
    make_string_python_safe,
    parse_format_string,
    parse_yield,
)


def test_default_line_item_string():
    assert _default_line_item(str) == ""


def test_default_line_item_string_non_default():
    assert _default_line_item(str, "default") == "default"


def test_default_line_item_int():
    assert _default_line_item(int) == 0


def test_default_line_item_int_non_default():
    assert _default_line_item(int, 10) == 10


def test_default_line_item_float():
    assert _default_line_item(float) == 0.0


def test_default_line_item_floa_non_default():
    assert _default_line_item(float, 0.1) == 0.1


def test_default_line_item_other():
    """Test that None is returned when not using string, int, or float."""
    assert _default_line_item(list, None) is None


def test_get_line_item_string():
    assert get_line_item("line      \n", (0, 2)) == "li"


def test_get_line_item_empty_string():
    assert get_line_item("") == ""


def test_get_line_item_empty_string_default():
    assert get_line_item("", default="default") == "default"


def test_get_line_item_int():
    assert get_line_item("12      \n", (0, 2), int) == 12


def test_get_line_item_empty_int():
    assert get_line_item("", cast_type=int) == 0


def test_get_line_item_empty_int_default():
    assert get_line_item("", cast_type=int, default=12) == 12


def test_get_line_item_float():
    assert get_line_item("12      \n", (0, 2), float) == 12


def test_get_line_item_empty_float():
    assert get_line_item("", cast_type=float) == 0.0


def test_get_line_item_empty_float_default():
    assert get_line_item("", cast_type=float, default=12.0) == 12.0


def test_get_line_item_casting_exceptions_thrown():
    with pytest.raises(ValueError):
        get_line_item("12.0", cast_type=int, catch_casting_exceptions=False)


def test_get_line_item_casting_exceptions_caught():
    assert (
        get_line_item(
            "12.0",
            cast_type=int,
            catch_casting_exceptions=True,
            default="default",
        )
        == "default"
    )


def test_parse_format_string():
    format_string = "IIIIIIPPPPPPPPPMMDDYYYYHHmmRRRRRRR"

    expected_result = {
        "I": (0, 6),
        "P": (6, 15),
        "M": (15, 17),
        "D": (17, 19),
        "Y": (19, 23),
        "H": (23, 25),
        "m": (25, 27),
        "R": (27, 34),
    }

    assert parse_format_string(format_string) == expected_result


def test_dict_elements_to_datetime():
    """Test that a dictionary of elements can be converted to a datetime
    object."""
    elements = {
        "month": 1,
        "day": 2,
        "year": 3,
        "hour": 4,
        "minute": 5,
        "test": "test",
    }

    expected_result = {
        "test": "test",  # ensure non-datetime elements are not modified
        "date_time": datetime.datetime(
            year=3,
            month=1,
            day=2,
            hour=4,
            minute=5,
        ),
    }

    assert dict_elements_to_datetime(elements) == expected_result


def test_parse_block_header_line_with_molblock(sample_mol_block_lines, sample_molecule_metadata):
    expected_result = sample_molecule_metadata
    expected_result["registry_number"] = "RRRRRR"
    # This is the only thing that is different from the
    # sample_molecule_metadata as the large regno is parsed seperately.
    # comment and name are held elsewhere, delete them from the dict
    del expected_result["comment"]
    del expected_result["molecule_name"]
    test_result = _parse_block_header_line(
        sample_mol_block_lines[1],
        CTF_DEFAULT_MOLBLOCK_HEADER_FORMAT_STRING,
        CTF_DEFAULT_LETTER_TO_FIELD_MAPPING,
    )

    assert test_result == expected_result


def test_parse_block_header_line_with_molblock_defaults():
    expected_result = dict(
        user_initials="",
        program_name="",
        date_time=None,
        dimensional_codes="",
        scaling_factor_1=0,
        scaling_factor_2=0.0,
        energy=0.0,
        registry_number="",
    )

    test_result = _parse_block_header_line(
        "",
        CTF_DEFAULT_MOLBLOCK_HEADER_FORMAT_STRING,
        CTF_DEFAULT_LETTER_TO_FIELD_MAPPING,
    )

    assert test_result == expected_result


def test_parse_block_header_line_with_rxnblock(sample_rxn_block_lines, sample_rxn_block_metadata):
    expected_result = sample_rxn_block_metadata
    for absent_key in [
        "comment",
        "reaction_name",
        "reactant_count",
        "product_count",
    ]:
        del expected_result[absent_key]  # not present in this line

    test_result = _parse_block_header_line(
        sample_rxn_block_lines[2],
        CTF_RXNBLOCK_HEADER_FORMAT_STRING,
        CTF_DEFAULT_LETTER_TO_FIELD_MAPPING,
    )

    assert test_result == expected_result


def test_parse_block_header_line_with_rxnblock_defaults():
    expected_result = dict(
        user_initials="",
        program_name="",
        date_time=None,
        registry_number="",
    )

    test_result = _parse_block_header_line(
        "",
        CTF_RXNBLOCK_HEADER_FORMAT_STRING,
        CTF_DEFAULT_LETTER_TO_FIELD_MAPPING,
    )

    assert test_result == expected_result


def test_get_whole_line_item():
    """Test that a line item can be retrieved from a line."""
    line = "line      \n"
    assert get_whole_line_item(line) == "line"


@pytest.mark.parametrize(
    "test_string,expected_result",
    [
        ("", None),
        (" ", None),
        (" \n\r1_:;'/test`string\n __", "_1_test_string"),
    ],
)
def test_make_string_python_safe(test_string, expected_result):
    test_result = make_string_python_safe(test_string)
    assert test_result == expected_result


def test_parse_yield():
    test_strings = [
        "17.0-17.0",
        "17",
        "17-17",
        "17.0 -- 17.0",
        "17 - 17",
        "17 -17",
        "17- 17",
        "17.0",
        "17  17",
        "17;17",
        "17:17",
        "17,17",
        "16-18",
    ]

    expected_result = 17.0

    for test_string in test_strings:
        assert parse_yield(test_string) == expected_result


def test_parse_yield_none():
    """Test that None is returned if the yield cannot be parsed."""
    test_strings = ["some text", "-1"]

    for test_string in test_strings:
        assert parse_yield(test_string) is None
