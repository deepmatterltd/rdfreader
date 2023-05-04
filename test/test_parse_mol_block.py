from rdfreader.parse.molblock import _parse_large_regno, get_mol_block_metadata


def test_parse_large_regno():
    assert _parse_large_regno(["", "M  REG     0123456    \n", ""]) == "0123456"


def test_get_mol_block_metadata(sample_mol_block, sample_molecule_metadata):
    """
    Components of this function are tested seperately, only need to test the
    once for integration.
    """
    expected_result = sample_molecule_metadata
    test_result = get_mol_block_metadata(sample_mol_block)

    assert test_result == expected_result
