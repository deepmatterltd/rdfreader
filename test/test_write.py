from tempfile import NamedTemporaryFile

from rdfreader.write import write_rdf


def test_write_rdf(first_sample_rxn):

    with NamedTemporaryFile("w+", suffix="rdf") as f:
        rxn_blocks = [first_sample_rxn] * 3

        write_rdf(f, rxn_blocks, [0, 1, 2])

        f.seek(0)

        rdf_text = f.read()

        # $RXN should occur once for each reaction block
        assert rdf_text.count("$RXN") == len(rxn_blocks)

        # first line should start with $RDFILE
        assert rdf_text.startswith("$RDFILE")
