from datetime import datetime
from io import TextIOWrapper


def write_rdf(
    f: TextIOWrapper,
    rxn_blocks: list[str],
    rxn_ids: list[str] = None,
):
    """Write a RDF file from a list of reaction blocks.

    Parameters
    ----------
    f : TextIOWrapper
        The file to write to.
    rxn_blocks : list[str]
        The reaction blocks to write.
    rxn_ids : list[str], optional
        The reaction IDs to use. Defaults to None. If None, sequential 5 digit
        numbers are used.
    """
    # current date and time as DD/MM/YY HH:MM
    datm = datetime.now().strftime("%d/%m/%y %H:%M")

    rdf_header = f"$RDFILE 1\n$DATM {datm}\n"
    rxn_header = "$RFMT $RIREG {}\n"

    if rxn_ids is None:
        # generate sequential 5 digit numbers
        rxn_ids = [f"{i:05d}" for i in range(1, len(rxn_blocks) + 1)]

    f.write(rdf_header)

    for rxn_id, rxn_block in zip(rxn_ids, rxn_blocks):
        f.write(rxn_header.format(rxn_id))
        f.write(rxn_block)
