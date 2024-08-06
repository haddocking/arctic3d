from io import StringIO
from pathlib import Path

import pytest

from arctic3d.modules.blast import blast_remote, parse_xml

from . import golden_data


@pytest.fixture
def fasta_file():
    return Path(golden_data, "1crn.fasta")


@pytest.fixture
def xml_file():
    return Path(golden_data, "blast.xml")


def test_blast_remote(mocker, fasta_file, xml_file):

    # Mock the remote blast call by passing the xml file as stringIo as return value
    with open(xml_file, "r") as f:
        xml = f.read()
    xml_string_io = StringIO(xml)

    mock_qbplast = mocker.patch("Bio.Blast.NCBIWWW.qblast")
    mock_qbplast.return_value = xml_string_io

    accession_id = blast_remote(fasta_file)

    assert accession_id == "P01541"


@pytest.mark.skip(reason="Not implemented")
def test_blast_local():
    pass


def test_parse_xml(xml_file):
    accession_id = parse_xml(xml_file)
    assert accession_id == "P01541"
