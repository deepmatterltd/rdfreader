from rdkit.Chem import Mol, MolFromMolBlock, MolToSmiles

from rdfreader.chem.mol import Molecule


def assert_molecule_from_mol_block(mol: Molecule, mol_block: str):
    """
    Assert that a molecule object created from a mol block string has the
    correct properties.
    """
    assert mol.rd_mol is not None
    assert mol.mol_block == mol_block


def test_molecule_from_mol_block(sample_mol_block):
    """
    Test that the Molecule.from_mol_block method works, creating a Molecule
    object from a mol block string and correctly creates class attributes.
    """
    mol = Molecule.from_mol_block(sample_mol_block)
    assert_molecule_from_mol_block(mol, sample_mol_block)


def test_molecule_init_from_mol_block(sample_mol_block):
    """
    Test we can init a molecule object directly by passing a mol block string.
    """
    mol = Molecule(sample_mol_block)
    assert_molecule_from_mol_block(mol, sample_mol_block)


def test_create_empty_molecule():
    """
    Test that a Molecule object can be created without a mol block string.
    """
    mol = Molecule()
    assert mol.mol_block is None


def test_molecule_to_rdkit_mol(sample_molecule, sample_mol_block):
    """
    Test the Molecule.rd_mol property.
    """
    rd_mol: Mol = sample_molecule.rd_mol
    # can't directly compare Mol objects, so we'll just check that it is of
    # the right type and is not None
    assert rd_mol is not None
    assert isinstance(rd_mol, Mol)


def test_molecule_to_smiles(sample_molecule, sample_mol_block):
    """
    Test the output of the Molecule.smiles property. Verifies the smiles
    matches that output by RDKit.
    """
    rd_mol: Mol = MolFromMolBlock(sample_mol_block)
    rd_smiles: str = MolToSmiles(rd_mol)
    assert sample_molecule.smiles == rd_smiles


def test_molecule_metadata(sample_molecule, sample_molecule_metadata):
    """
    Test the Molecule.metadata property.
    """
    assert sample_molecule.metadata == sample_molecule_metadata


def test_molecule_from_smiles():
    """
    Test the Molecule.from_smiles method.
    """
    smiles = "OCO"
    mol = Molecule.from_smiles(smiles)
    assert mol.rd_mol is not None
    assert mol.smiles == smiles
