import pytest

try:
    from openbabel.openbabel import OBMolAtomIter
except ModuleNotFoundError:
    from openbabel import OBMolAtomIter
from iktos.structure_interactions.mol_utils import read_obmol


def test_read_obmol_1():
    smiles = 'CCCO'
    obmol = read_obmol(smiles, as_string=True, fmt='smi')
    assert obmol is not None


def test_read_obmol_2():
    file = 'tests/data/template.sdf'
    obmol = read_obmol(file, as_string=False, fmt='sdf')
    assert obmol is not None
    # assert obmol.GetResidue(0)
    # assert obmol.GetResidue(0).GetChain() == 'A'


def test_read_obmol_3():
    file = 'tests/data/template_bis.sdf'
    obmol = read_obmol(file, as_string=False, fmt='sdf')
    assert obmol is not None


def test_read_obmol_4():
    file = 'tests/data/single_sdf_complete_new_format.sdf'
    obmol = read_obmol(file, as_string=False, fmt='sdf')
    assert obmol is not None


def test_read_obmol_5():
    file = 'tests/data/single_sdf_bug_obabel.sdf'
    obmol = read_obmol(file, as_string=False, fmt='sdf')
    assert obmol is not None
    # assert obmol.NumResidues() == 1
    # assert obmol.GetResidue(0).GetName() == 'UNL'


def test_read_obmol_6():
    # This SDF used to fail with an AttributeError related to a missing residue
    # The residue was correctly created by 'read_mol' but then disappeared
    # New strategy: 'read_mol' does not try to create a residue if there is none,
    # the rest of the code is made to handle with/without residue on the ligand side
    file = 'tests/data/ligand_bug_residue.sdf'
    obmol = read_obmol(file, as_string=False, fmt='sdf')
    assert obmol.NumResidues() == 0


def test_read_obmol_7():
    # Follow-up to previous test and bug...
    # To reproduce the bug, this SDF raises an error when we create the residue
    # and try to access it, eventhough we can see it being created:
    file = 'tests/data/ligand_bug_residue.sdf'
    obmol = read_obmol(file, as_string=False, fmt='sdf')
    assert obmol.NumResidues() == 0
    obres = obmol.NewResidue()
    obres.SetName('UNL')
    for obatom in OBMolAtomIter(obmol):
        obres.AddAtom(obatom)
    assert obmol.NumResidues() == 1
    with pytest.raises(AttributeError):
        obmol.GetAtom(1).GetResidue().GetName()


def test_read_obmol_8():
    # Follow-up to previous test and bug...
    # Whereas with this SDF, there is no problem:
    file = 'tests/data/lig_3S3M.sdf'
    obmol = read_obmol(file, as_string=False, fmt='sdf')
    assert obmol.NumResidues() == 0
    obres = obmol.NewResidue()
    obres.SetName('UNL')
    for obatom in OBMolAtomIter(obmol):
        obres.AddAtom(obatom)
    assert obmol.NumResidues() == 1
    assert obmol.GetAtom(1).GetResidue().GetName() == 'UNL'


def test_read_obmol_9():
    # Follow-up to previous test and bug...
    # In this case, we don't even need to create the missing residue:
    file = 'tests/data/lig_3S3M.sdf'
    obmol = read_obmol(file, as_string=False, fmt='sdf')
    assert obmol.NumResidues() == 0
    assert obmol.GetAtom(1).GetResidue().GetName() == 'UNL'
    assert obmol.NumResidues() == 1  # residue created out of the blue!


def test_read_obmol_10():
    smiles = 'toto'  # invalid -> raise error
    with pytest.raises(ValueError):
        read_obmol(smiles, as_string=True, fmt='smi')


def test_read_obmol_11():
    with pytest.raises(ValueError):
        read_obmol(
            'tests/data/lig_3S3M.sdf',
            as_string=True,
            fmt='mol2',  # wrong format
        )
