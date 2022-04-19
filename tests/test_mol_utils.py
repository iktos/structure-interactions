from iktos.structure_interactions.mol_utils import read_obmol


def test_read_obmol_1():
    smiles = 'CCCO'
    assert read_obmol(smiles, as_string=True, fmt='smi')


def test_read_obmol_2():
    file = 'tests/data/template.sdf'
    obmol = read_obmol(file, as_string=False, fmt='sdf')
    assert obmol
    assert obmol.GetResidue(0)
    assert obmol.GetResidue(0).GetChain() == 'A'


def test_read_obmol_3():
    file = 'tests/data/template_bis.sdf'
    assert read_obmol(file, as_string=False, fmt='sdf')


def test_read_obmol_4():
    file = 'tests/data/single_sdf_complete_new_format.sdf'
    assert read_obmol(file, as_string=False, fmt='sdf')


def test_read_obmol_5():
    file = 'tests/data/single_sdf_bug_obabel.sdf'
    obmol = read_obmol(file, as_string=False, fmt='sdf')
    assert obmol.NumResidues() == 1
    assert obmol.GetResidue(0).GetName() == 'UNL'
