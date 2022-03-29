from iktos.structure_interactions.InteractionProfiler import InteractionProfiler
from iktos.structure_interactions.utils import contacts_to_dict


def test_analyse_complex_3S3M():
    protein_path = 'tests/data/prot_3S3M.pdb'
    ligand_path = 'tests/data/lig_3S3M.sdf'

    plip = InteractionProfiler()
    contacts_ref = plip.analyse_complex(
        rec_coords=protein_path,
        lig_coords=ligand_path,
        as_string=False,
        refine=True,
        lig_format='sdf',
    )
    contacts_ref = contacts_to_dict(contacts_ref)

    assert len(contacts_ref["Hydrophobic"]) == 2
    assert len(contacts_ref["Pi_Stacking"]) == 3
    assert len(contacts_ref["Pi_Hydrophobic"]) == 1
    assert len(contacts_ref["Metal_Complex"]) == 2
    assert len(contacts_ref["Water_Bridge"]) == 1


def test_analyse_complex_5N9T():

    protein_path = 'tests/data/prot_5N9T.pdb'
    ligand_path = 'tests/data/lig_5N9T.sdf'

    plip = InteractionProfiler()
    contacts_ref = plip.analyse_complex(
        rec_coords=protein_path,
        lig_coords=ligand_path,
        as_string=False,
        refine=True,
        lig_format='sdf',
    )
    contacts_ref = contacts_to_dict(contacts_ref)

    assert len(contacts_ref["Hydrophobic"]) == 7
    assert len(contacts_ref["Pi_Hydrophobic"]) == 2
    assert len(contacts_ref["H_Bond"]) == 5
    assert len(contacts_ref["H_Bond_Weak"]) == 2
    assert len(contacts_ref["Multipolar"]) == 3
