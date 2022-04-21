from iktos.structure_interactions.analyse_complex import analyse_complex
from iktos.structure_interactions.utils import contacts_to_dict


def test_analyse_complex_3s3m():
    protein_path = 'tests/data/prot_3S3M.pdb'
    ligand_path = 'tests/data/lig_3S3M.sdf'

    contacts_ref = analyse_complex(
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


def test_analyse_complex_5n9t():

    protein_path = 'tests/data/prot_5N9T.pdb'
    ligand_path = 'tests/data/lig_5N9T.sdf'

    contacts_ref = analyse_complex(
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
