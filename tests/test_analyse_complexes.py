from iktos.structure_interactions.analyse_complex import analyse_complexes
from iktos.structure_interactions.utils import contacts_to_dict


def test_analyse_complexes():
    protein_path = 'tests/data/prot_6NO9.pdb'
    ligand_path = 'tests/data/lig_6NO9.sdf'

    contacts_raw = analyse_complexes(
        rec_coords=protein_path,
        lig_coords=[ligand_path],
        as_string=False,
        refine=True,
        lig_format='sdf',
    )

    dicts_contacts = [
        contacts_to_dict(contact) for contact in contacts_raw
    ]

    assert len(dicts_contacts[0]["H_Bond"]) == 3
    assert len(dicts_contacts[0]["H_Bond_Weak"]) == 1
    assert len(dicts_contacts[0]["Hydrophobic"]) == 2
    assert len(dicts_contacts[0]["Pi_Hydrophobic"]) == 2
    assert len(dicts_contacts[0]["Salt_Bridge"]) == 1
    assert len(dicts_contacts[0]["Water_Bridge"]) == 1
    assert len(dicts_contacts[0]) == 6
