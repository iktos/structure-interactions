from iktos.structure_interactions import analyse_interactions_intra, convert_to_dict_intra
from iktos.structure_interactions.InteractionParameters import InteractionParameters


def test_analyse_intra_ligand():
    interactions = analyse_interactions_intra(
        'tests/data/ligand_with_hbond.sdf',
        fmt='sdf',
        as_string=False,
        parameters=InteractionParameters(hbond_don_angle_min=130, hbond_acc_angle_min=80),
    )
    contacts = convert_to_dict_intra(interactions)
    expected_contacts = {
        'H_Bond': [
            {
                'partner_1': [14, 15],
                'partner_2': [10],
                'distance_ah': 1.7941293152947477,
                'distance_ad': 2.651371154704675,
                'angle_dha': 139.35255939604295,
                'type': 'strong'
            }
        ]
    }
    assert contacts == expected_contacts


def test_analyse_intra_peptide():
    interactions = analyse_interactions_intra(
        'tests/data/peptide.pdb',
        is_small_molecule=False,
        as_string=False,
    )
    contacts = convert_to_dict_intra(interactions)
    assert len(contacts) == 2
    assert len(contacts['H_Bond']) == 11
    assert len(contacts['Hydrophobic']) == 4
    

def test_analyse_intra_dna():
    interactions = analyse_interactions_intra(
        'tests/data/small_dna.pdb',
        is_small_molecule=False,
        as_string=False,
    )
    contacts = convert_to_dict_intra(interactions)
    assert len(contacts) == 4
    assert len(contacts['H_Bond']) == 14
    assert len(contacts['Hydrophobic']) == 1
    assert len(contacts['H_Bond_Weak']) == 1
    assert len(contacts['Pi_Stacking']) == 15
