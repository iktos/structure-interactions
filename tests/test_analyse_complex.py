from iktos.structure_interactions import analyse_complex, analyse_complexes, contacts_to_dict


def test_analyse_complex_3s3m():
    protein_path = 'tests/data/prot_3S3M.pdb'
    ligand_path = 'tests/data/lig_3S3M.sdf'
    contacts_raw = analyse_complex(
        rec_coords=protein_path,
        lig_coords=ligand_path,
        as_string=False,
        refine=True,
        lig_format='sdf',
    )
    dict_contacts = contacts_to_dict(contacts_raw)
    assert len(dict_contacts["Hydrophobic"]) == 2
    assert len(dict_contacts["Pi_Stacking"]) == 3
    assert len(dict_contacts["Pi_Hydrophobic"]) == 1
    assert len(dict_contacts["Metal_Complex"]) == 2
    assert len(dict_contacts["Water_Bridge"]) == 1
    assert len(dict_contacts) == 5


def test_analyse_complex_5n9t():
    protein_path = 'tests/data/prot_5N9T.pdb'
    ligand_path = 'tests/data/lig_5N9T.sdf'
    contacts_raw = analyse_complex(
        rec_coords=protein_path,
        lig_coords=ligand_path,
        as_string=False,
        refine=True,
        lig_format='sdf',
    )
    dict_contacts = contacts_to_dict(contacts_raw)
    assert len(dict_contacts["Hydrophobic"]) == 7
    assert len(dict_contacts["Pi_Hydrophobic"]) == 2
    assert len(dict_contacts["H_Bond"]) == 5
    assert len(dict_contacts["H_Bond_Weak"]) == 2
    assert len(dict_contacts["Multipolar"]) == 3
    assert len(dict_contacts) == 5


def test_analyse_complex_6nw6():
    protein_path = 'tests/data/prot_6NW6.pdb'
    ligand_path = 'tests/data/lig_6NW6.sdf'
    contacts_raw = analyse_complex(
        rec_coords=protein_path,
        lig_coords=ligand_path,
        as_string=False,
        refine=True,
        lig_format='sdf',
    )
    dict_contacts = contacts_to_dict(contacts_raw)
    assert len(dict_contacts["Hydrophobic"]) == 1
    assert len(dict_contacts["Pi_Hydrophobic"]) == 1
    assert len(dict_contacts["H_Bond"]) == 5
    assert len(dict_contacts["H_Bond_Weak"]) == 3
    assert len(dict_contacts["Pi_Amide"]) == 1
    assert len(dict_contacts) == 5


def test_analyse_complex_6nw6_refine():
    protein_path = 'tests/data/prot_6NW6.pdb'
    ligand_path = 'tests/data/lig_6NW6.sdf'
    contacts_raw = analyse_complex(
        rec_coords=protein_path,
        lig_coords=ligand_path,
        as_string=False,
        refine=False,
        lig_format='sdf',
    )
    dict_contacts = contacts_to_dict(contacts_raw)
    contacts_raw = analyse_complex(
        rec_coords=protein_path,
        lig_coords=ligand_path,
        as_string=False,
        refine=True,
        lig_format='sdf',
    )
    dict_contacts_refine = contacts_to_dict(contacts_raw)
    # check that refine removed 0 contact
    assert dict_contacts == dict_contacts_refine


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
    assert len(dicts_contacts[0]["Hydrophobic"]) == 3
    assert len(dicts_contacts[0]["Pi_Hydrophobic"]) == 3
    assert len(dicts_contacts[0]["Salt_Bridge"]) == 1
    assert len(dicts_contacts[0]["Water_Bridge"]) == 1
    assert len(dicts_contacts[0]) == 6
