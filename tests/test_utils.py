import json

from deepdiff import DeepDiff

from iktos.structure_interactions import analyse_interactions_inter, convert_to_dict_inter


def test_contacts_to_dict_1P5E():
    # test that the conversion to dict gives exactly
    # what is expected in SA - if this test breaks, you should consider
    # carefully the modifications that are being made to the package
    # as they will probably compromise the compatibility with SA
    protein_path = 'tests/data/prot_1P5E.pdb'
    ligand_path = 'tests/data/lig_1P5E.sdf'

    contacts_ref = analyse_interactions_inter(
        rec_coords=protein_path,
        lig_coords=ligand_path,
        as_string=False,
        refine=True,
        lig_format='sdf',
    )
    contacts_ref = convert_to_dict_inter(contacts_ref)

    with open('tests/data/contacts_1P5E.json', 'r') as f:
        expected_contacts = json.load(f)
    assert not DeepDiff(
        contacts_ref,
        expected_contacts,
        ignore_order=True,
        ignore_numeric_type_changes=True,
    )


def test_contacts_to_dict_3S3M():
    # test that the conversion to dict gives exactly
    # what is expected in SA - if this test breaks, you should consider
    # carefully the modifications that are being made to the package
    # as they will probably compromise the compatibility with SA
    protein_path = 'tests/data/prot_3S3M.pdb'
    ligand_path = 'tests/data/lig_3S3M.sdf'

    contacts_ref = analyse_interactions_inter(
        rec_coords=protein_path,
        lig_coords=ligand_path,
        as_string=False,
        refine=True,
        lig_format='sdf',
    )
    contacts_ref = convert_to_dict_inter(contacts_ref)

    with open('tests/data/contacts_3S3M.json', 'r') as f:
        expected_contacts = json.load(f)
    assert not DeepDiff(
        contacts_ref,
        expected_contacts,
        ignore_order=True,
        ignore_numeric_type_changes=True,
    )


def test_contacts_to_dict_5N9T():
    # test that the conversion to dict gives exactly
    # what is expected in SA - if this test breaks, you should consider
    # carefully the modifications that are being made to the package
    # as they will probably compromise the compatibility with SA
    protein_path = 'tests/data/prot_5N9T.pdb'
    ligand_path = 'tests/data/lig_5N9T.sdf'

    contacts_ref = analyse_interactions_inter(
        rec_coords=protein_path,
        lig_coords=ligand_path,
        as_string=False,
        refine=True,
        lig_format='sdf',
    )
    contacts_ref = convert_to_dict_inter(contacts_ref)

    with open('tests/data/contacts_5N9T.json', 'r') as f:
        expected_contacts = json.load(f)
    assert not DeepDiff(
        contacts_ref,
        expected_contacts,
        ignore_order=True,
        ignore_numeric_type_changes=True,
    )


def test_contacts_to_dict_5UIT():
    # test that the conversion to dict gives exactly
    # what is expected in SA - if this test breaks, you should consider
    # carefully the modifications that are being made to the package
    # as they will probably compromise the compatibility with SA
    protein_path = 'tests/data/prot_5UIT.pdb'
    ligand_path = 'tests/data/lig_5UIT.sdf'

    contacts_ref = analyse_interactions_inter(
        rec_coords=protein_path,
        lig_coords=ligand_path,
        as_string=False,
        refine=True,
        lig_format='sdf',
    )
    contacts_ref = convert_to_dict_inter(contacts_ref)

    with open('tests/data/contacts_5UIT.json', 'r') as f:
        expected_contacts = json.load(f)
    assert not DeepDiff(
        contacts_ref,
        expected_contacts,
        ignore_order=True,
        ignore_numeric_type_changes=True,
    )


def test_contacts_to_dict_6NO9():
    # test that the conversion to dict gives exactly
    # what is expected in SA - if this test breaks, you should consider
    # carefully the modifications that are being made to the package
    # as they will probably compromise the compatibility with SA
    protein_path = 'tests/data/prot_6NO9.pdb'
    ligand_path = 'tests/data/lig_6NO9.sdf'

    contacts_ref = analyse_interactions_inter(
        rec_coords=protein_path,
        lig_coords=ligand_path,
        as_string=False,
        refine=True,
        lig_format='sdf',
    )
    contacts_ref = convert_to_dict_inter(contacts_ref)

    with open('tests/data/contacts_6NO9.json', 'r') as f:
        expected_contacts = json.load(f)
    assert not DeepDiff(
        contacts_ref,
        expected_contacts,
        ignore_order=True,
        ignore_numeric_type_changes=True,
    )


def test_contacts_to_dict_6NW6():
    # test that the conversion to dict gives exactly
    # what is expected in SA - if this test breaks, you should consider
    # carefully the modifications that are being made to the package
    # as they will probably compromise the compatibility with SA
    protein_path = 'tests/data/prot_6NW6.pdb'
    ligand_path = 'tests/data/lig_6NW6.sdf'

    contacts_ref = analyse_interactions_inter(
        rec_coords=protein_path,
        lig_coords=ligand_path,
        as_string=False,
        refine=True,
        lig_format='sdf',
    )
    contacts_ref = convert_to_dict_inter(contacts_ref)

    with open('tests/data/contacts_6NW6.json', 'r') as f:
        expected_contacts = json.load(f)
    assert not DeepDiff(
        contacts_ref,
        expected_contacts,
        ignore_order=True,
        ignore_numeric_type_changes=True,
    )
