from deepdiff import DeepDiff

from iktos.structure_interactions import (
    analyse_interactions_inter,
    analyse_interactions_inter_multi,
    convert_to_dict_inter,
)
from iktos.structure_interactions.InteractionParameters import InteractionParameters


def test_analyse_inter_3s3m():
    protein_path = "tests/data/prot_3S3M.pdb"
    ligand_path = "tests/data/lig_3S3M.sdf"
    contacts_raw = analyse_interactions_inter(
        rec_coords=protein_path,
        lig_coords=ligand_path,
        as_string=False,
        refine=True,
        lig_format="sdf",
    )
    dict_contacts = convert_to_dict_inter(contacts_raw)
    assert len(dict_contacts["Hydrophobic"]) == 2
    assert len(dict_contacts["Pi_Stacking"]) == 3
    assert len(dict_contacts["Pi_Hydrophobic"]) == 1
    assert len(dict_contacts["Metal_Complex"]) == 2
    assert len(dict_contacts["Water_Bridge"]) == 1
    assert len(dict_contacts) == 5


def test_analyse_inter_5n9t():
    protein_path = "tests/data/prot_5N9T.pdb"
    ligand_path = "tests/data/lig_5N9T.sdf"
    contacts_raw = analyse_interactions_inter(
        rec_coords=protein_path,
        lig_coords=ligand_path,
        as_string=False,
        refine=True,
        lig_format="sdf",
    )
    dict_contacts = convert_to_dict_inter(contacts_raw)
    assert len(dict_contacts["Hydrophobic"]) == 7
    assert len(dict_contacts["Pi_Hydrophobic"]) == 2
    assert len(dict_contacts["H_Bond"]) == 6
    assert len(dict_contacts["H_Bond_Weak"]) == 3
    assert len(dict_contacts["Multipolar"]) == 3
    assert len(dict_contacts) == 6


def test_analyse_inter_6nw6():
    protein_path = "tests/data/prot_6NW6.pdb"
    ligand_path = "tests/data/lig_6NW6.sdf"
    contacts_raw = analyse_interactions_inter(
        rec_coords=protein_path,
        lig_coords=ligand_path,
        as_string=False,
        refine=True,
        lig_format="sdf",
    )
    dict_contacts = convert_to_dict_inter(contacts_raw)
    assert len(dict_contacts["Hydrophobic"]) == 1
    assert len(dict_contacts["Pi_Hydrophobic"]) == 1
    assert len(dict_contacts["H_Bond"]) == 5
    assert len(dict_contacts["H_Bond_Weak"]) == 4
    assert len(dict_contacts["Pi_Amide"]) == 1
    assert len(dict_contacts) == 5


def test_analyse_inter_6nw6_refine():
    protein_path = "tests/data/prot_6NW6.pdb"
    ligand_path = "tests/data/lig_6NW6.sdf"
    contacts_raw = analyse_interactions_inter(
        rec_coords=protein_path,
        lig_coords=ligand_path,
        as_string=False,
        refine=False,
        lig_format="sdf",
    )
    dict_contacts = convert_to_dict_inter(contacts_raw)
    contacts_raw = analyse_interactions_inter(
        rec_coords=protein_path,
        lig_coords=ligand_path,
        as_string=False,
        refine=True,
        lig_format="sdf",
    )
    dict_contacts_refine = convert_to_dict_inter(contacts_raw)
    # check that refine removed 0 contact
    assert not DeepDiff(dict_contacts, dict_contacts_refine, ignore_order=True)


def test_analyse_inter_multi():
    protein_path = "tests/data/prot_6NO9.pdb"
    ligand_path = "tests/data/lig_6NO9.sdf"
    contacts_raw = analyse_interactions_inter_multi(
        rec_coords=protein_path,
        lig_coords=[ligand_path],
        as_string=False,
        refine=True,
        lig_format="sdf",
    )
    dicts_contacts = [convert_to_dict_inter(contact) for contact in contacts_raw]
    assert len(dicts_contacts[0]["H_Bond"]) == 3
    assert len(dicts_contacts[0]["H_Bond_Weak"]) == 1
    assert len(dicts_contacts[0]["Hydrophobic"]) == 3
    assert len(dicts_contacts[0]["Pi_Hydrophobic"]) == 3
    assert len(dicts_contacts[0]["Salt_Bridge"]) == 1
    assert len(dicts_contacts[0]["Water_Bridge"]) == 1
    assert len(dicts_contacts[0]) == 6


def test_analyse_inter_3mz9():
    # no H-bond detected because Hs are not oriented correctly
    protein_path = "tests/data/prot_3MZ9.pdb"
    ligand_path = "tests/data/lig_3MZ9.sdf"
    contacts_raw = analyse_interactions_inter(
        rec_coords=protein_path,
        lig_coords=ligand_path,
        as_string=False,
        refine=False,
        lig_format="sdf",
    )
    dict_contacts = convert_to_dict_inter(contacts_raw)
    assert "H_Bond" not in dict_contacts


def test_analyse_inter_3mz9_rotation():
    # if we allow Hs to rotate, then we detect 2 H-bonds
    protein_path = "tests/data/prot_3MZ9.pdb"
    ligand_path = "tests/data/lig_3MZ9.sdf"
    contacts_raw = analyse_interactions_inter(
        rec_coords=protein_path,
        lig_coords=ligand_path,
        as_string=False,
        refine=False,
        lig_format="sdf",
        parameters=InteractionParameters(allow_h_rotation=True),
    )
    dict_contacts = convert_to_dict_inter(contacts_raw)
    assert len(dict_contacts["H_Bond"]) == 2


def test_analyse_inter_prot_prot_interaction():
    spike_path = "tests/data/6m0j_spike.pdb"
    ace_path = "tests/data/6m0j_ace.pdb"
    contacts_raw = analyse_interactions_inter(
        rec_coords=spike_path,
        lig_coords=ace_path,
        as_string=False,
        lig_format="pdb",
    )
    dict_contacts = convert_to_dict_inter(contacts_raw)
    assert all(
        [
            "res_l" in info_contact.keys()
            for _, list_contacts in dict_contacts.items()
            for info_contact in list_contacts
        ]
    )
