import json

from iktos.structure_interactions.InteractionProfiler import InteractionProfiler
from iktos.structure_interactions.utils import contacts_to_dict


def test_contacts_to_dict():
    # test that the conversion to dict gives exactly
    # what is expected in SA - if this test breaks, you should consider
    # carefully the modifications that are being made to the package
    # as they will probably compromise the compatibility with SA
    protein_path = 'tests/data/prot_5UIT.pdb'
    ligand_path = 'tests/data/lig_5UIT.sdf'

    plip = InteractionProfiler()
    contacts_ref = plip.analyse_complex(
        rec_coords=protein_path,
        lig_coords=ligand_path,
        as_string=False,
        refine=True,
        lig_format='sdf',
    )
    contacts_ref = contacts_to_dict(contacts_ref)

    with open('tests/data/contacts_5UIT.json', 'r') as f:
        expected_contacts = json.load(f)
    assert contacts_ref == expected_contacts
