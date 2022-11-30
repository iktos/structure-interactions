from iktos.structure_interactions import convert_to_dict_inter
from iktos.structure_interactions.Ligand import Ligand
from iktos.structure_interactions.Receptor import Receptor
from iktos.structure_interactions.detection import (
    find_metal_complexes,
    Metal_Complex,
)


def test_find_metal_complexes():
    ligand = Ligand('tests/data/lig_3S3M.sdf', 'sdf', as_string=False)
    receptor = Receptor('tests/data/prot_3S3M.pdb', as_string=False)
    receptor.detect_binding_site([ligand], distance=10.0)
    receptor.identify_functional_groups()
    complexes = find_metal_complexes(
        receptor.metals + ligand.metals,
        receptor.metal_binders + ligand.metal_binders + receptor.water_metal_binders,
        distance_max=3.0,
    )
    assert len(complexes) == 2
    expected_contacts = {
        'Metal_Complex': [
            {
                'at_m': [10098],
                'at_name_m': ['MG'],
                'res_p': 'MG|A|596',
                'at_l': [3, 5],
                'at_name_l': ['O3', 'O5'],
                'at_p': [1987, 2865, 10130, 10508],
                'at_name_p': ['OD1', 'OD1', 'O', 'O'],
                'num_partners': 6,
                'complex_num': 1,
            },
            {
                'at_m': [10099],
                'at_name_m': ['MG'],
                'res_p': 'MG|A|597',
                'at_l': [4, 5],
                'at_name_l': ['O4', 'O5'],
                'at_p': [1988, 3405, 3406, 10775],
                'at_name_p': ['OD2', 'OE1', 'OE2', 'O'],
                'num_partners': 6,
                'complex_num': 2,
            },
        ]
    }
    assert convert_to_dict_inter(complexes) == expected_contacts
