"""
Copyright (C) Iktos - All Rights Reserved
Unauthorized copying of this file, via any medium is strictly prohibited
Proprietary and confidential
"""

import json
import os
import re

from iktos.structure_interactions.visualization.pymol import (
    prepare_session_inter,
    prepare_session_inter_multistate,
    prepare_session_intra,
)


with open("tests/data/prot_5UIT.pdb", "r") as f:
    PROTEIN_BLOCK = f.read()
with open("tests/data/lig_5UIT.sdf", "r") as f:
    LIGAND_BLOCK = f.read()
with open("tests/data/contacts_5UIT.json") as f:
    CONTACTS = json.load(f)
with open("tests/data/poses_in_5UIT.sdf", "r") as f:
    POSES = re.split('MOL', f.read())[1:]
WEIGHTS = {
    'Hydrophobic': {'MET|A|31|CE': 1.0},
}


def test_prepare_session_inter():
    output_file_path = 'pymol.pse'
    prepare_session_inter(
        protein_pdb_block=PROTEIN_BLOCK,
        ligand_sdf_block=LIGAND_BLOCK,
        contacts=CONTACTS,
        extra_sdf_blocks=POSES,
        weights=WEIGHTS,
        output_file_path=output_file_path,
        color_bg="white",
        color_protein="cbaw",
        sphere_scale=2.0,
    )
    assert os.path.exists(output_file_path)
    os.remove(output_file_path)


def test_prepare_session_inter_multistate():
    output_file_path = 'pymol_multistate.pse'
    prepare_session_inter_multistate(
        protein_pdb_blocks=[PROTEIN_BLOCK] * 2,
        ligand_sdf_blocks=[LIGAND_BLOCK] * 2,
        contacts=[CONTACTS] * 2,
        extra_sdf_blocks=[POSES],
        weights=WEIGHTS,
        output_file_path=output_file_path,
        color_bg="white",
        color_protein="cbaw",
    )
    assert os.path.exists(output_file_path)
    os.remove(output_file_path)


def test_prepare_session_intra():
    output_file_path = 'pymol_intra.pse'
    contacts_intra = {
        'H_Bond': [
            {
                'partner_1': [173, 174],
                'partner_2': [131],
                'type': 'strong',
            },
        ],
    }
    prepare_session_intra(
        coords_block=PROTEIN_BLOCK,
        fmt='pdb',
        contacts=contacts_intra,
        output_file_path=output_file_path,
        color_bg="white",
        color_molecule="cbaw",
    )
    assert os.path.exists(output_file_path)
    os.remove(output_file_path)
