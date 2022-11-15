"""
Copyright (C) Iktos - All Rights Reserved
Unauthorized copying of this file, via any medium is strictly prohibited
Proprietary and confidential
"""

import json, re

from iktos.structure_interactions.visualization.pymol import (
    prepare_session, prepare_session_multistate
)


def parse_multi_mol2(multi_mol2_block: str):
    blocks = []
    pattern = "@<TRIPOS>MOLECULE"
    mol2_content_split = re.split(pattern, multi_mol2_block)[1:]
    for mol2_block in mol2_content_split:
        mol2_block = pattern + mol2_block
        mol2_block_clean = ''.join(
            [line for line in mol2_block.splitlines(True) if not line.startswith("#")]
        )
        blocks.append(mol2_block_clean)
    return blocks

def test_pymol():
    with open("tests/files/protein.pdb", "r") as f:
        protein_xray = f.read()
    with open("tests/files/ligand.sdf", "r") as f:
        ligand_xray = f.read()
    with open("tests/files/contacts_xray.json") as f:
        contacts_xray = json.load(f)
    with open("tests/files/multi_mol2.mol2", "r") as f:
    	docking_poses = parse_multi_mol2(f.read())
    docking_poses = docking_poses[0:8]
    weights = {
	'Hydrophobic': {'MET|A|31|CE': 1.0,
	'LYS|A|52|CB': 1.0,
  	'VAL|A|85|CG1': 1.0,
  	'LEU|A|157|CD1': 1.0},
 	'Pi_Hydrophobic': {'TYR|A|101|CD1+CD2+CE1+CE2+CG+CZ': 1.0,
  	'VAL|A|39|CG1': 1.0},
 	'H_Bond_Weak': {'GLU|A|33|O': 0.5, 'TYR|A|103|CA+HA': 0.5},
 	'H_Bond': {'VAL|A|102|O': 1.0,
  	'ALA|A|154|O': 1.0,
  	'MET|A|104|H+N': 1.0,
  	'SER|A|167|HG+OG': 1.0}
    }
    prepare_session(
        protein_pdb_block=protein_xray,
        ligand_sdf_block=ligand_xray,
        contacts=contacts_xray,
        extra_mol2_blocks=docking_poses,
        weights=weights,
        output_file_path="tests/files/pymol.pse",
        color_bg="white",
        color_protein="cbaw",
        sphere_scale=2.0,
    )
    prepare_session_multistate(
        protein_pdb_blocks=[protein_xray] * 8,
        ligand_sdf_blocks=[ligand_xray] * 8,
        contacts=[contacts_xray] * 8,
        extra_mol2_blocks=[docking_poses],
        weights=weights,
        output_file_path="tests/files/pymol_multistate.pse",
        color_bg="white",
        color_protein="cbaw",
    )
