"""
Copyright (C) Iktos - All Rights Reserved
Unauthorized copying of this file, via any medium is strictly prohibited
Proprietary and confidential
"""

from iktos.structure_interactions.Ligand import Ligand


def test_identify_functional_groups_5N9T():
    file = 'tests/data/lig_5N9T.sdf'
    ligand = Ligand(file, 'sdf', as_string=False)
    assert len(ligand.rings) == 4
    assert len(ligand.hydrophobics) == 19
    assert len(ligand.h_bond_acceptors) == 5
    assert len(ligand.h_bond_donors) == 23
    assert len(ligand.charged_atoms) == 1
    assert len(ligand.halogens) == 3
    assert len(ligand.x_bond_acceptors) == 5
    assert len(ligand.pi_carbons) == 1
    assert len(ligand.metal_binders) == 5
    assert len(ligand.metals) == 0


def test_identify_functional_groups_3S3M():
    file = 'tests/data/lig_3S3M.sdf'
    ligand = Ligand(file, 'sdf', as_string=False)
    assert len(ligand.rings) == 2
    assert len(ligand.hydrophobics) == 9
    assert len(ligand.h_bond_acceptors) == 5
    assert len(ligand.h_bond_donors) == 14
    assert len(ligand.charged_atoms) == 0
    assert len(ligand.halogens) == 2
    assert len(ligand.x_bond_acceptors) == 5
    assert len(ligand.pi_carbons) == 2
    assert len(ligand.metal_binders) == 5
    assert len(ligand.metals) == 0


def test_identify_functional_groups_old_bug():
    # This SDF used to fail with an AttributeError related to a missing residue
    # The residue was correctly created by 'read_mol' but then disappeared
    # New strategy: 'read_mol' does not try to create a residue if there is none,
    # the rest of the code is made to handle with/without residue on the ligand side
    file = 'tests/data/ligand_bug_residue.sdf'
    ligand = Ligand(file, 'sdf', as_string=False)
    assert len(ligand.rings) == 1
    assert len(ligand.hydrophobics) == 4
    assert len(ligand.h_bond_acceptors) == 4
    assert len(ligand.h_bond_donors) == 9
    assert len(ligand.charged_atoms) == 0
    assert len(ligand.halogens) == 0
    assert len(ligand.x_bond_acceptors) == 5
    assert len(ligand.pi_carbons) == 0  # acid excluded
    assert len(ligand.metal_binders) == 5
    assert len(ligand.metals) == 0
