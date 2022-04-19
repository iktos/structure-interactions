"""
Copyright (C) Iktos - All Rights Reserved
Unauthorized copying of this file, via any medium is strictly prohibited
Proprietary and confidential
"""

from iktos.structure_interactions.Ligand import Ligand
from iktos.structure_interactions.mol_utils import read_obmol


def test_identify_functional_groups_5N9T():
    obmol = read_obmol(
        'tests/data/lig_5N9T.sdf',
        as_string=False,
        fmt='sdf',
    )
    ligand = Ligand(obmol)
    ligand.identify_functional_groups()
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
    obmol = read_obmol(
        'tests/data/lig_3S3M.sdf',
        as_string=False,
        fmt='sdf',
    )
    ligand = Ligand(obmol)
    ligand.identify_functional_groups()
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
