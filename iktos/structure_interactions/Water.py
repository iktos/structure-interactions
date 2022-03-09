"""
Copyright (C) Iktos - All Rights Reserved
Unauthorized copying of this file, via any medium is strictly prohibited
Proprietary and confidential
"""

from __future__ import absolute_import

from iktos.logger import getLogger

from .Mol import Mol


logger = getLogger(__name__)


class Water(Mol):
    """Class to store water atoms and their properties"""

    def __init__(self, obmol_rec):
        Mol.__init__(self, mol_type='water')

    def identify_functional_groups(self, obatoms_bs):
        logger.debug('Identifying functional groups in water')
        self.atoms = obatoms_bs

        # Find H-bond donors and acceptors
        self.h_bond_acceptors = self.find_h_bond_acceptors(self.atoms)
        self.h_bond_donors = self.find_h_bond_donors(self.atoms)

        # Find halogen-bond acceptors
        self.x_bond_acceptors = self.find_x_bond_acceptors(self.atoms)

        # Find metal binding atoms
        self.metal_binders = self.find_metal_binders(self.atoms)
