from __future__ import absolute_import

try:
    from iktos.logger import getLogger
except ImportError:
    from logging import getLogger

from .atom_typing import (
    find_h_bond_acceptors,
    find_h_bond_donors,
    find_metal_binders,
    find_x_bond_acceptors,
)

logger = getLogger(__name__)


class Water:
    """Class to store water atoms and their properties"""

    def __init__(self):
        super().__init__()

    def identify_functional_groups(self, obatoms_bs):
        logger.debug('Identifying functional groups in water')
        self.atoms = obatoms_bs

        # Find H-bond donors and acceptors
        self.h_bond_acceptors = find_h_bond_acceptors(self.atoms)
        self.h_bond_donors = find_h_bond_donors(self.atoms)

        # Find halogen-bond acceptors
        self.x_bond_acceptors = find_x_bond_acceptors(self.atoms)

        # Find metal binding atoms
        self.metal_binders = find_metal_binders(self.atoms, 'water')
