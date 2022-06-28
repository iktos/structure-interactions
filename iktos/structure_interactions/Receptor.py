from __future__ import absolute_import

try:
    from iktos.logger import getLogger
except ImportError:
    from logging import getLogger

from .atom_typing import (
    find_charged_atoms,
    find_h_bond_acceptors,
    find_h_bond_donors,
    find_halogens,
    find_hydrophobics,
    find_metals,
    find_metal_binders,
    find_pi_carbons,
    find_rings,
    find_x_bond_acceptors,
)

logger = getLogger(__name__)


class Receptor:
    """Class to store binding site atoms and their properties"""

    def __init__(self, obmol_rec):
        super().__init__()

        # Find rings
        self.rings = find_rings(obmol_rec)
        self.num_rings = len(self.rings)

    def identify_functional_groups(self, obatoms_bs):
        logger.debug('Identifying functional groups in receptor')
        self.atoms = obatoms_bs  # receptor atoms in the binding site except water

        # Find hydrophobic atoms
        self.hydrophobics = find_hydrophobics(self.atoms)

        # Find H-bond donors and acceptors
        self.h_bond_acceptors = find_h_bond_acceptors(self.atoms)
        self.h_bond_donors = find_h_bond_donors(self.atoms)

        # Find charged atoms
        self.charged_atoms = find_charged_atoms(self.atoms)

        # Find metals and metal binders (atoms with lone pair)
        self.metals = find_metals(self.atoms, 'receptor')
        self.metal_binders = find_metal_binders(self.atoms, 'receptor')

        # Find halogens (there should not be any but just in case)
        self.halogens = find_halogens(self.atoms)

        # Find halogen-bond acceptors
        self.x_bond_acceptors = find_x_bond_acceptors(self.atoms)

        # Find pi-groups (for pi-stacking + multipolar interactions with F)
        self.pi_carbons = find_pi_carbons(self.atoms)
