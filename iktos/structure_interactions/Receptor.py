from __future__ import absolute_import

try:
    from iktos.logger import getLogger
except ImportError:
    from logging import getLogger

from .Mol import Mol

logger = getLogger(__name__)


class Receptor(Mol):
    """Class to store binding site atoms and their properties"""

    def __init__(self, obmol_rec):
        Mol.__init__(self, mol_type='receptor')

        # Find rings
        self.rings = self.find_rings(obmol_rec)
        self.num_rings = len(self.rings)

    def identify_functional_groups(self, obatoms_bs):
        logger.debug('Identifying functional groups in receptor')
        self.atoms = obatoms_bs

        # Find hydrophobic atoms
        self.hydrophobics = self.find_hydrophobics(self.atoms)

        # Find H-bond donors and acceptors
        self.h_bond_acceptors = self.find_h_bond_acceptors(self.atoms)
        self.h_bond_donors = self.find_h_bond_donors(self.atoms)

        # Find charged atoms
        self.charged_atoms = self.find_charged_atoms(self.atoms)

        # Find metals and metal binders (atoms with lone pair)
        self.metals = self.find_metals(self.atoms)
        self.metal_binders = self.find_metal_binders(self.atoms)

        # Find halogens (there should not be any but just in case)
        self.halogens = self.find_halogens(self.atoms)

        # Find halogen-bond acceptors
        self.x_bond_acceptors = self.find_x_bond_acceptors(self.atoms)

        # Find pi-groups (for pi-stacking + multipolar interactions with F)
        self.pi_carbons = self.find_pi_carbons(self.atoms)
