from __future__ import absolute_import

try:
    from iktos.logger import getLogger
except ImportError:
    from logging import getLogger

try:
    from openbabel.openbabel import OBMolAtomIter  # openbabel 3
except ModuleNotFoundError:
    from openbabel import OBMolAtomIter  # openbabel 2 (warning in structure-utils)

from .math_utils import get_centroid, get_euclidean_distance_3d
from .Mol import Mol
from .mol_utils import get_coords

logger = getLogger(__name__)


class Ligand(Mol):
    """Class to store ligand atoms and their properties"""

    def __init__(self, obmol_lig):
        Mol.__init__(self, mol_type='ligand')
        self.obmol = obmol_lig
        self.atoms = [a for a in OBMolAtomIter(self.obmol)]
        logger.debug(f'Found {len(self.atoms)} atoms in ligand')

        # Find rings
        self.rings = self.find_rings(obmol_lig)
        self.num_rings = len(self.rings)

    def identify_functional_groups(self):
        logger.debug('Identifying functional groups in ligand')

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

        # Find halogens
        self.halogens = self.find_halogens(self.atoms)

        # Find halogen-bond acceptors
        self.x_bond_acceptors = self.find_x_bond_acceptors(self.atoms)

        # Find pi-groups (for pi-stacking + multipolar interactions with F)
        self.pi_carbons = self.find_pi_carbons(self.atoms)

        # Define centroid for identification of binding site
        self.centroid = get_centroid([get_coords(a) for a in self.atoms])
        self.max_dist_to_center = max(
            [
                get_euclidean_distance_3d(self.centroid, get_coords(a))
                for a in self.atoms
            ]
        )
