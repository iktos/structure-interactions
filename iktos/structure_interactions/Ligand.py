from __future__ import absolute_import

try:
    from iktos.logger import getLogger
except ImportError:
    from logging import getLogger

try:
    from openbabel.openbabel import OBMolAtomIter  # openbabel 3
except ModuleNotFoundError:
    from openbabel import OBMolAtomIter  # openbabel 2 (warning in structure-utils)

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
from .math_utils import get_centroid, get_euclidean_distance_3d
from .mol_utils import get_coords

logger = getLogger(__name__)


class Ligand:
    """Class to store ligand atoms and their properties"""

    def __init__(self, obmol_lig):
        super().__init__()

        self.obmol = obmol_lig
        self.atoms = [a for a in OBMolAtomIter(self.obmol)]
        logger.debug(f'Found {len(self.atoms)} atoms in ligand')

        # Find rings
        self.rings = find_rings(obmol_lig)
        self.num_rings = len(self.rings)

    def identify_functional_groups(self):
        logger.debug('Identifying functional groups in ligand')

        # Find hydrophobic atoms
        self.hydrophobics = find_hydrophobics(self.atoms)

        # Find H-bond donors and acceptors
        self.h_bond_acceptors = find_h_bond_acceptors(self.atoms)
        self.h_bond_donors = find_h_bond_donors(self.atoms)

        # Find charged atoms
        self.charged_atoms = find_charged_atoms(self.atoms)

        # Find metals and metal binders (atoms with lone pair)
        self.metals = find_metals(self.atoms, 'ligand')
        self.metal_binders = find_metal_binders(self.atoms, 'ligand')

        # Find halogens
        self.halogens = find_halogens(self.atoms)

        # Find halogen-bond acceptors
        self.x_bond_acceptors = find_x_bond_acceptors(self.atoms)

        # Find pi-groups (for pi-stacking + multipolar interactions with F)
        self.pi_carbons = find_pi_carbons(self.atoms)

        # Define centroid for identification of binding site
        self.centroid = get_centroid([get_coords(a) for a in self.atoms])
        self.max_dist_to_center = max(
            [
                get_euclidean_distance_3d(self.centroid, get_coords(a))
                for a in self.atoms
            ]
        )
