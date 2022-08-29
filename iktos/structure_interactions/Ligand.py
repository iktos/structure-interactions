from __future__ import absolute_import

try:
    from iktos.logger import getLogger
except ImportError:
    from logging import getLogger

try:
    from openbabel.openbabel import OBMolAtomIter  # openbabel 3
except ModuleNotFoundError:
    from openbabel import (
        OBMolAtomIter,
    )  # openbabel 2 (warning in structure-utils)

from .atom_typing import (
    find_charged_atoms,
    find_h_bond_acceptors,
    find_h_bond_donors,
    find_halogens,
    find_hydrophobics,
    find_metal_binders,
    find_metals,
    find_pi_carbons,
    find_rings,
    find_x_bond_acceptors,
)
from .mol_utils import get_all_coordinates, map_atom_ids, read_obmol


logger = getLogger(__name__)


class Ligand:
    """Class to store ligand atoms and their properties"""

    def __init__(self, lig_coords: str, lig_format: str, as_string: bool):
        if lig_format != 'sdf':
            logger.warning(
                'It is recommended to use SDF blocks/files for the ligand; '
                'for other formats, make sure that formal (not partial) atomic charges '
                'are explicitely defined (otherwise negative charges like C(=O)[O-] '
                'will be missed)'
            )

        # Read and parse ligand file/string
        obmol = read_obmol(
            lig_coords, fmt=lig_format, title='ligand', as_string=as_string
        )

        # Map atom Id to their Idx value (1-based, needed by Pymol)
        self.obmol = map_atom_ids(obmol, None)

        # Store some info as public attributes
        self.coordinates = get_all_coordinates(self.obmol)
        logger.debug(f'Found {self.obmol.NumAtoms()} atoms in ligand')

        # Detect functional groups
        atoms = [a for a in OBMolAtomIter(self.obmol)]
        self.rings = find_rings(obmol)
        self.num_rings = len(self.rings)
        self.hydrophobics = find_hydrophobics(atoms)
        self.h_bond_acceptors = find_h_bond_acceptors(atoms)
        self.h_bond_donors = find_h_bond_donors(atoms)
        self.charged_atoms = find_charged_atoms(atoms)
        self.metals = find_metals(atoms, 'ligand')
        self.metal_binders = find_metal_binders(atoms, 'ligand')
        self.halogens = find_halogens(atoms)
        self.x_bond_acceptors = find_x_bond_acceptors(atoms)
        self.pi_carbons = find_pi_carbons(atoms)
