from typing import List, Optional

import numpy as np
from sklearn.neighbors import BallTree

try:
    from iktos.logger import getLogger
except ImportError:
    from logging import getLogger

try:
    from openbabel.openbabel import OBAtom, OBMolAtomIter  # openbabel 3
except ModuleNotFoundError:
    from openbabel import (
        OBAtom,
        OBMolAtomIter,
    )  # openbabel 2 (warning in structure-utils)

from .atom_typing import (
    find_charged_atoms,
    find_h_bond_acceptors,
    find_h_bond_donors,
    find_hydrophobics,
    find_metal_binders,
    find_metals,
    find_pi_carbons,
    find_rings,
    find_x_bond_acceptors,
)
from .Ligand import Ligand
from .mol_utils import get_all_coordinates, map_atom_ids, read_obmol
from .utils import parse_pdb

logger = getLogger(__name__)


class Receptor:
    """Class to store binding site atoms and their properties"""

    def __init__(self, rec_pdb_coords: str, as_string: bool):
        # Read and parse receptor file/string
        rec_coords_block, mapping = parse_pdb(rec_pdb_coords, as_string=as_string)
        obmol = read_obmol(
            rec_coords_block, fmt="pdb", title="receptor", as_string=True
        )

        # Map atom IDs back to their original value
        self.obmol = map_atom_ids(obmol, mapping)

        # Store some info as public attributes
        self.atoms = [a for a in OBMolAtomIter(self.obmol)]
        self.coordinates = get_all_coordinates(self.obmol)
        self.bs_atoms: Optional[List[OBAtom]] = None

        # Find rings (do it here because we need the whole OBMol
        self.rings = find_rings(obmol)
        self.num_rings = len(self.rings)

    def detect_binding_site(self, ligands: List[Ligand], distance: float = 7.5) -> None:
        """Identifies atoms that are part of the binding site.

        Args:
            ligands: list of ligands.
            distance: radius of the binding site (default: 7.5).
        """
        # Detect receptor atoms close to the ligand (binding site) and identify
        # functional groups for these atoms
        logger.debug("Selecting binding site residues and atoms")
        # Use Ball Tree to store coords and to efficiently compute
        # the distance between receptor and ligand's atoms
        tree = BallTree(self.coordinates, leaf_size=5)
        # Get the index of receptor atoms near the ligand
        ligands_coordinates = np.concatenate([lig.coordinates for lig in ligands])
        indexes = np.unique(
            np.concatenate(tree.query_radius(ligands_coordinates, r=distance))
        )
        self.bs_atoms = list(np.array(self.atoms)[indexes])
        logger.debug(f"Selected {len(self.bs_atoms)} atoms as binding site")

    def identify_functional_groups(self):
        """Identifies functional groups in the receptor.

        By default, the method will use all the atoms of the receptor. To limit
        the identification to atoms that are close the the ligand (or any ligand
        if a list if given), call `detect_binding_site` first."""
        if self.bs_atoms is not None:
            obatoms = self.bs_atoms
        else:
            obatoms = self.atoms

        # Separate water from non-water atoms
        protein_atoms, water_atoms = [], []
        for obatom in obatoms:
            if not obatom.GetResidue().GetResidueProperty(9):
                protein_atoms.append(obatom)
            else:
                water_atoms.append(obatom)

        # Identify functional groups for non-water residues
        logger.debug("Looking at non-water residues on the receptor side")
        self.hydrophobics = find_hydrophobics(protein_atoms)
        self.h_bond_acceptors = find_h_bond_acceptors(protein_atoms)
        self.h_bond_donors = find_h_bond_donors(protein_atoms)
        self.charged_atoms = find_charged_atoms(protein_atoms)
        self.metals = find_metals(protein_atoms, "receptor")
        self.metal_binders = find_metal_binders(protein_atoms, "receptor")
        # self.halogens = find_halogens(protein_atoms) (proteins normally don't contain halogens)
        self.x_bond_acceptors = find_x_bond_acceptors(protein_atoms)
        self.pi_carbons = find_pi_carbons(protein_atoms)

        # Now look at water molecules
        logger.debug("Looking at water molecules on the receptor side")
        self.water_h_bond_acceptors = find_h_bond_acceptors(water_atoms)
        self.water_h_bond_donors = find_h_bond_donors(water_atoms)
        self.water_x_bond_acceptors = find_x_bond_acceptors(water_atoms)
        self.water_metal_binders = find_metal_binders(water_atoms, "water")
