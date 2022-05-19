from __future__ import absolute_import

from typing import Sequence, Union

from sklearn.neighbors import BallTree

try:
    from iktos.logger import getLogger
except ImportError:
    from logging import getLogger

try:
    from openbabel.openbabel import OBResidueAtomIter, OBResidueIter  # openbabel 3
except ModuleNotFoundError:
    from openbabel import (
        OBResidueIter,
        OBResidueAtomIter,
    )  # openbabel 2 (warning in structure-utils)

import numpy as np

from .detection import (
    H_Bond,
    Halogen_Bond,
    Hydrophobic,
    Metal_Complex,
    Multipolar,
    Pi_Amide,
    Pi_Cation,
    Pi_Hydrophobic,
    Pi_Stacking,
    Salt_Bridge,
    Water_Bridge,
    find_h_bonds,
    find_hydrophobics,
    find_metal_complexes,
    find_multpipolar_interactions,
    find_pi_amides,
    find_pi_cations,
    find_pi_hydrophobics,
    find_pi_stackings,
    find_salt_bridges,
    find_water_bridges,
    find_x_bonds,
)
from .InteractionParameters import InteractionParameters
from .Ligand import Ligand
from .mol_utils import map_atom_ids, read_obmol
from .Receptor import Receptor
from .refinement import (
    refine_h_bonds,
    refine_hydrophobics,
    refine_pi_cations,
    refine_pi_hydrophobics,
    refine_water_bridges,
)
from .utils import parse_pdb
from .Water import Water

logger = getLogger(__name__)


class InteractionProfiler:
    """
    Analyse protein-ligand interactions using an in-house
    version of PLIP (Protein-Ligand Interaction Profiler)
    """

    def load_receptor(self, rec_coords: str, as_string: bool) -> bool:
        """Loads the receptor and initialises water object."""
        # Read and parse receptor file/string
        rec_coords_block, mapping = parse_pdb(rec_coords, as_string=as_string)
        obmol_rec = read_obmol(
            rec_coords_block, fmt='pdb', title='receptor', as_string=True
        )
        if obmol_rec is None:
            return False

        # Map atom Id back to their original value
        self.obmol_rec = map_atom_ids(obmol_rec, mapping)

        # Initialise water and receptor objects
        logger.info('Initialising water object')
        self.wat = Water()
        logger.info('Initialising receptor object')
        self.rec = Receptor(self.obmol_rec)
        return True

    def load_ligand(
        self, lig_coords: str, lig_format: str, as_string: bool, bs_dist: float
    ) -> bool:
        """Loads the ligand and initialises ligand and binding site objects."""
        if lig_format != 'sdf':
            logger.warning(
                'It is recommended to use SDF blocks/files for the ligand; '
                'for other formats, make sure that formal (not partial) atomic charges '
                'are explicitely defined (otherwise negative charges like C(=O)[O-] '
                'will be missed)'
            )

        # Read and parse ligand file/string
        obmol_lig = read_obmol(
            lig_coords, fmt=lig_format, title='ligand', as_string=as_string
        )
        if obmol_lig is None:
            return False

        # Map atom Id to their Idx value (1-based, needed by Pymol)
        obmol_lig = map_atom_ids(obmol_lig, None)

        # Initialise ligand object
        logger.info('Initialising ligand object')
        self.lig = Ligand(obmol_lig)
        self.lig.identify_functional_groups()

        # Finalise initialisation of receptor and water objects
        logger.info('Selecting binding site residues and atoms')
        bs_atoms_refined, water_atoms = self._extract_binding_site(bs_dist=bs_dist)
        self.rec.identify_functional_groups(bs_atoms_refined)
        self.wat.identify_functional_groups(water_atoms)
        return True

    def analyse_interactions(self, refine: bool, parameters: InteractionParameters):
        """Performs protein-ligand interaction analysis
        for 1 complex (after initialisation of the objects)."""
        # Setup PLIP config

        # Analyse contacts
        self._detect_interactions(parameters=parameters)

        # Refine contacts
        if refine:
            self._refine_interactions()

        # Parse results
        if refine:
            return self.interactions
        return self.interactions_all

    def _extract_binding_site(self, bs_dist: float):
        """
        Select atoms that belong to the binding site
        """
        obres_rec = [obres for obres in OBResidueIter(self.obmol_rec)]
        obatoms_rec = []
        water_rec = []
        for residue in obres_rec:
            if not residue.GetResidueProperty(9):
                obatoms_rec += [atom for atom in OBResidueAtomIter(residue)]
            else:
                water_rec += [atom for atom in OBResidueAtomIter(residue)]
        all_atoms = np.array(obatoms_rec + water_rec)
        # Extracting residue atoms coordinates:
        coords_atom_rec = np.array(
            [
                (obatom_r.GetX(), obatom_r.GetY(), obatom_r.GetZ())
                for obatom_r in all_atoms
            ]
        )
        # Extracting ligand atoms coordinates:
        coords_l = np.array(
            [
                (obatom_l.GetX(), obatom_l.GetY(), obatom_l.GetZ())
                for obatom_l in self.lig.atoms
            ]
        )
        # Use Ball Tree to store coords and to efficiently compute the distance between receptor and ligand's atoms
        tree = BallTree(coords_atom_rec, leaf_size=5)
        # Get the index of atoms near the ligand
        ind = np.unique(np.concatenate(tree.query_radius(coords_l, r=bs_dist)))
        ind_atom_receptor, *_ = np.where(ind < len(obatoms_rec))
        ind_water, *_ = np.where(ind >= len(obatoms_rec))
        bs_atoms_refined = all_atoms[ind[ind_atom_receptor]]
        water_atoms = all_atoms[ind[ind_water]]
        logger.debug(f'Selected {len(bs_atoms_refined)} atoms as binding site')
        return bs_atoms_refined, water_atoms

    def _detect_interactions(self, parameters: InteractionParameters):
        """
        Find all receptor-ligand interactions
        Call functions stroed in detection.py
        """
        logger.info('Analysing interactions')
        self.hydrophobics_all = find_hydrophobics(
            self.rec.hydrophobics,
            self.lig.hydrophobics,
            parameters.min_dist,
            parameters.hydrophobic_dist_max,
        )
        logger.debug(f'Found {len(self.hydrophobics_all)} hydrophobic interaction(s)')

        self.pi_stackings_all = find_pi_stackings(
            self.rec.rings,
            self.lig.rings,
            parameters.min_dist,
            parameters.pistacking_dist_max_t,
            parameters.pistacking_dist_max_f,
            parameters.pistacking_dist_max_p,
            parameters.pistacking_ang_dev,
            parameters.pistacking_offset_max,
        )
        logger.debug(f'Found {len(self.pi_stackings_all)} pi-stacking interaction(s)')

        self.pi_amides_all = find_pi_amides(
            self.rec.pi_carbons,
            self.lig.rings,
            parameters.min_dist,
            parameters.piother_dist_max,
            parameters.piother_offset_max,
            parameters.pistacking_ang_dev,
        ) + find_pi_amides(
            self.rec.rings,
            self.lig.pi_carbons,
            parameters.min_dist,
            parameters.piother_dist_max,
            parameters.piother_offset_max,
            parameters.pistacking_ang_dev,
        )
        logger.debug(f'Found {len(self.pi_amides_all)} pi-amide like interaction(s)')

        self.pi_cations_all = find_pi_cations(
            self.rec.rings,
            self.lig.charged_atoms,
            parameters.min_dist,
            parameters.piother_dist_max,
            parameters.piother_offset_max,
            True,
        ) + find_pi_cations(
            self.lig.rings,
            self.rec.charged_atoms,
            parameters.min_dist,
            parameters.piother_dist_max,
            parameters.piother_offset_max,
            False,
        )
        logger.debug(f'Found {len(self.pi_cations_all)} pi-cation interaction(s)')

        self.pi_hydrophobics_all = find_pi_hydrophobics(
            self.rec.rings,
            self.lig.hydrophobics,
            parameters.min_dist,
            parameters.piother_dist_max,
            parameters.piother_offset_max,
            True,
        ) + find_pi_hydrophobics(
            self.lig.rings,
            self.rec.hydrophobics,
            parameters.min_dist,
            parameters.piother_dist_max,
            parameters.piother_offset_max,
            False,
        )
        logger.debug(
            f'Found {len(self.pi_hydrophobics_all)} pi-hydrophobic interaction(s)'
        )

        self.h_bonds_all = (
            find_h_bonds(
                self.rec.h_bond_acceptors,
                self.lig.h_bond_donors,
                parameters.min_dist,
                parameters.hbond_dist_max,
                parameters.hbond_don_angle_min,
                parameters.hbond_acc_angle_min,
                False,
            )
            + find_h_bonds(
                self.wat.h_bond_acceptors,
                self.lig.h_bond_donors,
                parameters.min_dist,
                parameters.hbond_dist_max,
                parameters.hbond_don_angle_min,
                parameters.hbond_acc_angle_min,
                False,
            )
            + find_h_bonds(
                self.lig.h_bond_acceptors,
                self.rec.h_bond_donors,
                parameters.min_dist,
                parameters.hbond_dist_max,
                parameters.hbond_don_angle_min,
                parameters.hbond_acc_angle_min,
                True,
            )
            + find_h_bonds(
                self.lig.h_bond_acceptors,
                self.wat.h_bond_donors,
                parameters.min_dist,
                parameters.hbond_dist_max,
                parameters.hbond_don_angle_min,
                parameters.hbond_acc_angle_min,
                True,
            )
        )
        logger.debug(f'Found {len(self.h_bonds_all)} H-bond(s)')

        self.x_bonds_all = find_x_bonds(
            self.rec.x_bond_acceptors,
            self.lig.halogens,
            parameters.min_dist,
            parameters.xbond_dist_max,
            parameters.xbond_don_angle_min,
            parameters.xbond_acc_angle_min,
        ) + find_x_bonds(
            self.lig.x_bond_acceptors,
            self.rec.halogens,
            parameters.min_dist,
            parameters.xbond_dist_max,
            parameters.xbond_don_angle_min,
            parameters.xbond_acc_angle_min,
        )
        logger.debug(f'Found {len(self.x_bonds_all)} halogen bond(s)')

        self.mulipolar_all = find_multpipolar_interactions(
            self.rec.pi_carbons,
            self.lig.halogens,
            parameters.min_dist,
            parameters.multipolar_dist_max,
            parameters.multipolar_don_angle_min,
            parameters.multipolar_norm_angle_max,
        )
        logger.debug(f'Found {len(self.mulipolar_all)} multipolar interaction(s)')

        self.salt_bridges_all = find_salt_bridges(
            self.rec.charged_atoms,
            self.lig.charged_atoms,
            parameters.min_dist,
            parameters.saltbridge_dist_max,
        )
        logger.debug(f'Found {len(self.salt_bridges_all)} salt bridge(s)')

        self.water_bridges_all = find_water_bridges(
            self.rec.h_bond_acceptors,
            self.lig.h_bond_donors,
            self.wat.metal_binders,
            parameters.water_bridge_mindist,
            parameters.water_bridge_maxdist,
            parameters.water_bridge_omega_min,
            parameters.water_bridge_omega_max,
            parameters.hbond_acc_angle_min,
            parameters.hbond_don_angle_min,
            False,
        ) + find_water_bridges(
            self.lig.h_bond_acceptors,
            self.rec.h_bond_donors,
            self.wat.metal_binders,
            parameters.water_bridge_mindist,
            parameters.water_bridge_maxdist,
            parameters.water_bridge_omega_min,
            parameters.water_bridge_omega_max,
            parameters.hbond_acc_angle_min,
            parameters.hbond_don_angle_min,
            True,
        )
        logger.debug(f'Found {len(self.water_bridges_all)} water bridge(s)')

        self.metal_complexes = find_metal_complexes(
            self.rec.metals,
            self.rec.metal_binders,
            self.lig.metals,
            self.lig.metal_binders,
            self.wat.metal_binders,
            parameters.metal_dist_max,
        )
        logger.debug(f'Found {len(self.metal_complexes)} metal complex(es)')

        self.interactions_all: Sequence[
            Union[
                H_Bond,
                Hydrophobic,
                Metal_Complex,
                Multipolar,
                Pi_Amide,
                Pi_Cation,
                Pi_Hydrophobic,
                Pi_Stacking,
                Salt_Bridge,
                Water_Bridge,
                Halogen_Bond,
            ]
        ] = (
            self.hydrophobics_all  # type: ignore [operator]
            + self.pi_stackings_all  # type: ignore [operator]
            + self.pi_amides_all  # type: ignore [operator]
            + self.pi_cations_all  # type: ignore [operator]
            + self.pi_hydrophobics_all  # type: ignore [operator]
            + self.h_bonds_all  # type: ignore [operator]
            + self.x_bonds_all  # type: ignore [operator]
            + self.mulipolar_all  # type: ignore [operator]
            + self.salt_bridges_all  # type: ignore [operator]
            + self.metal_complexes  # type: ignore [operator]
            + self.water_bridges_all  # type: ignore [operator]
        )

    def _refine_interactions(self):
        """
        Find all receptor-ligand interactions
        """
        logger.info('Refining selection')
        self.hydrophobics = refine_hydrophobics(
            self.hydrophobics_all, self.pi_stackings_all, self.pi_hydrophobics_all
        )
        logger.debug(f'Kept {len(self.hydrophobics)} hydrophobic interaction(s)')

        self.pi_cations = refine_pi_cations(self.pi_cations_all, self.pi_stackings_all)
        logger.debug(f'Kept {len(self.pi_cations)} pi-cation interaction(s)')

        self.pi_hydrophobics = refine_pi_hydrophobics(self.pi_hydrophobics_all)
        logger.debug(f'Kept {len(self.pi_hydrophobics)} pi-hydrophobic interaction(s)')

        self.h_bonds = refine_h_bonds(
            self.h_bonds_all,
            self.salt_bridges_all,
            self.water_bridges_all,
            self.metal_complexes,
        )
        logger.debug(f'Kept {len(self.h_bonds)} H-bond(s)')

        self.water_bridges = refine_water_bridges(
            self.water_bridges_all, self.metal_complexes
        )
        logger.debug(f'Kept {len(self.water_bridges)} water bridge(s)')

        self.interactions = (
            self.hydrophobics
            + self.pi_stackings_all
            + self.pi_amides_all
            + self.pi_cations
            + self.pi_hydrophobics
            + self.h_bonds
            + self.x_bonds_all
            + self.mulipolar_all
            + self.salt_bridges_all
            + self.metal_complexes
            + self.water_bridges
        )
