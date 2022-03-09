from __future__ import absolute_import

import importlib
from itertools import product
from typing import Dict, List

from iktos.logger import getLogger

try:
    from openbabel.openbabel import OBResidueAtomIter, OBResidueIter  # openbabel 3
except ModuleNotFoundError:
    from openbabel import (
        OBResidueIter,
        OBResidueAtomIter,
    )  # openbabel 2 (warning in structure-utils)

import numpy as np
from iktos.structure_utils.file_manager.utils import file_to_blocks

from . import constants
from .detection import (
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
from .Ligand import Ligand
from .math_utils import get_centroid, get_euclidean_distance_3d
from .Receptor import Receptor
from .refinement import (
    refine_h_bonds,
    refine_hydrophobics,
    refine_pi_cations,
    refine_pi_hydrophobics,
    refine_water_bridges,
)
from .utils import get_coords, map_atom_ids, parse_pdb, read_mol
from .Water import Water

logger = getLogger(__name__)


class InteractionProfiler:
    """
    Analyse protein-ligand interactions using an in-house
    version of PLIP (Protein-Ligand Interaction Profiler)
    """

    def __init__(self, plip_config_version='default'):
        self.plip_config_version = plip_config_version

    def analyse_complex(
        self,
        rec_coords: str,
        lig_coords: str,
        lig_format: str = 'sdf',
        as_string: bool = True,
        refine: bool = True,
    ):
        """Performs protein-ligand interaction analysis for 1 complex.

        Note: The current version of the code fails to detect negatively charged groups
        (e.g. COO-) if ligand coords are given in MOL2 format.

        Args:
            rec_coords: name of receptor coords file or coords block (PDB format)
            lig_coords: name of ligand multi coords file or coords block
            lig_format: format of lig_coords (recommended format is 'sdf')
            as_string: whether to read coords in files or from blocks
            refine: if True, cleanup double-counted contacts (e.g. hydrophobics + pi-stacking)

        Returns:
            contacts in a formatted dict
        """
        status = self._load_receptor(rec_coords, as_string=as_string)
        if not status:
            return False
        status = self._load_ligand(
            lig_coords, lig_format=lig_format, as_string=as_string
        )
        if not status:
            return False
        return self._analyse_interactions(refine=refine)

    def analyse_complexes(
        self,
        rec_coords: str,
        lig_coords: [List, str],
        lig_format: str = 'sdf',
        as_string: bool = True,
        refine: bool = True,
    ):
        """Performs protein-ligand interaction analysis
        for multiple complexes that share the same recpetor.

        Note: The current version of the code fails to detect negatively charged groups
        (e.g. COO-) if ligand coords are given in MOL2 format.

        Args:
            rec_coords: name of receptor coords file or coords block (PDB format)
            lig_coords: name of ligand multi coords file or list of coords blocks
            lig_format: format of lig_coords (recommended format is 'sdf')
            as_string: whether to read coords in files or from blocks
            refine: if True, cleanup double-counted contacts (e.g. hydrophobics + pi-stacking)

        Returns:
            list of contacts in a formatted dict
        """
        status = self._load_receptor(rec_coords, as_string=as_string)
        if not status:
            return False
        contacts = []
        if as_string:
            lig_blocks = lig_coords
        else:
            lig_blocks = file_to_blocks(lig_coords, fmt=lig_format)
        for lig in lig_blocks:
            status = self._load_ligand(lig, lig_format=lig_format, as_string=as_string)
            if not status:
                contacts.append({})
            else:
                contacts.append(self._analyse_interactions(refine=refine))
        return contacts

    def _plip_alternate_config(self):
        constants.HBOND_DIST_MAX = 3.5
        constants.HBOND_DON_ANGLE_MIN = 150
        constants.XBOND_DIST_MAX = 3.5

    def _load_receptor(self, rec_coords: str, as_string: bool) -> bool:
        """Loads the receptor and initialises water object."""
        # Read and parse receptor file/string
        rec_coords_block, mapping = parse_pdb(rec_coords, as_string=as_string)
        obmol_rec = read_mol(
            rec_coords_block, fmt='pdb', title='receptor', as_string=True
        )
        if obmol_rec is None:
            return False

        # Map atom Id back to their original value
        self.obmol_rec = map_atom_ids(obmol_rec, mapping)

        # Initialise water and receptor objects
        logger.info('Initialising water object')
        self.wat = Water(self.obmol_rec)
        logger.info('Initialising receptor object')
        self.rec = Receptor(self.obmol_rec)
        return True

    def _load_ligand(self, lig_coords: str, lig_format: str, as_string: bool) -> bool:
        """Loads the ligand and initialises ligand and binding site objects."""
        if lig_format != 'sdf':
            logger.warning(
                'It is recommended to use SDF blocks/files for the ligand, '
                'some interactions might be missed otherwise'
            )

        # Read and parse ligand file/string
        obmol_lig = read_mol(
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
        bs_atoms_refined, water_atoms = self._extract_binding_site()
        self.rec.identify_functional_groups(bs_atoms_refined)
        self.wat.identify_functional_groups(water_atoms)
        return True

    def _analyse_interactions(self, refine: bool):
        """Performs protein-ligand interaction analysis
        for 1 complex (after initialisation of the objects)."""
        # Setup PLIP config
        importlib.reload(constants)
        if self.plip_config_version == 'custom':
            self._plip_alternate_config()

        # Analyse contacts
        self._detect_interactions()

        # Refine contacts
        if refine:
            self._refine_interactions()

        # Parse results
        if refine:
            contacts = self._parse_interactions(self.interactions)
            return contacts
        contacts = self._parse_interactions(self.interactions_all)
        return contacts

    def _extract_binding_site(self):
        """
        Select atoms that belong to the binding site
        """
        # Select binding site residues
        cutoff = self.lig.max_dist_to_center + constants.BS_DIST
        obres_rec = [
            obres
            for obres in OBResidueIter(self.obmol_rec)
            # if not obres.GetResidueProperty(9)
        ]
        bs_residues = []
        for obres in obres_rec:
            residue_centroid = get_centroid(
                [get_coords(a) for a in OBResidueAtomIter(obres)]
            )
            distance = get_euclidean_distance_3d(residue_centroid, self.lig.centroid)
            if distance < cutoff:
                bs_residues.append(obres)
        bs_atoms = [
            a
            for r in bs_residues
            for a in OBResidueAtomIter(r)
            if not r.GetResidueProperty(9)
        ]
        water_atoms = [
            a
            for r in bs_residues
            for a in OBResidueAtomIter(r)
            if r.GetResidueProperty(9)
        ]
        logger.debug(f'Found {len(bs_atoms)} atoms near ligand')

        # Extracting residue atoms coordinates:
        coords_r = np.array([get_coords(obatom_r) for obatom_r in bs_atoms])
        # Extracting ligand atoms coordinates:
        coords_l = np.array([get_coords(obatom_l) for obatom_l in self.lig.atoms])
        # We want to calculate the distance for each unique pair of (coords_r, coords_l)
        product_rl = np.array(list(product(coords_r, coords_l)))
        # Calculating the distances and keeping only the indexes for those being below our criterion:
        idx = np.argwhere(
            get_euclidean_distance_3d(product_rl[:, 0], product_rl[:, 1]).reshape(
                coords_r.shape[0], coords_l.shape[0]
            )
            < constants.BS_DIST
        )

        # Selecting the first atoms for each coords_r which respect the distance criterion:
        bs_atoms_refined = [bs_atoms[idx[0, 0]]]
        for (i_p, _), (i, _) in zip(idx[:-1], idx[1:]):
            if i_p != i:
                bs_atoms_refined.append(bs_atoms[i])

        logger.debug(f'Selected {len(bs_atoms_refined)} atoms as binding site')
        return bs_atoms_refined, water_atoms

    def _detect_interactions(self):
        """
        Find all receptor-ligand interactions
        Call functions stroed in detection.py
        """
        logger.info('Analysing interactions')
        self.hydrophobics_all = find_hydrophobics(
            self.rec.hydrophobics, self.lig.hydrophobics
        )
        logger.debug(f'Found {len(self.hydrophobics_all)} hydrophobic interaction(s)')

        self.pi_stackings_all = find_pi_stackings(self.rec.rings, self.lig.rings)
        logger.debug(f'Found {len(self.pi_stackings_all)} pi-stacking interaction(s)')

        self.pi_amides_all = find_pi_amides(
            self.rec.pi_carbons, self.lig.rings
        ) + find_pi_amides(self.rec.rings, self.lig.pi_carbons)
        logger.debug(f'Found {len(self.pi_amides_all)} pi-amide like interaction(s)')

        self.pi_cations_all = find_pi_cations(
            self.rec.rings, self.lig.charged_atoms, True
        ) + find_pi_cations(self.lig.rings, self.rec.charged_atoms, False)
        logger.debug(f'Found {len(self.pi_cations_all)} pi-cation interaction(s)')

        self.pi_hydrophobics_all = find_pi_hydrophobics(
            self.rec.rings, self.lig.hydrophobics, True
        ) + find_pi_hydrophobics(self.lig.rings, self.rec.hydrophobics, False)
        logger.debug(
            f'Found {len(self.pi_hydrophobics_all)} pi-hydrophobic interaction(s)'
        )

        self.h_bonds_all = (
            find_h_bonds(self.rec.h_bond_acceptors, self.lig.h_bond_donors, False)
            + find_h_bonds(self.wat.h_bond_acceptors, self.lig.h_bond_donors, False)
            + find_h_bonds(self.lig.h_bond_acceptors, self.rec.h_bond_donors, True)
            + find_h_bonds(self.lig.h_bond_acceptors, self.wat.h_bond_donors, True)
        )
        logger.debug(f'Found {len(self.h_bonds_all)} H-bond(s)')

        self.x_bonds_all = find_x_bonds(
            self.rec.x_bond_acceptors, self.lig.halogens
        ) + find_x_bonds(self.lig.x_bond_acceptors, self.rec.halogens)
        logger.debug(f'Found {len(self.x_bonds_all)} halogen bond(s)')

        self.mulipolar_all = find_multpipolar_interactions(
            self.rec.pi_carbons, self.lig.halogens
        )
        logger.debug(f'Found {len(self.mulipolar_all)} multipolar interaction(s)')

        self.salt_bridges_all = find_salt_bridges(
            self.rec.charged_atoms, self.lig.charged_atoms
        )
        logger.debug(f'Found {len(self.salt_bridges_all)} salt bridge(s)')

        self.water_bridges_all = find_water_bridges(
            self.rec.h_bond_acceptors,
            self.lig.h_bond_donors,
            self.wat.metal_binders,
            False,
        ) + find_water_bridges(
            self.lig.h_bond_acceptors,
            self.rec.h_bond_donors,
            self.wat.metal_binders,
            True,
        )
        logger.debug(f'Found {len(self.water_bridges_all)} water bridge(s)')

        self.metal_complexes = find_metal_complexes(
            self.rec.metals,
            self.rec.metal_binders,
            self.lig.metals,
            self.lig.metal_binders,
            self.wat.metal_binders,
        )
        logger.debug(f'Found {len(self.metal_complexes)} metal complex(es)')

        self.interactions_all = (
            self.hydrophobics_all
            + self.pi_stackings_all
            + self.pi_amides_all
            + self.pi_cations_all
            + self.pi_hydrophobics_all
            + self.h_bonds_all
            + self.x_bonds_all
            + self.mulipolar_all
            + self.salt_bridges_all
            + self.metal_complexes
            + self.water_bridges_all
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

    def _parse_interactions(self, plip_list: List) -> Dict:
        """
        Parse contacts and return them as a list of dicts
        formatted for further analysis/drawing
        """
        plip_dict = {}
        for contact in plip_list:
            contact_name = type(contact).__name__
            if contact_name == 'H_Bond':
                if contact.type == 'weak':
                    contact_name += '_Weak'
            if contact_name not in plip_dict:
                plip_dict[contact_name] = []
            parsed_contact = {}
            for field in contact._fields:
                if field == 'receptor':
                    atom_list = getattr(contact, field)
                    # Drop Hs
                    # atom_list = [a for a in atom_list if a.atomic_num != 1]
                    parsed_contact['at_p'] = [a.atom_id for a in atom_list]
                    parsed_contact['at_name_p'] = [a.atom_name for a in atom_list]
                    if contact_name != 'Metal_Complex':
                        parsed_contact['res_p'] = atom_list[0].residue_id
                elif field == 'ligand':
                    atom_list = getattr(contact, field)
                    parsed_contact['at_l'] = [a.atom_id for a in atom_list]
                    parsed_contact['at_name_l'] = [a.atom_name for a in atom_list]
                elif field == 'water':
                    atom_list = getattr(contact, field)
                    parsed_contact['at_owat'] = [a.atom_id for a in atom_list]
                    parsed_contact['res_w'] = atom_list[0].residue_id
                elif field == 'metal':
                    atom_list = getattr(contact, field)
                    parsed_contact['at_m'] = [a.atom_id for a in atom_list]
                    parsed_contact['at_name_m'] = [a.atom_name for a in atom_list]
                    parsed_contact['res_p'] = atom_list[0].residue_id
                else:
                    parsed_contact[field] = getattr(contact, field)
            plip_dict[contact_name].append(parsed_contact)
        return plip_dict
