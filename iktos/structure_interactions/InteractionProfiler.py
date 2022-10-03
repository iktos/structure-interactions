from __future__ import absolute_import

from typing import Any, Sequence

try:
    from iktos.logger import getLogger
except ImportError:
    from logging import getLogger

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
from .InteractionParameters import InteractionParameters
from .Ligand import Ligand
from .Receptor import Receptor
from .refinement import (
    refine_h_bonds,
    refine_hydrophobics,
    refine_pi_cations,
    refine_pi_hydrophobics,
    refine_water_bridges,
)


logger = getLogger(__name__)


class InteractionProfiler:
    """Class to analyse protein-ligand interactions using an in-house
    version of PLIP (Protein-Ligand Interaction Profiler).
    """

    def __init__(
        self,
        receptor: Receptor,
        ligand: Ligand,
        parameters: InteractionParameters = InteractionParameters(),
    ):
        """Inintialises the profiler with ligand and/or receptor.

        TODO: as is, the class is specific to intermolecular interactions,
        we want to refacto later to detect intramolecular interaction as well.
        """
        self.parameters = parameters

        if receptor is None and ligand is None:
            raise ValueError('You need to give a ligand or a receptor, or both')

        if receptor is None:
            logger.info('Analysing intramolecular interactions on the ligand')
            self.lig: Ligand = ligand

        elif ligand is None:
            logger.info('Analysing intramolecular interactions on the receptor')
            self.rec: Receptor = receptor
            self.rec.identify_functional_groups(self.rec.atoms)

        else:
            logger.info('Analysing intermolecular interactions')
            self.rec = receptor
            self.lig = ligand

    def detect_intermolecular_interactions(self) -> Sequence[Any]:
        """Finds all receptor-ligand interactions.

        Assumes that a ligand + a receptor have been loaded.
        Calls functions stored in detection.py
        """
        logger.debug('Analysing interactions')
        self.hydrophobics_all = find_hydrophobics(
            hydrophobics_1=self.rec.hydrophobics,
            hydrophobics_2=self.lig.hydrophobics,
            distance_min=self.parameters.min_dist,
            distance_max=self.parameters.hydrophobic_dist_max,
        )
        logger.debug(f'Found {len(self.hydrophobics_all)} hydrophobic interaction(s)')

        self.pi_stackings_all = find_pi_stackings(
            rings_1=self.rec.rings,
            rings_2=self.lig.rings,
            distance_min=self.parameters.min_dist,
            distance_max_t=self.parameters.pistacking_dist_max_t,
            distance_max_f=self.parameters.pistacking_dist_max_f,
            distance_max_p=self.parameters.pistacking_dist_max_p,
            angle_dev=self.parameters.pistacking_ang_dev,
            offset_max=self.parameters.pistacking_offset_max,
        )
        logger.debug(f'Found {len(self.pi_stackings_all)} pi-stacking interaction(s)')

        self.pi_amides_all = find_pi_amides(
            rings=self.lig.rings,
            pi_carbons=self.rec.pi_carbons,
            distance_min=self.parameters.min_dist,
            distance_max=self.parameters.piother_dist_max,
            offset_max=self.parameters.piother_offset_max,
            angle_dev=self.parameters.pistacking_ang_dev,
        ) + find_pi_amides(
            rings=self.rec.rings,
            pi_carbons=self.lig.pi_carbons,
            distance_min=self.parameters.min_dist,
            distance_max=self.parameters.piother_dist_max,
            offset_max=self.parameters.piother_offset_max,
            angle_dev=self.parameters.pistacking_ang_dev,
        )
        logger.debug(f'Found {len(self.pi_amides_all)} pi-amide like interaction(s)')

        self.pi_cations_all = find_pi_cations(
            rings=self.rec.rings,
            charged_atoms=self.lig.charged_atoms,
            distance_min=self.parameters.min_dist,
            distance_max=self.parameters.piother_dist_max,
            offset_max=self.parameters.piother_offset_max,
        ) + find_pi_cations(
            rings=self.lig.rings,
            charged_atoms=self.rec.charged_atoms,
            distance_min=self.parameters.min_dist,
            distance_max=self.parameters.piother_dist_max,
            offset_max=self.parameters.piother_offset_max,
        )
        logger.debug(f'Found {len(self.pi_cations_all)} pi-cation interaction(s)')

        self.pi_hydrophobics_all = find_pi_hydrophobics(
            rings=self.rec.rings,
            hydrophobics=self.lig.hydrophobics,
            distance_min=self.parameters.min_dist,
            distance_max=self.parameters.piother_dist_max,
            offset_max=self.parameters.piother_offset_max,
        ) + find_pi_hydrophobics(
            rings=self.lig.rings,
            hydrophobics=self.rec.hydrophobics,
            distance_min=self.parameters.min_dist,
            distance_max=self.parameters.piother_dist_max,
            offset_max=self.parameters.piother_offset_max,
        )
        logger.debug(
            f'Found {len(self.pi_hydrophobics_all)} pi-hydrophobic interaction(s)'
        )

        self.h_bonds_all = find_h_bonds(
            acceptors=self.rec.h_bond_acceptors + self.rec.water_h_bond_acceptors,
            donor_pairs=self.lig.h_bond_donors,
            distance_min=self.parameters.min_dist,
            distance_max=self.parameters.hbond_dist_max,
            donor_angle_min=self.parameters.hbond_don_angle_min,
            acceptor_angle_min=self.parameters.hbond_acc_angle_min,
        ) + find_h_bonds(
            acceptors=self.lig.h_bond_acceptors,
            donor_pairs=self.rec.h_bond_donors + self.rec.water_h_bond_donors,
            distance_min=self.parameters.min_dist,
            distance_max=self.parameters.hbond_dist_max,
            donor_angle_min=self.parameters.hbond_don_angle_min,
            acceptor_angle_min=self.parameters.hbond_acc_angle_min,
        )
        logger.debug(f'Found {len(self.h_bonds_all)} H-bond(s)')

        self.x_bonds_all = find_x_bonds(
            acceptors=self.rec.x_bond_acceptors + self.rec.water_x_bond_acceptors,
            donor_pairs=self.lig.halogens,
            distance_min=self.parameters.min_dist,
            distance_max=self.parameters.xbond_dist_max,
            donor_angle_min=self.parameters.xbond_don_angle_min,
            acceptor_angle_min=self.parameters.xbond_acc_angle_min,
        )
        logger.debug(f'Found {len(self.x_bonds_all)} halogen bond(s)')

        self.mulipolar_all = find_multpipolar_interactions(
            acceptors=self.rec.pi_carbons,
            donor_pairs=self.lig.halogens,
            distance_min=self.parameters.min_dist,
            distance_max=self.parameters.multipolar_dist_max,
            donor_angle_min=self.parameters.multipolar_don_angle_min,
            norm_angle_max=self.parameters.multipolar_norm_angle_max,
        )
        logger.debug(f'Found {len(self.mulipolar_all)} multipolar interaction(s)')

        self.salt_bridges_all = find_salt_bridges(
            charged_atoms_1=self.rec.charged_atoms,
            charged_atoms_2=self.lig.charged_atoms,
            distance_min=self.parameters.min_dist,
            distance_max=self.parameters.saltbridge_dist_max,
        )
        logger.debug(f'Found {len(self.salt_bridges_all)} salt bridge(s)')

        self.water_bridges_all = find_water_bridges(
            acceptors=self.rec.h_bond_acceptors,
            donor_pairs=self.lig.h_bond_donors,
            waters=self.rec.water_metal_binders,
            distance_min=self.parameters.water_bridge_mindist,
            distance_max=self.parameters.water_bridge_maxdist,
            omega_min=self.parameters.water_bridge_omega_min,
            omega_max=self.parameters.water_bridge_omega_max,
            hbond_acceptor_angle_min=self.parameters.hbond_acc_angle_min,
            hbond_donor_angle_min=self.parameters.hbond_don_angle_min,
        ) + find_water_bridges(
            acceptors=self.lig.h_bond_acceptors,
            donor_pairs=self.rec.h_bond_donors,
            waters=self.rec.water_metal_binders,
            distance_min=self.parameters.water_bridge_mindist,
            distance_max=self.parameters.water_bridge_maxdist,
            omega_min=self.parameters.water_bridge_omega_min,
            omega_max=self.parameters.water_bridge_omega_max,
            hbond_acceptor_angle_min=self.parameters.hbond_acc_angle_min,
            hbond_donor_angle_min=self.parameters.hbond_don_angle_min,
        )
        logger.debug(f'Found {len(self.water_bridges_all)} water bridge(s)')

        self.metal_complexes = find_metal_complexes(
            metals=self.rec.metals + self.lig.metals,
            metal_binders=self.rec.metal_binders
            + self.rec.water_metal_binders
            + self.lig.metal_binders,
            distance_max=self.parameters.metal_dist_max,
        )
        logger.debug(f'Found {len(self.metal_complexes)} metal complex(es)')

        interactions_all: Sequence[Any] = (
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
        return interactions_all

    def refine_intermolecular_interactions(self) -> Sequence[Any]:
        """Refines receptor-ligand interactions.

        Assumes that `detect_intermolecular_interactions` has already benn called.
        """
        logger.debug('Refining selection')
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

        interactions: Sequence[Any] = (
            self.hydrophobics  # type: ignore [operator]
            + self.pi_stackings_all  # type: ignore [operator]
            + self.pi_amides_all  # type: ignore [operator]
            + self.pi_cations  # type: ignore [operator]
            + self.pi_hydrophobics  # type: ignore [operator]
            + self.h_bonds  # type: ignore [operator]
            + self.x_bonds_all  # type: ignore [operator]
            + self.mulipolar_all  # type: ignore [operator]
            + self.salt_bridges_all  # type: ignore [operator]
            + self.metal_complexes  # type: ignore [operator]
            + self.water_bridges  # type: ignore [operator]
        )
        return interactions
