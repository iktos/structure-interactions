from __future__ import absolute_import

from collections import defaultdict
from itertools import product
from typing import List, NamedTuple, Union

try:
    from iktos.logger import getLogger
except ImportError:
    from logging import getLogger

import numpy as np

from . import constants
from .Atom import Atom
from .atom_typing import (
    ChargedGroup,
    HBondAcceptor,
    HBondDonor,
    HydrophobicAtom,
    Metal,
    MetalBinder,
    PiCarbon,
    Ring,
    XBondAcceptor,
    XBondDonor,
)
from .math_utils import (
    get_euclidean_distance_3d,
    get_vector,
    get_vector_angle,
    project_on_plane,
)

logger = getLogger(__name__)


class Hydrophobic(NamedTuple):
    ligand: List[Atom]
    receptor: List[Atom]
    distance: float


def find_hydrophobics(
    hydrophobics_rec: List[HydrophobicAtom],
    hydrophobics_lig: List[HydrophobicAtom],
    distance_min: float,
    distance_max: float,
) -> List[Hydrophobic]:
    """Detects hydrophobic interactions between
    hydrophobic atoms (C, S, Cl or F), excluding interactions between aromatic atoms.

    Definition: pairs of atoms within distance_max
    """

    pairings = []
    for rec, lig in product(hydrophobics_rec, hydrophobics_lig):
        if rec.atom_list[0].is_aromatic and lig.atom_list[0].is_aromatic:
            continue
        dist = get_euclidean_distance_3d(
            rec.atom_list[0].coords, lig.atom_list[0].coords
        )
        if not distance_min <= dist <= distance_max:
            continue
        contact = Hydrophobic(
            ligand=lig.atom_list, receptor=rec.atom_list, distance=dist
        )
        pairings.append(contact)
    return pairings


class Pi_Stacking(NamedTuple):
    ligand: List[Atom]
    receptor: List[Atom]
    distance: float
    angle: float
    offset: float
    type: str


def find_pi_stackings(
    groups_rec: List[Ring],
    groups_lig: List[Ring],
    distance_min: float,
    distance_max_t: float,
    distance_max_f: float,
    distance_max_p: float,
    angle_dev: float,
    offset_max: float,
) -> List[Pi_Stacking]:
    """Detects pi-stacking interactions between aromatic rings.

    Definiton: pairs of rings within distance_max,
        either // or |-- or |/, and offset <= offset_max
    """

    pairings = []
    for rec, lig in product(groups_rec, groups_lig):
        dist = get_euclidean_distance_3d(rec.center, lig.center)
        if not distance_min <= dist <= distance_max_t:
            continue

        # Take the smallest of two angles, depending on direction of normal
        angle_tmp = get_vector_angle(rec.normal, lig.normal)
        angle = min(angle_tmp, 180 - angle_tmp)
        if 0 <= angle <= angle_dev and dist <= distance_max_p:
            type = 'P'
        elif angle >= 90 - angle_dev:
            type = 'T'
        elif dist <= distance_max_f:
            type = 'F'
        else:
            continue

        # Project each ring center onto the other ring and calculate offset
        proj1 = project_on_plane(lig.normal, lig.center, rec.center)
        proj2 = project_on_plane(rec.normal, rec.center, lig.center)
        offset = min(
            get_euclidean_distance_3d(proj1, lig.center),
            get_euclidean_distance_3d(proj2, rec.center),
        )
        if not offset <= offset_max:
            logger.debug(
                f'Pi-staking ignored due to large offset {offset}, '
                f'angle {angle}, distance {dist}'
            )
            continue
        contact = Pi_Stacking(
            ligand=lig.atom_list,
            receptor=rec.atom_list,
            distance=dist,
            angle=angle,
            offset=offset,
            type=type,
        )
        pairings.append(contact)
    return pairings


class Pi_Amide(NamedTuple):
    ligand: List[Atom]
    receptor: List[Atom]
    distance: float
    angle: float
    offset: float
    type: str


def find_pi_amides(
    groups_rec: Union[List[Ring], List[PiCarbon]],
    groups_lig: Union[List[Ring], List[PiCarbon]],
    distance_min: float,
    distance_max: float,
    offset_max: float,
    angle_dev: float,
) -> List[Pi_Amide]:
    """Detects pi-stacking interactions between pi systems
    (aromatic ring + amide/guanidinium/carbamate).

    Definiton: pairs of (ring, pi-group) within distance_max,
        parallel (angle between normals approx 0), and offset <= offset_max
    """

    pairings = []
    for rec, lig in product(groups_rec, groups_lig):
        dist = get_euclidean_distance_3d(rec.center, lig.center)
        if not distance_min <= dist <= distance_max:
            continue

        # Take the smallest of two angles, depending on direction of normal
        angle_tmp = get_vector_angle(rec.normal, lig.normal)
        angle = min(angle_tmp, 180 - angle_tmp)
        if not 0 <= angle <= angle_dev:
            logger.debug(
                f'Pi-amide ignored due to large angle {angle}, distance {dist}'
            )
            continue

        # Project each ring center onto the other ring and calculate offset
        proj1 = project_on_plane(lig.normal, lig.center, rec.center)
        proj2 = project_on_plane(rec.normal, rec.center, lig.center)
        offset = min(
            get_euclidean_distance_3d(proj1, lig.center),
            get_euclidean_distance_3d(proj2, rec.center),
        )
        if not offset <= offset_max:
            logger.debug(
                f'Pi-amide ignored due to large offset {offset}, '
                f'angle {angle}, distance {dist}'
            )
            continue
        contact = Pi_Amide(
            ligand=lig.atom_list,
            receptor=rec.atom_list,
            distance=dist,
            angle=angle,
            offset=offset,
            type='P',
        )
        pairings.append(contact)
    return pairings


class Pi_Cation(NamedTuple):
    ligand: List[Atom]
    receptor: List[Atom]
    distance: float
    offset: float


def find_pi_cations(
    rings: List[Ring],
    charged_atoms: List[ChargedGroup],
    distance_min: float,
    distance_max: float,
    offset_max: float,
    pi_on_prot: bool = False,
) -> List[Pi_Cation]:
    """Detects pi-cation interactions between aromatic rings
    and positively charged groups.

    Definition: pairs of (ring, charged atom)
        within piother_dist_max and offset < offset_max
    """

    pairings = []
    for r, c in product(rings, charged_atoms):
        if c.charge != 'positive':
            continue
        dist = get_euclidean_distance_3d(r.center, c.center)
        if not distance_min <= dist <= distance_max:
            continue

        # Project the center of charge onto the ring and measure distance to ring center
        proj = project_on_plane(r.normal, r.center, c.center)
        offset = get_euclidean_distance_3d(proj, r.center)
        if not offset <= offset_max:
            logger.debug(
                f'Pi-cation ignored due to large offset {offset}, ' f'distance {dist}'
            )
            continue
        if pi_on_prot:
            contact = Pi_Cation(
                ligand=c.atom_list, receptor=r.atom_list, distance=dist, offset=offset
            )
        else:
            contact = Pi_Cation(
                ligand=r.atom_list, receptor=c.atom_list, distance=dist, offset=offset
            )
        pairings.append(contact)
    return pairings


class Pi_Hydrophobic(NamedTuple):
    ligand: List[Atom]
    receptor: List[Atom]
    distance: float
    offset: float


def find_pi_hydrophobics(
    rings: List[Ring],
    hydrophobics: List[HydrophobicAtom],
    distance_min: float,
    distance_max: float,
    offset_max: float,
    pi_on_prot: bool = False,
) -> List[Pi_Hydrophobic]:
    """Detects pi-hydrophobic interactions between aromatic rings
    and hydrophobic atoms (C, S, F, Cl).

    Definition: pairs of (ring, SP3 hydrophobic atom)
        within piother_dist_max and offset < offset_max
    """

    pairings = []
    for r, h in product(rings, hydrophobics):
        if h.atom_list[0].hybridisation != 3:
            continue
        dist = get_euclidean_distance_3d(r.center, h.atom_list[0].coords)
        if not distance_min <= dist <= distance_max:
            continue

        # Project the hydrophobic atom onto ring and measure distance to ring center
        proj = project_on_plane(r.normal, r.center, h.atom_list[0].coords)
        offset = get_euclidean_distance_3d(proj, r.center)
        if not offset <= offset_max:
            logger.debug(
                f'Pi-hydrophobic ignored due to large offset {offset}, '
                f'distance {dist}'
            )
            continue
        if pi_on_prot:
            contact = Pi_Hydrophobic(
                ligand=h.atom_list, receptor=r.atom_list, distance=dist, offset=offset
            )
        else:
            contact = Pi_Hydrophobic(
                ligand=r.atom_list, receptor=h.atom_list, distance=dist, offset=offset
            )
        pairings.append(contact)
    return pairings


class H_Bond(NamedTuple):
    ligand: List[Atom]
    receptor: List[Atom]
    distance_ah: float
    distance_ad: float
    angle_dha: float
    type: str


def find_h_bonds(
    acceptors: List[HBondAcceptor],
    donor_pairs: List[HBondDonor],
    distance_min: float,
    distance_max: float,
    donor_angle_min: float,
    acceptor_angle_min: float,
    donor_on_prot: bool = True,
) -> List[H_Bond]:
    """Detects H-bonds between acceptors and donor pairs.

    Definition: pairs of (H-bond acceptor and H-bond donor) within
        hbond_dist_max, DHA angle >= donor_angle_min
        and each YAD angle >= acceptor_angle_min (Y are A's neighbours,
        this check is to ensure the H-bond is on the lone pair side of A)

    Note: a.atom_list = [A], d.atom_list = [D, H]
    """

    pairings = []
    for a, d in product(acceptors, donor_pairs):
        dist_ad = get_euclidean_distance_3d(
            a.atom_list[0].coords, d.atom_list[0].coords
        )
        if not distance_min <= dist_ad <= distance_max:
            continue
        dist_ah = get_euclidean_distance_3d(
            a.atom_list[0].coords, d.atom_list[1].coords
        )
        vec_hd = get_vector(d.atom_list[1].coords, d.atom_list[0].coords)
        vec_ha = get_vector(d.atom_list[1].coords, a.atom_list[0].coords)
        angle_dha = get_vector_angle(vec_hd, vec_ha)
        if not angle_dha >= donor_angle_min:
            logger.debug(
                f'H-bond ignored due to small D-H--A angle {angle_dha}, '
                f'distance {dist_ad}'
            )
            continue
        angle_ok = True
        for y in a.neighbours:
            vec_ay = get_vector(a.atom_list[0].coords, y.coords)
            vec_ad = get_vector(a.atom_list[0].coords, d.atom_list[0].coords)
            angle_yad = get_vector_angle(vec_ay, vec_ad)
            if not angle_yad >= acceptor_angle_min:
                angle_ok = False
                logger.debug(
                    f'H-bond ignored due to small Y-A--D angle {angle_yad}, '
                    f'distance {dist_ad}'
                )
                break
        if not angle_ok:
            continue
        if donor_on_prot:
            contact = H_Bond(
                ligand=a.atom_list,
                receptor=d.atom_list,
                distance_ah=dist_ah,
                distance_ad=dist_ad,
                angle_dha=angle_dha,
                type=d.type,
            )
        else:
            contact = H_Bond(
                ligand=d.atom_list,
                receptor=a.atom_list,
                distance_ah=dist_ah,
                distance_ad=dist_ad,
                angle_dha=angle_dha,
                type=d.type,
            )
        pairings.append(contact)
    return pairings


class Halogen_Bond(NamedTuple):
    ligand: List[Atom]
    receptor: List[Atom]
    distance_ax: float
    angle_axd: float


def find_x_bonds(
    acceptors: List[XBondAcceptor],
    donor_pairs: List[XBondDonor],
    distance_min: float,
    distance_max: float,
    donor_angle_min: float,
    acceptor_angle_min: float,
) -> List[Halogen_Bond]:
    """Detects halogen bonds between acceptors and donor pairs (excluding F).

    Definition: pairs of (acceptor and C-X pair) within
        xbond_dist_max, DXA angle >= donor_angle_min
        and each YAX angle >= acceptor_angle_min (Y are A's neighbours,
        this check is to ensure the X-bond is on the lone pair side of A)

    https://www.pnas.org/content/101/48/16789 fig 3 for acceptor_angle_min

    Note: a.atom_list = [A], d.atom_list = [D, X]
    """

    pairings = []
    for a, d in product(acceptors, donor_pairs):
        # Exclude F
        if d.atom_list[1].atomic_num == 9:
            continue
        dist = get_euclidean_distance_3d(a.atom_list[0].coords, d.atom_list[1].coords)
        if not distance_min <= dist <= distance_max:
            continue
        vec_xa = get_vector(d.atom_list[1].coords, a.atom_list[0].coords)
        vec_xd = get_vector(d.atom_list[1].coords, d.atom_list[0].coords)
        angle_axd = get_vector_angle(vec_xa, vec_xd)
        if not angle_axd >= donor_angle_min:
            logger.debug(
                f'X-bond ignored due to small A--X-D angle {angle_axd}, '
                f'distance {dist}'
            )
            continue
        angle_ok = True
        for y in a.neighbours:
            vec_ay = get_vector(a.atom_list[0].coords, y.coords)
            vec_ax = get_vector(a.atom_list[0].coords, d.atom_list[1].coords)
            angle_yax = get_vector_angle(vec_ay, vec_ax)
            if not angle_yax >= acceptor_angle_min:
                angle_ok = False
                logger.debug(
                    f'X-bond ignored due to small Y-A--X angle {angle_yax}, '
                    f'distance {dist}'
                )
                break
        if not angle_ok:
            continue
        contact = Halogen_Bond(
            ligand=d.atom_list,
            receptor=a.atom_list,
            distance_ax=dist,
            angle_axd=angle_axd,
        )
        pairings.append(contact)
    return pairings


class Multipolar(NamedTuple):
    ligand: List[Atom]
    receptor: List[Atom]
    distance_ax: float
    angle_axd: float
    angle_xay: float


def find_multpipolar_interactions(
    acceptors: List[PiCarbon],
    donor_pairs: List[XBondDonor],
    distance_min: float,
    distance_max: float,
    donor_angle_min: float,
    norm_angle_max: float,
) -> List[Multipolar]:
    """Detects orthogonal multipolar interactions between F/Cl
    and polarised Csp2 (e.g. amide).

    Definition: pairs of (acceptor and C-X pair) within
        multipolar_dist_max, DXA angle >= multipolar_don_angle_min
        and angle between normal to amide plan and AX = 0 +/- norm_angle_max

    Note: a.atom_list = [C], d.atom_list = [D, X]
    """

    pairings = []
    for a, d in product(acceptors, donor_pairs):
        # Exclude Br and I
        if d.atom_list[1].atomic_num not in [9, 17]:
            continue
        dist = get_euclidean_distance_3d(a.atom_list[0].coords, d.atom_list[1].coords)
        if not distance_min <= dist <= distance_max:
            continue
        vec_xa = get_vector(d.atom_list[1].coords, a.atom_list[0].coords)
        vec_xd = get_vector(d.atom_list[1].coords, d.atom_list[0].coords)
        angle_axd = get_vector_angle(vec_xa, vec_xd)
        if not angle_axd >= donor_angle_min:
            logger.debug(
                f'Multipolar ignored due to small A--X-D angle {angle_axd}, '
                f'distance {dist}'
            )
            continue

        # Take the smallest of two angles, depending on direction of normal
        angle_tmp = get_vector_angle(a.normal, vec_xa)
        angle_xay = min(angle_tmp, 180 - angle_tmp)
        if not 0 <= angle_xay <= norm_angle_max:
            logger.debug(
                f'Multipolar ignored due to large X--A-Y angle {angle_xay}, '
                f'distance {dist}'
            )
            continue

        contact = Multipolar(
            ligand=d.atom_list,
            receptor=a.atom_list,
            distance_ax=dist,
            angle_axd=angle_axd,
            angle_xay=angle_xay,
        )
        pairings.append(contact)
    return pairings


class Salt_Bridge(NamedTuple):
    ligand: List[Atom]
    receptor: List[Atom]
    distance: float


def find_salt_bridges(
    charged_atoms_rec: List[ChargedGroup],
    charged_atoms_lig: List[ChargedGroup],
    distance_min: float,
    distance_max: float,
) -> List[Salt_Bridge]:
    """Detects salt bridges, i.e. interaction between positively charged
    and negatively charged groups.

    Definition: pairs of charged groups/atoms within saltbridge_dist_max
    """

    pairings = []
    for group_rec, group_lig in product(charged_atoms_rec, charged_atoms_lig):
        if group_rec.charge == group_lig.charge:
            continue
        dist = get_euclidean_distance_3d(group_rec.center, group_lig.center)
        if not distance_min <= dist <= distance_max:
            continue
        contact = Salt_Bridge(
            ligand=group_lig.atom_list, receptor=group_rec.atom_list, distance=dist
        )
        pairings.append(contact)
    return pairings


class Water_Bridge(NamedTuple):
    ligand: List[Atom]
    receptor: List[Atom]
    water: List[Atom]
    distance_aw: float
    distance_dw: float
    angle_dhw: float
    angle_awh: float


def find_water_bridges(
    acceptors: List[HBondAcceptor],
    donor_pairs: List[HBondDonor],
    waters: List,
    distance_min: float,
    distance_max: float,
    omega_min: float,
    omega_max: float,
    hbond_acceptor_angle_min: float,
    hbond_donor_angle_min: float,
    don_on_prot: bool = True,
) -> List[Water_Bridge]:
    """Detects water-bridged hydrogen bonds between ligand and protein

    Definition: pairs of (H-bond acceptor and H-bond donor pair)
        within distance_max of a water molecule,
        with appropriate AOH and OHD angles,
        and each YAH angle >= hbond_acceptor_angle_min (Y are A's neighbours,
        this check is to ensure the H-bond is on the lone pair side of A)

    Note: a.atom_list = [A], d.atom_list = [D, H]
    """
    pairings = []
    for a, d in product(acceptors, donor_pairs):
        if not d.type == 'strong':
            continue
        for w in waters:
            dist_aw = get_euclidean_distance_3d(
                a.atom_list[0].coords, w.atom_list[0].coords
            )
            if not distance_min <= dist_aw <= distance_max:
                continue
            dist_dw = get_euclidean_distance_3d(
                d.atom_list[0].coords, w.atom_list[0].coords
            )
            if not distance_min <= dist_dw <= distance_max:
                continue
            angle_ok = True
            for y in a.neighbours:
                vec_ay = get_vector(a.atom_list[0].coords, y.coords)
                vec_ao = get_vector(a.atom_list[0].coords, w.atom_list[0].coords)
                angle_yao = get_vector_angle(
                    vec_ay, vec_ao
                )  # angle with O_wat instead of H_wat
                if not angle_yao >= hbond_acceptor_angle_min:
                    angle_ok = False
                    logger.debug(
                        f'Water bridge ignored due to small Y-A--OW angle {angle_yao}, '
                        f'distance {dist_aw} and {dist_dw}'
                    )
                    break
            if not angle_ok:
                continue
            vec_hd = get_vector(d.atom_list[1].coords, d.atom_list[0].coords)
            vec_hw = get_vector(d.atom_list[1].coords, w.atom_list[0].coords)
            angle_dhw = get_vector_angle(vec_hd, vec_hw)
            if not angle_dhw >= hbond_donor_angle_min:
                logger.debug(
                    f'Water bridge ignored due to small D-H--OW angle {angle_dhw}, '
                    f'distance {dist_aw} and {dist_dw}'
                )
                continue
            vec_wa = get_vector(w.atom_list[0].coords, a.atom_list[0].coords)
            vec_wh = get_vector(w.atom_list[0].coords, d.atom_list[1].coords)
            angle_awh = get_vector_angle(vec_wa, vec_wh)
            if not (omega_min <= angle_awh <= omega_max):
                logger.debug(
                    f'Water bridge ignored due to invalid A--OW--H angle {angle_awh}, '
                    f'distance {dist_aw} and {dist_dw}'
                )
                continue
            if don_on_prot:
                contact = Water_Bridge(
                    ligand=a.atom_list,
                    receptor=d.atom_list,
                    water=w.atom_list,
                    distance_aw=dist_aw,
                    distance_dw=dist_dw,
                    angle_dhw=angle_dhw,
                    angle_awh=angle_awh,
                )
            else:
                contact = Water_Bridge(
                    ligand=d.atom_list,
                    receptor=a.atom_list,
                    water=w.atom_list,
                    distance_aw=dist_aw,
                    distance_dw=dist_dw,
                    angle_dhw=angle_dhw,
                    angle_awh=angle_awh,
                )
            pairings.append(contact)
    return pairings


class Metal_Complex(NamedTuple):
    metal: List[Atom]
    ligand: List[Atom]
    receptor: List[Atom]
    coordination_num: int
    rms: float
    geometry: str
    num_partners: int
    complex_num: int


def find_metal_complexes(  # noqa: C901
    metals_rec: List[Metal],
    metal_binders_rec: List[MetalBinder],
    metals_lig: List[Metal],
    metal_binders_lig: List[MetalBinder],
    metal_binders_wat: List[MetalBinder],
    distance_max: float,
) -> List[Metal_Complex]:
    """Detects metal-atom interactions between any metal (ligand or receptor side)
    and any appropriate group (ligand or receptor side), as well as water.

    Definition: set of L/R--M pairs within distance_max and with a predefined geometry

    Returns:
        1 set per detected complex, with M / L / R atoms in separate lists

    TODO: refacto
    """

    pairings = []
    metals = metals_lig + metals_rec
    metal_binders = metal_binders_lig + metal_binders_rec + metal_binders_wat

    pairings_dict = {}  # type: dict
    for m, b in product(metals, metal_binders):
        dist = get_euclidean_distance_3d(m.atom_list[0].coords, b.atom_list[0].coords)
        if not dist <= distance_max:
            continue
        m_id = m.atom_list[0].unique_id
        if m_id not in pairings_dict:
            pairings_dict[m_id] = []
        pairings_dict[m_id].append((m, b, dist))

    for cnum, m_id in enumerate(pairings_dict):
        logger.debug(f'Looking at metal complex {cnum + 1}')
        rms = 0.0
        excluded = []  # type: list
        contact_pairs = pairings_dict[m_id]

        # Cannot specify geometry if only one binder
        num_binders = len(contact_pairs)
        if num_binders == 1:
            final_geom = 'NA'
            final_coo = 1
            excluded = []
            rms = 0.0

        vectors_dict = defaultdict(list)
        for contact_pair in contact_pairs:
            m, b, dist = contact_pair
            b_idx = b.atom_list[0].unique_id
            vectors_dict[b_idx].append(
                get_vector(m.atom_list[0].coords, b.atom_list[0].coords)
            )

        angles_dict = {}
        for b_idx in vectors_dict:
            cur_vector = vectors_dict[b_idx]
            other_vectors = []
            for t in vectors_dict:
                if not t == b_idx:
                    [other_vectors.append(x) for x in vectors_dict[t]]  # type: ignore
            angles = [
                get_vector_angle(pair[0], pair[1])
                for pair in product(cur_vector, other_vectors)
            ]
            angles_dict[b_idx] = angles

        # Record fit information for each geometry tested
        all_total = []

        class gdata(NamedTuple):
            geometry: str
            rms: float
            coordination: int
            excluded: List
            diff_binders: int

        if num_binders > 1:
            # Start with highest coordination number and loop over
            # possible coordination numbers; score each combination
            for coo in sorted(constants.METAL_COMPLEX_COO, reverse=True):
                geometries = constants.METAL_COMPLEX_COO[coo]
                for geometry in geometries:

                    # Set ideal angles for geometry, from each perspective
                    signature = constants.METAL_COMPLEX_ANG[geometry]
                    geometry_total = 0
                    geometry_scores = []  # all scores for 1 geom
                    used_up_binders = []
                    not_used = []

                    # How many more observed binders are there?
                    coo_diff = num_binders - coo

                    # Find best match for each subsignature
                    for subsignature in signature:
                        # Ideal angles from one perspective
                        best_binder = None
                        # There's one best-matching binder for each subsignature
                        best_score = 999.0
                        for k, b_idx in enumerate(angles_dict):
                            if b_idx not in used_up_binders:
                                # Observed angles from perspective of 1 binder
                                observed_angles = angles_dict[b_idx]
                                single_binder_scores = []
                                used_up_observed_angles = []
                                for i, ideal_angle in enumerate(subsignature):
                                    # For each angle in the signature,
                                    # find the best-matching observed angle
                                    best_match = None
                                    best_match_diff = 999
                                    for j, observed_angle in enumerate(observed_angles):
                                        if j not in used_up_observed_angles:
                                            diff = abs(ideal_angle - observed_angle)
                                            if diff < best_match_diff:
                                                best_match_diff = diff
                                                best_match = j
                                    if best_match is not None:
                                        used_up_observed_angles.append(best_match)
                                        single_binder_scores.append(best_match_diff)
                                # Calculate RMS for binder angles
                                score = (
                                    sum([x**2 for x in single_binder_scores]) ** 0.5
                                )
                                if score < best_score:
                                    best_score = score
                                    best_binder = b_idx
                        used_up_binders.append(best_binder)
                        geometry_scores.append(best_score)
                        # Total score = mean of RMS values
                        geometry_total = np.mean(geometry_scores)
                    # Record binders not used for excluding
                    # them when deciding for a final geometry
                    [
                        not_used.append(b_id)  # type: ignore
                        for b_id in angles_dict
                        if b_id not in used_up_binders
                    ]
                    all_total.append(
                        gdata(
                            geometry=geometry,
                            rms=geometry_total,
                            coordination=coo,
                            excluded=not_used,
                            diff_binders=coo_diff,
                        )
                    )
            # Choose the complex geometry that best fits the contacts observed
            # Start with the geometry closer to ideal (num of partners close
            # to num of ideal partners), then check that the next
            # best solution is not too close in RMS (diff > 0.5)
            all_total = sorted(all_total, key=lambda x: abs(x.diff_binders))
            for i, total in enumerate(all_total):
                next_total = all_total[i + 1]
                this_rms, next_rms = total.rms, next_total.rms
                diff_to_next = next_rms - this_rms
                if diff_to_next > 0.5:
                    final_geom = total.geometry
                    final_coo = total.coordination
                    rms = total.rms
                    excluded = total.excluded
                    break
                elif next_total.rms < 3.5:
                    final_geom = next_total.geometry
                    final_coo = next_total.coordination
                    rms = next_total.rms
                    excluded = next_total.excluded
                    break
                elif i == len(all_total) - 2:
                    final_geom, final_coo = 'NA', 'NA'  # type: ignore
                    rms, excluded = float('nan'), []
                    break

        # Ignore complex if metal binded to water/receptor only
        if set([x[1].location for x in contact_pairs]) == {'water', 'receptor'}:
            continue
        logger.debug(
            f'--> metal ion {m_id} complexed with {final_geom} geometry '
            f'(coordination number {final_coo}/{num_binders} observed)'
        )

        # Save data for complex
        receptor_binders = []
        ligand_binders = []
        for contact_pair in contact_pairs:
            m, b, dist = contact_pair
            b_id = b.atom_list[0].unique_id
            if b_id in excluded:
                continue
            if b.location != 'ligand':
                receptor_binders.append(b.atom_list[0])
            else:
                ligand_binders.append(b.atom_list[0])
        contact = Metal_Complex(
            metal=m.atom_list,
            ligand=ligand_binders,
            receptor=receptor_binders,
            coordination_num=final_coo,
            rms=rms,
            num_partners=num_binders,
            complex_num=cnum + 1,
            geometry=final_geom,
        )
        pairings.append(contact)
    return pairings
