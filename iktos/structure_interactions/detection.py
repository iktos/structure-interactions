from __future__ import absolute_import
import numpy as np
from collections import defaultdict, namedtuple
from itertools import product
from typing import List, NamedTuple

from iktos.logger import getLogger

from . import constants
from .math_utils import (
    get_vector,
    get_vector_angle,
    get_euclidean_distance_3d,
    project_on_plane,
)


logger = getLogger(__name__)


def find_hydrophobics(
    hydrophobics_rec: List[NamedTuple], hydrophobics_lig: List[NamedTuple]
) -> List[NamedTuple]:
    """Detects hydrophobic interactions between
    hydrophobic atoms (C, S, Cl or F), excluding interactions between aromatic atoms.

    Definition: pairs of atoms within HYDROPHOBIC_DIST_MAX
    """
    # data = namedtuple('Hydrophobic', 'ligand receptor distance')

    class data(NamedTuple):
        name: str
        receptor: str
        distance: float

    pairings = []
    for rec, lig in product(hydrophobics_rec, hydrophobics_lig):
        if rec.atom_list[0].is_aromatic and lig.atom_list[0].is_aromatic:
            continue
        dist = get_euclidean_distance_3d(
            rec.atom_list[0].coords, lig.atom_list[0].coords
        )
        if not constants.MIN_DIST <= dist <= constants.HYDROPHOBIC_DIST_MAX:
            continue
        contact = data(ligand=lig.atom_list, receptor=rec.atom_list, distance=dist)
        pairings.append(contact)
    return pairings


def find_pi_stackings(
    groups_rec: List[NamedTuple], groups_lig: List[NamedTuple]
) -> List[NamedTuple]:
    """Detects pi-stacking interactions between aromatic rings.

    Definiton: pairs of rings within PISTACKING_DIST_MAX,
        either // or |-- or |/, and offset <= PISTACKING_OFFSET_MAX
    """
    data = namedtuple('Pi_Stacking', 'ligand receptor distance angle offset type')
    pairings = []
    for rec, lig in product(groups_rec, groups_lig):
        dist = get_euclidean_distance_3d(rec.center, lig.center)
        if not constants.MIN_DIST <= dist <= constants.PISTACKING_DIST_MAX_T:
            continue

        # Take the smallest of two angles, depending on direction of normal
        angle_tmp = get_vector_angle(rec.normal, lig.normal)
        angle = min(angle_tmp, 180 - angle_tmp)
        if (
            0 <= angle <= constants.PISTACKING_ANG_DEV
            and dist <= constants.PISTACKING_DIST_MAX_P
        ):
            type = 'P'
        elif angle >= 90 - constants.PISTACKING_ANG_DEV:
            type = 'T'
        elif dist <= constants.PISTACKING_DIST_MAX_F:
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
        if not offset <= constants.PISTACKING_OFFSET_MAX:
            logger.debug(
                f'Pi-staking ignored due to large offset {offset}, '
                f'angle {angle}, distance {dist}'
            )
            continue
        contact = data(
            ligand=lig.atom_list,
            receptor=rec.atom_list,
            distance=dist,
            angle=angle,
            offset=offset,
            type=type,
        )
        pairings.append(contact)
    return pairings


def find_pi_amides(
    groups_rec: List[NamedTuple], groups_lig: List[NamedTuple]
) -> List[NamedTuple]:
    """Detects pi-stacking interactions between pi systems
    (aromatic ring + amide/guanidinium/carbamate).

    Definiton: pairs of (ring, pi-group) within PIOTHER_DIST_MAX,
        parallel (angle between normals approx 0), and offset <= PIOTHER_OFFSET_MAX
    """
    data = namedtuple('Pi_Amide', 'ligand receptor distance angle offset type')
    pairings = []
    for rec, lig in product(groups_rec, groups_lig):
        dist = get_euclidean_distance_3d(rec.center, lig.center)
        if not constants.MIN_DIST <= dist <= constants.PIOTHER_DIST_MAX:
            continue

        # Take the smallest of two angles, depending on direction of normal
        angle_tmp = get_vector_angle(rec.normal, lig.normal)
        angle = min(angle_tmp, 180 - angle_tmp)
        if not 0 <= angle <= constants.PISTACKING_ANG_DEV:
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
        if not offset <= constants.PIOTHER_OFFSET_MAX:
            logger.debug(
                f'Pi-amide ignored due to large offset {offset}, '
                f'angle {angle}, distance {dist}'
            )
            continue
        contact = data(
            ligand=lig.atom_list,
            receptor=rec.atom_list,
            distance=dist,
            angle=angle,
            offset=offset,
            type='P',
        )
        pairings.append(contact)
    return pairings


def find_pi_cations(
    rings: List[NamedTuple], charged_atoms: List[NamedTuple], pi_on_prot: bool = False
) -> List[NamedTuple]:
    """Detects pi-cation interactions between aromatic rings
    and positively charged groups.

    Definition: pairs of (ring, charged atom)
        within PIOTHER_DIST_MAX and offset < PIOTHER_OFFSET_MAX
    """
    data = namedtuple('Pi_Cation', 'ligand receptor distance offset')
    pairings = []
    for r, c in product(rings, charged_atoms):
        if c.charge != 'positive':
            continue
        dist = get_euclidean_distance_3d(r.center, c.center)
        if not constants.MIN_DIST <= dist <= constants.PIOTHER_DIST_MAX:
            continue

        # Project the center of charge onto the ring and measure distance to ring center
        proj = project_on_plane(r.normal, r.center, c.center)
        offset = get_euclidean_distance_3d(proj, r.center)
        if not offset <= constants.PIOTHER_OFFSET_MAX:
            logger.debug(
                f'Pi-cation ignored due to large offset {offset}, ' f'distance {dist}'
            )
            continue
        if pi_on_prot:
            contact = data(
                ligand=c.atom_list, receptor=r.atom_list, distance=dist, offset=offset
            )
        else:
            contact = data(
                ligand=r.atom_list, receptor=c.atom_list, distance=dist, offset=offset
            )
        pairings.append(contact)
    return pairings


def find_pi_hydrophobics(
    rings: List[NamedTuple], hydrophobics: List[NamedTuple], pi_on_prot: bool = False
) -> List[NamedTuple]:
    """Detects pi-hydrophobic interactions between aromatic rings
    and hydrophobic atoms (C, S, F, Cl).

    Definition: pairs of (ring, SP3 hydrophobic atom)
        within PIOTHER_DIST_MAX and offset < PIOTHER_OFFSET_MAX
    """
    data = namedtuple('Pi_Hydrophobic', 'ligand receptor distance offset')
    pairings = []
    for r, h in product(rings, hydrophobics):
        if h.atom_list[0].hybridisation != 3:
            continue
        dist = get_euclidean_distance_3d(r.center, h.atom_list[0].coords)
        if not constants.MIN_DIST <= dist <= constants.PIOTHER_DIST_MAX:
            continue

        # Project the hydrophobic atom onto ring and measure distance to ring center
        proj = project_on_plane(r.normal, r.center, h.atom_list[0].coords)
        offset = get_euclidean_distance_3d(proj, r.center)
        if not offset <= constants.PIOTHER_OFFSET_MAX:
            logger.debug(
                f'Pi-hydrophobic ignored due to large offset {offset}, '
                f'distance {dist}'
            )
            continue
        if pi_on_prot:
            contact = data(
                ligand=h.atom_list, receptor=r.atom_list, distance=dist, offset=offset
            )
        else:
            contact = data(
                ligand=r.atom_list, receptor=h.atom_list, distance=dist, offset=offset
            )
        pairings.append(contact)
    return pairings


def find_h_bonds(
    acceptors: List[NamedTuple], donor_pairs: List[NamedTuple], don_on_prot: bool = True
) -> List[NamedTuple]:
    """Detects H-bonds between acceptors and donor pairs.

    Definition: pairs of (H-bond acceptor and H-bond donor) within
        HBOND_DIST_MAX, DHA angle >= HBOND_DON_ANGLE_MIN
        and each YAD angle >= HBOND_ACC_ANGLE_MIN (Y are A's neighbours,
        this check is to ensure the H-bond is on the lone pair side of A)

    Note: a.atom_list = [A], d.atom_list = [D, H]
    """
    data = namedtuple(
        'H_Bond', 'ligand receptor distance_ah distance_ad angle_dha type'
    )
    pairings = []
    for a, d in product(acceptors, donor_pairs):
        dist_ad = get_euclidean_distance_3d(
            a.atom_list[0].coords, d.atom_list[0].coords
        )
        if not constants.MIN_DIST <= dist_ad <= constants.HBOND_DIST_MAX:
            continue
        dist_ah = get_euclidean_distance_3d(
            a.atom_list[0].coords, d.atom_list[1].coords
        )
        vec_hd = get_vector(d.atom_list[1].coords, d.atom_list[0].coords)
        vec_ha = get_vector(d.atom_list[1].coords, a.atom_list[0].coords)
        angle_dha = get_vector_angle(vec_hd, vec_ha)
        if not angle_dha >= constants.HBOND_DON_ANGLE_MIN:
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
            if not angle_yad >= constants.HBOND_ACC_ANGLE_MIN:
                angle_ok = False
                logger.debug(
                    f'H-bond ignored due to small Y-A--D angle {angle_yad}, '
                    f'distance {dist_ad}'
                )
                break
        if not angle_ok:
            continue
        if don_on_prot:
            l, r = a, d
        else:
            l, r = d, a
        contact = data(
            ligand=l.atom_list,
            receptor=r.atom_list,
            distance_ah=dist_ah,
            distance_ad=dist_ad,
            angle_dha=angle_dha,
            type=d.type,
        )
        pairings.append(contact)
    return pairings


def find_x_bonds(
    acceptors: List[NamedTuple], donor_pairs: List[NamedTuple]
) -> List[NamedTuple]:
    """Detects halogen bonds between acceptors and donor pairs (excluding F).

    Definition: pairs of (acceptor and C-X pair) within
        HBOND_DIST_MAX, DXA angle >= XBOND_DON_ANGLE_MIN
        and each YAX angle >= XBOND_ACC_ANGLE_MIN (Y are A's neighbours,
        this check is to ensure the X-bond is on the lone pair side of A)

    https://www.pnas.org/content/101/48/16789 fig 3 for XBOND_ACC_ANGLE_MIN

    Note: a.atom_list = [A], d.atom_list = [D, X]
    """
    data = namedtuple('Halogen_Bond', 'ligand receptor distance_ax angle_axd')
    pairings = []
    for a, d in product(acceptors, donor_pairs):
        # Exclude F
        if d.atom_list[1].atomic_num == 9:
            continue
        dist = get_euclidean_distance_3d(a.atom_list[0].coords, d.atom_list[1].coords)
        if not constants.MIN_DIST <= dist <= constants.XBOND_DIST_MAX:
            continue
        vec_xa = get_vector(d.atom_list[1].coords, a.atom_list[0].coords)
        vec_xd = get_vector(d.atom_list[1].coords, d.atom_list[0].coords)
        angle_axd = get_vector_angle(vec_xa, vec_xd)
        if not angle_axd >= constants.XBOND_DON_ANGLE_MIN:
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
            if not angle_yax >= constants.XBOND_ACC_ANGLE_MIN:
                angle_ok = False
                logger.debug(
                    f'X-bond ignored due to small Y-A--X angle {angle_yax}, '
                    f'distance {dist}'
                )
                break
        if not angle_ok:
            continue
        contact = data(
            ligand=d.atom_list,
            receptor=a.atom_list,
            distance_ax=dist,
            angle_axd=angle_axd,
        )
        pairings.append(contact)
    return pairings


def find_multpipolar_interactions(
    acceptors: List[NamedTuple], donor_pairs: List[NamedTuple]
) -> List[NamedTuple]:
    """Detects orthogonal multipolar interactions between F/Cl
    and polarised Csp2 (e.g. amide).

    Definition: pairs of (acceptor and C-X pair) within
        MULTIPOLAR_DIST_MAX, DXA angle >= MULTIPOLAR_DON_ANGLE_MIN
        and angle between normal to amide plan and AX = 0 +/- MULTIPOLAR_NORM_ANGLE_MAX

    Note: a.atom_list = [C], d.atom_list = [D, X]
    """
    data = namedtuple('Multipolar', 'ligand receptor distance_ax angle_axd angle_xay')
    pairings = []
    for a, d in product(acceptors, donor_pairs):
        # Exclude Br and I
        if d.atom_list[1].atomic_num not in [9, 17]:
            continue
        dist = get_euclidean_distance_3d(a.atom_list[0].coords, d.atom_list[1].coords)
        if not constants.MIN_DIST <= dist <= constants.MULTIPOLAR_DIST_MAX:
            continue
        vec_xa = get_vector(d.atom_list[1].coords, a.atom_list[0].coords)
        vec_xd = get_vector(d.atom_list[1].coords, d.atom_list[0].coords)
        angle_axd = get_vector_angle(vec_xa, vec_xd)
        if not angle_axd >= constants.MULTIPOLAR_DON_ANGLE_MIN:
            logger.debug(
                f'Multipolar ignored due to small A--X-D angle {angle_axd}, '
                f'distance {dist}'
            )
            continue

        # Take the smallest of two angles, depending on direction of normal
        angle_tmp = get_vector_angle(a.normal, vec_xa)
        angle_xay = min(angle_tmp, 180 - angle_tmp)
        if not 0 <= angle_xay <= constants.MULTIPOLAR_NORM_ANGLE_MAX:
            logger.debug(
                f'Multipolar ignored due to large X--A-Y angle {angle_xay}, '
                f'distance {dist}'
            )
            continue

        contact = data(
            ligand=d.atom_list,
            receptor=a.atom_list,
            distance_ax=dist,
            angle_axd=angle_axd,
            angle_xay=angle_xay,
        )
        pairings.append(contact)
    return pairings


def find_salt_bridges(
    charged_atoms_rec: List[NamedTuple], charged_atoms_lig: List[NamedTuple]
) -> List[NamedTuple]:
    """Detects salt bridges, i.e. interaction between positively charged
    and negatively charged groups.

    Definition: pairs of charged groups/atoms within SALTBRIDGE_DIST_MAX
    """
    data = namedtuple('Salt_Bridge', 'ligand receptor distance')
    pairings = []
    for group_rec, group_lig in product(charged_atoms_rec, charged_atoms_lig):
        if group_rec.charge == group_lig.charge:
            continue
        dist = get_euclidean_distance_3d(group_rec.center, group_lig.center)
        if not constants.MIN_DIST <= dist <= constants.SALTBRIDGE_DIST_MAX:
            continue
        contact = data(
            ligand=group_lig.atom_list, receptor=group_rec.atom_list, distance=dist
        )
        pairings.append(contact)
    return pairings


def find_water_bridges(
    acceptors: List[NamedTuple],
    donor_pairs: List[NamedTuple],
    waters: List[NamedTuple],
    don_on_prot: bool = True,
) -> List[NamedTuple]:
    """Detects water-bridged hydrogen bonds between ligand and protein

    Definition: pairs of (H-bond acceptor and H-bond donor pair)
        within WATER_BRIDGE_MAXDIST of a water molecule,
        with appropriate AOH and OHD angles,
        and each YAH angle >= HBOND_ACC_ANGLE_MIN (Y are A's neighbours,
        this check is to ensure the H-bond is on the lone pair side of A)

    Note: a.atom_list = [A], d.atom_list = [D, H]
    """
    data = namedtuple(
        'Water_Bridge',
        'ligand receptor water distance_aw distance_dw angle_dhw angle_awh',
    )
    pairings = []
    for a, d in product(acceptors, donor_pairs):
        if not d.type == 'strong':
            continue
        for w in waters:
            dist_aw = get_euclidean_distance_3d(
                a.atom_list[0].coords, w.atom_list[0].coords
            )
            if (
                not constants.WATER_BRIDGE_MINDIST
                <= dist_aw
                <= constants.WATER_BRIDGE_MAXDIST
            ):
                continue
            dist_dw = get_euclidean_distance_3d(
                d.atom_list[0].coords, w.atom_list[0].coords
            )
            if (
                not constants.WATER_BRIDGE_MINDIST
                <= dist_dw
                <= constants.WATER_BRIDGE_MAXDIST
            ):
                continue
            angle_ok = True
            for y in a.neighbours:
                vec_ay = get_vector(a.atom_list[0].coords, y.coords)
                vec_ao = get_vector(a.atom_list[0].coords, w.atom_list[0].coords)
                angle_yao = get_vector_angle(
                    vec_ay, vec_ao
                )  # angle with O_wat instead of H_wat
                if not angle_yao >= constants.HBOND_ACC_ANGLE_MIN:
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
            if not angle_dhw >= constants.HBOND_DON_ANGLE_MIN:
                logger.debug(
                    f'Water bridge ignored due to small D-H--OW angle {angle_dhw}, '
                    f'distance {dist_aw} and {dist_dw}'
                )
                continue
            vec_wa = get_vector(w.atom_list[0].coords, a.atom_list[0].coords)
            vec_wh = get_vector(w.atom_list[0].coords, d.atom_list[1].coords)
            angle_awh = get_vector_angle(vec_wa, vec_wh)
            if not (
                constants.WATER_BRIDGE_OMEGA_MIN
                <= angle_awh
                <= constants.WATER_BRIDGE_OMEGA_MAX
            ):
                logger.debug(
                    f'Water bridge ignored due to invalid A--OW--H angle {angle_awh}, '
                    f'distance {dist_aw} and {dist_dw}'
                )
                continue
            if don_on_prot:
                l, r = a, d
            else:
                l, r = d, a
            contact = data(
                ligand=l.atom_list,
                receptor=r.atom_list,
                water=w.atom_list,
                distance_aw=dist_aw,
                distance_dw=dist_dw,
                angle_dhw=angle_dhw,
                angle_awh=angle_awh,
            )
            pairings.append(contact)
    return pairings


def find_metal_complexes(  # noqa: C901
    metals_rec: List[NamedTuple],
    metal_binders_rec: List[NamedTuple],
    metals_lig: List[NamedTuple],
    metal_binders_lig: List[NamedTuple],
    metal_binders_wat: List[NamedTuple],
) -> List[NamedTuple]:
    """Detects metal-atom interactions between any metal (ligand or receptor side)
    and any appropriate group (ligand or receptor side), as well as water.

    Definition: set of L/R--M pairs within METAL_DIST_MAX and with a predefined geometry

    Returns:
        1 set per detected complex, with M / L / R atoms in separate lists

    TODO: refacto
    """
    data = namedtuple(
        'Metal_Complex',
        'metal ligand receptor coordination_num rms geometry num_partners complex_num',
    )
    pairings = []
    metals = metals_lig + metals_rec
    metal_binders = metal_binders_lig + metal_binders_rec + metal_binders_wat

    pairings_dict = {}
    for m, b in product(metals, metal_binders):
        dist = get_euclidean_distance_3d(m.atom_list[0].coords, b.atom_list[0].coords)
        if not dist <= constants.METAL_DIST_MAX:
            continue
        m_id = m.atom_list[0].unique_id
        if m_id not in pairings_dict:
            pairings_dict[m_id] = []
        pairings_dict[m_id].append((m, b, dist))

    for cnum, m_id in enumerate(pairings_dict):
        logger.debug(f'Looking at metal complex {cnum + 1}')
        rms = 0.0
        excluded = []
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
                    [other_vectors.append(x) for x in vectors_dict[t]]
            angles = [
                get_vector_angle(pair[0], pair[1])
                for pair in product(cur_vector, other_vectors)
            ]
            angles_dict[b_idx] = angles

        # Record fit information for each geometry tested
        all_total = []
        gdata = namedtuple('gdata', 'geometry rms coordination excluded diff_binders')

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
                        best_score = 999
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
                                    sum([x ** 2 for x in single_binder_scores]) ** 0.5
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
                        not_used.append(b_id)
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
                    final_geom, final_coo = 'NA', 'NA'
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
        contact = data(
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
