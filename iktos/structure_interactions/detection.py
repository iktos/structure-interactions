from __future__ import absolute_import

from itertools import product
from typing import List, NamedTuple

try:
    from iktos.logger import getLogger
except ImportError:
    from logging import getLogger

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
    partner_1: List[Atom]
    partner_2: List[Atom]
    distance: float


def find_hydrophobics(
    hydrophobics_1: List[HydrophobicAtom],
    hydrophobics_2: List[HydrophobicAtom],
    distance_min: float = 2.0,
    distance_max: float = 4.0,
) -> List[Hydrophobic]:
    """Detects hydrophobic interactions.

    Considers hydrophobic atoms (C, S, Cl or F) but excludes interactions
    between aromatic atoms.

    Args:
        hydrophobics_1: list of hydrophobic atoms.
        hydrophobics_2: list of hydrophobic atoms.
        distance_min: distance min (default: 2.0 Ang).
        distance_max: distance max (default: 4.0 Ang).

    Returns:
        list of detected interactions.
    """
    pairings: List[Hydrophobic] = []
    for partner_1, partner_2 in product(hydrophobics_1, hydrophobics_2):
        # Ignore interactions between aromatic atoms
        if partner_1.atom_list[0].is_aromatic and partner_2.atom_list[0].is_aromatic:
            continue
        # Check distance
        dist = get_euclidean_distance_3d(
            partner_1.atom_list[0].coords, partner_2.atom_list[0].coords
        )
        if not distance_min <= dist <= distance_max:
            continue
        # If the same contact has already been seen (only relevant for intra), pass
        if any(
            [
                partner_2.atom_list == pairing.partner_1
                and partner_1.atom_list == pairing.partner_2
                for pairing in pairings
            ]
        ):
            continue
        # If atoms are bound (up to 3 bonds apart, relevant for intra), pass
        ids_1 = [atom.unique_id for atom in partner_1.neighbours_radius_2]
        ids_2 = [atom.unique_id for atom in partner_2.neighbours_radius_2]
        if any([i in ids_2 for i in ids_1]):
            continue
        # Save contact
        contact = Hydrophobic(
            partner_1=partner_1.atom_list, partner_2=partner_2.atom_list, distance=dist
        )
        pairings.append(contact)
    return pairings


class Pi_Stacking(NamedTuple):
    partner_1: List[Atom]
    partner_2: List[Atom]
    distance: float
    angle: float
    offset: float
    type: str


def find_pi_stackings(
    rings_1: List[Ring],
    rings_2: List[Ring],
    distance_min: float = 2.0,
    distance_max_t: float = 5.5,
    distance_max_f: float = 4.75,
    distance_max_p: float = 4.0,
    angle_dev: float = 30,
    offset_max: float = 2.5,
) -> List[Pi_Stacking]:
    """Detects pi-stacking interactions between aromatic rings.

    Args:
        rings_1: list of rings.
        rings_2: list of rings.
        distance_min: distance min (default: 2.0 Ang).
        distance_max_t: distance max for the detection of T-shaped pi-stacking
            (default: 5.5 Ang).
        distance_max_f: distance max for the detection of face-to-face pi-stacking
            (intermediate between T-shaped and parallel, default: 4.75 Ang).
        distance_max_p: distance max for the detection of parallel pi-stacking
            (default: 4.0).
        angle_dev: deviation to optimal angle (the optimal angle differs for P/T/F
            pi-stacking, default: 30 deg).
        offset_max: max offset between the centers (default: 2.5 Ang).

    Returns:
        list of detected interactions.
    """
    pairings: List[Pi_Stacking] = []
    for ring_1, ring_2 in product(rings_1, rings_2):
        # Check distance
        dist = get_euclidean_distance_3d(ring_1.center, ring_2.center)
        if not distance_min <= dist <= distance_max_t:
            continue
        # Take the smallest of two angles, depending on direction of normal
        angle_tmp = get_vector_angle(ring_1.normal, ring_2.normal)
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
        proj1 = project_on_plane(ring_1.normal, ring_1.center, ring_2.center)
        proj2 = project_on_plane(ring_2.normal, ring_2.center, ring_1.center)
        offset = min(
            get_euclidean_distance_3d(proj1, ring_1.center),
            get_euclidean_distance_3d(proj2, ring_2.center),
        )
        if not offset <= offset_max:
            logger.debug(
                f'Pi-staking ignored due to large offset {round(offset, 1)}, '
                f'angle {round(angle, 1)}, distance {round(dist, 1)}'
            )
            continue
        # If the same contact has already been seen (only relevant for intra), pass
        if any(
            [
                ring_2.atom_list == pairing.partner_1
                and ring_1.atom_list == pairing.partner_2
                for pairing in pairings
            ]
        ):
            continue
        # If the rings share any atom (only relevant for intra), pass
        ids_1 = [atom.unique_id for atom in ring_1.atom_list]
        ids_2 = [atom.unique_id for atom in ring_2.atom_list]
        if any([i in ids_2 for i in ids_1]):
            continue
        # Save contact
        contact = Pi_Stacking(
            partner_1=ring_1.atom_list,
            partner_2=ring_2.atom_list,
            distance=dist,
            angle=angle,
            offset=offset,
            type=type,
        )
        pairings.append(contact)
    return pairings


class Pi_Amide(NamedTuple):
    ring: List[Atom]
    pi_carbon: List[Atom]
    distance: float
    angle: float
    offset: float


def find_pi_amides(
    rings: List[Ring],
    pi_carbons: List[PiCarbon],
    distance_min: float = 2.0,
    distance_max: float = 4.0,
    offset_max: float = 2.0,
    angle_dev: float = 30.0,
) -> List[Pi_Amide]:
    """Detects pi-stacking interactions between pi systems
    (aromatic ring + amide/guanidinium/carbamate).

    Args:
        rings: list of rings.
        pi_carbons: list of pi-carbons.
        distance_min: distance min (default: 2.0 Ang).
        distance_max: distance max (default: 4.0 Ang).
        offset_max: max offset between the centers (default: 2.0 Ang).
        angle_dev: deviation to the optimal angle of 180 deg (default: 30 deg).

    Returns:
        list of detected interactions.
    """
    pairings: List[Pi_Amide] = []
    for ring, pi_carbon in product(rings, pi_carbons):
        # Check distance
        dist = get_euclidean_distance_3d(ring.center, pi_carbon.center)
        if not distance_min <= dist <= distance_max:
            continue
        # Take the smallest of two angles, depending on direction of normal
        angle_tmp = get_vector_angle(ring.normal, pi_carbon.normal)
        angle = min(angle_tmp, 180 - angle_tmp)
        if not 0 <= angle <= angle_dev:
            logger.debug(
                f'Pi-amide ignored due to large angle {round(angle, 1)}, '
                f'distance {round(dist, 1)}'
            )
            continue
        # Project each ring center onto the other ring and calculate offset
        proj1 = project_on_plane(pi_carbon.normal, pi_carbon.center, ring.center)
        proj2 = project_on_plane(ring.normal, ring.center, pi_carbon.center)
        offset = min(
            get_euclidean_distance_3d(proj1, pi_carbon.center),
            get_euclidean_distance_3d(proj2, ring.center),
        )
        if not offset <= offset_max:
            logger.debug(
                f'Pi-amide ignored due to large offset {round(offset, 1)}, '
                f'angle {round(angle, 1)}, distance {round(dist, 1)}'
            )
            continue
        # Save contact
        contact = Pi_Amide(
            ring=ring.atom_list,
            pi_carbon=pi_carbon.atom_list,
            distance=dist,
            angle=angle,
            offset=offset,
        )
        pairings.append(contact)
    return pairings


class Pi_Cation(NamedTuple):
    ring: List[Atom]
    cation: List[Atom]
    distance: float
    offset: float


def find_pi_cations(
    rings: List[Ring],
    charged_atoms: List[ChargedGroup],
    distance_min: float = 2.0,
    distance_max: float = 4.0,
    offset_max: float = 2.0,
) -> List[Pi_Cation]:
    """Detects pi-cation interactions between aromatic rings
    and positively charged groups.

    Args:
        rings: list of rings.
        charged_atoms: list of charged groups.
        distance_min: distance min (default: 2.0 Ang).
        distance_max: distance max (default: 4.0 Ang).
        offset_max: max offset between the centers (default: 2.0 Ang).

    Returns:
        list of detected interactions.
    """
    pairings: List[Pi_Cation] = []
    for r, c in product(rings, charged_atoms):
        if c.charge != 'positive':
            continue
        # Check distance
        dist = get_euclidean_distance_3d(r.center, c.center)
        if not distance_min <= dist <= distance_max:
            continue
        # Project the center of charge onto the ring and measure distance to ring center
        proj = project_on_plane(r.normal, r.center, c.center)
        offset = get_euclidean_distance_3d(proj, r.center)
        if not offset <= offset_max:
            logger.debug(
                f'Pi-cation ignored due to large offset {round(offset, 1)}, '
                f'distance {round(dist, 1)}'
            )
            continue
        # Save contact
        contact = Pi_Cation(
            ring=r.atom_list, cation=c.atom_list, distance=dist, offset=offset
        )
        pairings.append(contact)
    return pairings


class Pi_Hydrophobic(NamedTuple):
    ring: List[Atom]
    hydrophobic: List[Atom]
    distance: float
    offset: float


def find_pi_hydrophobics(
    rings: List[Ring],
    hydrophobics: List[HydrophobicAtom],
    distance_min: float = 2.0,
    distance_max: float = 4.0,
    offset_max: float = 2.0,
) -> List[Pi_Hydrophobic]:
    """Detects pi-hydrophobic interactions between aromatic rings
    and hydrophobic atoms (C, S, F, Cl).

    Args:
        rings: list of rings.
        hydrophobics: list of hydrophobic atoms.
        distance_min: distance min (default: 2.0 Ang).
        distance_max: distance max (default: 4.0 Ang).
        offset_max: max offset between the centers (default: 2.0 Ang).

    Returns:
        list of detected interactions.
    """
    pairings: List[Pi_Hydrophobic] = []
    for r, h in product(rings, hydrophobics):
        if h.atom_list[0].hybridisation != 3:
            continue
        # Check distance
        dist = get_euclidean_distance_3d(r.center, h.atom_list[0].coords)
        if not distance_min <= dist <= distance_max:
            continue
        # Project the hydrophobic atom onto ring and measure distance to ring center
        proj = project_on_plane(r.normal, r.center, h.atom_list[0].coords)
        offset = get_euclidean_distance_3d(proj, r.center)
        if not offset <= offset_max:
            logger.debug(
                f'Pi-hydrophobic ignored due to large offset {round(offset, 1)}, '
                f'distance {round(dist, 1)}'
            )
            continue
        # Save contact
        contact = Pi_Hydrophobic(
            hydrophobic=h.atom_list, ring=r.atom_list, distance=dist, offset=offset
        )
        pairings.append(contact)
    return pairings


class H_Bond(NamedTuple):
    donor: List[Atom]
    acceptor: List[Atom]
    distance_ah: float
    distance_ad: float
    angle_dha: float
    type: str


def find_h_bonds(
    acceptors: List[HBondAcceptor],
    donor_pairs: List[HBondDonor],
    distance_min: float = 2.0,
    distance_max: float = 4.0,
    donor_angle_min: float = 140.0,
    acceptor_angle_min: float = 100.0,
) -> List[H_Bond]:
    """Detects H-bonds between acceptors and donor pairs.

    Definition: pairs of (H-bond acceptor and H-bond donor) within dist max,
    DHA angle >= donor_angle_min and each YAD angle >= acceptor_angle_min
    (Y are A's neighbours, this check is to ensure the H-bond is on the lone
    pair side of A).

    Args:
        acceptors: list of acceptors.
        donor_pairs: list of donors.
        distance_min: distance min (default: 2.0 Ang).
        distance_max: distance max (default: 4.0 Ang).
        donor_angle_min: min angle D-H-A (default: 140 deg).
        acceptor_angle_min: min angle Y-A-D (default: 100 deg).

    Returns:
        list of detected interactions.
    """
    # Note: a.atom_list = [A], d.atom_list = [D, H]
    pairings: List[H_Bond] = []
    for a, d in product(acceptors, donor_pairs):
        # Check distance
        dist_ad = get_euclidean_distance_3d(
            a.atom_list[0].coords, d.atom_list[0].coords
        )
        if not distance_min <= dist_ad <= distance_max:
            continue
        # Check angle around H
        vec_hd = get_vector(d.atom_list[1].coords, d.atom_list[0].coords)
        vec_ha = get_vector(d.atom_list[1].coords, a.atom_list[0].coords)
        angle_dha = get_vector_angle(vec_hd, vec_ha)
        if not angle_dha >= donor_angle_min:
            logger.debug(
                f'H-bond ignored due to small D-H--A angle {round(angle_dha, 1)}, '
                f'distance {round(dist_ad, 1)}'
            )
            continue
        # Check acceptor angle
        angle_ok = True
        for y in a.neighbours:
            vec_ay = get_vector(a.atom_list[0].coords, y.coords)
            vec_ad = get_vector(a.atom_list[0].coords, d.atom_list[0].coords)
            angle_yad = get_vector_angle(vec_ay, vec_ad)
            if not angle_yad >= acceptor_angle_min:
                angle_ok = False
                logger.debug(
                    f'H-bond ignored due to small Y-A--D angle {round(angle_yad, 1)}, '
                    f'distance {round(dist_ad, 1)}'
                )
                break
        if not angle_ok:
            continue
        # Save contact
        dist_ah = get_euclidean_distance_3d(
            a.atom_list[0].coords, d.atom_list[1].coords
        )
        contact = H_Bond(
            acceptor=a.atom_list,
            donor=d.atom_list,
            distance_ah=dist_ah,
            distance_ad=dist_ad,
            angle_dha=angle_dha,
            type=d.type,
        )
        pairings.append(contact)
    return pairings


class Halogen_Bond(NamedTuple):
    acceptor: List[Atom]
    donor: List[Atom]
    distance_ax: float
    angle_axd: float


def find_x_bonds(
    acceptors: List[XBondAcceptor],
    donor_pairs: List[XBondDonor],
    distance_min: float = 2.0,
    distance_max: float = 4.0,
    donor_angle_min: float = 140.0,
    acceptor_angle_min: float = 90.0,
) -> List[Halogen_Bond]:
    """Detects halogen bonds between acceptors and donor pairs (excluding F).

    Definition: pairs of (acceptor and C-X pair) within dist max, DXA angle
    >= donor_angle_min and each YAX angle >= acceptor_angle_min (Y are A's neighbours,
    this check is to ensure the X-bond is on the lone pair side of A).
    https://www.pnas.org/content/101/48/16789 fig 3 for acceptor_angle_min

    Args:
        acceptors: list of acceptors.
        donor_pairs: list of donors.
        distance_min: distance min (default: 2.0 Ang).
        distance_max: distance max (default: 4.0 Ang).
        donor_angle_min: min angle D-X-A (default: 140 deg).
        acceptor_angle_min: min angle Y-A-D (default: 90 deg).

    Returns:
        list of detected interactions.
    """
    # Note: a.atom_list = [A], d.atom_list = [D, X]
    pairings: List[Halogen_Bond] = []
    for a, d in product(acceptors, donor_pairs):
        # Exclude F
        if d.atom_list[1].atomic_num == 9:
            continue
        # Check distance
        dist = get_euclidean_distance_3d(a.atom_list[0].coords, d.atom_list[1].coords)
        if not distance_min <= dist <= distance_max:
            continue
        # Check angle around X
        vec_xa = get_vector(d.atom_list[1].coords, a.atom_list[0].coords)
        vec_xd = get_vector(d.atom_list[1].coords, d.atom_list[0].coords)
        angle_axd = get_vector_angle(vec_xa, vec_xd)
        if not angle_axd >= donor_angle_min:
            logger.debug(
                f'X-bond ignored due to small A--X-D angle {round(angle_axd, 1)}, '
                f'distance {round(dist, 1)}'
            )
            continue
        # Check acceptor angles
        angle_ok = True
        for y in a.neighbours:
            vec_ay = get_vector(a.atom_list[0].coords, y.coords)
            vec_ax = get_vector(a.atom_list[0].coords, d.atom_list[1].coords)
            angle_yax = get_vector_angle(vec_ay, vec_ax)
            if not angle_yax >= acceptor_angle_min:
                angle_ok = False
                logger.debug(
                    f'X-bond ignored due to small Y-A--X angle {round(angle_yax, 1)}, '
                    f'distance {round(dist, 1)}'
                )
                break
        if not angle_ok:
            continue
        # Save contact
        contact = Halogen_Bond(
            donor=d.atom_list,
            acceptor=a.atom_list,
            distance_ax=dist,
            angle_axd=angle_axd,
        )
        pairings.append(contact)
    return pairings


class Multipolar(NamedTuple):
    acceptor: List[Atom]
    donor: List[Atom]
    distance_ax: float
    angle_axd: float
    angle_xay: float


def find_multpipolar_interactions(
    acceptors: List[PiCarbon],
    donor_pairs: List[XBondDonor],
    distance_min: float = 2.0,
    distance_max: float = 4.0,
    donor_angle_min: float = 90.0,
    norm_angle_max: float = 40.0,
) -> List[Multipolar]:
    """Detects orthogonal multipolar interactions between F/Cl
    and polarised Csp2 (e.g. amide).

    Definition: pairs of (acceptor and C-X pair) within distmax, DXA angle
    >= multipolar_don_angle_min and angle between normal to amide plan
    and AX = 0 +/- norm_angle_max.

    Args:
        acceptors: list of pi carbons.
        donor_pairs: list of X-bond donors.
        distance_min: distance min (default: 2.0 Ang).
        distance_max: distance max (default: 4.0 Ang).
        donor_angle_min: min angle C-X-A (default: 90 deg).
        norm_angle_min: min angle C-X--amide (default: 40 deg).

    Returns:
        list of detected interactions.
    """
    # Note: a.atom_list = [C], d.atom_list = [D, X]
    pairings: List[Multipolar] = []
    for a, d in product(acceptors, donor_pairs):
        # Exclude Br and I
        if d.atom_list[1].atomic_num not in [9, 17]:
            continue
        # Check distance
        dist = get_euclidean_distance_3d(a.atom_list[0].coords, d.atom_list[1].coords)
        if not distance_min <= dist <= distance_max:
            continue
        # Check angle around X
        vec_xa = get_vector(d.atom_list[1].coords, a.atom_list[0].coords)
        vec_xd = get_vector(d.atom_list[1].coords, d.atom_list[0].coords)
        angle_axd = get_vector_angle(vec_xa, vec_xd)
        if not angle_axd >= donor_angle_min:
            logger.debug(
                f'Multipolar ignored due to small A--X-D angle {round(angle_axd, 1)}, '
                f'distance {round(dist, 1)}'
            )
            continue
        # Take the smallest of two angles, depending on direction of normal
        angle_tmp = get_vector_angle(a.normal, vec_xa)
        angle_xay = min(angle_tmp, 180 - angle_tmp)
        if not 0 <= angle_xay <= norm_angle_max:
            logger.debug(
                f'Multipolar ignored due to large X--A-Y angle {round(angle_xay, 1)}, '
                f'distance {round(dist, 1)}'
            )
            continue
        # Save contact
        contact = Multipolar(
            donor=d.atom_list,
            acceptor=a.atom_list,
            distance_ax=dist,
            angle_axd=angle_axd,
            angle_xay=angle_xay,
        )
        pairings.append(contact)
    return pairings


class Salt_Bridge(NamedTuple):
    partner_1: List[Atom]
    partner_2: List[Atom]
    distance: float


def find_salt_bridges(
    charged_atoms_1: List[ChargedGroup],
    charged_atoms_2: List[ChargedGroup],
    distance_min: float = 2.0,
    distance_max: float = 5.5,
) -> List[Salt_Bridge]:
    """Detects salt bridges, i.e. interaction between charged groups.

    Args:
        charged_atoms_1: list of charged groups.
        charged_atoms_2: list of charged groups.
        distance_min: distance min (default: 2.0 Ang).
        distance_max: distance max (default: 5.5 Ang).

    Returns:
        list of detected interactions.
    """
    pairings: List[Salt_Bridge] = []
    for group_1, group_2 in product(charged_atoms_1, charged_atoms_2):
        if group_1.charge == group_2.charge:
            continue
        # Check distance
        dist = get_euclidean_distance_3d(group_1.center, group_2.center)
        if not distance_min <= dist <= distance_max:
            continue
        # If the same contact has already been seen (only relevant for intra), pass
        if any(
            [
                group_2.atom_list == pairing.partner_1
                and group_1.atom_list == pairing.partner_2
                for pairing in pairings
            ]
        ):
            continue
        # Save contact
        contact = Salt_Bridge(
            partner_1=group_1.atom_list, partner_2=group_2.atom_list, distance=dist
        )
        pairings.append(contact)
    return pairings


class Water_Bridge(NamedTuple):
    donor: List[Atom]
    acceptor: List[Atom]
    water: List[Atom]
    distance_aw: float
    distance_dw: float
    angle_dhw: float
    angle_awh: float


def find_water_bridges(
    acceptors: List[HBondAcceptor],
    donor_pairs: List[HBondDonor],
    waters: List[MetalBinder],
    distance_min: float = 2.0,
    distance_max: float = 4.1,
    omega_min: float = 71,
    omega_max: float = 140,
    hbond_acceptor_angle_min: float = 100,
    hbond_donor_angle_min: float = 140,
) -> List[Water_Bridge]:
    """Detects water-bridged hydrogen bonds between ligand and protein

    Definition: pairs of (H-bond acceptor and H-bond donor pair) within dist max
    of a water molecule, with appropriate AOH and OHD angles, and each YAH angle
    >= hbond_acceptor_angle_min (Y are A's neighbours, this check is to ensure
    the H-bond is on the lone pair side of A).

    Args:
        acceptors: list of acceptors.
        donor_pairs: list of donors.
        distance_min: distance min (default: 2.0 Ang).
        distance_max: distance max between A and OW and between D and OW (default: 4.1 Ang).
        omega_min: min angle between A, OW and D hydrogen (default: 71 deg).
        omega_max: max angle between A, OW and D hydrogen (default: 140 deg).
        hbond_donor_angle_min: min angle D-H-A (default: 140 deg).
        hbond_acceptor_angle_min: min angle Y-A-D (default: 100 deg).

    Returns:
        list of detected interactions.
    """
    # Note: a.atom_list = [A], d.atom_list = [D, H]
    pairings: List[Water_Bridge] = []
    for a, d in product(acceptors, donor_pairs):
        if not d.type == 'strong':
            continue
        for w in waters:
            # Check distances
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
            # Check angle around acceptor
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
                        f'Water bridge ignored due to small Y-A--OW angle {round(angle_yao, 1)}, '
                        f'distance {round(dist_aw, 1)} and {round(dist_dw, 1)}'
                    )
                    break
            if not angle_ok:
                continue
            # Check angle around H
            vec_hd = get_vector(d.atom_list[1].coords, d.atom_list[0].coords)
            vec_hw = get_vector(d.atom_list[1].coords, w.atom_list[0].coords)
            angle_dhw = get_vector_angle(vec_hd, vec_hw)
            if not angle_dhw >= hbond_donor_angle_min:
                logger.debug(
                    f'Water bridge ignored due to small D-H--OW angle {round(angle_dhw, 1)}, '
                    f'distance {round(dist_aw, 1)} and {round(dist_dw, 1)}'
                )
                continue
            # Check angle around OW
            vec_wa = get_vector(w.atom_list[0].coords, a.atom_list[0].coords)
            vec_wh = get_vector(w.atom_list[0].coords, d.atom_list[1].coords)
            angle_awh = get_vector_angle(vec_wa, vec_wh)
            if not (omega_min <= angle_awh <= omega_max):
                logger.debug(
                    f'Water bridge ignored due to invalid A--OW--H angle {round(angle_awh, 1)}, '
                    f'distance {round(dist_aw, 1)} and {round(dist_dw, 1)}'
                )
                continue
            # Save contacts
            contact = Water_Bridge(
                acceptor=a.atom_list,
                donor=d.atom_list,
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
    num_partners: int
    complex_num: int


def find_metal_complexes(
    metals: List[Metal],
    metal_binders: List[MetalBinder],
    distance_max: float = 3.0,
) -> List[Metal_Complex]:
    """Detects metal-atom interactions between any metal (ligand or receptor side)
    and any appropriate group (ligand or receptor side), as well as water.

    Definition: set of L/R--M pairs within distance_max and with a predefined geometry

    Warnings:
        Ignores complexes where the metal is not connected to any ligand atom
            or any receptor atom.
        This function is only used to detect intermolecular interactions, i.e.
            with binders on the ligand and receptor sides.

    Args:
        metals: list of metal atoms.
        metal_binders: list of metal binders.
        distance_max: all binders within that distance will be considered as bound
            to the metal (default: 3.0 Ang).

    Returns:
        list of Metal_Complex with M / L / R atoms in separate lists.
    """
    pairings: List[Metal_Complex] = []
    pairings_dict = {}  # type: dict
    # List possible pairs
    for metal, binder in product(metals, metal_binders):
        dist = get_euclidean_distance_3d(
            metal.atom_list[0].coords,
            binder.atom_list[0].coords,
        )
        if not dist <= distance_max:
            continue
        metal_id = metal.atom_list[0].unique_id
        if metal_id not in pairings_dict:
            pairings_dict[metal_id] = []
        pairings_dict[metal_id].append((metal, binder))
    # Save relevant contacts
    for i, (metal_id, contact_pairs) in enumerate(pairings_dict.items()):
        logger.debug(f'Looking at metal: {metal_id}')
        # If the list of binders doesn't include receptor and ligand, discard
        locations = set([x[1].location for x in contact_pairs])
        if 'ligand' not in locations:
            logger.warning(
                f'--> ignoring metal {metal_id} because it is not connected to the ligand'
            )
            continue
        elif 'receptor' not in locations:
            logger.warning(
                f'--> ignoring metal {metal_id} because it is not connected to the receptor'
            )
            continue
        # Save contact - note: this contact explicitely says ligand and receptor,
        # the function is written to detect intermolecular complexes
        num_binders = len(contact_pairs)
        logger.debug(
            f'--> metal ion {metal_id} complexed by {num_binders} metal binders'
        )
        receptor_binders, ligand_binders = [], []
        for metal, binder in contact_pairs:
            if binder.location != 'ligand':
                receptor_binders.append(binder.atom_list[0])
            else:
                ligand_binders.append(binder.atom_list[0])
        contact = Metal_Complex(
            metal=metal.atom_list,
            ligand=ligand_binders,
            receptor=receptor_binders,
            num_partners=num_binders,
            complex_num=i + 1,
        )
        pairings.append(contact)
    return pairings
