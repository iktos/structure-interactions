from __future__ import absolute_import

from itertools import combinations, product
from typing import List

from .detection import (
    H_Bond,
    Hydrophobic,
    Metal_Complex,
    Pi_Cation,
    Pi_Hydrophobic,
    Pi_Stacking,
    Salt_Bridge,
    Water_Bridge,
)

try:
    from iktos.logger import getLogger
except ImportError:
    from logging import getLogger

logger = getLogger(__name__)


def refine_hydrophobics(
    hydrophobics_all: List[Hydrophobic],
    pi_stackings_all: List[Pi_Stacking],
    pi_hydrophobics_all: List[Pi_Hydrophobic],
) -> List[Hydrophobic]:
    """Refines hydrophobic interactions to avoid double counts.

    Removes hydrophobic contacts that either involve atoms that are also involved
    in a pi-stacking or pi-hydrophobic interaction, or involve an atom (ligand or receptor)
    that is also involved in another hydrophobic interaction (in that case, the contact with
    the shortest distance will be kept).
    """
    # Remove hydrophobic contacts for atoms that also interact via stacking
    filter_out: List[Hydrophobic] = []
    for h, s in product(hydrophobics_all, pi_stackings_all):
        h_partners = h.partner_1 + h.partner_2
        s_partners = s.partner_1 + s.partner_2
        h_rec_ids = {i.unique_id for i in h_partners if "receptor" in i.unique_id}
        s_rec_ids = {i.unique_id for i in s_partners if "receptor" in i.unique_id}
        if h_rec_ids.isdisjoint(s_rec_ids):
            continue
        h_lig_ids = {i.unique_id for i in h_partners if "ligand" in i.unique_id}
        s_lig_ids = {i.unique_id for i in s_partners if "ligand" in i.unique_id}
        if h_lig_ids.isdisjoint(s_lig_ids):
            continue
        logger.debug(
            f"Removing hydrophobic: {(h.partner_1[0].unique_id, h.partner_2[0].unique_id)}"
        )
        filter_out.append(h)
    selection1 = [contact for contact in hydrophobics_all if contact not in filter_out]

    # Remove hydrophobic contacts for atoms that also interact via pi-hydrophobic
    filter_out = []
    for h, p in product(selection1, pi_hydrophobics_all):
        h_partners = h.partner_1 + h.partner_2
        p_partners = p.ring + p.hydrophobic
        h_rec_ids = {i.unique_id for i in h_partners if "receptor" in i.unique_id}
        p_rec_ids = {i.unique_id for i in p_partners if "receptor" in i.unique_id}
        if h_rec_ids.isdisjoint(p_rec_ids):
            continue
        h_lig_ids = {i.unique_id for i in h_partners if "ligand" in i.unique_id}
        p_lig_ids = {i.unique_id for i in p_partners if "ligand" in i.unique_id}
        if h_lig_ids.isdisjoint(p_lig_ids):
            continue
        logger.debug(
            f"Removing hydrophobic: {(h.partner_1[0].unique_id, h.partner_2[0].unique_id)}"
        )
        filter_out.append(h)
    selection2 = [contact for contact in selection1 if contact not in filter_out]

    # Filter out contacts that involve an atom involved in another hydrophobic
    filter_out = []
    for h1, h2 in combinations(selection2, 2):
        if h1 == h2:
            raise ValueError(
                "Your list contains exact duplicates, which is not supported "
                "(not needed normally)"
            )
        if h1 in filter_out or h2 in filter_out:
            continue
        h1_partners = h1.partner_1 + h1.partner_2
        h2_partners = h2.partner_1 + h2.partner_2
        h1_lig_ids = {i.unique_id for i in h1_partners if "ligand" in i.unique_id}
        h2_lig_ids = {i.unique_id for i in h2_partners if "ligand" in i.unique_id}
        h1_rec_ids = {i.unique_id for i in h1_partners if "receptor" in i.unique_id}
        h2_rec_ids = {i.unique_id for i in h2_partners if "receptor" in i.unique_id}
        if h1_lig_ids == h2_lig_ids or h1_rec_ids == h2_rec_ids:
            h = h1 if h1.distance > h2.distance else h2
            logger.debug(
                f"Removing hydrophobic: {(h.partner_1[0].unique_id, h.partner_2[0].unique_id)}"
            )
            filter_out.append(h)
    return [contact for contact in selection2 if contact not in filter_out]


def refine_pi_cations(
    pi_cations_all: List[Pi_Cation],
    pi_stackings_all: List[Pi_Stacking],
) -> List[Pi_Cation]:
    """Refines pi-cation interactions to avoid double counts.

    Removes pi-cations that involve atoms that are also involved in a pi-stacking
    interaction (e.g. histidine rings also positively charged).
    """
    filter_out: List[Pi_Cation] = []
    for c, s in product(pi_cations_all, pi_stackings_all):
        c_partners = c.ring + c.cation
        s_partners = s.partner_1 + s.partner_2
        c_rec_ids = {i.unique_id for i in c_partners if "receptor" in i.unique_id}
        s_rec_ids = {i.unique_id for i in s_partners if "receptor" in i.unique_id}
        if c_rec_ids.isdisjoint(s_rec_ids):
            continue
        c_lig_ids = {i.unique_id for i in c_partners if "ligand" in i.unique_id}
        s_lig_ids = {i.unique_id for i in s_partners if "ligand" in i.unique_id}
        if c_lig_ids.isdisjoint(s_lig_ids):
            continue
        logger.debug(
            f"Removing pi-cation: {(c.ring[0].unique_id, c.cation[0].unique_id)}"
        )
        filter_out.append(c)
    return [contact for contact in pi_cations_all if contact not in filter_out]


def refine_pi_hydrophobics(
    pi_hydrohobics_all: List[Pi_Hydrophobic],
) -> List[Pi_Hydrophobic]:
    """Refines pi-hydrophobic interactions to avoid double counts.

    When 2 pi-hydrophobics involve the same group of atoms (ligand or receptor)
    interacting with atoms of the same residue, this function will keep the contact
    with the smallest offset.
    """
    filter_out: List[Pi_Hydrophobic] = []
    for p1, p2 in combinations(pi_hydrohobics_all, 2):
        if p1 == p2:
            raise ValueError(
                "Your list contains exact duplicates, which is not supported "
                "(not needed normally)"
            )
        if p1 in filter_out or p2 in filter_out:
            continue
        p1_partners = p1.ring + p1.hydrophobic
        p2_partners = p2.ring + p2.hydrophobic
        p1_lig_ids = {i.unique_id for i in p1_partners if "ligand" in i.unique_id}
        p2_lig_ids = {i.unique_id for i in p2_partners if "ligand" in i.unique_id}
        p1_rec_ids = {i.unique_id for i in p1_partners if "receptor" in i.unique_id}
        p2_rec_ids = {i.unique_id for i in p2_partners if "receptor" in i.unique_id}
        if p1_lig_ids == p2_lig_ids:
            id_1 = [i.residue_id for i in p1_partners if "receptor" in i.unique_id][0]
            id_2 = [i.residue_id for i in p2_partners if "receptor" in i.unique_id][0]
            if id_1 == id_2:
                p = p1 if p1.offset > p2.offset else p2
            else:
                continue
        elif p1_rec_ids == p2_rec_ids:
            id_1 = [i.residue_id for i in p1_partners if "ligand" in i.unique_id][0]
            id_2 = [i.residue_id for i in p2_partners if "ligand" in i.unique_id][0]
            if id_1 == id_2:
                p = p1 if p1.offset > p2.offset else p2
            else:
                continue
        else:
            continue
        logger.debug(
            f"Removing pi-hydrophobic: {(p.ring[0].unique_id, p.hydrophobic[0].unique_id)}"
        )
        filter_out.append(p)
    return [contact for contact in pi_hydrohobics_all if contact not in filter_out]


def refine_h_bonds(
    h_bonds_all: List[H_Bond],
    salt_bridges_all: List[Salt_Bridge],
    water_bridges_all: List[Water_Bridge],
    metal_complexes: List[Metal_Complex],
) -> List[H_Bond]:
    """Refines H-bonds to avoid double counts.

    If a ligand H-bond donor and acceptor are also involved in a separate
    salt/water bridge, it favors the bridge; if an H-bond involves a water coordinated
    to a metal, it discards the H-bond; if the same h-bond donor is involved in more
    than 1 H-bond, it keeps the one with the smallest D--A distance or the largest
    A-H-D angle if the distances are similar.
    """
    # Remove contact if donor AND acceptor form a salt bridge
    filter_out: List[H_Bond] = []
    for h, b in product(h_bonds_all, salt_bridges_all):
        h_partners = h.acceptor + h.donor
        b_partners = b.partner_1 + b.partner_2
        h_rec_ids = {i.unique_id for i in h_partners if "receptor" in i.unique_id}
        b_rec_ids = {i.unique_id for i in b_partners if "receptor" in i.unique_id}
        if h_rec_ids.isdisjoint(b_rec_ids):
            continue
        h_lig_ids = {i.unique_id for i in h_partners if "ligand" in i.unique_id}
        b_lig_ids = {i.unique_id for i in b_partners if "ligand" in i.unique_id}
        if h_lig_ids.isdisjoint(b_lig_ids):
            continue
        logger.debug(
            "Removing H-bond (duplicate): "
            f"{(h.acceptor[0].unique_id, h.donor[0].unique_id)}"
        )
        filter_out.append(h)
    selection1 = [contact for contact in h_bonds_all if contact not in filter_out]

    # Remove contact if H-bond actually part of a water bridge
    filter_out = []
    for h, w in product(selection1, water_bridges_all):
        h_partners = h.acceptor + h.donor
        w_partners = w.acceptor + w.donor
        h_lig_ids = {i.unique_id for i in h_partners if "ligand" in i.unique_id}
        h_rec_ids = {i.unique_id for i in h_partners if "receptor" in i.unique_id}
        w_lig_ids = {i.unique_id for i in w_partners if "ligand" in i.unique_id}
        w_wat_ids = {i.unique_id for i in w.water}
        if h_rec_ids.isdisjoint(w_wat_ids) or h_lig_ids != w_lig_ids:
            continue
        logger.debug(
            "Removing H-bond (water bridge): "
            f"{(h.acceptor[0].unique_id, h.donor[0].unique_id)}"
        )
        filter_out.append(h)
    selection2 = [h for h in selection1 if h not in filter_out]

    # Remove contact if ligand atom part of a metal complex and H-bond with a binder
    # of that same complex
    filter_out = []
    for h, c in product(selection2, metal_complexes):
        h_partners = h.acceptor + h.donor
        h_rec_ids = {i.unique_id for i in h_partners if "receptor" in i.unique_id}
        c_rec_ids = {i.unique_id for i in c.receptor}
        if h_rec_ids.isdisjoint(c_rec_ids):
            continue
        logger.debug(
            "Removing H-bond (metal complex): "
            f"{(h.acceptor[0].unique_id, h.donor[0].unique_id)}"
        )
        filter_out.append(h)
    selection3 = [h for h in selection2 if h not in filter_out]

    # Allow only one H-bond per donor (compare strong with strong and weak with weak)
    filter_out = []
    for h1, h2 in combinations(selection3, 2):
        if h1 == h2:
            raise ValueError(
                "Your list contains exact duplicates, which is not supported "
                "(not needed normally)"
            )
        if h1 in filter_out or h2 in filter_out:
            continue
        if h1.type != h2.type:
            continue
        h1_partners = h1.acceptor + h1.donor
        h2_partners = h2.acceptor + h2.donor
        h1_lig_ids = {i.unique_id for i in h1_partners if "ligand" in i.unique_id}
        h2_lig_ids = {i.unique_id for i in h2_partners if "ligand" in i.unique_id}
        if len(h1_lig_ids) == 2 and h1_lig_ids != h2_lig_ids:
            continue
        h1_rec_ids = {i.unique_id for i in h1_partners if "receptor" in i.unique_id}
        h2_rec_ids = {i.unique_id for i in h2_partners if "receptor" in i.unique_id}
        if len(h2_rec_ids) == 2 and h1_rec_ids != h2_rec_ids:
            continue
        # Check the distance first, if one is much shorter than the other, keep this one,
        # otherwise, keep the one with the largest donor angle
        if abs(h1.distance_ad - h2.distance_ad) > 0.1:
            h = h1 if h1.distance_ad > h2.distance_ad else h2
        else:
            h = h1 if h1.angle_dha < h2.angle_dha else h2
        logger.debug(
            "Removing H-bond (donor with multiple H-bonds): "
            f"{(h.acceptor[0].unique_id, h.donor[0].unique_id)}"
        )
        filter_out.append(h)
    return [contact for contact in selection3 if contact not in filter_out]


def refine_water_bridges(
    water_bridges_all: List[Water_Bridge], metal_complexes: List[Metal_Complex]
) -> List[Water_Bridge]:
    """Refines water bridges to avoid double counts.

    If a water bridge involves a water also involved in a metal complex, this function
    will drop it (this is because water molecules coordinated to metals tend to form
    lots of water bridges and H-bonds just because they are close, but angles can be unrealistic).
    """
    filter_out = []
    for w, c in product(water_bridges_all, metal_complexes):
        if w in filter_out:
            continue
        w_wat_id = w.water[0].unique_id
        c_lig_ids = {i.unique_id for i in c.ligand}
        c_rec_ids = {i.unique_id for i in c.receptor}
        if w_wat_id in c_lig_ids or w_wat_id in c_rec_ids:
            logger.debug(
                f"Removing water bridge: "
                f"{(w.acceptor[0].unique_id, w.donor[0].unique_id)}"
            )
            filter_out.append(w)
    selection = [w for w in water_bridges_all if w not in filter_out]
    return selection


def drop_duplicated_h_bonds(
    h_bonds_all: List[H_Bond],
) -> List[H_Bond]:
    """Drops duplicated H-bonds (same donor and same acceptor but different Hs).

    This is needed when we allow Hs to rotate for groups like R-NH2 or R-NH3+.
    """
    filter_out = []
    for h1, h2 in combinations(h_bonds_all, 2):
        if h1 == h2:
            raise ValueError(
                "Your list contains exact duplicates, which is not supported "
                "(not needed normally)"
            )
        if h1 in filter_out or h2 in filter_out:
            continue
        h1_partners = h1.acceptor + h1.donor[:1]
        h2_partners = h2.acceptor + h2.donor[:1]
        h1_lig_ids = {i.unique_id for i in h1_partners if "ligand" in i.unique_id}
        h2_lig_ids = {i.unique_id for i in h2_partners if "ligand" in i.unique_id}
        if h1_lig_ids != h2_lig_ids:
            continue
        h1_rec_ids = {i.unique_id for i in h1_partners if "receptor" in i.unique_id}
        h2_rec_ids = {i.unique_id for i in h2_partners if "receptor" in i.unique_id}
        if h1_rec_ids != h2_rec_ids:
            continue
        # Same donor AND same acceptor (but a different H) -> keep the H-bond
        # with the shortest A---H distance or largest D-H---A angle if distances are similar
        if abs(h1.distance_ah - h2.distance_ah) > 0.1:
            h = h1 if h1.distance_ah > h2.distance_ad else h2
        else:
            h = h1 if h1.angle_dha < h2.angle_dha else h2
        logger.debug(
            "Removing H-bond (duplicate): "
            f"{(h.acceptor[0].unique_id, h.donor[0].unique_id)}"
        )
        filter_out.append(h)
    return [contact for contact in h_bonds_all if contact not in filter_out]


def drop_duplicated_water_bridges(
    water_bridges_all: List[Water_Bridge],
) -> List[Water_Bridge]:
    """Drops duplicated water bridges (same donor / acceptor / water but different Hs).

    This is needed when we allow Hs to rotate for groups like R-NH2 or R-NH3+.
    """
    filter_out = []
    for w1, w2 in combinations(water_bridges_all, 2):
        if w1 == w2:
            raise ValueError(
                "Your list contains exact duplicates, which is not supported "
                "(not needed normally)"
            )
        if w1 in filter_out or w2 in filter_out:
            continue
        w1_partners = w1.acceptor + w1.donor[:1] + w1.water
        w2_partners = w2.acceptor + w2.donor[:1] + w2.water
        w1_lig_ids = {i.unique_id for i in w1_partners if "ligand" in i.unique_id}
        w2_lig_ids = {i.unique_id for i in w2_partners if "ligand" in i.unique_id}
        if w1_lig_ids != w2_lig_ids:
            continue
        w1_rec_ids = {i.unique_id for i in w1_partners if "receptor" in i.unique_id}
        w2_rec_ids = {i.unique_id for i in w2_partners if "receptor" in i.unique_id}
        if w1_rec_ids != w2_rec_ids:
            continue
        # Same donor AND same acceptor (but a different H) -> keep the water bridge
        # with the largest D-H---W angle
        w = w1 if w1.angle_dhw < w2.angle_dhw else w2
        logger.debug(
            "Removing water bridge (duplicate): "
            f"{(w.acceptor[0].unique_id, w.donor[0].unique_id)}"
        )
        filter_out.append(w)
    return [contact for contact in water_bridges_all if contact not in filter_out]
