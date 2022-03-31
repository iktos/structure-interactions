from __future__ import absolute_import

from itertools import combinations, product
from logging import getLogger

logger = getLogger(__name__)


def filter_contact_set(set1, set2):
    """
    Filter out contacts that are in set1 (probe) and set2 (ref)
    e.g. L and P atoms involved in hydrophobic + pi-stacking
    Loop over contact pairs in product(set1, set2) and compare
    both receptor and ligand ids (note: this function does not work
    with water and metal bridges)
    Return: list of contacts to filter out from set1
    """
    filter_out = []
    for pair1, pair2 in product(set1, set2):
        if pair1 in filter_out:
            continue
        pair1_rec_ids = [i.unique_id for i in pair1.receptor]
        pair2_rec_ids = [i.unique_id for i in pair2.receptor]
        pair1_lig_ids = [i.unique_id for i in pair1.ligand]
        pair2_lig_ids = [i.unique_id for i in pair2.ligand]
        if not set(pair1_lig_ids).isdisjoint(pair2_lig_ids) and not set(
            pair1_rec_ids
        ).isdisjoint(pair2_rec_ids):
            filter_out.append(pair1)
    logger.debug(
        f'Removing the following contacts: '
        f'{[(p.receptor[0].unique_id, p.ligand[0].unique_id) for p in filter_out]}'
    )
    return filter_out


def select_best_contact(set1, mode='hydrophobic'):
    """
    Filter out contacts that involve an atom already involved
    in another contact of the same type -> keep the 'best' contact
    e.g. atom involved in 2 different hydrophobics
    If mode='hydrophobic' -> keep shortest bond
    If mode='pi_hydrophobic' -> keep smallest offset
    If mode='h_bond' -> keep largest donor angle
    For H-bonds and water bridges, focus on donor pairs
    Return: list of contacts to remove from set1
    """
    filter_out = []
    for pair1, pair2 in combinations(set1, 2):
        if pair1 in filter_out:
            continue
        pair1_lig_ids = [i.unique_id for i in pair1.ligand]
        pair2_lig_ids = [i.unique_id for i in pair2.ligand]
        pair1_rec_ids = [i.unique_id for i in pair1.receptor]
        pair2_rec_ids = [i.unique_id for i in pair2.receptor]
        if not set(pair1_lig_ids).isdisjoint(pair2_lig_ids) or not set(
            pair1_rec_ids
        ).isdisjoint(pair2_rec_ids):
            if mode == 'hydrophobic':
                filter_out += [pair1] if pair1.distance > pair2.distance else [pair2]
            elif mode == 'pi_hydrophobic':
                filter_out += [pair1] if pair1.offset > pair2.offset else [pair2]
            elif mode == 'h_bond' and len(pair1_lig_ids) == len(pair2_lig_ids) == 2:
                filter_out += [pair1] if pair1.angle_dha > pair2.angle_dha else [pair2]
    logger.debug(
        f'Removing the following contacts: '
        f'{[(p.receptor[0].unique_id, p.ligand[0].unique_id) for p in filter_out]}'
    )
    return filter_out


def refine_hydrophobics(hydrophobics_all, pi_stackings_all, pi_hydrophobics_all):
    """
    Refine hydrophobic interactions to avoid double counts
    """
    # Remove hydrophobic contacts for atoms that
    # also interact via stacking or pi-hydrophobic
    filter_out = filter_contact_set(
        hydrophobics_all, pi_stackings_all + pi_hydrophobics_all
    )
    selection1 = [h for h in hydrophobics_all if h not in filter_out]

    # Filter out contacts that involve an atom involved in another hydrophobic
    filter_out = select_best_contact(selection1, mode='hydrophobic')
    selection2 = [h for h in selection1 if h not in filter_out]
    return selection2


def refine_pi_cations(pi_cations_all, pi_stackings_all):
    """
    Refine pi-cation interactions to avoid double counts
    If pi-cation and pi-stacking reported for the same L-P residues,
    keep the pi-stacking (e.g. histidine rings also positively charged)
    """
    filter_out = filter_contact_set(pi_cations_all, pi_stackings_all)
    selection = [c for c in pi_cations_all if c not in filter_out]
    return selection


def refine_pi_hydrophobics(pi_hydrohobics_all):
    """
    Refine pi-hydrophobic interactions to avoid double counts
    """
    filter_out = select_best_contact(pi_hydrohobics_all, mode='pi_hydrophobic')
    selection = [h for h in pi_hydrohobics_all if h not in filter_out]
    return selection


def refine_h_bonds(h_bonds_all, salt_bridges_all, water_bridges_all, metal_complexes):
    """
    Refine H-bonds to avoid double counts
    If ligand H-bond donor or acceptor also involved in a separate
    salt/water bridge, favor the bridge
    If H-bond with a water coordinated to a metal, discard H-bond
    """
    # Remove contact if donor AND acceptor form a salt bridge
    filter_out = filter_contact_set(h_bonds_all, salt_bridges_all)

    # Remove contact if H-bond actually part of a water bridge
    for h, w in product(h_bonds_all, water_bridges_all):
        if h in filter_out:
            continue
        h_lig_ids = [i.unique_id for i in h.ligand]
        h_rec_ids = [i.unique_id for i in h.receptor]
        w_lig_ids = [i.unique_id for i in w.ligand]
        w_wat_ids = [i.unique_id for i in w.water]
        if h_lig_ids == w_lig_ids and not set(h_rec_ids).isdisjoint(w_wat_ids):
            logger.debug(
                f'Removing H-bond: '
                f'{(h.receptor[0].unique_id, h.ligand[0].unique_id)}'
            )
            filter_out.append(h)
    selection = [h for h in h_bonds_all if h not in filter_out]

    # Remove contact if ligand atom part of a metal complex
    # and H-bond with a binder of that same complex
    for h, c in product(h_bonds_all, metal_complexes):
        if h in filter_out:
            continue
        h_rec_ids = [i.unique_id for i in h.receptor]
        c_rec_ids = [i.unique_id for i in c.receptor]
        if not set(h_rec_ids).isdisjoint(c_rec_ids):
            logger.debug(
                f'Removing H-bond: '
                f'{(h.receptor[0].unique_id, h.ligand[0].unique_id)}'
            )
            filter_out.append(h)
    selection = [h for h in h_bonds_all if h not in filter_out]
    h_bonds = selection

    # Allow only one H-bond per donor
    filter_out = select_best_contact(h_bonds, mode='h_bond')
    selection = [h for h in h_bonds if h not in filter_out]
    return selection


def refine_water_bridges(water_bridges_all, metal_complexes):
    """
    Refine water bridges to avoid double counts
    If water bridge with a water involved in a metal complex, keep the metal complex
    (water molecules coordinated to metals tend to form lots of water bridges
    and H-bonds just because they are close, but angles can be unrealistic)
    """
    filter_out = []
    for w, c in product(water_bridges_all, metal_complexes):
        if w in filter_out:
            continue
        w_wat_id = w.water[0].unique_id
        c_lig_ids = [i.unique_id for i in c.ligand]
        c_rec_ids = [i.unique_id for i in c.receptor]
        if w_wat_id in c_lig_ids + c_rec_ids:
            logger.debug(
                f'Removing water bridge: '
                f'{(w.receptor[0].unique_id, w.ligand[0].unique_id)}'
            )
            filter_out.append(w)
    selection = [w for w in water_bridges_all if w not in filter_out]
    return selection
