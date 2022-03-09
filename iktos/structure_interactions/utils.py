from __future__ import absolute_import

from copy import deepcopy
from itertools import product
from typing import List, Optional, Tuple

import numpy as np
from iktos.logger import getLogger
from iktos.structure_utils.toolkits.obabel import read_mol as read_obmol
from iktos.structure_utils.toolkits.obabel import rename_residue

try:
    from openbabel.openbabel import (  # openbabel 3
        OBAtom,
        OBAtomAtomIter,
        OBMol,
        OBMolAtomIter,
        obErrorLog,
    )
except ModuleNotFoundError:
    from openbabel import (
        OBAtom,
        obErrorLog,
        OBMol,
        OBMolAtomIter,
        OBAtomAtomIter,
    )  # openbabel 2 (warning in structure-utils)

from . import constants
from .math_utils import get_vector, get_vector_angle

logger = getLogger(__name__)


def get_coords(obatom):
    return np.array([obatom.GetX(), obatom.GetY(), obatom.GetZ()])


def map_atom_ids(obmol, mapping):
    """Maps atom Id back to their original value.

    Returns:
        updated OBMol object

    Note: this function changes the internal atom Id, not the id
    that would be printed out in the corresponding pdb string/file
    mapping: dict[internal_atom_id] = orignal_atom_id
    """
    mapping_copy = deepcopy(mapping)
    for obatom in OBMolAtomIter(obmol):
        internal_atom_id = obatom.GetIdx()
        if mapping is not None:
            original_atom_id = mapping_copy[internal_atom_id]
            mapping_copy.pop(internal_atom_id)
        else:
            original_atom_id = internal_atom_id
        obatom.SetId(original_atom_id)
    return obmol


def ring_is_planar(ring, r_atoms):
    """Given a set of ring atoms, checks whether the ring
    is sufficiently planar to be considered aromatic."""
    normals = []
    for a in r_atoms:
        obneighs = OBAtomAtomIter(a)
        n_coords = [
            get_coords(obneigh) for obneigh in obneighs if ring.IsMember(obneigh)
        ]
        vec1, vec2 = (
            get_vector(get_coords(a), n_coords[0]),
            get_vector(get_coords(a), n_coords[1]),
        )
        normals.append(np.cross(vec1, vec2))
    # Given all normals of ring atoms and their neighbours,
    # the angle between any has to be 5.0 deg or less
    for n1, n2 in product(normals, repeat=2):
        arom_angle = get_vector_angle(n1, n2)
        if all(
            [
                arom_angle > constants.AROMATIC_PLANARITY,
                arom_angle < 180.0 - constants.AROMATIC_PLANARITY,
            ]
        ):
            return False
    return True


def ring_is_aromatic(r_atoms: List[OBAtom]) -> bool:
    """Given a set of ring atoms, checks whether the ring
    contains only SP2 atoms, in which case it can be considered as aromatic.
    """
    for a in r_atoms:
        if a.GetHyb() != 2:
            return False
    return True


def atom_has_lone_pair(obatom: OBAtom) -> bool:
    """Checks if input atom has a free lone pair
    (-> potential metal binder and halogen-bond acceptor)
    Considers N, O, S, P as binders if not positively charged
    except for N (and P) SP2 with 3 neighbours
    """
    if obatom.GetFormalCharge() < 0:
        return True
    if obatom.GetFormalCharge() > 0:
        return False
    if obatom.GetAtomicNum() in [8, 16]:
        obneighs = list(OBAtomAtomIter(obatom))
        if len(obneighs) <= 2:
            return True
    if obatom.GetAtomicNum() in [7, 15]:
        obneighs = list(OBAtomAtomIter(obatom))
        if len(obneighs) > 3:
            return False
        if obatom.GetHyb() == 2 and len(obneighs) == 3:
            return False
        return True
    return False


def identify_charged_group(obatom: OBAtom) -> Tuple[str, str, List[OBAtom]]:
    """Identifies functional group a charged atom belongs to.

    Returns:
        fonctional group name, charge and list of obatoms involved
    """
    # Check the neighbours, if some of them are charged -> log and ignore
    obneighs = list(OBAtomAtomIter(obatom))
    if not all(a.GetFormalCharge() == 0 for a in obneighs):
        logger.debug('Found a group with multiple charges (e.g. NO2, NO) -> ignoring')
        return 'unknown', 'unknown', [obatom]
    if obatom.GetAtomicNum() == 8:
        return _identify_group_with_oxygen(obatom)
    if obatom.GetAtomicNum() == 7:
        return _identify_group_with_nitrogen(obatom)
    if obatom.GetAtomicNum() == 6:
        if obatom.GetFormalCharge() < 0:
            logger.debug('Found a carbanion; please check your input')
            return 'carbanion', 'negative', [obatom]
        logger.debug('Found a carbocation; please check your input')
        return 'carbocation', 'positive', [obatom]
    if obatom.GetAtomicNum() == 15:  # phosphorus
        if obatom.GetFormalCharge() > 0:
            return 'phopsphonium', 'positive', [obatom]
        return 'phospho-anion', 'negative', [obatom]
    if obatom.GetAtomicNum() == 16:  # sulfur
        if obatom.GetFormalCharge() > 0:
            return 'sulfonium', 'positive', [obatom]
        return 'sulfo-anion', 'negative', [obatom]
    logger.debug('Failed to identify charged group')
    return 'unknown', 'unknown', [obatom]


def _identify_group_with_oxygen(obatom: OBAtom) -> Tuple[str, str, List[OBAtom]]:
    """Identifies functional group a charged oxygen belongs to.

    Returns:
        fonctional group name, charge and list of obatoms involved
    """
    if obatom.IsCarboxylOxygen():
        obcarbon = [obneigh for obneigh in OBAtomAtomIter(obatom)][0]
        obneighs = get_o_neighbours(obcarbon)
        return 'carboxylate', 'negative', [obcarbon] + obneighs
    if obatom.IsPhosphateOxygen():  # R-PO3
        obphosphorus = [obneigh for obneigh in OBAtomAtomIter(obatom)][0]
        obneighs = get_o_neighbours(obphosphorus)
        return 'phosphate', 'negative', [obphosphorus] + obneighs
    if obatom.IsSulfateOxygen():  # R-SO3
        obsulfur = list(OBAtomAtomIter(obatom))[0]
        obneighs = get_o_neighbours(obsulfur)
        return 'sulfate', 'negative', [obsulfur] + obneighs
    if obatom.MatchesSMARTS('[$([#8-])]'):
        return 'alkoxide', 'negative', [obatom]
    return 'oxycation', 'positive', [obatom]


def _identify_group_with_nitrogen(obatom: OBAtom) -> Tuple[str, str, List[OBAtom]]:
    """Identifies functional group a charged nitrogen belongs to.

    Returns:
        fonctional group name, charge and list of obatoms involved
    """
    if obatom.MatchesSMARTS('[$([#7+X4])]'):
        return 'ammonium', 'positive', [obatom]
    if obatom.MatchesSMARTS('[$([#7+X3]=[#6X3])]'):
        obneighs = list(OBAtomAtomIter(obatom))
        atom_list = [obatom]
        for obneigh in obneighs:
            obnitrogens = get_n_neighbours(obneigh)
            # Remove current atom from the list of its neighbor's neighbors
            # Remove nitrogens that are not SP2 or are positively charged
            obnitrogens = [
                a
                for a in obnitrogens
                if a != obatom and a.GetHyb() == 2 and a.GetFormalCharge() <= 0
            ]
            if obnitrogens:
                atom_list.append(obneigh)
                atom_list += obnitrogens
                logger.debug('Found a delocalised iminium')
        return 'iminium', 'positive', atom_list
    if obatom.MatchesSMARTS('[$([#7+])]'):
        return 'azacation', 'positive', [obatom]
    return 'azanide', 'negative', [obatom]


def get_o_neighbours(obatom: OBAtom) -> List[OBAtom]:
    return [a for a in OBAtomAtomIter(obatom) if a.GetAtomicNum() == 8]


def get_n_neighbours(obatom: OBAtom) -> List[OBAtom]:
    return [a for a in OBAtomAtomIter(obatom) if a.GetAtomicNum() == 7]


def get_h_neighbours(obatom: OBAtom) -> List[OBAtom]:
    return [a for a in OBAtomAtomIter(obatom) if a.GetAtomicNum() == 1]


def read_mol(
    coords: str, fmt: str, as_string: bool, title: str = ''
) -> Optional[OBMol]:
    """Reads a coords file/string (PDB, MOL2, SDF) and returns an OBMol object."""
    obErrorLog.StopLogging()  # to avoid huge logs
    if as_string:
        obmol = read_obmol(coords, fmt=fmt)
    else:
        obmol = read_obmol(coords, fmt=fmt, from_file=True)
    if obmol is None:
        return None
    if obmol.NumResidues() == 0:
        # logger.warning('No residue perceived for input mol -> creating a `UNL` residue')
        obmol = rename_residue(obmol, 'UNL')
    obmol.PerceiveBondOrders()  # assign multiple bonds
    correct_implicit_count(obmol)
    correct_hybridisation(obmol)
    obmol.SetTitle(title)
    return obmol


def correct_implicit_count(obmol) -> None:
    """Sets to 0 the number of implicit Hs."""
    obatoms = list(OBMolAtomIter(obmol))
    for obatom in obatoms:
        try:
            obatom.SetImplicitHCount(0)
        except AttributeError:
            obatom.SetImplicitValence(0)


def correct_hybridisation(obmol) -> None:
    """Corrects the perception of hybridisation for neutral C and N atoms:
        - first, loops over carbons
          C with 2 neighbours -> sp
          C with 3 neighbours -> sp2
          C with 4 neighbours -> sp3
        - then loops over nitrogens
          N with 1 neighbour -> sp
          N with 2 neighbours -> sp2
          N with 3 neighbours, including 1 Csp2 -> sp2 (delocalisation)
          N with 3 neighbours, including 0 Csp2 -> sp3

    Note: this is a quick fix to correct the original perception
    by OpenBabel. Not all atoms and charge states need to be considered,
    only the ones for which the original perception leads to an incorrect
    property assignment by PLIP (e.g. for aromaticity or lone pairs).
    It assumes that the input mol contains 0 implicit Hs.
    """
    obatoms = [
        obatom
        for obatom in OBMolAtomIter(obmol)
        if obatom.GetFormalCharge() == 0 and obatom.GetAtomicNum() == 6
    ]
    for obatom in obatoms:
        hyb = obatom.GetHyb()
        idx = obatom.GetIdx()
        obneighs = list(OBAtomAtomIter(obatom))
        if len(obneighs) == 2 and hyb != 1:
            logger.debug(f'Changing hybridisation perceived for C{idx} to SP')
            obatom.SetHyb(1)
        if len(obneighs) == 3 and hyb != 2:
            logger.debug(f'Changing hybridisation perceived for C{idx} to SP2')
            obatom.SetHyb(2)
        if len(obneighs) == 4 and hyb != 3:
            logger.debug(f'Changing hybridisation perceived for C{idx} to SP3')
            obatom.SetHyb(3)
    obatoms = [
        obatom
        for obatom in OBMolAtomIter(obmol)
        if obatom.GetFormalCharge() == 0 and obatom.GetAtomicNum() == 7
    ]
    for obatom in obatoms:
        hyb = obatom.GetHyb()
        idx = obatom.GetIdx()
        obneighs = list(OBAtomAtomIter(obatom))
        if len(obneighs) == 1 and hyb != 1:
            logger.debug(f'Changing hybridisation perceived for N{idx} to SP')
            obatom.SetHyb(1)
        if len(obneighs) == 2 and hyb != 2:
            logger.debug(f'Changing hybridisation perceived for N{idx} to SP2')
            obatom.SetHyb(2)
        if len(obneighs) == 3:
            obneighs_csp2 = [
                obneigh
                for obneigh in obneighs
                if obneigh.GetAtomicNum() == 6 and obneigh.GetHyb() == 2
            ]
            if len(obneighs_csp2) == 0 and hyb != 3:
                logger.debug(f'Changing hybridisation perceived for N{idx} to SP3')
                obatom.SetHyb(3)
            if len(obneighs_csp2) != 0 and hyb != 2:
                logger.debug(f'Changing hybridisation perceived for N{idx} to SP2')
                obatom.SetHyb(2)
    # Set some tags on mol to force a new perception
    # of properties related to hydridisation
    if obmol.HasAtomTypesPerceived():
        try:
            obmol.SetAtomTypesPerceived(False)
        except TypeError:
            logger.warning(
                'Failed to set `HasAtomTypesPerceived` to False, '
                'be careful with hybridisation perceived by OpenBabel'
            )


def parse_pdb(pdb, as_string=True):
    """Parses a PDB file/string and extract additional information.

    Note: OpenBabel starts numbering at 1 and counts TER lines,
          so we need to keep track of original ids
    """
    if as_string:
        pdb_lines = pdb.splitlines(True)
    else:
        with open(pdb) as f:
            pdb_lines = f.readlines()
    i = 1
    mapping = {}
    pdb_block = ''
    for line in pdb_lines:
        if line.startswith(('ATOM', 'HETATM')):
            atom_idx = int(line[6:11])
            location = line[16]
            if location != ' ':
                raise ValueError('Found unsupported alternate location value')
            mapping[i] = atom_idx
            pdb_block += line
            i += 1
    return pdb_block, mapping
