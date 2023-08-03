import math
from builtins import filter
from typing import List, NamedTuple

import numpy as np
import numpy.typing as npt

try:
    from iktos.logger import getLogger
except ImportError:
    from logging import getLogger

try:
    from openbabel.openbabel import (  # openbabel 3
        OBAtom,
        OBMol,
        OBAtomAtomIter,
        OBMolAtomIter,
        OBResidueAtomIter,
    )
except ModuleNotFoundError:
    from openbabel import (
        OBAtom,
        OBAtomAtomIter,
        OBMol,
        OBMolAtomIter,
        OBResidueAtomIter,
    )  # openbabel 2 (warning in structure-utils)

from . import constants
from .Atom import Atom
from .math_utils import get_centroid, get_vector, normalize_vector, rotate
from .mol_utils import (
    atom_has_lone_pair,
    get_atom_coordinates,
    get_h_neighbours,
    get_n_neighbours,
    get_o_neighbours,
    identify_charged_group,
    ring_is_aromatic,
    ring_is_planar,
)

logger = getLogger(__name__)


class Ring(NamedTuple):
    """Class to store information about aromatic rings.

    Attributes:
        atom_list: list of atoms in the ring.
        normal: vector defining the normal to the ring.
        center: coordinates of the center of the ring.
    """

    atom_list: List[Atom]
    normal: npt.NDArray
    center: npt.NDArray


def find_rings(obmol: OBMol) -> List[Ring]:
    """Finds all aromatic rings.

    Aromatic rings include 5- and 6-membered rings detected by OpenBabel as aromatic,
    or sufficiently planar.
    """
    aromatic_amino = ['TYR', 'TRP', 'HIS', 'HID', 'HIE', 'HIP', 'HSD', 'HSE', 'PHE']
    all_rings = obmol.GetSSSR()
    selection = []
    obmol_atoms = list(OBMolAtomIter(obmol))
    for obring in all_rings:
        ring_atoms = [a for a in obmol_atoms if obring.IsMember(a)]

        # Residue name is only relevant for the proteins
        residue_name = ''
        if ring_atoms[0].GetResidue() is not None:
            residue_name = ring_atoms[0].GetResidue().GetName()
        if len(ring_atoms) not in [5, 6]:
            continue  # ignore rings of unusual size
        if (
            residue_name in aromatic_amino
            or obring.IsAromatic()
            or ring_is_aromatic(ring_atoms)
            or ring_is_planar(obring, ring_atoms)
        ):
            # Tag atoms in ring as aromatic (for later)
            for a in ring_atoms:
                a.SetAromatic()
            atoms = [Atom(a) for a in ring_atoms]
            coords = [get_atom_coordinates(ring_atoms[i]) for i in [0, 2, 4]]
            ringv1 = get_vector(coords[0], coords[1])
            ringv2 = get_vector(coords[2], coords[0])
            normal = normalize_vector(np.cross(ringv1, ringv2))
            center = get_centroid([get_atom_coordinates(ra) for ra in ring_atoms])
            selection.append(Ring(atom_list=atoms, normal=normal, center=center))
    logger.debug(f'Found {len(selection)} aromatic ring(s)')
    return selection


class HydrophobicAtom(NamedTuple):
    """Class to store information about hydrophobic atoms.

    Attributes:
        atom_list: hydrophobic atom (in a list to be consistent with other classes).
        neighbours_radius_2: neighbours of hydrophobic atom up to radius 2 - this is used
            to discard interactios between bound atoms (relevant for intra).
    """

    atom_list: List[Atom]
    neighbours_radius_2: List[Atom]


def find_hydrophobics(obatoms: List[OBAtom]) -> List[HydrophobicAtom]:
    """Finds all hydrophobic atoms.

    Hydrophobic atoms include carbon and sulfur atoms which have only carbons/sulfurs/hydrogens
    as direct neighbors + Cl and F. Sulfur atoms that are not SP3 are excluded.
    """
    selection = []
    obatoms_filtered = filter(
        lambda obatom: (
            obatom.GetAtomicNum() in [6, 16, 9, 17]
            and set(
                [obneigh.GetAtomicNum() for obneigh in OBAtomAtomIter(obatom)]
            ).issubset({1, 6, 16})
        ),
        obatoms,
    )
    for obatom in obatoms_filtered:
        if obatom.GetAtomicNum() == 16 and obatom.GetHyb() != 3:
            logger.debug(
                f'Excluding S{obatom.GetId()} from list of hydrophobics (not SP3)'
            )
            continue
        # List neighbours radius 1
        obneighs_1 = [
            obneigh for obneigh in OBAtomAtomIter(obatom) if obneigh.GetAtomicNum() != 1
        ]
        obneighs_all = [a for a in obneighs_1]
        # Extend with neighbours radius 2
        for obneigh_1 in obneighs_1:
            obneighs_2 = [
                obneigh
                for obneigh in OBAtomAtomIter(obneigh_1)
                if obneigh.GetAtomicNum() != 1
                and obneigh != obatom
                and obneigh not in obneighs_all
            ]
            obneighs_all.extend(obneighs_2)
        selection.append(
            HydrophobicAtom(
                atom_list=[Atom(obatom)],
                neighbours_radius_2=[Atom(obneigh) for obneigh in obneighs_all],
            )
        )
    logger.debug(f'Found {len(selection)} hydrophobic atoms')
    return selection


class HBondAcceptor(NamedTuple):
    """Class to store information about H-bond acceptors.

    Attributes:
        atom_list: H-bond acceptor atom (in a list to be consistent with other classes).
        neighbours: list of neighbouring atoms (needed to check the angles around A
            and ensure that the H is on the lone pair side).
    """

    atom_list: List[Atom]
    neighbours: List[Atom]


def find_h_bond_acceptors(obatoms: List[OBAtom]) -> List[HBondAcceptor]:
    """Finds all H-bond acceptors.

    Halogens and sulfur are NOT considered as H-bond acceptors.

    Neighbours of the acceptor are listed as well and are used in the code
    to make shure H-bonds are detected on the lone pair side.
    """
    selection = []
    obatoms_filtered = filter(
        lambda obatom: (obatom.IsHbondAcceptor() or atom_has_lone_pair(obatom))
        and obatom.GetAtomicNum() not in [9, 17, 35, 53, 16],
        obatoms,
    )
    for obatom in obatoms_filtered:
        obneighs = list(OBAtomAtomIter(obatom))
        selection.append(
            HBondAcceptor(
                atom_list=[Atom(obatom)],
                neighbours=[Atom(a) for a in obneighs],
            )
        )
    logger.debug(f'Found {len(selection)} H-bond acceptors')
    return selection


class HBondDonor(NamedTuple):
    """Class to store information about H-bond donors.

    Attributes:
        atom_list: H-bond donor pair (D, H).
        type: weak or strong.
        alt_locations: alternative locations of H (only for e.g. R-OH, R-NH2, R-NH3+).
    """

    atom_list: List[Atom]
    type: str
    alt_locations: List[npt.NDArray]


def find_h_bond_donors(obatoms: List[OBAtom]) -> List[HBondDonor]:
    """Finds all strong and weak H-bond donors (i.e. with polarised C-H bonds).

    For strong rotatable H-bond donors (R-OH, R-NH2), alternative positions
    of H are added to the object `HBondDonor`.
    """
    selection = []

    # Loop over atoms flagged by OpenBabel as H-bond donor
    obatoms_filtered = filter(lambda obatom: obatom.IsHbondDonor(), obatoms)
    for obatom in obatoms_filtered:
        # Loop over Hs bonded to this H-bond donor
        obneighs_h = [
            obneigh for obneigh in OBAtomAtomIter(obatom) if obneigh.GetAtomicNum() == 1
        ]
        for i, obneigh in enumerate(obneighs_h):
            alt_locations = []

            # If the H-bond donor is e.g. OH or NH2, consider alternative locations
            if obatom.GetHyb() == 3 and obatom.GetHvyDegree() == 1:
                if obatom.GetResidue() is not None:
                    logger.debug(
                        'Identified a rotatable H-bond donor in residue '
                        f'{obatom.GetResidue().GetName()}|{obatom.GetResidue().GetNum()}'
                    )
                coords_h_initial = get_atom_coordinates(obneigh)
                obneighs_not_h = [
                    obneigh
                    for obneigh in OBAtomAtomIter(obatom)
                    if obneigh.GetAtomicNum() != 1
                ]
                coords_neigh = get_atom_coordinates(obneighs_not_h[0])
                coords_donor = get_atom_coordinates(obatom)
                for angle in np.arange(30, 360, 30):
                    coords_h_final = rotate(
                        coords_h_initial,
                        coords_donor,
                        coords_neigh,
                        math.radians(angle),
                    )
                    alt_locations.append(coords_h_final)
            selection.append(
                HBondDonor(
                    atom_list=[Atom(obatom), Atom(obneigh)],
                    type='strong',
                    alt_locations=alt_locations,
                )
            )
    logger.debug(f'Found {len(selection)} strong H-bond donor pairs')

    # Find polarised C-H (near O, or near N, or SP2 or SP) and S-H
    for obatom in obatoms:
        if obatom.GetAtomicNum() == 16:
            obhs = get_h_neighbours(obatom)
            for obh in obhs:
                selection.append(
                    HBondDonor(
                        atom_list=[Atom(obatom), Atom(obh)],
                        type='weak',
                        alt_locations=[],
                    )
                )
        elif obatom.GetAtomicNum() == 6:
            obhs = get_h_neighbours(obatom)
            if not obhs:
                continue
            obneighs = get_n_neighbours(obatom) + get_o_neighbours(obatom)
            obneighs_filtered = filter(lambda a: a.GetFormalCharge() <= 0, obneighs)
            if obatom.GetHyb() in [1, 2] or len(list(obneighs_filtered)) != 0:
                for obh in obhs:
                    selection.append(
                        HBondDonor(
                            atom_list=[Atom(obatom), Atom(obh)],
                            type='weak',
                            alt_locations=[],
                        )
                    )
    logger.debug(f'Found {len(selection)} H-bond donor pairs (incl. weak)')
    return selection


class XBondAcceptor(NamedTuple):
    """Class to store information about X-bond acceptors.

    Attributes:
        atom_list: X-bond acceptor atom (in a list to be consistent with other classes).
        neighbours: list of neighbouring atoms (needed to check the angles around A).
    """

    atom_list: List[Atom]
    neighbours: List[Atom]


def find_x_bond_acceptors(obatoms: List[OBAtom]) -> List[XBondAcceptor]:
    """Finds all X-bond acceptors (O, N, S with lone pair)."""
    selection = []
    obatoms_filtered = filter(
        lambda obatom: obatom.GetAtomicNum() in [7, 8, 16], obatoms
    )
    for obatom in obatoms_filtered:
        if atom_has_lone_pair(obatom):
            obneighs = [a for a in OBAtomAtomIter(obatom)]
            selection.append(
                XBondAcceptor(
                    atom_list=[Atom(obatom)], neighbours=[Atom(a) for a in obneighs]
                )
            )
    logger.debug(f'Found {len(selection)} halogen bond acceptor(s)')
    return selection


class XBondDonor(NamedTuple):
    """Class to store information about X-bond donors.

    Attributes:
        atom_list: X-bond donor pair (C, X).
        type: weak or strong.
    """

    atom_list: List[Atom]


def find_halogens(obatoms: List[OBAtom]) -> List[XBondDonor]:
    """Finds all X-bond donors.

    X-bond donors include C-X, with X = F, Cl, Br, I, and are used for the detection
    of X-bonds (except for F) and multipolar interactions (only for F).
    """
    selection = []
    obatoms_filtered = filter(
        lambda obatom: obatom.GetAtomicNum() in [9, 17, 35, 53], obatoms
    )
    for obatom in obatoms_filtered:
        obneighs = [a for a in OBAtomAtomIter(obatom)]
        if len(obneighs) == 1 and obneighs[0].GetAtomicNum() == 6:
            selection.append(XBondDonor(atom_list=[Atom(obneighs[0]), Atom(obatom)]))
    logger.debug(f'Found {len(selection)} halogen(s)')
    return selection


class PiCarbon(NamedTuple):
    """Class to store information about pi-carbons.

    Attributes:
        atom_list: central C atom (in a list, to make it consistent with the class Ring).
        normal: vector defining the normal to the pi-system.
        center: coordinates of the center of the pi-system.
    """

    atom_list: List[Atom]
    normal: npt.NDArray
    center: npt.NDArray


def find_pi_carbons(obatoms: List[OBAtom]) -> List[PiCarbon]:
    """Finds all pi-carbons, e.g. amide, guanidinium.

    Ester and acid are excluded for now, carbamate are included
    (at least 1 N amongst neighbours).
    """
    selection = []
    obatoms_filtered = filter(
        lambda obatom: obatom.GetAtomicNum() == 6
        and obatom.GetHyb() == 2
        and not obatom.IsAromatic(),
        obatoms,
    )
    for obatom in obatoms_filtered:
        obneighs = list(OBAtomAtomIter(obatom))
        obneighs_n = get_n_neighbours(obatom)
        obneighs_o = get_o_neighbours(obatom)
        if len(obneighs_n) + len(obneighs_o) >= 2 and len(obneighs_n) > 0:
            coords = [get_atom_coordinates(a) for a in obneighs]
            v1 = get_vector(coords[0], coords[1])
            v2 = get_vector(coords[2], coords[0])
            normal = normalize_vector(np.cross(v1, v2))
            center = get_atom_coordinates(obatom)
            selection.append(
                PiCarbon(
                    atom_list=[Atom(obatom)],
                    normal=normal,
                    center=center,
                )
            )
    logger.debug(f'Found {len(selection)} pi-carbon(s)')
    return selection


class Metal(NamedTuple):
    """Class to store information about metals.

    Attributes:
        atom_list: metal atom (in a list to be consistent with other classes).
        location: where the metal is located (ligand or receptor).
    """

    atom_list: List[Atom]
    location: str


def find_metals(obatoms: List[OBAtom], location: str) -> List[Metal]:
    """Finds all metal atoms."""
    selection = []
    obatoms_filtered = filter(
        lambda obatom: obatom.IsMetal() or obatom.GetType() in constants.METAL_IONS,
        obatoms,
    )
    for obatom in obatoms_filtered:
        selection.append(Metal(atom_list=[Atom(obatom)], location=location))
    logger.debug(f'Found {len(selection)} metal atom(s)')
    return selection


class MetalBinder(NamedTuple):
    """Class to store information about metal binders.

    Attributes:
        atom_list: atom with lone pair (in a list to be consistent with other classes).
        location: where the metal is located (ligand or receptor).
    """

    atom_list: List[Atom]
    location: str


def find_metal_binders(obatoms: List[OBAtom], location: str) -> List[MetalBinder]:
    """Finds atoms/groups that can be involved in binding metal ions.

    Metal binders include oxygens from carboxylate, phophoryl, phenolate, alcohol;
    nitrogen from imidazole; sulfur from thiolate.
    """
    selection = []
    obatoms_filtered = filter(lambda obatom: atom_has_lone_pair(obatom), obatoms)
    for obatom in obatoms_filtered:
        selection.append(MetalBinder(atom_list=[Atom(obatom)], location=location))
    logger.debug(f'Found {len(selection)} metal binder(s)')
    return selection


class ChargedGroup(NamedTuple):
    """Class to store information about charged groups.

    Attributes:
        atom_list: list of atoms that define the group.
        charge: positive or negative.
        center: coordinates of the center of the group.
        fgroup: name of the functional group.
    """

    atom_list: List[Atom]
    charge: str
    center: npt.NDArray
    fgroup: str


def find_charged_atoms(obatoms: List[OBAtom]) -> List[ChargedGroup]:
    """Finds  all charged atoms/groups.

    Ref: 'Cation-pi interactions in ligand recognition and catalysis'
         (Zacharias et al., 2002)
    """
    selection = []
    for obatom in obatoms:
        charge = None
        obres = obatom.GetResidue()
        if obres is not None:
            residue_name = obres.GetName().replace(' ', '')
            atom_name = obres.GetAtomID(obatom).replace(' ', '')
            # For protein residues, charged groups are tabulated
            # If the residue does not appear in the list of known charged residues,
            # then we check the charge as we do with ligands (note: a ligand can be
            # a peptide and appear in CHARGED_RESIDUES)
            if residue_name in constants.CHARGED_RESIDUES:
                for charged_group in constants.CHARGED_RESIDUES[residue_name]:
                    if atom_name == charged_group['on_atoms'][0]:
                        obatms = [
                            obatm
                            for obatm in OBResidueAtomIter(obres)
                            if obres.GetAtomID(obatm).replace(' ', '')
                            in charged_group['on_atoms']
                        ]
                        fgroup = charged_group['fgroup']
                        charge = charged_group['charge']
            elif obatom.GetFormalCharge() != 0:
                fgroup, charge, obatms = identify_charged_group(obatom)
        elif obatom.GetFormalCharge() != 0:
            fgroup, charge, obatms = identify_charged_group(obatom)
        if charge in ['positive', 'negative']:
            centroid = get_centroid([get_atom_coordinates(a) for a in obatms])
            atoms = [Atom(a) for a in obatms]
            selection.append(
                ChargedGroup(
                    atom_list=atoms, charge=charge, center=centroid, fgroup=fgroup
                )
            )
    logger.debug(f'Found {len(selection)} charged atom(s)')
    return selection
