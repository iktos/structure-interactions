from __future__ import absolute_import

from builtins import filter
from collections import namedtuple
from typing import List, NamedTuple

import numpy as np
from iktos.logger import getLogger
from iktos.structure_utils.pdb.constants import CHARGED_RESIDUES

try:
    from openbabel.openbabel import (  # openbabel 3
        OBAtom,
        OBAtomAtomIter,
        OBMolAtomIter,
        OBResidueAtomIter,
    )
except ModuleNotFoundError:
    from openbabel import (
        OBAtom,
        OBAtomAtomIter,
        OBMolAtomIter,
        OBResidueAtomIter,
    )  # openbabel 2 (warning in structure-utils)

from . import constants
from .Atom import Atom
from .math_utils import get_centroid, get_vector, normalize_vector
from .utils import (
    atom_has_lone_pair,
    get_coords,
    get_h_neighbours,
    get_n_neighbours,
    get_o_neighbours,
    identify_charged_group,
    ring_is_aromatic,
    ring_is_planar,
)

logger = getLogger(__name__)


class Mol:
    def __init__(self, mol_type):
        """
        Class with functions to identify atomic properties
        (e.g. hydrophobic, charged, etc)
        """
        self.mol_type = mol_type

    def identify_functional_groups(self):
        raise NotImplementedError

    def find_rings(self, obmol) -> List[NamedTuple]:
        """Finds aromatic rings, i.e. rings detected by OpenBabel
        as aromatic or sufficiently planar.

        Atom list contains all atoms in ring
        """
        data = namedtuple('aromatic_ring', 'atom_list normal center type')
        aromatic_amino = ['TYR', 'TRP', 'HIS', 'HID', 'HIE', 'HIP', 'HSD', 'HSE', 'PHE']
        all_rings = obmol.GetSSSR()
        selection = []
        obmol_atoms = list(OBMolAtomIter(obmol))
        for obring in all_rings:
            ring_atoms = [a for a in obmol_atoms if obring.IsMember(a)]
            obresidue = ring_atoms[0].GetResidue().GetName()
            if len(ring_atoms) not in [5, 6]:
                continue  # ignore rings of unusual size
            if (
                obresidue in aromatic_amino
                or obring.IsAromatic()
                or ring_is_aromatic(ring_atoms)
                or ring_is_planar(obring, ring_atoms)
            ):
                # Tag atoms in ring as aromatic (for later)
                for a in ring_atoms:
                    a.SetAromatic()
                atoms = [Atom(a) for a in ring_atoms]
                type = f'{len(ring_atoms)}-membered'
                coords = [get_coords(ring_atoms[i]) for i in [0, 2, 4]]
                ringv1 = get_vector(coords[0], coords[1])
                ringv2 = get_vector(coords[2], coords[0])
                normal = normalize_vector(np.cross(ringv1, ringv2))
                center = get_centroid([get_coords(ra) for ra in ring_atoms])
                selection.append(
                    data(atom_list=atoms, type=type, normal=normal, center=center)
                )
        logger.debug(f'Found {len(selection)} aromatic ring(s)')
        return selection

    def find_hydrophobics(self, obatoms: List[OBAtom]) -> List[NamedTuple]:
        """Finds all hydrophobic atoms, i.e. carbon and sulfur atoms
        which have only carbons/sulfurs/hydrogens as direct neighbors
        + Cl and F. Sulfur atoms that are not SP3 are excluded.

        Atom list contains [C|S|Cl|F]
        """
        data = namedtuple('hydrophobic', 'atom_list')
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
                logger.debug(f'Excluding S{obatom.GetId()} from list of hydrophobics')
                continue
            selection.append(data(atom_list=[Atom(obatom)]))
        logger.debug(f'Found {len(selection)} hydrophobic atoms')
        return selection

    def find_h_bond_acceptors(self, obatoms: List[OBAtom]) -> List[NamedTuple]:
        """Finds H-bond acceptors; excluding halogens and sulfur.
        Atom list contains [A]
        Stores information about A's neighbours to check angle about A
        """
        data = namedtuple('hbond_acceptor', 'atom_list neighbours')
        selection = []
        obatoms_filtered = filter(
            lambda obatom: (obatom.IsHbondAcceptor() or atom_has_lone_pair(obatom))
            and obatom.GetAtomicNum() not in [9, 17, 35, 53, 16],
            obatoms,
        )
        for obatom in obatoms_filtered:
            obneighs = list(OBAtomAtomIter(obatom))
            selection.append(
                data(atom_list=[Atom(obatom)], neighbours=[Atom(a) for a in obneighs])
            )
        logger.debug(f'Found {len(selection)} H-bond acceptors')
        return selection

    def find_h_bond_donors(self, obatoms: List[OBAtom]) -> List[NamedTuple]:
        """Finds strong and weak H-bond donors (i.e. with polarised C-H bonds).

        Atom list contains [D], [H]
        """
        data = namedtuple('hbond_donor', 'atom_list type')
        selection = []
        # Loop over atoms flagged by OpenBabel as H-bond donor
        obatoms_filtered = filter(lambda obatom: obatom.IsHbondDonor(), obatoms)
        for obatom in obatoms_filtered:
            for obneigh in OBAtomAtomIter(obatom):
                if obneigh.IsHbondDonorH():
                    selection.append(
                        data(atom_list=[Atom(obatom), Atom(obneigh)], type='strong')
                    )
        logger.debug(f'Found {len(selection)} strong H-bond donor pairs')
        # Find polarised C-H (near O, or near N, or SP2 or SP) and S-H
        for obatom in obatoms:
            if obatom.GetAtomicNum() == 16:
                obhs = get_h_neighbours(obatom)
                for obh in obhs:
                    selection.append(
                        data(atom_list=[Atom(obatom), Atom(obh)], type='weak')
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
                            data(atom_list=[Atom(obatom), Atom(obh)], type='weak')
                        )
        logger.debug(f'Found {len(selection)} H-bond donor pairs (incl. weak)')
        return selection

    def find_x_bond_acceptors(self, obatoms: List[OBAtom]) -> List[NamedTuple]:
        """Finds X-bond acceptors (look for O, N, S with lone pair).

        Atom list contains [A]
        Store information about A's neighbours to check angle later
        """
        data = namedtuple('xbond_acceptor', 'atom_list neighbours')
        selection = []
        obatoms_filtered = filter(
            lambda obatom: obatom.GetAtomicNum() in [7, 8, 16], obatoms
        )
        for obatom in obatoms_filtered:
            if atom_has_lone_pair(obatom):
                obneighs = [a for a in OBAtomAtomIter(obatom)]
                selection.append(
                    data(
                        atom_list=[Atom(obatom)], neighbours=[Atom(a) for a in obneighs]
                    )
                )
        logger.debug(f'Found {len(selection)} halogen bond acceptor(s)')
        return selection

    def find_halogens(self, obatoms: List[OBAtom]) -> List[NamedTuple]:
        """Finds halogens (C-X, with X = F, Cl, Br, I) for the detection of halogen
        bonds (except for F) and multipolar interactions (only for F).

        Atom list contains [C], [X]
        """
        data = namedtuple('halogen', 'atom_list')
        selection = []
        obatoms_filtered = filter(
            lambda obatom: obatom.GetAtomicNum() in [9, 17, 35, 53], obatoms
        )
        for obatom in obatoms_filtered:
            obneighs = [a for a in OBAtomAtomIter(obatom)]
            if len(obneighs) == 1 and obneighs[0].GetAtomicNum() == 6:
                selection.append(data(atom_list=[Atom(obneighs[0]), Atom(obatom)]))
        logger.debug(f'Found {len(selection)} halogen(s)')
        return selection

    def find_pi_carbons(self, obatoms: List[OBAtom]) -> List[NamedTuple]:
        """Finds pi carbons, e.g. amide, guanidinium.
        Ester and acid are excluded for now, carbamate are included
        (at least 1 N amongst neighbours).

        Atom list contains the central C atom
        Store the same info as for aromatic rings (center + normal)
        """
        data = namedtuple('pi_group', 'atom_list normal center type')
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
                coords = [get_coords(a) for a in obneighs]
                v1 = get_vector(coords[0], coords[1])
                v2 = get_vector(coords[2], coords[0])
                normal = normalize_vector(np.cross(v1, v2))
                center = get_coords(obatom)
                selection.append(
                    data(
                        atom_list=[Atom(obatom)],
                        type='amide/ester/guanidinium',
                        normal=normal,
                        center=center,
                    )
                )
        logger.debug(f'Found {len(selection)} pi-carbon(s)')
        return selection

    def find_metals(self, obatoms: List[OBAtom]) -> List[NamedTuple]:
        """Finds metals.

        Atom list contains [M]
        """
        data = namedtuple('metal', 'atom_list location')
        selection = []
        obatoms_filtered = filter(
            lambda obatom: obatom.IsMetal() or obatom.GetType() in constants.METAL_IONS,
            obatoms,
        )
        for obatom in obatoms_filtered:
            selection.append(data(atom_list=[Atom(obatom)], location=self.mol_type))
        logger.debug(f'Found {len(selection)} metal atom(s)')
        return selection

    def find_metal_binders(self, obatoms: List[OBAtom]) -> List[NamedTuple]:
        """Finds atoms/groups that can be involved in binding metal ions
        e.g. oxygen from carboxylate, phophoryl, phenolate, alcohol;
        nitrogen from imidazole; sulfur from thiolate.

        Atom list contains [D] (only 1 atom per binder)
        """
        data = namedtuple('metal_binder', 'atom_list location')
        selection = []
        obatoms_filtered = filter(lambda obatom: atom_has_lone_pair(obatom), obatoms)
        for obatom in obatoms_filtered:
            selection.append(data(atom_list=[Atom(obatom)], location=self.mol_type))
        logger.debug(f'Found {len(selection)} metal binder(s)')
        return selection

    def find_charged_atoms(self, obatoms: List[OBAtom]) -> List[NamedTuple]:
        """Finds charged atoms/groups.

        Ref: 'Cation-pi interactions in ligand recognition and catalysis'
             (Zacharias et al., 2002)

        Atom list contains atoms in charged group, in no particular order
        """
        data = namedtuple('charged', 'atom_list charge center fgroup')
        selection = []
        for obatom in obatoms:
            charge = None
            obres = obatom.GetResidue()
            residue_name = obres.GetName().replace(' ', '')
            atom_name = obres.GetAtomID(obatom).replace(' ', '')
            if residue_name in CHARGED_RESIDUES:
                for charged_group in CHARGED_RESIDUES[residue_name]:
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
            if charge in ['positive', 'negative']:
                centroid = get_centroid([get_coords(a) for a in obatms])
                atoms = [Atom(a) for a in obatms]
                selection.append(
                    data(atom_list=atoms, charge=charge, center=centroid, fgroup=fgroup)
                )
        logger.debug(f'Found {len(selection)} charged atom(s)')
        return selection
