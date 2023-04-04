from __future__ import absolute_import

try:
    from iktos.logger import getLogger
except ImportError:
    from logging import getLogger

from .mol_utils import get_atom_name, get_atom_coordinates

logger = getLogger(__name__)


class Atom:
    """
    Class to store atoms and their properties
    --> to avoid openbabel objects as much as possible
    (they are useful to find atoms/groups that have specific
    properties, e.g. hydrophobic, but are a pain when dealing with coordinates)
    """

    def __init__(self, obatom):
        self.mol_title = obatom.GetParent().GetTitle()
        self.atomic_num = obatom.GetAtomicNum()
        self.atom_id = obatom.GetId()

        # Get atom name: for ligands, concatenate element + ID
        # because the property does not exist in SDF files,
        # for protein/water, use `GetAtomID` (field in PDB files)
        if self.mol_title == 'ligand':
            self.atom_name = get_atom_name(obatom)
        else:
            self.atom_name = obatom.GetResidue().GetAtomID(obatom).strip(' ')
        self.is_aromatic = obatom.IsAromatic()
        self.hybridisation = obatom.GetHyb()
        self.heavy_degree = obatom.GetHvyDegree()

        # Get info relative to residue - for ligands, the residue might not exist
        # (especially if input format is SDF) so we assign generic descriptors
        # (creating an OBResidue has been tested and leads to various problems, segfault, etc)
        # Note: proteins SHOULD have the residues initialised by OpenBabel, the code will crash
        # if they don't, that's what we want
        # Note 2: ligands can also have the residue initialised by OpenBabel and in that case
        # we want to use this info (e.g. if input format is MOL2, useful for peptides)
        if self.mol_title == 'ligand' and obatom.GetResidue() is None:
            self.residue_name = 'UNL'
            self.residue_num = 1
            self.residue_chain = ' '
        else:
            self.residue_name = obatom.GetResidue().GetName()
            self.residue_num = obatom.GetResidue().GetNum()
            self.residue_chain = obatom.GetResidue().GetChain()
        self.residue_id = f'{self.residue_name}|{self.residue_chain}|{self.residue_num}'
        self.unique_id = f'{self.mol_title}|{self.residue_id}|{self.atom_id}'
        self.coords = get_atom_coordinates(obatom)
