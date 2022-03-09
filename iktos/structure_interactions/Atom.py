from __future__ import absolute_import

from iktos.logger import getLogger

from .utils import get_coords

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
        self.atom_name = obatom.GetResidue().GetAtomID(obatom).strip(' ')
        self.is_aromatic = obatom.IsAromatic()
        self.hybridisation = obatom.GetHyb()
        self.residue_name = obatom.GetResidue().GetName()
        self.residue_num = obatom.GetResidue().GetNum()
        self.residue_chain = obatom.GetResidue().GetChain()
        self.residue_id = f'{self.residue_name}|{self.residue_chain}|{self.residue_num}'
        self.unique_id = f'{self.mol_title}|{self.residue_id}|{self.atom_id}'
        self.coords = get_coords(obatom)
