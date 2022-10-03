from __future__ import absolute_import

from typing import Any, Dict, List, Sequence


# List of attributes in interaction classes that correspond
# to the partners involved in the interaction
PARTNERS = [
    'receptor',
    'ligand',
    'partner_1',
    'partner_2',
    'ring',
    'cation',
    'pi_carbon',
    'hydrophobic',
    'donor',
    'acceptor',
]


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


def contacts_to_dict(plip_list: Sequence) -> Dict[str, List[Dict[str, Any]]]:
    """Converts detected contacts to dict of contacts as expected
    by functions in SA.

    Args:
        plip_list: list of contacts detected by `InteractionProfiler`.
    """
    plip_dict: Dict[str, Any] = {}
    for contact in plip_list:
        contact_name = type(contact).__name__
        if contact_name == 'H_Bond':
            if contact.type == 'weak':
                contact_name += '_Weak'
        if contact_name not in plip_dict:
            plip_dict[contact_name] = []
        parsed_contact = {}
        for field in contact._fields:
            if field in PARTNERS:
                # Check if partner is receptor (protein) or ligand
                atom_list = getattr(contact, field)
                if atom_list[0].mol_title == 'receptor':
                    parsed_contact['at_p'] = [a.atom_id for a in atom_list]
                    parsed_contact['at_name_p'] = [a.atom_name for a in atom_list]
                    if contact_name != 'Metal_Complex':
                        parsed_contact['res_p'] = atom_list[0].residue_id
                else:
                    parsed_contact['at_l'] = [a.atom_id for a in atom_list]
                    parsed_contact['at_name_l'] = [a.atom_name for a in atom_list]
            elif field == 'water':
                atom_list = getattr(contact, field)
                parsed_contact['at_owat'] = [a.atom_id for a in atom_list]
                parsed_contact['res_w'] = atom_list[0].residue_id
            elif field == 'metal':
                atom_list = getattr(contact, field)
                parsed_contact['at_m'] = [a.atom_id for a in atom_list]
                parsed_contact['at_name_m'] = [a.atom_name for a in atom_list]
                parsed_contact['res_p'] = atom_list[0].residue_id
            else:
                parsed_contact[field] = getattr(contact, field)
        plip_dict[contact_name].append(parsed_contact)
    return plip_dict
