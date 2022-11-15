"""
Copyright (C) Iktos - All Rights Reserved
Unauthorized copying of this file, via any medium is strictly prohibited
Proprietary and confidential
"""

import errno
import os
import random
import tempfile
from typing import Any, Dict, List, Optional

import logging

from .constants import CONTACT_COLOR


logger = logging.getLogger(__name__)


try:
    from pymol import cmd as pymol_cmd, util, querying, CmdException  # noqa: F401
except ModuleNotFoundError:
    logger.warning('Module Pymol not found')


def get_colors() -> List:
    """Returns a list of available pymol colors.
    Function found on pymol wiki.

    Returns:
        list of colors (str, e.g. 'black')
    """
    pymol_color_list = []
    for tuplepair in querying.get_color_indices(''):
        pymol_color_list.append(tuplepair[0])
    pymol_color_list.sort()
    return pymol_color_list


def get_random_color() -> str:
    """Returns a random color from available pymol colors.
    Function found on pymol wiki.

    Returns:
        color (str, e.g. 'black') selected at random
    """
    return random.choice(get_colors())


def _show_binding_site(
    chain_id: str,
    label_protein: str,
    state: int,
    show_binding_site_as: str,
    residue_id: Optional[int] = None,
    atom_id: Optional[int] = None,
):
    """Shows residue as e.g. lines. Used to visualise residues involved
    in contacts or weights.

    The residue can be selected by its sequence number OR by giving the index
    of any atom of the residue, in which case the selection will be extended
    to the whole residue (with byres operator) - this is to support the special case
    of metal bridges.

    Args:
        chain_id: chain ID.
        label_protein: name of the protein object in Pymol.
        state: index of the Pymol state.
        show_binding_site_as: show binding site residues as (usually, lines).
        residue_id (optional): residue sequence number (the function takes either
            an atom index OR a residue index).
        atom_id (optional): atom index (the function takes either an atom index
            OR a residue index; note: if `residue_id` is defined, `atom_id` will be ignored).
    """
    if residue_id is not None:
        selection = f'resi {residue_id} and {label_protein} and state {state}'
    elif atom_id is not None:
        selection = f'byres id {atom_id} and {label_protein} and state {state}'
    else:
        raise ValueError('You need to provide a residue number or an atom ID')
    if chain_id != ' ':
        selection += f' and chain {chain_id}'
    pymol_cmd.show(show_binding_site_as, selection)


def draw_contacts(
    contacts: Dict[str, Any],
    label_protein: str,
    label_ligand: str,
    state: int = 1,
    show_binding_site_as: str = 'lines',
) -> None:
    """Draws all the contacts for one complex. Modifies an existing Pymol scene object.

    Args:
        contacts: dict with contacts found for one complex (from PLIP analysis)
        label_protein: name of the protein object in Pymol
        label_ligand: name of the ligand object in Pymol
        state: index of the Pymol state, useful for ensemble docking
        show_binding_site_as: show binding site residues as (default: `lines`)
    """
    labels = f'{label_ligand} '
    for contact_type in CONTACT_COLOR:
        labels += f'{contact_type}_{label_ligand} '
        if contact_type not in contacts:
            # To avoid bug with multistate sessions
            pymol_cmd.distance(
                f'{contact_type}_{label_ligand}',
                f'id 1 and {label_ligand} and state {state}',
                f'id 1 and {label_ligand} and state {state}',
            )
            continue
        for contact in contacts[contact_type]:
            _draw_contact_dash(
                contact,
                label_contact=f'{contact_type}_{label_ligand}',
                label_protein=label_protein,
                label_ligand=label_ligand,
                color=CONTACT_COLOR[contact_type],
                state=state,
            )
            chain_id, residue_id = contact['res_p'].split('|')[1:3]
            _show_binding_site(
                residue_id=int(residue_id),
                chain_id=chain_id,
                label_protein=label_protein,
                state=state,
                show_binding_site_as=show_binding_site_as,
            )
            if 'at_m' in contact:  # -> metal complex
                for atom_id in contact['at_p']:
                    _show_binding_site(
                        atom_id=atom_id,
                        chain_id=chain_id,
                        label_protein=label_protein,
                        state=state,
                        show_binding_site_as=show_binding_site_as,
                    )

    pymol_cmd.group(f'{label_ligand}_all', f'{labels}')
    pymol_cmd.disable(f'{label_ligand}_all')


def show_weights(
    weights: Dict[str, Dict[str, float]],
    label_protein: str,
    label_weights: str,
    state: int = 1,
    sphere_scale: float = 1.0,
    show_binding_site_as: str = 'lines',
):
    """Adds weighted spheres on protein atoms/centroids. Modifies an existing Pymol scene object.

    Args:
        weights: dict of weights (from ContactScorer)
        label_protein: name of the protein object in Pymol
            (object needs to exist in current session)
        label_weights: name of the weights object in Pymol
        state: index of the Pymol state, useful for ensemble docking
        sphere_scale: scale spheres used to show weights (default: 1.0)
        show_binding_site_as: show binding site residues as (default: `lines`)
    """
    labels = ''
    for contact_type, contacts in weights.items():
        if contact_type not in CONTACT_COLOR:
            logger.warning(f'Contact type {contact_type} not handled')
            continue
        for prot_id, weight in contacts.items():
            # To avoid bug with multistate sessions
            weight = max(weight, 0.1)
            # Add centroid and show as sphere
            # prot_id looks like 'VAL|A|296|H+N'
            chain_id, residue_id, at_names = prot_id.split('|')[1:4]
            if chain_id != ' ':
                selection = (
                    f'name {at_names} and resi {residue_id} '
                    f'and chain {chain_id} and {label_protein} '
                    f'and state {state}'
                )
            else:
                selection = (
                    f'name {at_names} and resi {residue_id} '
                    f'and {label_protein} and state {state}'
                )
            try:
                centroid = pymol_cmd.centerofmass(selection, state=state)
            except CmdException:
                raise ValueError(f'Invalid selection `{selection}`')
            label = f'{contact_type}_{residue_id}_{label_weights}'
            labels += label + ' '
            color = CONTACT_COLOR[contact_type]
            pymol_cmd.pseudoatom(
                object=label, pos=centroid, state=state, vdw=weight, color=color
            )
            pymol_cmd.show_as('spheres', label)
            # Show residues as e.g. lines
            _show_binding_site(
                residue_id=int(residue_id),
                chain_id=chain_id,
                label_protein=label_protein,
                state=state,
                show_binding_site_as=show_binding_site_as,
            )
    pymol_cmd.group(label_weights, f'{labels}')
    max_weight = 0.0
    for contact_type in weights:
        for contact in weights[contact_type]:
            max_weight = max(weights[contact_type][contact], max_weight)
    scale = sphere_scale / max_weight
    pymol_cmd.set('sphere_scale', scale, label_weights, state=state)
    pymol_cmd.disable(label_weights)


def initialise_session(start_from_scratch: bool = True, color_bg: str = 'grey'):
    """Initialises a Pymol session.

    Args:
        start_from_scratch: whether to start the session from scratch (delete all first)
        color_bg: background color (default: 'grey')
    """
    if start_from_scratch:
        logger.info('Starting a new Pymol session')
        pymol_cmd.delete('all')
    else:
        logger.info('Will continue with existing Pymol session (if any)')
    pymol_cmd.bg_color(color=color_bg)


def load_block(
    coords_block: str,
    label: str,
    format: str = 'mol2',
    color_molecule: Optional[str] = None,
    state: int = 1,
) -> None:
    """Loads a coords block to an existing Pymol session (if any, else creates a new session).

    Args:
        mol2_block: coords block
        format: format of coords block (allows: mol2, sdf, pdb)
            mol2 and sdf are loaded with discrete=1 to avoid connectivity issues
        label: name of molecule in Pymol
        color_molecule (optional): color molecule by atom
        state: index of the Pymol state, useful for ensemble docking
    """
    _id, path = tempfile.mkstemp(suffix=f'.{format}')
    try:
        with os.fdopen(_id, 'w') as tmp:
            tmp.write(coords_block)
        if format in ['mol2', 'sdf']:
            load_file(
                path,
                label=label,
                color_molecule=color_molecule,
                discrete=1,
                state=state,
            )
        else:
            load_file(
                path,
                label=label,
                color_molecule=color_molecule,
                discrete=-1,
                state=state,
            )
    finally:
        os.remove(path)


def load_file(
    filename: str,
    label: str,
    color_molecule: Optional[str] = None,
    discrete: int = -1,
    state: int = 1,
    wdir: str = './',
) -> None:
    """Loads a coord file to an existing Pymol session (if any, else creates a new session).

    Args:
        filename: name of a coords file (PDB, SDF, MOL2, etc)
        label: name of molecule in Pymol
        color_molecule (optional): color molecule by atom
        discrete: if = 0, load the file partially, assuming same connectivity as state 1
                  if = 1, load the entire file, states are therefore independent
                  default = -1, choose depending on the format (SDF -> 1, PDB and MOL2 -> 0)
        state: index of the Pymol state, useful for ensemble docking
        wdir: path to coords file directory (default: current directory)
    """
    logger.debug(f'Loading {filename}')
    filepath = os.path.join(wdir, filename)
    if not os.path.isfile(filepath):
        raise FileNotFoundError(errno.ENOENT, os.strerror(errno.ENOENT), filepath)
    pymol_cmd.load(filepath, label, state=state, discrete=discrete)
    if color_molecule:
        util = globals()['util']  # noqa: F811
        cba = getattr(util, color_molecule)
        cba(label)
    # Decrease size of spheres (metal ions/atoms)
    pymol_cmd.set('sphere_scale', 0.5, label, state=state)


def load_complex_from_blocks(
    protein_block: str,
    ligand_block: Optional[str],
    protein_format: str = 'pdb',
    ligand_format: str = 'mol2',
    contacts: Optional[Dict[str, Any]] = None,
    weights: Optional[Dict[str, Dict[str, float]]] = None,
    label_protein: str = 'Protein_0',
    label_ligand: str = 'Ligand_0',
    label_weights: str = 'Weights_0',
    color_protein: Optional[str] = None,
    color_ligand: Optional[str] = None,
    sphere_scale: float = 1.0,
    show_binding_site_as: str = 'lines',
    state: int = 1,
    wdir: str = './',
) -> None:
    """Loads a complex (pocket-ligand) to an existing Pymol session,
    with their contacts if given, starting from coords blocks.

    Args:
        protein_block: coords block of protein
        ligand_block (optional): coords block of ligand
        protein_format: format of protein block (default: pdb)
        ligand_format: format of ligand block (default: mol2)
        contacts (optional): dict of contacts (from PLIP analysis)
        weights: dict of contact weights (from ContactScorer)
        label_protein (optional): name for protein object in Pymol
        label_ligand (optional): name for ligand object in Pymol
        label_weights (optional): name of the weights object in Pymol
        color_protein: color of protein atoms (e.g. cbaw, default: standard coloring)
        color_ligand: color of ligand atoms (e.g. cbaw, default: standard coloring)
        sphere_scale: scale spheres used to show weights (default: 1.0)
        show_binding_site_as: show binding site residues as (default: `lines`)
        start_from_scratch: whether to start the session from scratch (delete all in the beginning)
    """
    logger.debug('Attempting to load a single protein-ligand complex')
    if protein_block:
        load_block(
            coords_block=protein_block,
            label=label_protein,
            color_molecule=color_protein,
            format=protein_format,
            state=state,
        )
        if weights:
            show_weights(
                weights=weights,
                label_protein=label_protein,
                label_weights=label_weights,
                state=state,
                sphere_scale=sphere_scale,
                show_binding_site_as=show_binding_site_as,
            )
        if ligand_block:
            load_block(
                coords_block=ligand_block,
                label=label_ligand,
                color_molecule=color_ligand,
                format=ligand_format,
                state=state,
            )
            if contacts is not None:
                logger.debug('Found contacts for protein-ligand complex')
                if not isinstance(contacts, Dict):
                    raise ValueError('Bad input: `contacts` should be a dict')
                draw_contacts(
                    contacts=contacts,
                    label_protein=label_protein,
                    label_ligand=label_ligand,
                    show_binding_site_as=show_binding_site_as,
                    state=state,
                )
    return None


def load_complex_from_files(
    protein_file_path: str,
    ligand_file_path: Optional[str],
    contacts: Optional[Dict[str, Any]] = None,
    weights: Optional[Dict[str, Dict[str, float]]] = None,
    label_protein: str = 'Protein_0',
    label_ligand: str = 'Ligand_0',
    label_weights: str = 'Weights_0',
    color_protein: Optional[str] = None,
    color_ligand: Optional[str] = None,
    sphere_scale: float = 1.0,
    show_binding_site_as: str = 'lines',
    state: int = 1,
    wdir: str = './',
) -> None:
    """Loads a complex (pocket-ligand) to an existing Pymol session, with their contacts if given.

    Args:
        protein_file_path: path to the reference protein file (e.g. PDB)
        ligand_file_path (optional): path to the reference ligand file (PDB, MOL2, SDF, etc)
        contacts (optional): dict of contacts (from PLIP analysis)
        weights: dict of contact weights (from ContactScorer)
        label_protein (optional): name for protein object in Pymol
        label_ligand (optional): name for ligand object in Pymol
        label_weights (optional): name of the weights object in Pymol
        color_protein: color of protein atoms (e.g. cbaw, default: standard coloring)
        color_ligand: color of ligand atoms (e.g. cbaw, default: standard coloring)
        sphere_scale: scale spheres used to show weights (default: 1.0)
        show_binding_site_as: show binding site residues as (default: `lines`)
        start_from_scratch: whether to start the session from scratch (delete all in the beginning)
        wdir: path to coords file directory (default: current directory)
    """
    logger.debug('Attempting to load a single protein-ligand complex')
    if contacts is not None and not isinstance(contacts, Dict):
        raise ValueError('Bad input: `contacts` should be a dict')
    load_file(
        protein_file_path, label_protein, color_molecule=color_protein, state=state
    )
    if ligand_file_path:
        load_file(
            ligand_file_path, label_ligand, color_molecule=color_ligand, state=state
        )
    if weights:
        show_weights(
            weights=weights,
            label_protein=label_protein,
            label_weights=label_weights,
            state=state,
            sphere_scale=sphere_scale,
            show_binding_site_as=show_binding_site_as,
        )
    if contacts is not None:
        logger.debug('Found contacts for protein-ligand complex')
        draw_contacts(
            contacts=contacts,
            label_protein=label_protein,
            label_ligand=label_ligand,
            show_binding_site_as=show_binding_site_as,
            state=state,
        )


def load_ligand_from_block(
    label_protein: str,
    mol2_block: str,
    contacts: Dict = None,
    label_ligand: str = 'Mol_',
    show_binding_site_as: str = 'lines',
    state: int = 1,
) -> None:
    """Loads a single MOL2 block (e.g. docking pose) to an existing Pymol session,
    with its contacts if given.

    Args:
        label_protein: name of protein object in Pymol (object needs to exist in current session)
        mol2_block: MOL2 block (e.g. docking pose)
        contacts (optional): dict of contacts (from PLIP analysis)
        label_ligand (optional): name for new object in Pymol
        show_binding_site_as: show binding site residues as (default: `lines`)
        state: index of the Pymol state, useful for ensemble docking
    """
    logger.debug('Attempting to load a single ligand/pose')
    if mol2_block is None:
        logger.debug('MOL2 block is None -> pass')
        return
    load_block(mol2_block, label=label_ligand, state=state)
    if contacts is not None:
        draw_contacts(
            contacts=contacts,
            label_protein=label_protein,
            label_ligand=label_ligand,
            show_binding_site_as=show_binding_site_as,
            state=state,
        )


def _add_centroids(
    contact: Dict[str, Any],
    label_ligand: str,
    label_protein: str,
    label_contact: str,
    label_centroid_ligand: str,
    label_centroid_protein: str,
    state: int,
) -> None:
    """Creates centroid objects on atoms involved in the contact
    (one centroid on the ligand and one on the protein). Uses atom IDs
    to select ligand atoms and (atom names, residue seq number, chain ID)
    to select protein atoms.

    Note: since PLIP analysis is done on SDF blocks for the ligand,
    atom names in the dict of contacts do not match that in MOL2 blocks,
    so we need to stick to atom IDs for the ligand.

    Args:
        contact: dict with one contact (from PLIP analysis)
        label_ligand: name of the ligand object in Pymol
        label_protein: name of the protein object in Pymol
        label_contact: name of the contact object in Pymol
        label_centroid_ligand: name of the ligand centroid objects in Pymol
        label_centroid_protein: name of the protein centroid objects in Pymol
        state: index of the Pymol state
    """
    # Select atoms on ligand side and add centroid
    sel1 = 'id ' + ' or id '.join([str(i) for i in contact['at_l']])
    sel1 = f'({sel1}) and {label_ligand} and state {state}'
    try:
        centroid = pymol_cmd.centerofmass(sel1, state=state)
        pymol_cmd.pseudoatom(label_centroid_ligand, pos=centroid, state=state)
    except CmdException:
        raise ValueError(f'Invalid selection `{sel1}`')

    # Select atoms on protein side (using residue seq nbs
    # and chain ID if available) and add centroid
    chain_id, residue_id = contact['res_p'].split('|')[1:3]
    sel2 = 'name ' + ' or name '.join(contact['at_name_p'])
    if chain_id != ' ':
        sel2 = f'resi {residue_id} and chain {chain_id} and ({sel2})'
    else:
        sel2 = f'resi {residue_id} and ({sel2})'
    sel2 = f'{sel2} and {label_protein} and state {state}'
    try:
        centroid = pymol_cmd.centerofmass(sel2, state=state)
        pymol_cmd.pseudoatom(label_centroid_protein, pos=centroid, state=state)
    except CmdException:
        raise ValueError(f'Invalid selection `{sel2}`')


def _draw_simple_contact(
    contact: Dict[str, Any],
    label_ligand: str,
    label_protein: str,
    label_contact: str,
    state: int,
) -> None:
    """Creates centroid objects on the ligand and the protein
    and draws a line between them.

    Args:
        contact: dict with one contact (from PLIP analysis)
        label_ligand: name of the ligand object in Pymol
        label_protein: name of the protein object in Pymol
        label_contact: name of the contact object in Pymol
        state: index of the Pymol state
    """
    _add_centroids(
        contact=contact,
        label_ligand=label_ligand,
        label_protein=label_protein,
        label_contact=label_contact,
        label_centroid_ligand='PiS1',
        label_centroid_protein='PiS2',
        state=state,
    )
    pymol_cmd.distance(label_contact, 'PiS1', 'PiS2')


def _draw_water_bridge(
    contact: Dict[str, Any],
    label_ligand: str,
    label_protein: str,
    label_contact: str,
    state: int,
) -> None:
    """Creates centroid objects on the ligand and the protein
    and draws lines between the ligand centroid, OW and the protein centroid.

    Args:
        contact: dict with one contact (from PLIP analysis)
        label_ligand: name of the ligand object in Pymol
        label_protein: name of the protein object in Pymol
        label_contact: name of the contact object in Pymol
        state: index of the Pymol state
    """
    _add_centroids(
        contact=contact,
        label_ligand=label_ligand,
        label_protein=label_protein,
        label_contact=label_contact,
        label_centroid_ligand='PiS1',
        label_centroid_protein='PiS2',
        state=state,
    )
    chain_id, residue_id = contact['res_w'].split('|')[1:3]
    if chain_id != ' ':
        sel3 = f"name OW and resi {residue_id} and chain {chain_id}"
    else:
        sel3 = f"name OW and resi {residue_id}"
    sel3 = f'{sel3} and {label_protein} and state {state}'
    try:
        pymol_cmd.distance(label_contact, 'PiS1', sel3)
        pymol_cmd.distance(label_contact, 'PiS2', sel3)
    except CmdException:
        raise ValueError(f'Invalid selection `{sel3}`')


def _draw_metal_bridge(
    contact: Dict[str, Any],
    label_ligand: str,
    label_protein: str,
    label_contact: str,
    state: int,
) -> None:
    """Draws a metal bridge, i.e. lines between each ligand and protein atom
    involved in the contact and the metal.

    Note: uses atom IDs for now - would be better to use atom names for the protein
    (not currently available in the dict of contacts).

    Args:
        contact: dict with one contact (from PLIP analysis)
        label_ligand: name of the ligand object in Pymol
        label_protein: name of the protein object in Pymol
        label_contact: name of the contact object in Pymol
        state: index of the Pymol state
    """
    sel1 = 'id ' + ' or id '.join([str(i) for i in contact['at_l']])
    sel1 = f'({sel1}) and {label_ligand} and state {state}'
    sel2 = 'id ' + ' or id '.join([str(i) for i in contact['at_p']])
    sel2 = f'({sel2}) and {label_protein} and state {state}'
    sel3 = f"id {contact['at_m'][0]}"
    sel3 = f'({sel3}) and {label_protein} and state {state}'
    try:
        pymol_cmd.distance(label_contact, sel1, sel3)
    except CmdException:
        raise ValueError(f'Invalid selection `{sel1}` or `{sel3}`')
    try:
        pymol_cmd.distance(label_contact, sel2, sel3)
    except CmdException:
        raise ValueError(f'Invalid selection `{sel2}` or `{sel3}`')


def _draw_contact_dash(
    contact: Dict[str, Any],
    label_ligand: str,
    label_protein: str,
    label_contact: str,
    color: str,
    state: int,
) -> None:
    """Draws a single contact between protein and ligand using a dotted line
    representation. Modifies an existing Pymol scene object.

    For M-L interactions, draws each M-L contact;
    otherwise, draws contact between centroids.
    For water bridges, draws P-WAT-L.

    Args:
        contact: dict with one contact (from PLIP analysis)
        label_ligand: name of the ligand object in Pymol
        label_protein: name of the protein object in Pymol
        label_contact: name of the contact object in Pymol
        color: color of the contact object in Pymol
        state: index of the Pymol state, useful for ensemble docking
    """
    if not contact['at_l'] or not contact['at_name_p']:
        logger.warning(f'Invalid contact dict `{contact}`')
        return None
    if 'at_m' in contact:  # -> metal complex
        _draw_metal_bridge(
            contact=contact,
            label_ligand=label_ligand,
            label_protein=label_protein,
            label_contact=label_contact,
            state=state,
        )
    elif 'at_owat' in contact:  # -> water bridge
        _draw_water_bridge(
            contact=contact,
            label_ligand=label_ligand,
            label_protein=label_protein,
            label_contact=label_contact,
            state=state,
        )
    else:
        _draw_simple_contact(
            contact=contact,
            label_ligand=label_ligand,
            label_protein=label_protein,
            label_contact=label_contact,
            state=state,
        )
    pymol_cmd.set('dash_color', color, label_contact, state=state)
    pymol_cmd.set('dash_gap', 0.2, label_contact, state=state)
    pymol_cmd.set('dash_radius', 0.08, label_contact, state=state)
    pymol_cmd.delete('PiS1')
    pymol_cmd.delete('PiS2')
    pymol_cmd.delete('PiS3')
