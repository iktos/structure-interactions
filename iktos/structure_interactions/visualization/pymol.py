"""
Copyright (C) Iktos - All Rights Reserved
Unauthorized copying of this file, via any medium is strictly prohibited
Proprietary and confidential
"""

from typing import Any, Dict, List, Optional

try:
    from iktos.logger import getLogger
except ImportError:
    from logging import getLogger


from .utils import (
    draw_contacts_intra,
    initialise_session,
    load_block,
    load_complex_from_blocks,
    load_ligand_from_block,
)

LOGGER = getLogger(__name__)


try:
    from pymol import cmd as pymol_cmd
except ModuleNotFoundError:
    pass


def prepare_session_inter(
    protein_pdb_block: str,
    ligand_sdf_block: Optional[str] = None,
    contacts: Optional[Dict[str, Any]] = None,
    weights: Optional[Dict[str, Dict[str, float]]] = None,
    extra_sdf_blocks: Optional[List[str]] = None,
    extra_contacts: Optional[List[Dict[str, Any]]] = None,
    label_protein: str = 'Protein_0',
    label_ligand: str = 'Ligand_0',
    label_weights: str = 'Weights_0',
    labels_extra: Optional[List[str]] = None,
    color_protein: Optional[str] = None,
    color_ligand: Optional[str] = None,
    color_bg: str = 'grey',
    sphere_scale: float = 0.8,
    show_binding_site_as: str = 'lines',
    hide_non_polar_hs: bool = True,
    output_file_path: Optional[str] = None,
    start_from_scratch: bool = True,
) -> None:
    """Prepares a standard Pymol session, i.e. with all the molecules loaded in the same state.

    Takes one protein-ligand complex and a list of extra ligands (e.g. docking poses).
    If corresponding contacts are given, includes them in the session.
    If weights are given, adds spheres on protein atoms involved in important contacts; radius <-> weight.

    Args:
        protein_pdb_block: PDB block of reference protein
        ligand_sdf_block (optional): SDF block of reference ligand
        contacts (optional): dict with contacts for reference protein-ligand complex (from PLIP analysis)
        weights: dict of contact weights (from ContactScorer)
        extra_sdf_blocks (optional): list of sdf blocks (e.g. poses from docking)
        extra_contacts (optional): list of dicts with contacts for extra protein-ligand complexes (from PLIP analysis)
        label_protein (optional): name of protein object in Pymol
        label_ligand (optional): name of ligand object in Pymol
        label_weights (optional): name of the weights object in Pymol
        labels_extra (optional): list of names for extra ligand objects in Pymol
        color_protein: color of protein atoms (e.g. cbaw, default: standard coloring)
        color_ligand: color of ligand atoms (e.g. cbaw, default: standard coloring)
        color_bg: background color (default: `grey`)
        sphere_scale: scale spheres used to show weights (default: 1.0)
        show_binding_site_as: show binding site residues as (default: `lines`)
        hide_non_polar_hs: whether to hide non-polar Hs (default: True)
        output_file_path: path and name of final Pymol session (.pse, default: None -> session not saved)
        start_from_scratch: whether to start the session from scratch (delete all in the beginning)
    """
    initialise_session(start_from_scratch=start_from_scratch, color_bg=color_bg)
    load_complex_from_blocks(
        protein_block=protein_pdb_block,
        ligand_block=ligand_sdf_block,
        protein_format='pdb',
        ligand_format='sdf',
        contacts=contacts,
        weights=weights,
        label_protein=label_protein,
        label_ligand=label_ligand,
        label_weights=label_weights,
        color_protein=color_protein,
        color_ligand=color_ligand,
        sphere_scale=sphere_scale,
        show_binding_site_as=show_binding_site_as,
    )
    if extra_sdf_blocks:
        if not labels_extra:  # can be None or []
            labels_extra = [f'Mol_{x}' for x in range(len(extra_sdf_blocks))]
        for i, sdf_block in enumerate(extra_sdf_blocks):
            if extra_contacts is None:
                single_contacts_dict = None
            else:
                single_contacts_dict = extra_contacts[i]
            load_ligand_from_block(
                label_protein=label_protein,
                coords_block=sdf_block,
                contacts=single_contacts_dict,
                label_ligand=labels_extra[i],
                show_binding_site_as=show_binding_site_as,
            )
    if hide_non_polar_hs:
        pymol_cmd.hide('everything', 'h. and (e. c extend 1)')
    if output_file_path:
        LOGGER.info(f'Saving Pymol session in `{output_file_path}`')
        pymol_cmd.save(output_file_path)


def prepare_session_inter_multistate(
    protein_pdb_blocks: List[str],
    ligand_sdf_blocks: Optional[List[str]] = None,
    contacts: Optional[List[Dict[str, Any]]] = None,
    weights: Optional[Dict[str, Dict[str, float]]] = None,
    extra_sdf_blocks: Optional[List[List[str]]] = None,
    extra_contacts: Optional[List[List[Dict[str, Any]]]] = None,
    label_protein: str = 'Protein_0',
    label_ligand: str = 'Ligand_0',
    label_weights: str = 'Weights_0',
    labels_extra: Optional[List[str]] = None,
    color_protein: Optional[str] = None,
    color_ligand: Optional[str] = None,
    color_bg: str = 'grey',
    sphere_scale: float = 0.8,
    show_binding_site_as: str = 'lines',
    hide_non_polar_hs: bool = True,
    output_file_path: Optional[str] = None,
    start_from_scratch: bool = True,
) -> None:
    """Prepares a multistate Pymol session, i.e. with each molecule loaded in a different state (movie).

    Takes a list of protein-ligand complexes (e.g. pockets for ensemble docking) and a list of lists
    with e.g. docking poses (each sublist should be the same size as filenames_protein).
    If corresponding contacts are given, includes them in the session.
    If weights are given, adds spheres on protein atoms involved in important contacts; radius <-> weight.

    Args:
        protein_pdb_blocks: PDB blocks of reference protein
        ligand_sdf_blocks (optional): SDF blocks of reference ligand
        contacts (optional): list of dicts with contacts for reference protein-ligand complexes (from PLIP analysis)
        weights: dict of contact weights (from ContactScorer)
        extra_sdf_blocks (optional): list of lists of sdf blocks (e.g. poses from ensemble docking)
        extra_contacts (optional): list of lists of dicts with contacts for extra protein-ligand complexes (from PLIP analysis)
        label_protein (optional): name of protein object in Pymol
        label_ligand (optional): name of ligand object in Pymol
        label_weights (optional): name of the weights object in Pymol
        labels_extra (optional): list of names for extra ligand objects in Pymol (1 name per ligand, not pose)
        color_protein: color of protein atoms (e.g. cbaw, default: standard coloring)
        color_ligand: color of ligand atoms (e.g. cbaw, default: standard coloring)
        color_bg: background color (default: `grey`)
        sphere_scale: scale spheres used to show weights (default: 1.0)
        show_binding_site_as: show binding site residues as (default: `lines`)
        hide_non_polar_hs: whether to hide non-polar Hs (default: True)
        output_file_path: path and name of final Pymol session (.pse, default: None -> session not saved)
        start_from_scratch: whether to start the session from scratch (delete all in the beginning)
    """
    initialise_session(start_from_scratch=start_from_scratch, color_bg=color_bg)

    for j, protein_pdb_block in enumerate(protein_pdb_blocks):
        if ligand_sdf_blocks is None:
            ligand_sdf_block = None
        else:
            ligand_sdf_block = ligand_sdf_blocks[j]
        if contacts is None:
            single_contacts_dict = None
        else:
            single_contacts_dict = contacts[j]
        load_complex_from_blocks(
            protein_block=protein_pdb_block,
            ligand_block=ligand_sdf_block,
            protein_format='pdb',
            ligand_format='sdf',
            contacts=single_contacts_dict,
            weights=weights,
            label_protein=label_protein,
            label_ligand=label_ligand,
            label_weights=label_weights,
            color_protein=color_protein,
            color_ligand=color_ligand,
            sphere_scale=sphere_scale,
            show_binding_site_as=show_binding_site_as,
            state=j + 1,
        )
    if extra_sdf_blocks:
        if not labels_extra:  # can be None or []
            labels_extra = [f'Mol_{x}' for x in range(len(extra_sdf_blocks))]
        for i, sdf_blocks in enumerate(extra_sdf_blocks):
            for j, sdf_block in enumerate(sdf_blocks):
                state = j + 1
                if extra_contacts is None:
                    single_contacts_dict = None
                else:
                    single_contacts_dict = extra_contacts[i][j]
                load_ligand_from_block(
                    coords_block=sdf_block,
                    contacts=single_contacts_dict,
                    label_protein=label_protein,
                    label_ligand=labels_extra[i],
                    show_binding_site_as=show_binding_site_as,
                    state=state,
                )
    if hide_non_polar_hs:
        pymol_cmd.hide('everything', 'h. and (e. c extend 1)')
    if output_file_path:
        LOGGER.info(f'Saving Pymol session in `{output_file_path}`')
        pymol_cmd.save(output_file_path)


def prepare_session_intra(
    coords_block: str,
    fmt: str,
    contacts: Dict[str, Any],
    color_molecule: Optional[str] = None,
    color_bg: str = 'grey',
    hide_non_polar_hs: bool = True,
    output_file_path: Optional[str] = None,
    start_from_scratch: bool = True,
) -> None:
    """Prepares a standard Pymol session with intramolecular contacts.

    Args:
        coords_block: coords block.
        fmt: format of coords block.
        contacts: dict with intramolecular contacts.
        color_molecule: color of atoms (e.g. cbaw, default: standard coloring).
        color_bg: background color (default: `grey`).
        hide_non_polar_hs: whether to hide non-polar Hs (default: True).
        output_file_path: name of final Pymol session (.pse, default: None -> session not saved).
        start_from_scratch: whether to start the session from scratch (delete all in the beginning).
    """
    initialise_session(start_from_scratch=start_from_scratch, color_bg=color_bg)
    load_block(
        coords_block=coords_block,
        fmt=fmt,
        label='Mol',
        color_molecule=color_molecule,
    )
    draw_contacts_intra(
        contacts=contacts,
        label='Mol',
    )
    pymol_cmd.show('lines', '*')
    if hide_non_polar_hs:
        pymol_cmd.hide('everything', 'h. and (e. c extend 1)')
    if output_file_path:
        LOGGER.info(f'Saving Pymol session in `{output_file_path}`')
        pymol_cmd.save(output_file_path)
