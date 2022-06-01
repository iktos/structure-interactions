from typing import List, Sequence

from .InteractionParameters import InteractionParameters
from .InteractionProfiler import InteractionProfiler


def analyse_complex(
    rec_coords: str,
    lig_coords: str,
    lig_format: str = 'sdf',
    as_string: bool = True,
    refine: bool = True,
    parameters: InteractionParameters = InteractionParameters(),
) -> Sequence:
    """Performs protein-ligand interaction analysis for 1 complex.

    Warnings:
        The current version of the code fails to detect negatively charged groups
            (e.g. COO-) if ligand coords are given in MOL2 format.

    Args:
        rec_coords: name of receptor coords file or coords block (PDB format).
        lig_coords: name of ligand multi coords file or coords block.
        lig_format: format of lig_coords (recommended format is 'sdf').
        as_string: whether to read coords in files or from blocks.
        refine: if True, cleanup double-counted contacts (e.g. hydrophobics + pi-stacking).
        parameters: cutoffs to use for the detection of interactions.

    Returns:
        contact object.

    Raises:
        FileNotFoudError: if an input file does not exist (only when `as_string`=False).
        ValueError: if one of the input coords block or file is invalid
            (i.e. obmol is None or number of atoms or bonds == 0).
    """
    plip = InteractionProfiler()
    plip.load_receptor(rec_coords, as_string=as_string)
    plip.load_ligand(
        lig_coords,
        lig_format=lig_format,
        as_string=as_string,
        bs_dist=parameters.bs_dist,
    )
    return plip.analyse_interactions(refine=refine, parameters=parameters)


def analyse_complexes(
    rec_coords: str,
    lig_coords: List[str],
    lig_format: str = 'sdf',
    as_string: bool = True,
    refine: bool = True,
    parameters: InteractionParameters = InteractionParameters(),
) -> List[Sequence]:
    """Performs protein-ligand interaction analysis for multiple complexes
    that share the same receptor.

    Warnings:
        The current version of the code fails to detect negatively charged groups
            (e.g. COO-) if ligand coords are given in MOL2 format.

    Args:
        rec_coords: name of receptor coords file or coords block (PDB format).
        lig_coords: list of file paths or list of coords blocks.
        lig_format: format of lig_coords (recommended format is 'sdf').
        as_string: whether to read coords in files or from blocks.
        refine: if True, cleanup double-counted contacts (e.g. hydrophobics + pi-stacking).
        parameters: cutoffs to use for the detection of interactions.

    Returns:
        list of contact objects.

    Raises:
        FileNotFoudError: if an input file does not exist (only when `as_string`=False).
        ValueError: if one of the input coords block or file is invalid
            (i.e. obmol is None or number of atoms or bonds == 0).
    """
    plip = InteractionProfiler()
    plip.load_receptor(rec_coords, as_string=as_string)
    contacts = []  # type: List
    for lig in lig_coords:
        plip.load_ligand(
            lig,
            lig_format=lig_format,
            as_string=as_string,
            bs_dist=parameters.bs_dist,
        )
        contacts.append(plip.analyse_interactions(refine=refine, parameters=parameters))
    return contacts
