from typing import List, Sequence

from .InteractionParameters import InteractionParameters
from .InteractionProfiler import InteractionProfiler
from .Ligand import Ligand
from .Receptor import Receptor


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
    # Load ligand
    ligand = Ligand(lig_coords, lig_format, as_string=as_string)

    # Load receptor and use the ligand to define the binding site
    receptor = Receptor(rec_coords, as_string=as_string)
    receptor.detect_binding_site([ligand], parameters.bs_dist)
    receptor.identify_functional_groups()

    # Init interaction profiler
    plip = InteractionProfiler(receptor, ligand, parameters)
    if refine:
        plip.detect_intermolecular_interactions()
        return plip.refine_intermolecular_interactions()
    return plip.detect_intermolecular_interactions()


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
    # Load all the ligands
    ligands = [Ligand(coords, lig_format, as_string=as_string) for coords in lig_coords]

    # Load receptor and use the ligands to define the binding site
    receptor = Receptor(rec_coords, as_string=as_string)
    receptor.detect_binding_site(ligands, parameters.bs_dist)
    receptor.identify_functional_groups()

    # Init interaction profiler
    contacts = []  # type: List
    for ligand in ligands:
        plip = InteractionProfiler(receptor, ligand, parameters)
        if refine:
            plip.detect_intermolecular_interactions()
            interactions = plip.refine_intermolecular_interactions()
        else:
            interactions = plip.detect_intermolecular_interactions()
        contacts.append(interactions)
    return contacts
