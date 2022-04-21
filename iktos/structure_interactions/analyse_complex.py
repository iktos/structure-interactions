from typing import List, Optional, Sequence

from .InteractionParameters import InteractionParameters
from .InteractionProfiler import InteractionProfiler


def analyse_complex(
    rec_coords: str,
    lig_coords: str,
    lig_format: str = 'sdf',
    as_string: bool = True,
    refine: bool = True,
    parameters: InteractionParameters = InteractionParameters(),
) -> Optional[Sequence]:
    """Performs protein-ligand interaction analysis for 1 complex.

    Note: The current version of the code fails to detect negatively charged groups
    (e.g. COO-) if ligand coords are given in MOL2 format.

    Args:
        rec_coords: name of receptor coords file or coords block (PDB format)
        lig_coords: name of ligand multi coords file or coords block
        lig_format: format of lig_coords (recommended format is 'sdf')
        as_string: whether to read coords in files or from blocks
        refine: if True, cleanup double-counted contacts (e.g. hydrophobics + pi-stacking)
        parameters: parameters of class InteractionParameters

    Returns:
        contacts in a formatted dict
    """
    plip = InteractionProfiler()
    status = plip.load_receptor(rec_coords, as_string=as_string)
    if not status:
        return None
    status = plip.load_ligand(
        lig_coords,
        lig_format=lig_format,
        as_string=as_string,
        bs_dist=parameters.bs_dist,
    )
    if not status:
        return None
    return plip.analyse_interactions(refine=refine, parameters=parameters)


def analyse_complexes(
    self,
    rec_coords: str,
    lig_coords: List[str],
    lig_format: str = 'sdf',
    as_string: bool = True,
    refine: bool = True,
    parameters: InteractionParameters = InteractionParameters(),
) -> Optional[Sequence]:
    """Performs protein-ligand interaction analysis
    for multiple complexes that share the same recpetor.

    Note: The current version of the code fails to detect negatively charged groups
    (e.g. COO-) if ligand coords are given in MOL2 format.

    Args:
        rec_coords: name of receptor coords file or coords block (PDB format)
        lig_coords: list of file paths or list of coords blocks
        lig_format: format of lig_coords (recommended format is 'sdf')
        as_string: whether to read coords in files or from blocks
        refine: if True, cleanup double-counted contacts (e.g. hydrophobics + pi-stacking)
        parameters: parameters of class InteractionParameters

    Returns:
        list of contacts in a formatted dict
    """

    plip = InteractionProfiler()
    status = plip.load_receptor(rec_coords, as_string=as_string)
    if not status:
        return None
    contacts = []  # type: List
    for lig in lig_coords:
        status = self._load_ligand(
            lig,
            lig_format=lig_format,
            as_string=as_string,
            bs_dist=parameters.bs_dist,
        )
        if not status:
            contacts.append({})
        else:
            contacts.append(
                self._analyse_interactions(refine=refine, parameters=parameters)
            )
    return contacts
