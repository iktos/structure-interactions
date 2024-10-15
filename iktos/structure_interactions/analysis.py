from typing import List, Optional, Sequence, Union

try:
    from iktos.logger import getLogger
except ImportError:
    from logging import getLogger

from .detection import (
    find_h_bonds,
    find_hydrophobics,
    find_metal_complexes,
    find_multpipolar_interactions,
    find_pi_amides,
    find_pi_cations,
    find_pi_hydrophobics,
    find_pi_stackings,
    find_salt_bridges,
    find_water_bridges,
    find_x_bonds,
)
from .InteractionParameters import InteractionParameters
from .Ligand import Ligand
from .Receptor import Receptor
from .refinement import (
    drop_duplicated_h_bonds,
    drop_duplicated_water_bridges,
    refine_h_bonds,
    refine_hydrophobics,
    refine_pi_cations,
    refine_pi_hydrophobics,
    refine_water_bridges,
)

LOGGER = getLogger(__name__)


def _detect_interactions_inter(
    rec: Receptor,
    lig: Ligand,
    parameters: InteractionParameters,
    refine: bool,
) -> Sequence:
    """Detects and refines all intermolecular interactions."""
    LOGGER.debug("Analysing intermolecular interactions")
    hydrophobics_all = find_hydrophobics(
        hydrophobics_1=rec.hydrophobics,
        hydrophobics_2=lig.hydrophobics,
        distance_min=parameters.min_dist,
        distance_max=parameters.hydrophobic_dist_max,
    )
    LOGGER.debug(f"Found {len(hydrophobics_all)} hydrophobic interaction(s)")

    pi_stackings_all = find_pi_stackings(
        rings_1=rec.rings,
        rings_2=lig.rings,
        distance_min=parameters.min_dist,
        distance_max_t=parameters.pistacking_dist_max_t,
        distance_max_f=parameters.pistacking_dist_max_f,
        distance_max_p=parameters.pistacking_dist_max_p,
        angle_dev=parameters.pistacking_ang_dev,
        offset_max=parameters.pistacking_offset_max,
    )
    LOGGER.debug(f"Found {len(pi_stackings_all)} pi-stacking interaction(s)")

    pi_amides_all = find_pi_amides(
        rings=lig.rings,
        pi_carbons=rec.pi_carbons,
        distance_min=parameters.min_dist,
        distance_max=parameters.piother_dist_max,
        offset_max=parameters.piother_offset_max,
        angle_dev=parameters.pistacking_ang_dev,
    ) + find_pi_amides(
        rings=rec.rings,
        pi_carbons=lig.pi_carbons,
        distance_min=parameters.min_dist,
        distance_max=parameters.piother_dist_max,
        offset_max=parameters.piother_offset_max,
        angle_dev=parameters.pistacking_ang_dev,
    )
    LOGGER.debug(f"Found {len(pi_amides_all)} pi-amide like interaction(s)")

    pi_cations_all = find_pi_cations(
        rings=rec.rings,
        charged_atoms=lig.charged_atoms,
        distance_min=parameters.min_dist,
        distance_max=parameters.pication_dist_max,
        offset_max=parameters.pication_offset_max,
    ) + find_pi_cations(
        rings=lig.rings,
        charged_atoms=rec.charged_atoms,
        distance_min=parameters.min_dist,
        distance_max=parameters.pication_dist_max,
        offset_max=parameters.pication_offset_max,
    )
    LOGGER.debug(f"Found {len(pi_cations_all)} pi-cation interaction(s)")

    pi_hydrophobics_all = find_pi_hydrophobics(
        rings=rec.rings,
        hydrophobics=lig.hydrophobics,
        distance_min=parameters.min_dist,
        distance_max=parameters.piother_dist_max,
        offset_max=parameters.piother_offset_max,
    ) + find_pi_hydrophobics(
        rings=lig.rings,
        hydrophobics=rec.hydrophobics,
        distance_min=parameters.min_dist,
        distance_max=parameters.piother_dist_max,
        offset_max=parameters.piother_offset_max,
    )
    LOGGER.debug(f"Found {len(pi_hydrophobics_all)} pi-hydrophobic interaction(s)")

    h_bonds_all = find_h_bonds(
        acceptors=rec.h_bond_acceptors + rec.water_h_bond_acceptors,
        donor_pairs=lig.h_bond_donors,
        distance_min=parameters.min_dist,
        distance_max=parameters.hbond_dist_max,
        donor_angle_min=parameters.hbond_don_angle_min,
        acceptor_angle_min=parameters.hbond_acc_angle_min,
        allow_h_rotation=parameters.allow_h_rotation,
    ) + find_h_bonds(
        acceptors=lig.h_bond_acceptors,
        donor_pairs=rec.h_bond_donors + rec.water_h_bond_donors,
        distance_min=parameters.min_dist,
        distance_max=parameters.hbond_dist_max,
        donor_angle_min=parameters.hbond_don_angle_min,
        acceptor_angle_min=parameters.hbond_acc_angle_min,
        allow_h_rotation=parameters.allow_h_rotation,
    )
    if parameters.allow_h_rotation:
        # Drop duplicated H-bonds for e.g. R-NH2 (in cases like this,
        # both Hs can be detected in interaction with the same acceptor)
        h_bonds_all = drop_duplicated_h_bonds(h_bonds_all)
    LOGGER.debug(f"Found {len(h_bonds_all)} H-bond(s)")

    x_bonds_all = find_x_bonds(
        acceptors=rec.x_bond_acceptors + rec.water_x_bond_acceptors,
        donor_pairs=lig.halogens,
        distance_min=parameters.min_dist,
        distance_max=parameters.xbond_dist_max,
        donor_angle_min=parameters.xbond_don_angle_min,
        acceptor_angle_min=parameters.xbond_acc_angle_min,
    )
    LOGGER.debug(f"Found {len(x_bonds_all)} halogen bond(s)")

    mulipolar_all = find_multpipolar_interactions(
        acceptors=rec.pi_carbons,
        donor_pairs=lig.halogens,
        distance_min=parameters.min_dist,
        distance_max=parameters.multipolar_dist_max,
        donor_angle_min=parameters.multipolar_don_angle_min,
        norm_angle_max=parameters.multipolar_norm_angle_max,
    )
    LOGGER.debug(f"Found {len(mulipolar_all)} multipolar interaction(s)")

    salt_bridges_all = find_salt_bridges(
        charged_atoms_1=rec.charged_atoms,
        charged_atoms_2=lig.charged_atoms,
        distance_min=parameters.min_dist,
        distance_max=parameters.saltbridge_dist_max,
    )
    LOGGER.debug(f"Found {len(salt_bridges_all)} salt bridge(s)")

    water_bridges_all = find_water_bridges(
        acceptors=rec.h_bond_acceptors,
        donor_pairs=lig.h_bond_donors,
        waters=rec.water_metal_binders,
        distance_min=parameters.water_bridge_mindist,
        distance_max=parameters.water_bridge_maxdist,
        omega_min=parameters.water_bridge_omega_min,
        omega_max=parameters.water_bridge_omega_max,
        hbond_acceptor_angle_min=parameters.hbond_acc_angle_min,
        hbond_donor_angle_min=parameters.hbond_don_angle_min,
        allow_h_rotation=parameters.allow_h_rotation,
    ) + find_water_bridges(
        acceptors=lig.h_bond_acceptors,
        donor_pairs=rec.h_bond_donors,
        waters=rec.water_metal_binders,
        distance_min=parameters.water_bridge_mindist,
        distance_max=parameters.water_bridge_maxdist,
        omega_min=parameters.water_bridge_omega_min,
        omega_max=parameters.water_bridge_omega_max,
        hbond_acceptor_angle_min=parameters.hbond_acc_angle_min,
        hbond_donor_angle_min=parameters.hbond_don_angle_min,
        allow_h_rotation=parameters.allow_h_rotation,
    )
    if parameters.allow_h_rotation:
        # Drop duplicated water bridges for e.g. R-NH2 (in cases like this,
        # both Hs can be detected in interaction with the same acceptor)
        water_bridges_all = drop_duplicated_water_bridges(water_bridges_all)
    LOGGER.debug(f"Found {len(water_bridges_all)} water bridge(s)")

    metal_complexes = find_metal_complexes(
        metals=rec.metals + lig.metals,
        metal_binders=rec.metal_binders + rec.water_metal_binders + lig.metal_binders,
        distance_max=parameters.metal_dist_max,
    )
    LOGGER.debug(f"Found {len(metal_complexes)} metal complex(es)")

    if refine:
        LOGGER.debug("Refining selection")
        hydrophobics = refine_hydrophobics(
            hydrophobics_all, pi_stackings_all, pi_hydrophobics_all
        )
        LOGGER.debug(f"Kept {len(hydrophobics)} hydrophobic interaction(s)")

        pi_cations = refine_pi_cations(pi_cations_all, pi_stackings_all)
        LOGGER.debug(f"Kept {len(pi_cations)} pi-cation interaction(s)")

        pi_hydrophobics = refine_pi_hydrophobics(pi_hydrophobics_all)
        LOGGER.debug(f"Kept {len(pi_hydrophobics)} pi-hydrophobic interaction(s)")

        h_bonds = refine_h_bonds(
            h_bonds_all,
            salt_bridges_all,
            water_bridges_all,
            metal_complexes,
        )
        LOGGER.debug(f"Kept {len(h_bonds)} H-bond(s)")

        water_bridges = refine_water_bridges(water_bridges_all, metal_complexes)
        LOGGER.debug(f"Kept {len(water_bridges)} water bridge(s)")
        return (
            hydrophobics  # type: ignore [operator]
            + pi_stackings_all  # type: ignore [operator]
            + pi_amides_all  # type: ignore [operator]
            + pi_cations  # type: ignore [operator]
            + pi_hydrophobics  # type: ignore [operator]
            + h_bonds  # type: ignore [operator]
            + x_bonds_all  # type: ignore [operator]
            + mulipolar_all  # type: ignore [operator]
            + salt_bridges_all  # type: ignore [operator]
            + metal_complexes  # type: ignore [operator]
            + water_bridges  # type: ignore [operator]
        )

    return (
        hydrophobics_all  # type: ignore [operator]
        + pi_stackings_all  # type: ignore [operator]
        + pi_amides_all  # type: ignore [operator]
        + pi_cations_all  # type: ignore [operator]
        + pi_hydrophobics_all  # type: ignore [operator]
        + h_bonds_all  # type: ignore [operator]
        + x_bonds_all  # type: ignore [operator]
        + mulipolar_all  # type: ignore [operator]
        + salt_bridges_all  # type: ignore [operator]
        + metal_complexes  # type: ignore [operator]
        + water_bridges_all  # type: ignore [operator]
    )


def _detect_interactions_intra(
    mol: Union[Receptor, Ligand],
    parameters: InteractionParameters,
) -> Sequence:
    """Detects all intramolecular interactions.

    Note:
        Refine not implemented!
    """
    LOGGER.debug("Analysing intramolecular interactions")
    hydrophobics = find_hydrophobics(
        hydrophobics_1=mol.hydrophobics,
        hydrophobics_2=mol.hydrophobics,
        distance_min=parameters.min_dist,
        distance_max=parameters.hydrophobic_dist_max,
    )
    LOGGER.debug(f"Found {len(hydrophobics)} hydrophobic interaction(s)")

    pi_stackings = find_pi_stackings(
        rings_1=mol.rings,
        rings_2=mol.rings,
        distance_min=parameters.min_dist,
        distance_max_t=parameters.pistacking_dist_max_t,
        distance_max_f=parameters.pistacking_dist_max_f,
        distance_max_p=parameters.pistacking_dist_max_p,
        angle_dev=parameters.pistacking_ang_dev,
        offset_max=parameters.pistacking_offset_max,
    )
    LOGGER.debug(f"Found {len(pi_stackings)} pi-stacking interaction(s)")

    pi_amides = find_pi_amides(
        rings=mol.rings,
        pi_carbons=mol.pi_carbons,
        distance_min=parameters.min_dist,
        distance_max=parameters.piother_dist_max,
        offset_max=parameters.piother_offset_max,
        angle_dev=parameters.pistacking_ang_dev,
    )
    LOGGER.debug(f"Found {len(pi_amides)} pi-amide like interaction(s)")

    pi_hydroph = find_pi_hydrophobics(
        rings=mol.rings,
        hydrophobics=mol.hydrophobics,
        distance_min=parameters.min_dist,
        distance_max=parameters.piother_dist_max,
        offset_max=parameters.piother_offset_max,
    )
    LOGGER.debug(f"Found {len(pi_hydroph)} pi-hydrophobic interaction(s)")

    hbonds = find_h_bonds(
        acceptors=mol.h_bond_acceptors,
        donor_pairs=mol.h_bond_donors,
        distance_min=parameters.min_dist,
        distance_max=parameters.hbond_dist_max,
        donor_angle_min=parameters.hbond_don_angle_min,
        acceptor_angle_min=parameters.hbond_acc_angle_min,
        allow_h_rotation=parameters.allow_h_rotation,
    )
    if parameters.allow_h_rotation:
        # Drop duplicated H-bonds for e.g. R-NH2 (in cases like this,
        # both Hs can be detected in interaction with the same acceptor)
        hbonds = drop_duplicated_h_bonds(hbonds)
    LOGGER.debug(f"Found {len(hbonds)} H-bond(s)")

    salt_bridges = find_salt_bridges(
        charged_atoms_1=mol.charged_atoms,
        charged_atoms_2=mol.charged_atoms,
        distance_min=parameters.min_dist,
        distance_max=parameters.saltbridge_dist_max,
    )
    LOGGER.debug(f"Found {len(salt_bridges)} salt bridge(s)")

    return (
        hbonds  # type: ignore [operator]
        + hydrophobics  # type: ignore [operator]
        + pi_amides  # type: ignore [operator]
        + pi_hydroph  # type: ignore [operator]
        + pi_stackings  # type: ignore [operator]
        + salt_bridges  # type: ignore [operator]
    )


def analyse_interactions_inter(
    rec_coords: str,
    lig_coords: str,
    lig_format: str = "sdf",
    as_string: bool = True,
    refine: bool = True,
    parameters: InteractionParameters = InteractionParameters(),
) -> Sequence:
    """Analyses protein-ligand interactions for 1 complex.

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

    return _detect_interactions_inter(
        rec=receptor,
        lig=ligand,
        parameters=parameters,
        refine=refine,
    )


def analyse_interactions_inter_multi(
    rec_coords: str,
    lig_coords: List[str],
    lig_format: str = "sdf",
    as_string: bool = True,
    refine: bool = True,
    parameters: InteractionParameters = InteractionParameters(),
) -> List[Sequence]:
    """Analyses protein-ligand interactions for multiple complexes that share the same receptor.

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
        interactions = _detect_interactions_inter(
            rec=receptor,
            lig=ligand,
            parameters=parameters,
            refine=refine,
        )
        contacts.append(interactions)
    return contacts


def analyse_interactions_intra(
    coords: str,
    fmt: str = "sdf",
    is_small_molecule: bool = True,
    as_string: bool = True,
    parameters: InteractionParameters = InteractionParameters(),
    lig_coords: Optional[str] = None,
    lig_format: Optional[str] = None,
) -> Sequence:
    """Analyses intramolecular interactions.

    Notes:
       For now, the following interactions are detected: H-bonds, hydrophobic, pi-stacking,
           pi-amide, pi-hydrophobic, salt bridges; please extend as needed.
        The function can work with small molecules (default) or polymeric chains. The difference lies
            in how atom types are detected:
                * for a small molecule, the recommended format is SDF (MOL2 are supported but negatively
                  charged groups are often missed) and the atom typing is done based on the connectivity;
                * for polymeric chains (peptides or DNA), only PDB is supported and the atom typing
                  is partly done based on atom names and residue names (for the detection of rings
                  and charged groups).
                  It is possible to add a ligand to only detect intramolecular interactions in the binding site rather than in the whole protein.

    Args:
        coords: coords block or path to the coords file.
        fmt: format of the coords blocks (only relevant for small molecules, can be SDF or MOL2).
        is_small_molecule: whether to treat the input as a small molecule or a polymeric chain (default: True).
        as_string: whether to read coords in files or from blocks (default: True).
        parameters: cutoffs to use for the detection of interactions.
        lig_coords: optional, name of ligand multi coords file or coords block to detect only intramolecular interactions in the binding site.
        lig_format: optional, format of lig_coords (recommended format is 'sdf').

    Returns:
        list of contact objects.

    Raises:
        FileNotFoudError: if an input file does not exist (only when `as_string`=False).
        ValueError: if one of the input coords block or file is invalid
            (i.e. obmol is None or number of atoms or bonds == 0).
    """
    if is_small_molecule:
        mol: Union[Ligand, Receptor] = Ligand(coords, fmt, as_string=as_string)
    else:
        mol = Receptor(coords, as_string=as_string)
        if (lig_coords is not None) and (lig_format is not None):
            # Load ligand
            ligand = Ligand(lig_coords, lig_format, as_string=as_string)
            # Restrain search of intramolecular interactions only in the binding site of the protein, as defined by the previously loaded ligand
            mol.detect_binding_site([ligand], parameters.bs_dist)
        mol.identify_functional_groups()
    return _detect_interactions_intra(
        mol=mol,
        parameters=parameters,
    )
