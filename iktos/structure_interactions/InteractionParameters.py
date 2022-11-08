from dataclasses import dataclass


@dataclass
class InteractionParameters:
    """Thresholds for the detection of interactions.

    Attributes:
        bs_dist: max dist to include binding site residues (default=7.5)
        min_dist: min dist for all distance thresholds (default=2.0)
        hydrophobic_dist_max: max dist for detection of hydrophobic contacts (default=4.0)
        hbond_dist_max: max dist between d and a (hubbard & haider, 2001) + 0.5 a (default=4.0)
        hbond_don_angle_min: minimal angle for d-h-a (default=140)
        hbond_acc_angle_min: minimal angle for each h-a-y, where y are a's neighbours (custom) (default=100)
        xbond_dist_max: max dist between a and x (ref + 0.5) (default=4.0)
        xbond_don_angle_min: minimal angle for c-x-a (default=140)
        xbond_acc_angle_min: minimal angle for each x-a-y, where y are a's neighbours (custom) (default=90)
        multipolar_dist_max: max dist between a and x (default=4.0)
        multipolar_don_angle_min: minimal angle for c-x-a (default=90)
        multipolar_norm_angle_max: max angle for c-x--camide (default=40)
        pistacking_dist_max_t: max dist for T-shaped pi-stacking (mcgaughey, 1998) (default=5.5)
        pistacking_dist_max_f: between pistacking_dist_max_t and pistacking_dist_max_p (default=4.75)
        pistacking_dist_max_p: max dist for parallel pi-stacking (default=4.0)
        pistacking_ang_dev: max deviation to optimal angle (0 for //, 90 for |--, between -> face/edge-to-face) (default=30)
        pistacking_offset_max: max offset of the two r (corresponds to the radius of benzene + 1 a) (default=2.5)
        piother_dist_max: ref 1, table s2 (default=4.0)
        piother_offset_max: max offset for pi-hydrophobic, pi-cation and pi-amide interactions (custom) (default=2.0)
        saltbridge_dist_max: max dist between centers of charge for salt bridges (barlow and thornton, 1983) + 1.5 (default=5.5)
        water_bridge_mindist: min dist between o_wat and polar atom (jiang et al., 2005) - 0.1 (default=2.5)
        water_bridge_maxdist: max dist between o_wat and polar atom (jiang et al., 2005) + 0.5 (default=4.1)
        water_bridge_omega_min: min angle between a, o_wat and d hydrogen (jiang et al., 2005) - 9 (default=71)
        water_bridge_omega_max: max angle between a, o_wat and d hydrogen (jiang et al., 2005) (default=140)
        metal_dist_max: max dist between m and interacting atom (harding, 2001) (default=3.0)

    Note:
        Some distance thresholds were extended (max. 1.0 A) if too restrictive to account for low-quality structures.

    References:
        C. Bissantz, B. Kuhn, M. Stahl A Medicinal Chemist’s Guide to Molecular Interactions, 2010,
            J. Med. Chem., 53, 5061–5084
        R. Ferreira de Freitas, M.Schapira A systematic analysis of atomic protein–ligand interactions
            in the PDB, 2017, Med. Chem. Commun., 8, 1970-1981
        E. Nittinger, T. Inhester, S. Bietz, A. Meyder, K. T. Schomburg, G. Lange, R. Klein, M. Rarey
            Large-Scale Analysis of Hydrogen Bond Interaction Patterns in Protein−Ligand Interfaces,
            2017, J. Med. Chem., 60, 4245−4257
        N. K. Shinada, A. G. de Brevern, P. Schmidtke Halogens in Protein–Ligand Binding Mechanism:
            A Structural Perspective, 2019, J. Med. Chem., 62, 21, 9341–9356
        B. Kuhn, E. Gilberg, R. Taylorn J. Cole, O. Korb How Significant Are Unusual Protein−Ligand
            Interactions? Insights from Database Mining, 2019, J. Med. Chem, 62, 22, 10441–10455
    """

    bs_dist: float = 7.5
    min_dist: float = 2.0

    # hydrophobic interactions
    hydrophobic_dist_max: float = 4.0

    # h-bonds
    hbond_dist_max: float = 4.0
    hbond_don_angle_min: float = 140
    hbond_acc_angle_min: float = 90

    # halogen bonds (cl/br/i with e.g. [o]=c)
    xbond_dist_max: float = 4.0
    xbond_don_angle_min: float = 140
    xbond_acc_angle_min: float = 90

    # orthogonal multipolar interactions (f/cl with e.g. [c]=o)
    # see ppt anna vulpetti, cambridge, 2013, slide 11
    multipolar_dist_max: float = 4.0
    multipolar_don_angle_min: float = 90
    multipolar_norm_angle_max: float = 40

    # pi interactions (aromatic with aromatic/cation/hydrophobic/amide)
    pistacking_dist_max_t: float = 5.5
    pistacking_dist_max_f: float = 4.75
    pistacking_dist_max_p: float = 4.0
    pistacking_ang_dev: float = 30
    pistacking_offset_max: float = 2.5
    piother_dist_max: float = 4.0
    piother_offset_max: float = 2.0

    # other interactions
    saltbridge_dist_max: float = 5.5
    water_bridge_mindist: float = 2.5
    water_bridge_maxdist: float = 4.1
    water_bridge_omega_min: float = 71
    water_bridge_omega_max: float = 140
    metal_dist_max: float = 3.0
