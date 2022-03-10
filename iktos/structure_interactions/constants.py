# fmt: off

# Thresholds for detection
# Note: some distance thresholds were extended (max. 1.0 A) if too restrictive
# to account for low-quality structures

# References:
# 1/ C. Bissantz, B. Kuhn, M. Stahl A Medicinal Chemist’s Guide to Molecular Interactions, 2010, J. Med. Chem., 53, 5061–5084
# 2/ R. Ferreira de Freitas, M.Schapira A systematic analysis of atomic protein–ligand interactions in the PDB, 2017, Med. Chem. Commun., 8, 1970-1981
# 3/ E. Nittinger, T. Inhester, S. Bietz, A. Meyder, K. T. Schomburg, G. Lange, R. Klein, M. Rarey Large-Scale Analysis of Hydrogen Bond Interaction Patterns in Protein−Ligand Interfaces, 2017, J. Med. Chem., 60, 4245−4257
# 4/ N. K. Shinada, A. G. de Brevern, P. Schmidtke Halogens in Protein–Ligand Binding Mechanism: A Structural Perspective, 2019, J. Med. Chem., 62, 21, 9341–9356
# 5/ B. Kuhn, E. Gilberg, R. Taylorn J. Cole, O. Korb How Significant Are Unusual Protein−Ligand Interactions? Insights from Database Mining, 2019, J. Med. Chem, 62, 22, 10441–10455

BS_DIST = 7.5  # max dist to include binding site residues
AROMATIC_PLANARITY = 5.0  # max allowed deviation from planarity in aromatic rings
MIN_DIST = 0.5  # min dist for all distance thresholds

# Hydrophobic interactions
HYDROPHOBIC_DIST_MAX = 4.0  # max dist for detection of hydrophobic contacts

# H-bonds
HBOND_DIST_MAX = 4.0  # max dist between D and A (Hubbard & Haider, 2001) + 0.5 A
HBOND_DON_ANGLE_MIN = 140  # minimal angle for D-H-A
HBOND_ACC_ANGLE_MIN = 100  # minimal angle for each H-A-Y, where Y are A's neighbours (custom)

# Halogen bonds (Cl/Br/I with e.g. [O]=C)
XBOND_DIST_MAX = 4.0  # max dist between A and X (ref + 0.5)
XBOND_DON_ANGLE_MIN = 140  # minimal angle for C-X-A
XBOND_ACC_ANGLE_MIN = 90  # minimal angle for each X-A-Y, where Y are A's neighbours (custom)

# Orthogonal multipolar interactions (F/Cl with e.g. [C]=O)
# See ppt Anna Vulpetti, Cambridge, 2013, slide 11
MULTIPOLAR_DIST_MAX = 4.0  # max dist between A and X
MULTIPOLAR_DON_ANGLE_MIN = 90  # minimal angle for C-X-A
MULTIPOLAR_NORM_ANGLE_MAX = 40  # max angle for C-X--Camide

# Pi interactions (aromatic with aromatic/cation/hydrophobic/amide)
PISTACKING_DIST_MAX_T = 5.5  # max dist for T-shaped pi-stacking (McGaughey, 1998)
PISTACKING_DIST_MAX_F = 4.75  # between PISTACKING_DIST_MAX_T and PISTACKING_DIST_MAX_P
PISTACKING_DIST_MAX_P = 4.0  # max dist for parallel pi-stacking
PISTACKING_ANG_DEV = 30  # max deviation to optimal angle (0 for //, 90 for |--, between -> face/edge-to-face)
PISTACKING_OFFSET_MAX = 2.5  # max offset of the two R (corresponds to the radius of benzene + 1 A)
PIOTHER_DIST_MAX = 4.0  # ref 1, table S2
PIOTHER_OFFSET_MAX = 2.0  # max offset for pi-hydrophobic, pi-cation and pi-amide interactions (custom)

# Other interactions
SALTBRIDGE_DIST_MAX = 5.5  # max dist between centers of charge for salt bridges (Barlow and Thornton, 1983) + 1.5
WATER_BRIDGE_MINDIST = 2.5  # min dist between O_wat and polar atom (Jiang et al., 2005) - 0.1
WATER_BRIDGE_MAXDIST = 4.1  # max dist between O_wat and polar atom (Jiang et al., 2005) + 0.5
WATER_BRIDGE_OMEGA_MIN = 71  # min angle between A, O_wat and D hydrogen (Jiang et al., 2005) - 9
WATER_BRIDGE_OMEGA_MAX = 140  # max angle between A, O_wat and D hydrogen (Jiang et al., 2005)
METAL_DIST_MAX = 3.0  # max dist between M and interacting atom (Harding, 2001)

# Useful references
# X-bonds: Halogen bonds in biological molecules, Auffinger

# Coordination numbers and geometries for metal complexes detection
METAL_COMPLEX_COO = {
    2: ['linear'],
    3: ['trigonal.planar', 'trigonal.pyramidal'],
    4: ['tetrahedral', 'square.planar'],
    5: ['trigonal.bipyramidal', 'square.pyramidal'],
    6: ['octahedral'],
}

# Angle signatures for each geometry (as seen from each target atom)
METAL_COMPLEX_ANG = {
    'linear': [[180.0]] * 2,
    'trigonal.planar': [[120.0, 120.0]] * 3,
    'trigonal.pyramidal': [[109.5, 109.5]] * 3,
    'tetrahedral': [[109.5, 109.5, 109.5, 109.5]] * 4,
    'square.planar': [[90.0, 90.0, 90.0, 90.0]] * 4,
    'trigonal.bipyramidal': [[120.0, 120.0, 90.0, 90.0]] * 3 + [[90.0, 90.0, 90.0, 180.0]] * 2,
    'square.pyramidal': [[90.0, 90.0, 90.0, 180.0]] * 4 + [[90.0, 90.0, 90.0, 90.0]],
    'octahedral': [[90.0, 90.0, 90.0, 90.0, 180.0]] * 6,
}

# Metal cations
METAL_IONS = [
    'CA',
    'CO',
    'MG',
    'MN',
    'FE',
    'FE1',
    'FE2',
    'FE3',
    'FE4',
    'LI',
    'NA',
    'K',
    'RB',
    'SR',
    'CS',
    'BA',
    'CR',
    'NI',
    'CU',
    'ZN',
    'RU',
    'RU1',
    'RH',
    'RH1',
    'PD',
    'AG',
    'CD',
    'LA',
    'W',
    'W1',
    'OS',
    'IR',
    'PT',
    'PT1',
    'AU',
    'HG',
    'CE',
    'PR',
    'SM',
    'EU',
    'GD',
    'TB',
    'YB',
    'LU',
    'AL',
    'GA',
    'IN',
    'SB',
    'TL',
    'PB',
]

# fmt: on
