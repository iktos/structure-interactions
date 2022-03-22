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

AROMATIC_PLANARITY = 5.0  # max allowed deviation from planarity in aromatic rings
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
