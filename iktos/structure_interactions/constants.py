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

from typing import Any, Dict, List

AROMATIC_PLANARITY = 5.0  # max allowed deviation from planarity in aromatic rings

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

# Dict of charged AAs and cofactors, used by PLIP analysis
# Note: protonated HIS should be HIP, protonated ASP should be ASH, etc.
# but here we consider all of these AAs as potentially charged
CHARGED_RESIDUES: Dict[str, List[Dict[str, Any]]] = {
    'ARG': [{
        'fgroup': 'iminium',
        'charge': 'positive',
        'on_atoms': ['CZ', 'NE', 'NH1', 'NH2'],
    }],
    'ASP': [{
        'fgroup': 'carboxylate',
        'charge': 'negative',
        'on_atoms': ['CG', 'OD1', 'OD2'],
    }],
    'ASH': [{
        'fgroup': 'carboxylate',
        'charge': 'negative',
        'on_atoms': ['CG', 'OD1', 'OD2'],
    }],
    'GLU': [{
        'fgroup': 'carboxylate',
        'charge': 'negative',
        'on_atoms': ['CD', 'OE1', 'OE2'],
    }],
    'GLH': [{
        'fgroup': 'carboxylate',
        'charge': 'negative',
        'on_atoms': ['CD', 'OE1', 'OE2'],
    }],
    'LYS': [{'fgroup': 'ammonium', 'charge': 'positive', 'on_atoms': ['NZ']}],
    'LYN': [{'fgroup': 'ammonium', 'charge': 'positive', 'on_atoms': ['NZ']}],
    'HIS': [{
        'fgroup': 'iminium',
        'charge': 'positive',
        'on_atoms': ['ND1', 'NE2', 'CE1'],
    }],
    'HID': [{
        'fgroup': 'iminium',
        'charge': 'positive',
        'on_atoms': ['ND1', 'NE2', 'CE1'],
    }],
    'HIE': [{
        'fgroup': 'iminium',
        'charge': 'positive',
        'on_atoms': ['ND1', 'NE2', 'CE1'],
    }],
    'HIP': [{
        'fgroup': 'iminium',
        'charge': 'positive',
        'on_atoms': ['ND1', 'NE2', 'CE1'],
    }],
    'NAD': [{
        'fgroup': 'phosphate',
        'charge': 'negative',
        'on_atoms': ['PA', 'O1A', 'O2A'],
    }, {
        'fgroup': 'phosphate',
        'charge': 'negative',
        'on_atoms': ['PN', 'O1N', 'O2N'],
    }],
}

# fmt: on
