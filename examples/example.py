# Examples to detect contacts in complexes (PDB of protein + SDF of ligand),
# compute contact score and prepare Pymol sessions

from pathlib import Path

from iktos.structure_analysis.scorers import ContactScorer
from iktos.structure_analysis.visualization.pymol import (
    prepare_session,
    prepare_session_multistate,
)

from iktos.structure_interactions.InteractionProfiler import InteractionProfiler

# Note: PDB files need to be clean!
# Atom names should be UNIQUE (i.e. CA = backbone carbone, H = hydrogen bonded to backbone nitrogen, other Hs are e.g. HE21)
# Be very careful with custom naming, this might lead to errors
# (e.g. a H named HO34 will be interpreted as Holmium by OpenBabel, leading to the detection of metal interactions)

here = Path(__file__).parent
data_path = here / "data"

protein_pdb_list = [
    data_path / 'f_0_prot.pdb',
    data_path / 'f_146_prot.pdb',
    data_path / 'f_277_prot.pdb',
    data_path / 'f_332_prot.pdb',
]
ligand_sdf_list = [
    data_path / 'f_0_lig.sdf',
    data_path / 'f_146_lig.sdf',
    data_path / 'f_277_lig.sdf',
    data_path / 'f_332_lig.sdf',
]

# Analyse contacts in multiple complexes
plip = InteractionProfiler(plip_config_version='default')  # set plip_config here
contacts_ref = [
    plip.analyse_complex(
        str(protein_pdb_list[i]),
        str(ligand_sdf_list[i]),
        as_string=False,
        refine=True,
        lig_format='sdf',
    )
    for i in range(len(protein_pdb_list))
]

# Score contacts
# 1. Use the 1st complex to instantiate the scorer
print('** Example 1: score contacts using the 1st complex to initialise the score')
scorer = ContactScorer(contacts_ref[:1], mode='md')
scores = scorer.score_many(contacts_ref)
print('Contact scores =', scores)
scores = scorer.score_ensemble_docking(
    {f'{i}': [contacts] for i, contacts in enumerate(contacts_ref)}
)
print('Ensemble contact score =', scores)
print()

# 2. Use the reference complexes to instantiate the scorer
print(
    '** Example 2: score contacts using the 4 complexes to initialise the score, with mode `md`'
)
scorer = ContactScorer(contacts_ref, mode='md')
scores = scorer.score_many(contacts_ref)
print('Contact scores =', scores)
scores = scorer.score_ensemble_docking(
    {f'{i}': [contacts] for i, contacts in enumerate(contacts_ref)}
)
print('Ensemble contact score =', scores)
print()

# 3. Use custom weights
print(
    '** Example 3: score contacts using custom weights (here, we want to ignore some hydrophobic and weak H-bonds)'
)
scorer = ContactScorer(
    contacts_ref[:1],
    custom_weights={
        'Hydrophobic': {
            'MET|A|196|SD': 0.0,
            'PHE|A|198|CB': 0.0,
            'PHE|A|198|CE1': 0.0,
            'LYS|A|209|CG': 0.0,
            'HID|A|250|CB': 0.0,
            'TYR|A|303|CB': 0.0,
            'TYR|A|303|CD2': 0.0,
        },
        'H_Bond_Weak': {'ASP|A|84|OD2': 0.0, 'TYR|A|13|CD2+HD2': 0.0},
    },
)
scores = scorer.score_many(contacts_ref)
print('Contact scores =', scores)
scores = scorer.score_ensemble_docking(
    {f'{i}': [contacts] for i, contacts in enumerate(contacts_ref)}
)
print('Ensemble contact score =', scores)
print()


# Prepare Pymol session
# 1. Standard session with the 1st complex + contacts and weights
with open(protein_pdb_list[0]) as f:
    protein_pdb_block = f.read()
with open(ligand_sdf_list[0]) as f:
    ligand_sdf_block = f.read()
scorer = ContactScorer(contacts_ref[:1], mode='md')
prepare_session(
    protein_pdb_block=protein_pdb_block,
    ligand_sdf_block=ligand_sdf_block,
    contacts=contacts_ref[0],
    weights=scorer.weights,
    filename_session='edelris_f_0.pse',
    color_bg='white',
    color_protein='cbaw',  # show carbon atoms of the protein in white/grey
)

# 2. Multistate session with the 4 complexes loaded as different states + contacts + weights
protein_pdb_blocks = []
for file in protein_pdb_list:
    with open(file) as f:
        protein_pdb_blocks.append(f.read())
ligand_sdf_blocks = []
for file in ligand_sdf_list:
    with open(file) as f:
        ligand_sdf_blocks.append(f.read())
scorer = ContactScorer(contacts_ref, mode='md')
prepare_session_multistate(
    protein_pdb_blocks=protein_pdb_blocks,
    ligand_sdf_blocks=ligand_sdf_blocks,
    contacts=contacts_ref,
    weights=scorer.weights,
    filename_session='edelris_multistate.pse',
    color_bg='white',
    color_protein='cbaw',  # show carbon atoms of the protein in white/grey
)
