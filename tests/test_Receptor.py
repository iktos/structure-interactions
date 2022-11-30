from iktos.structure_interactions.Ligand import Ligand
from iktos.structure_interactions.Receptor import Receptor


def test_Receptor():
    # Load all the ligands
    with open('tests/data/poses_in_5UIT.sdf') as f:
        block = f.read()
    lig_coords = block.split('$$$$\n')[:-1]
    ligands = [Ligand(coords, 'sdf', as_string=True) for coords in lig_coords]

    # Load receptor and use the first ligand to define the binding site
    rec_coords = 'tests/data/prot_5UIT.pdb'
    receptor = Receptor(rec_coords, as_string=False)
    receptor.detect_binding_site(ligands[:1])
    receptor.identify_functional_groups()
    assert len(receptor.bs_atoms) == 502

    # Load receptor and use all the ligands to define the binding site
    receptor = Receptor(rec_coords, as_string=False)
    receptor.detect_binding_site(ligands)
    receptor.identify_functional_groups()
    assert len(receptor.bs_atoms) == 623
