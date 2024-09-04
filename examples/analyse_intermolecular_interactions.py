from iktos.structure_interactions import (
    analyse_interactions_inter,
    analyse_interactions_inter_multi,
    convert_to_dict_inter,
)

"""
NOTES:
 - Hs need to be explicit and in 3D both for the protein and the ligand
 - ligand files are assumed to contain only 1 molecule (multi-SDF files are not supported)
 - if cofactors are given with the protein, the PDB should contain the CONECT lines and explicit charges
   (only NAD has been fully integrated and tested, for other cofactors you should read the logs carefully)
 - atom names in the protein need to follow conventions: CA = backbone carbone, H = hydrogen bonded to backbone nitrogen,
   other Hs are e.g. HE21, etc
 - be extremely careful with custom naming, this might lead to errors (e.g. a H named HO34 will be interpreted
   as Holmium by OpenBabel, leading to the detection of metal interactions)
"""

# Example to analyse contacts in a single complex:
protein_path = "data/prot_5UIT.pdb"  # PDB file
ligand_path = "data/lig_5UIT.sdf"  # SDF file is recommended
contacts_raw = analyse_interactions_inter(
    rec_coords=protein_path,
    lig_coords=ligand_path,
    as_string=False,  # set to True if you want to give the content of the files instead
    refine=True,  # whether to drop redundant contacts (e.g. hydrophobic + pi-stacking)
    lig_format="sdf",
)
# The previous function returns a list of contact objects - to use it in SA,
# you need to convert them to a dict:
dict_contacts = convert_to_dict_inter(contacts_raw)
print(f"Dict of contacts: {dict_contacts}")

# Example to analyse contacts in multiple complexes (different ligands with the same protein,
# e.g. to analyse docking poses, this function is faster than looping over the complexes):
protein_path = "data/prot_5UIT.pdb"  # PDB file
ligand_path = "data/lig_5UIT.sdf"  # SDF file is recommended
contacts_raw_multi = analyse_interactions_inter_multi(
    rec_coords=protein_path,
    lig_coords=[ligand_path] * 5,
    as_string=False,  # set to True if you want to give the content of the files instead
    refine=True,  # whether to drop redundant contacts (e.g. hydrophobic + pi-stacking)
    lig_format="sdf",
)
# The previous function returns a list of list of contact objects - to use it in SA,
# you need to convert them to a list of dicts:
dict_contacts_multi = [
    convert_to_dict_inter(contacts_raw) for contacts_raw in contacts_raw_multi
]
print(f"List of dicts of contacts: {dict_contacts_multi}")
