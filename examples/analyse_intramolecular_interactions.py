from iktos.structure_interactions import (
    analyse_interactions_intra,
    convert_to_dict_intra,
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

# Example to analyse contacts in a peptide/protein:
protein_path = 'data/prot_5UIT.pdb'  # PDB file
contacts_raw = analyse_interactions_intra(
    protein_path,
    is_small_molecule=False,
    as_string=False,  # set to True if you want to give the content of the files instead
)
# The previous function returns a list of contact objects - to use it in SA,
# you need to convert them to a dict:
dict_contacts = convert_to_dict_intra(contacts_raw)
print(f'Dict of contacts: {dict_contacts}')
