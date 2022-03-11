from iktos.structure_analysis.scorers import ContactScorer

from iktos.structure_interactions.InteractionProfiler import InteractionProfiler


def test_analyse_complex():
    protein_path = 'tests/data/prot_3S3M.pdb'
    ligand_path = 'tests/data/lig_3S3M.sdf'

    plip = InteractionProfiler(plip_config_version='default')  # set plip_config here
    contacts_ref = [
        plip.analyse_complex(
            rec_coords=protein_path,
            lig_coords=ligand_path,
            as_string=False,
            refine=True,
            lig_format='sdf',
        )
    ]

    scorer = ContactScorer(contacts_ref=contacts_ref[:1])
    scores = scorer.score_many(contacts_ref)
    print('Contact scores =', scores)
    assert scores == [1.0]
