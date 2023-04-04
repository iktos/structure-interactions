from typing import List, NamedTuple, Optional

from iktos.structure_interactions.detection import (
    H_Bond,
    Hydrophobic,
    Metal_Complex,
    Pi_Cation,
    Pi_Hydrophobic,
    Pi_Stacking,
    Water_Bridge,
)
from iktos.structure_interactions.refinement import (
    drop_duplicated_h_bonds,
    drop_duplicated_water_bridges,
    refine_h_bonds,
    refine_hydrophobics,
    refine_pi_cations,
    refine_pi_hydrophobics,
    refine_water_bridges,
)


# Create a dummy class Atom to mimic (in a simplified way) that used in the package
class Atom(NamedTuple):
    unique_id: str
    residue_id: Optional[str] = None


def test_refine_hbonds_1():
    # weak and strong H-bonds -> keep both
    hbond_1 = H_Bond(
        donor=[Atom(unique_id='ligand|UNL|A|1|1'), Atom(unique_id='ligand|UNL|A|1|2')],
        acceptor=[Atom(unique_id='receptor|SER|A|1|OG')],
        distance_ah=0.,  # not used
        distance_ad=3.5,
        angle_dha=150,
        type='strong',
    )
    hbond_2 = H_Bond(
        donor=[Atom(unique_id='ligand|UNL|A|1|5'), Atom(unique_id='ligand|UNL|A|1|6')],
        acceptor=[Atom(unique_id='receptor|SER|A|1|OG')],
        distance_ah=0.,  # not used
        distance_ad=3.5,
        angle_dha=150,
        type='weak',
    )
    selection = refine_h_bonds([hbond_1, hbond_2], [], [], [])
    assert len(selection) == 2


def test_refine_hbonds_2():
    # SER in 2 H-bonds with ligand, as an acceptor in both cases -> keep both
    hbond_1 = H_Bond(
        donor=[Atom(unique_id='ligand|UNL|A|1|1'), Atom(unique_id='ligand|UNL|A|1|2')],
        acceptor=[Atom(unique_id='receptor|SER|A|1|OG')],
        distance_ah=0.,  # not used
        distance_ad=3.5,
        angle_dha=150,
        type='strong',
    )
    hbond_2 = H_Bond(
        donor=[Atom(unique_id='ligand|UNL|A|1|11'), Atom(unique_id='ligand|UNL|A|1|12')],
        acceptor=[Atom(unique_id='receptor|SER|A|1|OG')],
        distance_ah=0.,  # not used
        distance_ad=3.0,
        angle_dha=150,
        type='strong',
    )
    selection = refine_h_bonds([hbond_1, hbond_2], [], [], [])
    assert len(selection) == 2


def test_refine_hbonds_3():
    # SER in 2 H-bonds with ligand, as a donor in both cases -> keep the 'best' one
    hbond_1 = H_Bond(
        acceptor=[Atom(unique_id='ligand|UNL|A|1|1')],
        donor=[Atom(unique_id='receptor|SER|A|1|OG'), Atom(unique_id='receptor|SER|A|1|HG')],
        distance_ah=0.,  # not used
        distance_ad=3.5,
        angle_dha=150,
        type='strong',
    )
    hbond_2 = H_Bond(
        acceptor=[Atom(unique_id='ligand|UNL|A|1|11')],
        donor=[Atom(unique_id='receptor|SER|A|1|OG'), Atom(unique_id='receptor|SER|A|1|HG')],
        distance_ah=0.,  # not used in refinement
        distance_ad=3.0,
        angle_dha=150,
        type='strong',
    )
    selection = refine_h_bonds([hbond_1, hbond_2], [], [], [])
    assert len(selection) == 1
    assert selection == [hbond_2]


def test_refine_hbonds_4():
    # ligand in 2 H-bonds with receptor, as an acceptor in both cases -> keep both
    hbond_1 = H_Bond(
        acceptor=[Atom(unique_id='ligand|UNL|A|1|1')],
        donor=[Atom(unique_id='receptor|SER|A|1|OG'), Atom(unique_id='receptor|SER|A|1|HG')],
        distance_ah=0.,  # not used
        distance_ad=3.5,
        angle_dha=150,
        type='strong',
    )
    hbond_2 = H_Bond(
        acceptor=[Atom(unique_id='ligand|UNL|A|1|1')],
        donor=[Atom(unique_id='receptor|TYR|A|5|N'), Atom(unique_id='receptor|TYR|A|5|NH')],
        distance_ah=0.,  # not used
        distance_ad=3.4,
        angle_dha=160,
        type='strong',
    )
    selection = refine_h_bonds([hbond_1, hbond_2], [], [], [])
    assert len(selection) == 2


def test_refine_hbonds_5():
    # ligand in 2 H-bonds with receptor, as a donor in both cases -> keep the 'best' one
    hbond_1 = H_Bond(
        donor=[Atom(unique_id='ligand|UNL|A|1|1'), Atom(unique_id='ligand|UNL|A|1|2')],
        acceptor=[Atom(unique_id='receptor|SER|A|1|OG')],
        distance_ah=0.,  # not used
        distance_ad=3.5,
        angle_dha=150,
        type='strong',
    )
    hbond_2 = H_Bond(
        donor=[Atom(unique_id='ligand|UNL|A|1|1'), Atom(unique_id='ligand|UNL|A|1|2')],
        acceptor=[Atom(unique_id='receptor|TYR|A|5|N')],
        distance_ah=0.,  # not used
        distance_ad=3.4,
        angle_dha=160,
        type='strong',
    )
    selection = refine_h_bonds([hbond_1, hbond_2], [], [], [])
    assert len(selection) == 1
    assert selection == [hbond_2]


def test_refine_hydrophobics_1():
    # 2 hydrophobic interactions -> keep both
    hydroph_1 = Hydrophobic(
        partner_1=[Atom(unique_id='ligand|UNL|A|1|1')],
        partner_2=[Atom(unique_id='receptor|LEU|A|1|CD1')],
        distance=3.8,
    )
    hydroph_2 = Hydrophobic(
        partner_1=[Atom(unique_id='ligand|UNL|A|1|6')],
        partner_2=[Atom(unique_id='receptor|VAL|A|5|CG1')],
        distance=3.4,
    )
    selection = refine_hydrophobics([hydroph_1, hydroph_2], [], [])
    assert len(selection) == 2


def test_refine_hydrophobics_2():
    # 2 hydrophobic interactions with the same ligand atom -> keep the 'best' one
    hydroph_1 = Hydrophobic(
        partner_1=[Atom(unique_id='ligand|UNL|A|1|1')],
        partner_2=[Atom(unique_id='receptor|LEU|A|1|CD1')],
        distance=3.8,
    )
    hydroph_2 = Hydrophobic(
        partner_1=[Atom(unique_id='ligand|UNL|A|1|1')],
        partner_2=[Atom(unique_id='receptor|VAL|A|5|CG1')],
        distance=3.4,
    )
    selection = refine_hydrophobics([hydroph_1, hydroph_2], [], [])
    assert len(selection) == 1
    assert selection == [hydroph_2]


def test_refine_hydrophobics_3():
    # 2 hydrophobic interactions with the same receptor atom -> keep the 'best' one
    hydroph_1 = Hydrophobic(
        partner_1=[Atom(unique_id='ligand|UNL|A|1|1')],
        partner_2=[Atom(unique_id='receptor|LEU|A|1|CD1')],
        distance=3.8,
    )
    hydroph_2 = Hydrophobic(
        partner_1=[Atom(unique_id='ligand|UNL|A|1|6')],
        partner_2=[Atom(unique_id='receptor|LEU|A|1|CD1')],
        distance=3.4,
    )
    selection = refine_hydrophobics([hydroph_1, hydroph_2], [], [])
    assert len(selection) == 1
    assert selection == [hydroph_2]


def test_refine_hydrophobics_4():
    # hydrophobic interaction + pi-hydrophobic interaction -> drop
    hydroph_1 = Hydrophobic(
        partner_1=[Atom(unique_id='ligand|UNL|A|1|1')],
        partner_2=[Atom(unique_id='receptor|LEU|A|1|CD1')],
        distance=3.8,
    )
    hydroph_2 = Hydrophobic(
        partner_1=[Atom(unique_id='ligand|UNL|A|1|6')],
        partner_2=[Atom(unique_id='receptor|LEU|A|1|CD1')],
        distance=3.4,
    )
    pi_hydroph = Pi_Hydrophobic(
        ring=[Atom(unique_id='ligand|UNL|A|1|1')],
        hydrophobic=[Atom(unique_id='receptor|LEU|A|1|CD1')],
        distance=3.5,
        offset=0.,  # not used
    )
    selection = refine_hydrophobics([hydroph_1, hydroph_2], [], [pi_hydroph])
    assert len(selection) == 1
    assert selection == [hydroph_2]


def test_refine_hydrophobics_5():
    # hydrophobic interaction + pi-hydrophobic interaction -> drop
    hydroph = Hydrophobic(
        partner_1=[Atom(unique_id='ligand|UNL|A|1|1')],
        partner_2=[Atom(unique_id='receptor|PHE|A|1|CD1')],
        distance=3.8,
    )
    pi_stacking = Pi_Stacking(
        partner_1=[Atom(unique_id='ligand|UNL|A|1|1')],
        partner_2=[Atom(unique_id='receptor|PHE|A|1|CD1'), Atom(unique_id='receptor|PHE|A|1|CE1')],
        distance=3.5,
        offset=0.,  # not used
        angle=0.,  # not used
        type='',  # not used
    )
    selection = refine_hydrophobics([hydroph], [pi_stacking], [])
    assert len(selection) == 0


def test_refine_pi_cations_1():
    # 2 pi-cations -> keep both
    pi_cation_1 = Pi_Cation(
        cation=[Atom(unique_id='ligand|UNL|A|1|1')],
        ring=[Atom(unique_id='receptor|HIS|A|9|CE1')],
        distance=4.0,
        offset=0.,  # not used
    )
    pi_cation_2 = Pi_Cation(
        cation=[Atom(unique_id='ligand|UNL|A|1|1')],
        ring=[Atom(unique_id='receptor|PHE|A|9|CE1')],
        distance=4.0,
        offset=0.,  # not used
    )
    selection = refine_pi_cations([pi_cation_1, pi_cation_2], [])
    assert len(selection) == 2


def test_refine_pi_cations_2():
    # pi-cation + pi-stacking -> drop
    pi_cation = Pi_Cation(
        ring=[Atom(unique_id='ligand|UNL|A|1|1')],
        cation=[Atom(unique_id='receptor|HIS|A|9|CE1')],
        distance=4.0,
        offset=0.,  # not used
    )
    pi_stacking = Pi_Stacking(
        partner_1=[Atom(unique_id='ligand|UNL|A|1|1')],
        partner_2=[Atom(unique_id='receptor|HIS|A|9|CE1'), Atom(unique_id='receptor|HIS|A|9|ND1')],
        distance=4.0,
        offset=0.,  # not used
        angle=0.,  # not used
        type='',  # not used
    )
    selection = refine_pi_cations([pi_cation], [pi_stacking])
    assert len(selection) == 0


def test_refine_pi_hydrophobics_1():
    # 2 pi-hydrophobics -> keep both
    pi_hydroph_1 = Pi_Hydrophobic(
        ring=[Atom(unique_id='ligand|UNL|A|1|1', residue_id='UNL|A|1')],
        hydrophobic=[Atom(unique_id='receptor|LEU|A|1|CD1', residue_id='LEU|A|1')],
        distance=3.5,
        offset=10.,  # used here!
    )
    pi_hydroph_2 = Pi_Hydrophobic(
        ring=[Atom(unique_id='ligand|UNL|A|1|5', residue_id='UNL|A|1')],
        hydrophobic=[Atom(unique_id='receptor|LEU|A|3|CG', residue_id='LEU|A|3')],
        distance=4.0,
        offset=30.,  # used here!
    )
    selection = refine_pi_hydrophobics([pi_hydroph_1, pi_hydroph_2])
    assert len(selection) == 2


def test_refine_pi_hydrophobics_2():
    # 2 pi-hydrophobics  between the same ligand atom(s) and 2 different receptor atoms
    # of the same residue -> keep the 'best' one
    pi_hydroph_1 = Pi_Hydrophobic(
        ring=[Atom(unique_id='ligand|UNL|A|1|1', residue_id='UNL|A|1')],
        hydrophobic=[Atom(unique_id='receptor|LEU|A|1|CD1', residue_id='LEU|A|1')],
        distance=3.5,
        offset=10.,  # used here!
    )
    pi_hydroph_2 = Pi_Hydrophobic(
        ring=[Atom(unique_id='ligand|UNL|A|1|1', residue_id='UNL|A|1')],
        hydrophobic=[Atom(unique_id='receptor|LEU|A|1|CG', residue_id='LEU|A|1')],
        distance=4.0,
        offset=30.,  # used here!
    )
    selection = refine_pi_hydrophobics([pi_hydroph_1, pi_hydroph_2])
    assert len(selection) == 1
    assert selection == [pi_hydroph_1]


def test_refine_pi_hydrophobics_3():
    # 2 pi-hydrophobics between the same ligand atom(s) and 2 different receptor residues
    # -> keep both
    pi_hydroph_1 = Pi_Hydrophobic(
        ring=[Atom(unique_id='ligand|UNL|A|1|1', residue_id='UNL|A|1')],
        hydrophobic=[Atom(unique_id='receptor|LEU|A|1|CD1', residue_id='LEU|A|1')],
        distance=3.5,
        offset=10.,  # used here!
    )
    pi_hydroph_2 = Pi_Hydrophobic(
        ring=[Atom(unique_id='ligand|UNL|A|1|1', residue_id='UNL|A|1')],
        hydrophobic=[Atom(unique_id='receptor|LEU|A|5|CG', residue_id='LEU|A|5')],
        distance=4.0,
        offset=30.,  # used here!
    )
    selection = refine_pi_hydrophobics([pi_hydroph_1, pi_hydroph_2])
    assert len(selection) == 2


def test_refine_pi_hydrophobics_4():
    # 2 pi-hydrophobics between different ligand atoms and the same receptor atoms
    # -> keep the 'best' one
    pi_hydroph_1 = Pi_Hydrophobic(
        ring=[Atom(unique_id='ligand|UNL|A|1|1', residue_id='UNL|A|1')],
        hydrophobic=[Atom(unique_id='receptor|LEU|A|1|CD1', residue_id='LEU|A|1')],
        distance=3.5,
        offset=10.,  # used here!
    )
    pi_hydroph_2 = Pi_Hydrophobic(
        ring=[Atom(unique_id='ligand|UNL|A|1|5', residue_id='UNL|A|1')],
        hydrophobic=[Atom(unique_id='receptor|LEU|A|1|CD1', residue_id='LEU|A|1')],
        distance=4.0,
        offset=30.,  # used here!
    )
    selection = refine_pi_hydrophobics([pi_hydroph_1, pi_hydroph_2])
    assert len(selection) == 1
    assert selection == [pi_hydroph_1]


def test_refine_water_bridges_1():
    # water bridge + metal complex -> drop water bridge
    water_bridge = Water_Bridge(
        acceptor=[Atom(unique_id='ligand|UNL|A|1|1')],
        donor=[Atom(unique_id='receptor|LEU|A|1|O')],
        water=[Atom(unique_id='receptor|SOL|A|7|OW')],
        distance_aw=0.,  # not used
        distance_dw=0.,  # not used
        angle_dhw=0.,  # not used
        angle_awh=0.,  # not used
    )
    metal_complex = Metal_Complex(
        ligand=[Atom(unique_id='ligand|UNL|A|1|1')],
        receptor=[Atom(unique_id='receptor|SOL|A|7|OW')],
        metal=[],  # not used
        num_partners=0,  # not used
        complex_num=0,  # not used
    )
    selection = refine_water_bridges([water_bridge], [metal_complex])
    assert len(selection) == 0


def test_drop_duplicated_h_bonds():
    # duplicated H-bonds (differ by H only), keep the 2nd one (shorter A---H distance)
    hbond_1 = H_Bond(
        donor=[Atom(unique_id='ligand|UNL|A|1|1'), Atom(unique_id='ligand|UNL|A|1|2')],
        acceptor=[Atom(unique_id='receptor|SER|A|1|OG')],
        distance_ah=3.5,
        distance_ad=0.,  # not used
        angle_dha=150,
        type='strong',
    )
    hbond_2 = H_Bond(
        donor=[Atom(unique_id='ligand|UNL|A|1|1'), Atom(unique_id='ligand|UNL|A|1|3')],  # change H
        acceptor=[Atom(unique_id='receptor|SER|A|1|OG')],
        distance_ah=3.0,
        distance_ad=0.,  # not used
        angle_dha=160,
        type='strong',
    )
    hbonds_filtered = drop_duplicated_h_bonds([hbond_1, hbond_2])
    assert len(hbonds_filtered) == 1
    assert hbonds_filtered[0] == hbond_2


def test_drop_duplicated_water_bridges():
    # duplicated water bridges (differ by H only), keep the 2nd one (larger D-H---W angle)
    water_bridge_1 = Water_Bridge(
        acceptor=[Atom(unique_id='ligand|UNL|A|1|1')],
        donor=[Atom(unique_id='receptor|LYS|A|1|NZ'), Atom(unique_id='receptor|LYS|A|1|HZ1')],
        water=[Atom(unique_id='receptor|SOL|A|7|OW')],
        distance_aw=0.,  # not used
        distance_dw=0.,  # not used
        angle_dhw=150.,
        angle_awh=0.,  # not used
    )
    water_bridge_2 = Water_Bridge(
        acceptor=[Atom(unique_id='ligand|UNL|A|1|1')],
        donor=[Atom(unique_id='receptor|LYS|A|1|NZ'), Atom(unique_id='receptor|LYS|A|1|HZ2')],
        water=[Atom(unique_id='receptor|SOL|A|7|OW')],
        distance_aw=0.,  # not used
        distance_dw=0.,  # not used
        angle_dhw=160.,
        angle_awh=0.,  # not used
    )
    water_bridges_filtered = drop_duplicated_water_bridges([water_bridge_1, water_bridge_2])
    assert len(water_bridges_filtered) == 1
    assert water_bridges_filtered[0] == water_bridge_2
