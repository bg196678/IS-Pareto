import pytest
from pathlib import Path

from insiliciopt.species import (
    Species,
    TransitionState,
)


@pytest.fixture(scope="session")
def test_data_path():
    return Path(__file__).parent / "test_data"

@pytest.fixture
def test_fchk_path(tmp_path):
    fchk_path = tmp_path / "test.fchk"
    fchk_path.touch()
    return fchk_path

@pytest.fixture
def test_tab_path(tmp_path):
    tab_path = tmp_path / "test.tab"
    tab_path.touch()
    return tab_path

@pytest.fixture(scope="session")
def test_data_system_1_path(test_data_path):
    return test_data_path / "system_1"

@pytest.fixture(scope="session")
def test_data_system_1_gaussian(test_data_system_1_path):
    return test_data_system_1_path / "gaussian"

@pytest.fixture(scope="session")
def test_data_system_1_cosmo(test_data_system_1_path):
    return test_data_system_1_path / "cosmotherm"

@pytest.fixture(scope="session")
def test_data_system_1_kinetics(test_data_system_1_path):
    return test_data_system_1_path / "kinetics" / "combined_kinetics.xlsx"

@pytest.fixture(scope="session")
def test_data_system_1_gsolv(test_data_system_1_path):
    return test_data_system_1_path / "solvation" / "gsolv.xlsx"

@pytest.fixture
def construct_system_1(test_data_system_1_gaussian, test_data_system_1_cosmo):
    ##############
    # Species
    ##############
    substrate = Species(
        name="Substrate",
        mass=0.159,
        fchk_file_path=test_data_system_1_gaussian / "Substrate.fchk",
        tab_file_path=test_data_system_1_cosmo / "Substrate.tab",
        energy=-635.185999076,
    )

    nucleophilic = Species(
        name="Nucleophilic",
        mass=0.071,
        fchk_file_path=test_data_system_1_gaussian / "Nucleophilic.fchk",
        tab_file_path=test_data_system_1_cosmo / "Nucleophilic.tab",
        energy=-212.576801654,
    )

    its_1 = Species(
        name="ITS1",
        mass=0.159 + 0.071,
        fchk_file_path=test_data_system_1_gaussian / "ITS1.fchk",
        tab_file_path=test_data_system_1_cosmo / "ITS1.tab",
        energy=-847.760037351,
    )

    its_2 = Species(
        name="ITS2",
        mass=0.159 + 0.071,
        fchk_file_path=test_data_system_1_gaussian / "ITS2.fchk",
        tab_file_path=test_data_system_1_cosmo / "ITS2.tab",
        energy=-847.734323929,
    )

    its_3 = Species(
        name="ITS3",
        mass=0.159 + 0.071 + 0.071,
        fchk_file_path=test_data_system_1_gaussian / "ITS3.fchk",
        tab_file_path=test_data_system_1_cosmo / "ITS3.tab",
        energy=-959.9276068,
    )

    its_4 = Species(
        name="ITS4",
        mass=0.159 + 0.071 + 0.071,
        fchk_file_path=test_data_system_1_gaussian / "ITS4.fchk",
        tab_file_path=test_data_system_1_cosmo / "ITS4.tab",
        energy=-959.881772303,
    )

    product_1 = Species(
        name="Product1",
        mass=0.21008046,
        fchk_file_path=test_data_system_1_gaussian / "Product1.fchk",
        tab_file_path=test_data_system_1_cosmo / "Product1.tab",
        energy=-747.339396317,
    )

    product_2 = Species(
        name="Product2",
        mass=0.21008046,
        fchk_file_path=test_data_system_1_gaussian / "Product2.fchk",
        tab_file_path=test_data_system_1_cosmo / "Product2.tab",
        energy=-747.339639927,
    )

    product_3 = Species(
        name="Product3",
        mass=0.26114773,
        fchk_file_path=test_data_system_1_gaussian / "Product3.fchk",
        tab_file_path=test_data_system_1_cosmo / "Product3.tab",
        energy=-859.487278239,
    )

    leaving_group = Species(
        name="LeavingGroup",
        mass=0.020,
        fchk_file_path=test_data_system_1_gaussian / "Leaving_Group.fchk",
        tab_file_path=test_data_system_1_cosmo / "Leaving_Group.tab",
        energy=-100.453158015,
    )

    ##################
    # Transition State
    ##################
    ts_1_fwd = TransitionState(
        name="TS1_fwd",
        fchk_file_path=test_data_system_1_gaussian / "TS1_fwd.fchk",
        tab_file_path=test_data_system_1_cosmo / "TS1.tab",
        energy=-847.759104357,
    )

    ts_1_rev = TransitionState(
        name="TS1_rev",
        fchk_file_path=test_data_system_1_gaussian / "TS1_rev.fchk",
        tab_file_path=test_data_system_1_cosmo / "TS1.tab",
        energy=-847.759104357,
    )

    ts_2_fwd = TransitionState(
        name="TS2_fwd",
        fchk_file_path=test_data_system_1_gaussian / "TS2_fwd.fchk",
        tab_file_path=test_data_system_1_cosmo / "TS2.tab",
        energy=-847.748732583,
    )

    ts_2_rev = TransitionState(
        name="TS2_rev",
        fchk_file_path=test_data_system_1_gaussian / "TS2_rev.fchk",
        tab_file_path=test_data_system_1_cosmo / "TS2.tab",
        energy=-847.748732583,
    )

    ts_3_fwd = TransitionState(
        name="TS3_fwd",
        fchk_file_path=test_data_system_1_gaussian / "TS3_fwd.fchk",
        tab_file_path=test_data_system_1_cosmo / "TS3.tab",
        energy=-959.907069926,
    )

    ts_3_rev = TransitionState(
        name="TS3_rev",
        fchk_file_path=test_data_system_1_gaussian / "TS3_rev.fchk",
        tab_file_path=test_data_system_1_cosmo / "TS3.tab",
        energy=-959.907069926,
    )

    ts_4_fwd = TransitionState(
        name="TS4_fwd",
        fchk_file_path=test_data_system_1_gaussian / "TS4_fwd.fchk",
        tab_file_path=test_data_system_1_cosmo / "TS4.tab",
        energy=-959.894825835,
    )

    ts_4_rev = TransitionState(
        name="TS4_rev",
        fchk_file_path=test_data_system_1_gaussian / "TS4_rev.fchk",
        tab_file_path=test_data_system_1_cosmo / "TS4.tab",
        energy=-959.894825835,
    )

    ts_1_2_fwd = TransitionState(
        name="TS12_fwd",
        fchk_file_path=test_data_system_1_gaussian / "TS12_fwd.fchk",
        tab_file_path=test_data_system_1_cosmo / "TS12.tab",
        energy=-847.749041202,
    )

    ts_1_2_rev = TransitionState(
        name="TS12_rev",
        fchk_file_path=test_data_system_1_gaussian / "TS12_rev.fchk",
        tab_file_path=test_data_system_1_cosmo / "TS12.tab",
        energy=-847.749041202,
    )

    ts_2_2_fwd = TransitionState(
        name="TS22_fwd",
        fchk_file_path=test_data_system_1_gaussian / "TS22_fwd.fchk",
        tab_file_path=test_data_system_1_cosmo / "TS22.tab",
        energy=-847.729373816,
    )

    ts_2_2_rev = TransitionState(
        name="TS22_rev",
        fchk_file_path=test_data_system_1_gaussian / "TS22_rev.fchk",
        tab_file_path=test_data_system_1_cosmo / "TS22.tab",
        energy=-847.729373816,
    )

    ts_3_2_fwd = TransitionState(
        name="TS32_fwd",
        fchk_file_path=test_data_system_1_gaussian / "TS32_fwd.fchk",
        tab_file_path=test_data_system_1_cosmo / "TS32.tab",
        energy=-959.900414613,
    )

    ts_3_2_rev = TransitionState(
        name="TS32_rev",
        fchk_file_path=test_data_system_1_gaussian / "TS32_rev.fchk",
        tab_file_path=test_data_system_1_cosmo / "TS32.tab",
        energy=-959.900414613,
    )

    ts_4_2_fwd = TransitionState(
        name="TS42_fwd",
        fchk_file_path=test_data_system_1_gaussian / "TS42_fwd.fchk",
        tab_file_path=test_data_system_1_cosmo / "TS42.tab",
        energy=-959.874575137,
    )

    ts_4_2_rev = TransitionState(
        name="TS42_rev",
        fchk_file_path=test_data_system_1_gaussian / "TS42_rev.fchk",
        tab_file_path=test_data_system_1_cosmo / "TS42.tab",
        energy=-959.874575137,
    )

    ##############
    # Reactions
    ##############
    # R1_fwd: Substrate + Nucleophilic -> ITS1
    reaction_1_fwd = its_1 - substrate - nucleophilic
    reaction_1_fwd.name = "R1_fwd"
    reaction_1_fwd.transition_state = ts_1_fwd

    # R1_rev: ITS1 -> Substrate + Nucleophilic
    reaction_1_rev = substrate + nucleophilic - its_1
    reaction_1_rev.name = "R1_rev"
    reaction_1_rev.transition_state = ts_1_rev

    # R2_fwd: Substrate + Nucleophilic -> ITS2
    reaction_2_fwd = its_2 - substrate - nucleophilic
    reaction_2_fwd.name = "R2_fwd"
    reaction_2_fwd.transition_state = ts_2_fwd

    # R2_rev: ITS2 -> Substrate + Nucleophilic
    reaction_2_rev = substrate + nucleophilic - its_2
    reaction_2_rev.name = "R2_rev"
    reaction_2_rev.transition_state = ts_2_rev

    # R3_fwd: Product2 + Nucleophilic -> ITS3
    reaction_3_fwd = its_3 - product_2 - nucleophilic
    reaction_3_fwd.name = "R3_fwd"
    reaction_3_fwd.transition_state = ts_3_fwd

    # R3_rev: ITS3 -> Product2 + Nucleophilic
    reaction_3_rev = product_2 + nucleophilic - its_3
    reaction_3_rev.name = "R3_rev"
    reaction_3_rev.transition_state = ts_3_rev

    # R4_fwd: Product1 + Nucleophilic -> ITS4
    reaction_4_fwd = its_4 - product_1 - nucleophilic
    reaction_4_fwd.name = "R4_fwd"
    reaction_4_fwd.transition_state = ts_4_fwd

    # R4_rev: ITS4 -> Product1 + Nucleophilic
    reaction_4_rev = product_1 + nucleophilic - its_4
    reaction_4_rev.name = "R4_rev"
    reaction_4_rev.transition_state = ts_4_rev

    # R5_fwd: ITS1 -> Product1 + LeavingGroup
    reaction_5_fwd = product_1 + leaving_group - its_1
    reaction_5_fwd.name = "R5_fwd"
    reaction_5_fwd.transition_state = ts_1_2_fwd

    # R5_rev: Product1 + LeavingGroup -> ITS1
    reaction_5_rev = its_1 - product_1 - leaving_group
    reaction_5_rev.name = "R5_rev"
    reaction_5_rev.transition_state = ts_1_2_rev

    # R6_fwd: ITS2 -> Product2 + LeavingGroup
    reaction_6_fwd = product_2 + leaving_group - its_2
    reaction_6_fwd.name = "R6_fwd"
    reaction_6_fwd.transition_state = ts_2_2_fwd

    # R6_rev: Product2 + LeavingGroup -> ITS2
    reaction_6_rev = its_2 - product_2 - leaving_group
    reaction_6_rev.name = "R6_rev"
    reaction_6_rev.transition_state = ts_2_2_rev

    # R7_fwd: ITS3 -> Product3 + LeavingGroup
    reaction_7_fwd = product_3 + leaving_group - its_3
    reaction_7_fwd.name = "R7_fwd"
    reaction_7_fwd.transition_state = ts_3_2_fwd

    # R7_rev: Product3 + LeavingGroup -> ITS3
    reaction_7_rev = its_3 - product_3 - leaving_group
    reaction_7_rev.name = "R7_rev"
    reaction_7_rev.transition_state = ts_3_2_rev

    # R8_fwd: ITS4 -> Product3 + LeavingGroup
    reaction_8_fwd = product_3 + leaving_group - its_4
    reaction_8_fwd.name = "R8_fwd"
    reaction_8_fwd.transition_state = ts_4_2_fwd

    # R8_rev: Product3 + LeavingGroup -> ITS4
    reaction_8_rev = its_4 - product_3 - leaving_group
    reaction_8_rev.name = "R8_rev"
    reaction_8_rev.transition_state = ts_4_2_rev

    reactions = [
        reaction_1_fwd, reaction_1_rev,
        reaction_2_fwd, reaction_2_rev,
        reaction_3_fwd, reaction_3_rev,
        reaction_4_fwd, reaction_4_rev,
        reaction_5_fwd, reaction_5_rev,
        reaction_6_fwd, reaction_6_rev,
        reaction_7_fwd, reaction_7_rev,
        reaction_8_fwd, reaction_8_rev
    ]

    return reactions
