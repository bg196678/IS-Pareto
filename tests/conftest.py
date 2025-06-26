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
    return test_data_system_1_path / "kinetics" / "kinetics.xlsx"

@pytest.fixture(scope="session")
def test_data_system_1_gsolv(test_data_system_1_path):
    return test_data_system_1_path / "solvation" / "gsolv.xlsx"

@pytest.fixture(scope="session")
def test_data_system_1_experimental(test_data_system_1_path):
    return test_data_system_1_path / "experimental" / "experimental.xlsx"

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
        energy=-635.334856191927,
    )

    nucleophilic = Species(
        name="Nucleophilic",
        mass=0.071,
        fchk_file_path=test_data_system_1_gaussian / "Nucleophilic.fchk",
        tab_file_path=test_data_system_1_cosmo / "Nucleophilic.tab",
        energy=-212.567127673868,
    )

    its_1 = Species(
        name="ITS1",
        mass=0.159 + 0.071,
        fchk_file_path=test_data_system_1_gaussian / "ITS1.fchk",
        tab_file_path=test_data_system_1_cosmo / "ITS1.tab",
        energy=-847.894967749602,
    )

    its_2 = Species(
        name="ITS2",
        mass=0.159 + 0.071,
        fchk_file_path=test_data_system_1_gaussian / "ITS2.fchk",
        tab_file_path=test_data_system_1_cosmo / "ITS2.tab",
        energy=-847.87075447606,
    )

    its_3 = Species(
        name="ITS3",
        mass=0.159 + 0.071 + 0.071,
        fchk_file_path=test_data_system_1_gaussian / "ITS3.fchk",
        tab_file_path=test_data_system_1_cosmo / "ITS3.tab",
        energy=-960.014334953842,
    )

    its_4 = Species(
        name="ITS4",
        mass=0.159 + 0.071 + 0.071,
        fchk_file_path=test_data_system_1_gaussian / "ITS4.fchk",
        tab_file_path=test_data_system_1_cosmo / "ITS4.tab",
        energy=-959.989508763147,
    )

    product_1 = Species(
        name="Product1",
        mass=0.21008046,
        fchk_file_path=test_data_system_1_gaussian / "Product1.fchk",
        tab_file_path=test_data_system_1_cosmo / "Product1.tab",
        energy=-747.459040149113,
    )

    product_2 = Species(
        name="Product2",
        mass=0.21008046,
        fchk_file_path=test_data_system_1_gaussian / "Product2.fchk",
        tab_file_path=test_data_system_1_cosmo / "Product2.tab",
        energy=-747.459616153585,
    )

    product_3 = Species(
        name="Product3",
        mass=0.26114773,
        fchk_file_path=test_data_system_1_gaussian / "Product3.fchk",
        tab_file_path=test_data_system_1_cosmo / "Product3.tab",
        energy=-859.579422928392,
    )

    leaving_group = Species(
        name="LeavingGroup",
        mass=0.020,
        fchk_file_path=test_data_system_1_gaussian / "LeavingGroup.fchk",
        tab_file_path=test_data_system_1_cosmo / "LeavingGroup.tab",
        energy=-100.467619260434,
    )

    ##################
    # Transition State
    ##################
    ts_1_fwd = TransitionState(
        name="TS1_fwd",
        fchk_file_path=test_data_system_1_gaussian / "TS1.fchk",
        tab_file_path=test_data_system_1_cosmo / "TS1.tab",
        energy=-847.893645136721,
    )

    ts_1_rev = TransitionState(
        name="TS1_rev",
        fchk_file_path=test_data_system_1_gaussian / "TS1.fchk",
        tab_file_path=test_data_system_1_cosmo / "TS1.tab",
        energy=-847.893645136721,
    )

    ts_2_fwd = TransitionState(
        name="TS2_fwd",
        fchk_file_path=test_data_system_1_gaussian / "TS2.fchk",
        tab_file_path=test_data_system_1_cosmo / "TS2.tab",
        energy=-847.88361098324,
    )

    ts_2_rev = TransitionState(
        name="TS2_rev",
        fchk_file_path=test_data_system_1_gaussian / "TS2.fchk",
        tab_file_path=test_data_system_1_cosmo / "TS2.tab",
        energy=-847.88361098324,
    )

    ts_3_fwd = TransitionState(
        name="TS3_fwd",
        fchk_file_path=test_data_system_1_gaussian / "TS3.fchk",
        tab_file_path=test_data_system_1_cosmo / "TS3.tab",
        energy=-960.013678294337,
    )

    ts_3_rev = TransitionState(
        name="TS3_rev",
        fchk_file_path=test_data_system_1_gaussian / "TS3.fchk",
        tab_file_path=test_data_system_1_cosmo / "TS3.tab",
        energy=-960.013678294337,
    )

    ts_4_fwd = TransitionState(
        name="TS4_fwd",
        fchk_file_path=test_data_system_1_gaussian / "TS4.fchk",
        tab_file_path=test_data_system_1_cosmo / "TS4.tab",
        energy=-960.001694263016,
    )

    ts_4_rev = TransitionState(
        name="TS4_rev",
        fchk_file_path=test_data_system_1_gaussian / "TS4.fchk",
        tab_file_path=test_data_system_1_cosmo / "TS4.tab",
        energy=-960.001694263016,
    )

    ts_1_2_fwd = TransitionState(
        name="TS12_fwd",
        fchk_file_path=test_data_system_1_gaussian / "TS12.fchk",
        tab_file_path=test_data_system_1_cosmo / "TS12.tab",
        energy=-847.881008554995,
    )

    ts_1_2_rev = TransitionState(
        name="TS12_rev",
        fchk_file_path=test_data_system_1_gaussian / "TS12.fchk",
        tab_file_path=test_data_system_1_cosmo / "TS12.tab",
        energy=-847.881008554995,
    )

    ts_2_2_fwd = TransitionState(
        name="TS22_fwd",
        fchk_file_path=test_data_system_1_gaussian / "TS22.fchk",
        tab_file_path=test_data_system_1_cosmo / "TS22.tab",
        energy=-847.862395909096,
    )

    ts_2_2_rev = TransitionState(
        name="TS22_rev",
        fchk_file_path=test_data_system_1_gaussian / "TS22.fchk",
        tab_file_path=test_data_system_1_cosmo / "TS22.tab",
        energy=-847.862395909096,
    )

    ts_3_2_fwd = TransitionState(
        name="TS32_fwd",
        fchk_file_path=test_data_system_1_gaussian / "TS32.fchk",
        tab_file_path=test_data_system_1_cosmo / "TS32.tab",
        energy=-960.005129677018,
    )

    ts_3_2_rev = TransitionState(
        name="TS32_rev",
        fchk_file_path=test_data_system_1_gaussian / "TS32.fchk",
        tab_file_path=test_data_system_1_cosmo / "TS32.tab",
        energy=-960.005129677018,
    )

    ts_4_2_fwd = TransitionState(
        name="TS42_fwd",
        fchk_file_path=test_data_system_1_gaussian / "TS42.fchk",
        tab_file_path=test_data_system_1_cosmo / "TS42.tab",
        energy=-959.980523196275,
    )

    ts_4_2_rev = TransitionState(
        name="TS42_rev",
        fchk_file_path=test_data_system_1_gaussian / "TS42.fchk",
        tab_file_path=test_data_system_1_cosmo / "TS42.tab",
        energy=-959.980523196275,
    )

    ##############
    # Reactions
    ##############
    reaction_1_fwd = its_1 - substrate - nucleophilic
    reaction_1_fwd.name = "R1_fwd"
    reaction_1_fwd.transition_state = ts_1_fwd

    reaction_1_rev = substrate + nucleophilic - its_1
    reaction_1_rev.name = "R1_rev"
    reaction_1_rev.transition_state = ts_1_rev

    reaction_2_fwd = its_2 - substrate - nucleophilic
    reaction_2_fwd.name = "R2_fwd"
    reaction_2_fwd.transition_state = ts_2_fwd

    reaction_2_rev = substrate + nucleophilic - its_2
    reaction_2_rev.name = "R2_rev"
    reaction_2_rev.transition_state = ts_2_rev

    reaction_3_fwd = its_3 - product_2 - nucleophilic
    reaction_3_fwd.name = "R3_fwd"
    reaction_3_fwd.transition_state = ts_3_fwd

    reaction_3_rev = product_2 + nucleophilic - its_3
    reaction_3_rev.name = "R3_rev"
    reaction_3_rev.transition_state = ts_3_rev

    reaction_4_fwd = its_4 - product_1 - nucleophilic
    reaction_4_fwd.name = "R4_fwd"
    reaction_4_fwd.transition_state = ts_4_fwd

    reaction_4_rev = product_1 + nucleophilic - its_4
    reaction_4_rev.name = "R4_rev"
    reaction_4_rev.transition_state = ts_4_rev

    reaction_5_fwd = product_1 + leaving_group - its_1
    reaction_5_fwd.name = "R5_fwd"
    reaction_5_fwd.transition_state = ts_1_2_fwd

    reaction_5_rev = its_1 - product_1 - leaving_group
    reaction_5_rev.name = "R5_rev"
    reaction_5_rev.transition_state = ts_1_2_rev

    reaction_6_fwd = product_2 + leaving_group - its_2
    reaction_6_fwd.name = "R6_fwd"
    reaction_6_fwd.transition_state = ts_2_2_fwd

    reaction_6_rev = its_2 - product_2 - leaving_group
    reaction_6_rev.name = "R6_rev"
    reaction_6_rev.transition_state = ts_2_2_rev

    reaction_7_fwd = product_3 + leaving_group - its_3
    reaction_7_fwd.name = "R7_fwd"
    reaction_7_fwd.transition_state = ts_3_2_fwd

    reaction_7_rev = its_3 - product_3 - leaving_group
    reaction_7_rev.name = "R7_rev"
    reaction_7_rev.transition_state = ts_3_2_rev

    reaction_8_fwd = product_3 + leaving_group - its_4
    reaction_8_fwd.name = "R8_fwd"
    reaction_8_fwd.transition_state = ts_4_2_fwd

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
