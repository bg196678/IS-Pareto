import pytest
import numpy as np
import pandas as pd
from pathlib import Path

from insiliciopt.kinetics import Kinetics
from insiliciopt.reactor import (
    Reactor,
    ReactorConditions,
)
from insiliciopt.species import Species
from insiliciopt.solvation import Solvation

class PatchedKinetics(Kinetics):

    def __init__(self, reactions: list, excel_path: Path):
        super().__init__(reactions)
        self.kinetics = self._load_kinetics(excel_path)

    @staticmethod
    def _load_kinetics(excel_path: Path) -> pd.DataFrame:
        kinetics = pd.read_excel(excel_path)
        return kinetics

    def k(self, reaction, temperature) -> float:
        transition_state = reaction.transition_state
        name = transition_state.name
        return float(np.interp(
            temperature,
            self.kinetics["Temperature"].values.tolist(),
            self.kinetics[name].values.tolist(),
        ))


class TestSystem1PatchedReactor:

    @pytest.fixture
    def setup_kinetics(self, construct_system_1, test_data_system_1_kinetics):
        excel_path = test_data_system_1_kinetics.parent / "kinetics.xlsx"
        kinetics = PatchedKinetics(
            construct_system_1,
            excel_path,
        )
        return kinetics

    @pytest.fixture
    def setup_solvation(self, construct_system_1):
        solvation = Solvation(construct_system_1)
        return solvation

    @pytest.fixture
    def setup_reactants(
            self,
            test_data_system_1_gaussian,
            test_data_system_1_cosmo,
    ):
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
        return substrate, nucleophilic

    @pytest.fixture
    def setup_product(
            self,
            test_data_system_1_gaussian,
            test_data_system_1_cosmo,
    ):
        product_1 = Species(
            name="Product1",
            mass=0.21008046,
            fchk_file_path=test_data_system_1_gaussian / "Product1.fchk",
            tab_file_path=test_data_system_1_cosmo / "Product1.tab",
            energy=-747.339396317,
        )
        return product_1

    @pytest.mark.slow
    def test_reactor(
            self,
            construct_system_1,
            setup_kinetics,
            setup_solvation,
            setup_reactants,
            setup_product,
            test_data_system_1_experimental,
    ):
        experimental_excel_data = pd.read_excel(
            test_data_system_1_experimental)

        reactor = Reactor(
            reactions=construct_system_1,
            kinetics=setup_kinetics,
            solvation=setup_solvation,
        )

        substrate, nucleophilic = setup_reactants

        for _, row in experimental_excel_data.iterrows():
            time = row["tres/min"]
            ratio = row["2:1"]
            concentration = row["Conc 1/M"] * 1000
            temperature = row["Temp/°C"]
            expected_E = row["E-factor"]
            exptected_STY = row["STY/kg m-3 h-1"]

            concentrations = {
                substrate: concentration,
                nucleophilic: concentration * ratio,
            }

            reactor_conditions = ReactorConditions(
                temperature=temperature,
                concentrations=concentrations,
                products=[setup_product, ],
                time=time,
            )

            STY, E = reactor.simulate(reactor_conditions)

            assert np.isclose(
                E, expected_E, rtol=0.15, atol=0.3
            ), (
                f"(Conditions: '{row}') "
                f"at {temperature}K: Calculated={E:.4e}, "
                f"Excel={expected_E:.4e}"
            )
            print("=" * 20)
            print("STY:", STY)
            print("Excel STY:", exptected_STY)
            print("E:", E)
            print("Excel E:", expected_E)
            print("=" * 20)

            assert np.isclose(
                STY, exptected_STY, rtol=0.10, atol=1250
            ), (
                f"(Conditions: '{row}') "
                f"at {temperature}K: Calculated={STY:.4e}, "
                f"Excel={exptected_STY:.4e}"
            )


class TestSystem1Reactor:

    @pytest.fixture
    def setup_kinetics(self, construct_system_1):
        kinetics = Kinetics(
            construct_system_1
        )
        return kinetics

    @pytest.fixture
    def setup_solvation(self, construct_system_1):
        solvation = Solvation(construct_system_1)
        return solvation

    @pytest.fixture
    def setup_reactants(
            self,
            test_data_system_1_gaussian,
            test_data_system_1_cosmo,
    ):
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
        return substrate, nucleophilic

    @pytest.fixture
    def setup_product(
            self,
            test_data_system_1_gaussian,
            test_data_system_1_cosmo,
    ):
        product_1 = Species(
            name="Product1",
            mass=0.21008046,
            fchk_file_path=test_data_system_1_gaussian / "Product1.fchk",
            tab_file_path=test_data_system_1_cosmo / "Product1.tab",
            energy=-747.339396317,
        )
        return product_1

    @pytest.mark.slow
    def test_reactor(
            self,
            construct_system_1,
            setup_kinetics,
            setup_solvation,
            setup_reactants,
            setup_product,
            test_data_system_1_experimental,
    ):
        experimental_excel_data = pd.read_excel(test_data_system_1_experimental)

        reactor = Reactor(
            reactions=construct_system_1,
            kinetics=setup_kinetics,
            solvation=setup_solvation,
        )

        substrate, nucleophilic = setup_reactants

        for _, row in experimental_excel_data.iterrows():
            time = row["tres/min"]
            ratio = row["2:1"]
            concentration = row["Conc 1/M"]
            temperature = row["Temp/°C"]
            expected_E = row["E-factor"]
            exptected_STY = row["STY/kg m-3 h-1"]

            concentrations = {
                substrate: concentration,
                nucleophilic: concentration * ratio,
            }

            reactor_conditions = ReactorConditions(
                temperature=temperature,
                concentrations=concentrations,
                products=[setup_product,],
                time=time,
            )

            STY, E = reactor.simulate(reactor_conditions)

            assert np.isclose(
                E, expected_E, rtol=0.15, atol=0.3
            ), (
                f"(Conditions: '{row}') "
                f"at {temperature}K: Calculated={E:.4e}, "
                f"Excel={expected_E:.4e}"
            )
            print("="*20)
            print("STY:", STY)
            print("Excel STY:", exptected_STY)
            print("E:", E)
            print("Excel E:", expected_E)
            print("="*20)

            assert np.isclose(
                STY, exptected_STY, rtol=0.10, atol=1250
            ), (
                f"(Conditions: '{row}') "
                f"at {temperature}K: Calculated={STY:.4e}, "
                f"Excel={exptected_STY:.4e}"
            )




