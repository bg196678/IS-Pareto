import pytest
import numpy as np
import pandas as pd

from insiliciopt.kinetics import Kinetics
from insiliciopt.reactor import (
    Reactor,
    ReactorConditions,
)
from insiliciopt.species import Species
from insiliciopt.solvation import Solvation


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
    def setup_products(
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

        return product_1, product_2, product_3

    def test_reactor(
            self,
            construct_system_1,
            setup_kinetics,
            setup_solvation,
            setup_reactants,
            setup_products,
            test_data_system_1_experimental,
    ):
        experimental_excel_data = pd.read_excel(test_data_system_1_experimental)

        reactor = Reactor(
            reactions=construct_system_1,
            kinetics=setup_kinetics,
            solvation=setup_solvation,
        )

        substrate, nucleophilic = setup_reactants
        product_1, product_2, product_3 = setup_products

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
                products=[
                    product_1,
                    product_2,
                    product_3,
                ],
                time=time,
            )

            STY, E = reactor.simulate(reactor_conditions)

            # assert np.isclose(
            #     E, expected_E, rtol=0.15, atol=0.3
            # ), (
            #     f"(Conditions: '{row}') "
            #     f"at {temperature}K: Calculated={E:.4e}, "
            #     f"Excel={expected_E:.4e}"
            # )
            print("="*20)
            print("STY:", STY)
            print("Excel STY:", exptected_STY)
            print("E:", E)
            print("Excel E:", expected_E)
            print("="*20)

            # assert np.isclose(
            #     STY, exptected_STY, rtol=0.10, atol=1250
            # ), (
            #     f"(Conditions: '{row}') "
            #     f"at {temperature}K: Calculated={STY:.4e}, "
            #     f"Excel={exptected_STY:.4e}"
            # )




