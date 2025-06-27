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
            energy=-635.334856191927,
        )

        nucleophilic = Species(
            name="Nucleophilic",
            mass=0.071,
            fchk_file_path=test_data_system_1_gaussian / "Nucleophilic.fchk",
            tab_file_path=test_data_system_1_cosmo / "Nucleophilic.tab",
            energy=-212.567127673868,
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
            energy=-747.459040149113,
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
            test_data_system_1_experimental
        )

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
            energy=-635.334856191927,
        )

        nucleophilic = Species(
            name="Nucleophilic",
            mass=0.071,
            fchk_file_path=test_data_system_1_gaussian / "Nucleophilic.fchk",
            tab_file_path=test_data_system_1_cosmo / "Nucleophilic.tab",
            energy=-212.567127673868,
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
            energy=-747.459040149113,
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

        # Simulated Reactor Data
        E_reactor = []
        STY_reactor = []

        # Experimental Data
        E_experimental = []
        STY_experimental = []

        # Conditions
        temps = []
        ratios = []
        concs = []

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
                products=[setup_product,],
                time=time,
            )

            STY, E = reactor.simulate(reactor_conditions)

            E_reactor.append(E)
            STY_reactor.append(STY)
            E_experimental.append(expected_E)
            STY_experimental.append(exptected_STY)
            temps.append(temperature)
            ratios.append(ratio)
            concs.append(concentration)

            assert np.isclose(
                E, expected_E, rtol=0.15, atol=0.3
            ), (
                f"(Conditions: '{row}') "
                f"at {temperature}K: Calculated={E:.4e}, "
                f"Excel={expected_E:.4e}"
            )

            assert np.isclose(
                STY, exptected_STY, rtol=0.10, atol=1250
            ), (
                f"(Conditions: '{row}') "
                f"at {temperature}K: Calculated={STY:.4e}, "
                f"Excel={exptected_STY:.4e}"
            )

        import matplotlib.pyplot as plt

        E_reactor = np.array(E_reactor)
        STY_reactor = np.array(STY_reactor)
        E_experimental = np.array(E_experimental)
        STY_experimental = np.array(STY_experimental)

        E_abs_errors = np.abs(E_reactor - E_experimental)
        E_rel_errors = np.abs(
            E_reactor - E_experimental) / E_experimental * 100
        STY_abs_errors = np.abs(STY_reactor - STY_experimental)
        STY_rel_errors = np.abs(
            STY_reactor - STY_experimental
        ) / STY_experimental * 100

        E_correlation = np.corrcoef(E_reactor, E_experimental)[0, 1]
        STY_correlation = np.corrcoef(STY_reactor, STY_experimental)[0, 1]
        E_r2 = E_correlation ** 2
        STY_r2 = STY_correlation ** 2

        print("\n" + "=" * 50)
        print("ERROR STATISTICS")
        print("=" * 50)

        print("\nE-FACTOR ERRORS:")
        print(f"Mean Absolute Error: {np.mean(E_abs_errors):.4f}")
        print(f"Std Dev of Absolute Error: {np.std(E_abs_errors):.4f}")
        print(
            f"RMSE: {np.sqrt(np.mean((E_reactor - E_experimental) ** 2)):.4f}")
        print(f"Mean Relative Error: {np.mean(E_rel_errors):.2f}%")
        print(f"Max Absolute Error: {np.max(E_abs_errors):.4f}")
        print(f"Min Absolute Error: {np.min(E_abs_errors):.4f}")
        print(f"R²: {E_r2:.4f}")

        print("\nSTY ERRORS:")
        print(f"Mean Absolute Error: {np.mean(STY_abs_errors):.4f}")
        print(f"Std Dev of Absolute Error: {np.std(STY_abs_errors):.4f}")
        print(
            f"RMSE: {np.sqrt(np.mean((STY_reactor - STY_experimental) ** 2)):.4f}")
        print(f"Mean Relative Error: {np.mean(STY_rel_errors):.2f}%")
        print(f"Max Absolute Error: {np.max(STY_abs_errors):.4f}")
        print(f"Min Absolute Error: {np.min(STY_abs_errors):.4f}")
        print(f"R²: {STY_r2:.4f}")

        # E-factor: Model vs Experimental
        plt.figure(figsize=(8, 6))
        plt.scatter(E_experimental, E_reactor, alpha=0.6)
        plt.plot([E_experimental.min(), E_experimental.max()],
                 [E_experimental.min(), E_experimental.max()], 'r--')
        plt.xlabel('Experimental E-factor')
        plt.ylabel('Model E-factor')
        plt.title('E-factor: Model vs Experimental')
        plt.tight_layout()
        plt.savefig("e_factor_model_vs_experimental.png")
        plt.show()

        # E-factor Relative Error Distribution
        plt.figure(figsize=(8, 6))
        plt.hist(E_rel_errors, bins=20, edgecolor='black')
        plt.xlabel('Relative Error (%)')
        plt.ylabel('Frequency')
        plt.title('E-factor Relative Error Distribution')
        plt.tight_layout()
        plt.savefig("e_factor_relative_error_distribution.png")
        plt.show()

        # E-factor Absolute Error Box Plot
        plt.figure(figsize=(8, 6))
        plt.boxplot([E_abs_errors], labels=['E-factor'])
        plt.ylabel('Absolute Error')
        plt.title('E-factor Absolute Error Box Plot')
        plt.tight_layout()
        plt.savefig("e_factor_absolute_error_boxplot.png")
        plt.show()

        # STY: Model vs Experimental
        plt.figure(figsize=(8, 6))
        plt.scatter(STY_experimental, STY_reactor, alpha=0.6)
        plt.plot([STY_experimental.min(), STY_experimental.max()],
                 [STY_experimental.min(), STY_experimental.max()], 'r--')
        plt.xlabel('Experimental STY')
        plt.ylabel('Model STY')
        plt.title('STY: Model vs Experimental')
        plt.tight_layout()
        plt.savefig("sty_model_vs_experimental.png")
        plt.show()

        # STY Relative Error Distribution
        plt.figure(figsize=(8, 6))
        plt.hist(STY_rel_errors, bins=20, edgecolor='black')
        plt.xlabel('Relative Error (%)')
        plt.ylabel('Frequency')
        plt.title('STY Relative Error Distribution')
        plt.tight_layout()
        plt.savefig("sty_relative_error_distribution.png")
        plt.show()

        # STY Absolute Error Box Plot
        plt.figure(figsize=(8, 6))
        plt.boxplot([STY_abs_errors], labels=['STY'])
        plt.ylabel('Absolute Error')
        plt.title('STY Absolute Error Box Plot')
        plt.tight_layout()
        plt.savefig("sty_absolute_error_boxplot.png")
        plt.show()

        df_sim = pd.DataFrame({
            "Temp": temps,
            "Ratio_2_1": ratios,
            "Conc": concs,
            "STY": STY_reactor,
            "E_factor": E_reactor
        })

        df_exp = pd.DataFrame({
            "Temp": temps,
            "Ratio_2_1": ratios,
            "Conc": concs,
            "STY": STY_experimental,
            "E_factor": E_experimental
        })

        params = ["Temp", "Ratio_2_1", "Conc"]

        corr_exp = df_exp[params + ["STY", "E_factor"]].corr()
        corr_sim = df_sim[params + ["STY", "E_factor"]].corr()

        sty_corr = pd.DataFrame({
            "Experiment": corr_exp["STY"].drop("STY"),
            "Simulation": corr_sim["STY"].drop("STY")
        })

        e_corr = pd.DataFrame({
            "Experiment": corr_exp["E_factor"].drop("E_factor"),
            "Simulation": corr_sim["E_factor"].drop("E_factor")
        })

        plt.figure(figsize=(8, 5))
        sty_corr.plot(kind="bar", ax=plt.gca())
        plt.title("Correlation of Parameters with STY")
        plt.ylabel("Correlation Coefficient")
        plt.axhline(0, color="black", linewidth=0.8)
        plt.xticks(rotation=0)
        plt.tight_layout()
        plt.savefig("sty_correlation.png")
        plt.show()

        plt.figure(figsize=(8, 5))
        e_corr.plot(kind="bar", ax=plt.gca())
        plt.title("Correlation of Parameters with E-Factor")
        plt.ylabel("Correlation Coefficient")
        plt.axhline(0, color="black", linewidth=0.8)
        plt.xticks(rotation=0)
        plt.tight_layout()
        plt.savefig("e_factor_correlation.png")
        plt.show()




