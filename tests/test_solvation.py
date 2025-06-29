import pytest
import numpy as np
import pandas as pd

from quantumpareto.solvation import Solvation


class TestSolvation:
    @pytest.fixture
    def solvation_instance(self, construct_system_1):
        """Fixture to create a Solvation instance with system 1 reactions."""
        return Solvation(reactions=construct_system_1)

    def test_solvation_initialization(self, solvation_instance, construct_system_1):
        """Test if Solvation class initializes and extracts G values."""
        assert solvation_instance._g_values is not None
        all_species_in_reactions = set()
        for reaction in construct_system_1:
            all_species_in_reactions.update(reaction.reactants)
            all_species_in_reactions.add(reaction.transition_state)

        expected_species_keys = solvation_instance.species.union(
            solvation_instance.transition_states)


        assert len(solvation_instance._g_values) == len(expected_species_keys)
        for species in expected_species_keys:
            assert species in solvation_instance._g_values
            assert isinstance(solvation_instance._g_values[species], dict)
            assert len(solvation_instance._g_values[species]) > 0

    def test_parse_cosmo_therm_file(
            self, construct_system_1, solvation_instance
    ):
        """Test the _parse_cosmo_therm_file static method."""
        test_species = construct_system_1[0].reactants[0]

        g_values = solvation_instance._parse_cosmo_therm_file(test_species)
        assert isinstance(g_values, dict)
        for temp, g_val in g_values.items():
            assert isinstance(temp, float)
            assert isinstance(g_val, float)
        if g_values:
             assert len(g_values) > 0


    def test_g_method_interpolation(self, solvation_instance, construct_system_1):
        """Test the _g method for a species and temperature, including interpolation."""
        test_species = construct_system_1[0].reactants[0]

        known_temp = list(solvation_instance._g_values[test_species].keys())[0]
        expected_g = solvation_instance._g_values[test_species][known_temp]
        assert np.isclose(
            solvation_instance._g(test_species, known_temp), expected_g
        )
        temps = sorted(list(solvation_instance._g_values[test_species].keys()))
        g_vals = [
            solvation_instance._g_values[test_species][t] for t in temps
        ]
        interp_temp = (temps[0] + temps[1]) / 2.0
        expected_interp_g = g_vals[0] + (g_vals[1] - g_vals[0]) * \
                            ((interp_temp - temps[0]) / (temps[1] - temps[0]))
        assert np.isclose(
            solvation_instance._g(test_species, interp_temp), expected_interp_g
        )


    def test_correction_factor(self, solvation_instance, construct_system_1):
        """Test the correction_factor method for a reaction."""
        test_reaction = construct_system_1[0]
        temperature = 298.15

        g_reactants_sum = sum(
            [solvation_instance._g(r, temperature) for r in test_reaction.reactants]
        )
        g_ts = solvation_instance._g(test_reaction.transition_state, temperature)

        delta_g_kcal = g_ts - g_reactants_sum
        delta_g_si = delta_g_kcal * 4184

        expected_factor = np.exp(
            -delta_g_si / (solvation_instance.gas_constant * temperature)
        )

        calculated_factor = solvation_instance.correction_factor(test_reaction, temperature)
        assert np.isclose(calculated_factor, expected_factor)

    def test_correction_factor_multiple_reactants(self, construct_system_1, test_data_system_1_cosmo):
        """Test correction_factor with a reaction having multiple reactants."""
        reaction_1_fwd = construct_system_1[0]
        assert len(reaction_1_fwd.reactants) == 2

        solv = Solvation(reactions=[reaction_1_fwd])
        temperature = 298.15

        substrate = reaction_1_fwd.reactants[0]
        nucleophilic = reaction_1_fwd.reactants[1]
        ts1_fwd = reaction_1_fwd.transition_state

        g_substrate = solv._g(substrate, temperature)
        g_nucleophilic = solv._g(nucleophilic, temperature)
        g_ts1_fwd = solv._g(ts1_fwd, temperature)

        delta_g_kcal = g_ts1_fwd - (g_substrate + g_nucleophilic)
        delta_g_si = delta_g_kcal * 4184

        expected_factor = np.exp(-delta_g_si / (Solvation.gas_constant * temperature))
        calculated_factor = solv.correction_factor(reaction_1_fwd, temperature)

        assert np.isclose(calculated_factor, expected_factor)

class TestSystem1Solvation:

    def test_compare_gsolv(
            self,
            construct_system_1,
            test_data_system_1_gsolv,
    ):
        gsolv_excel_data = pd.read_excel(test_data_system_1_gsolv)

        for reaction in construct_system_1:
            solvation = Solvation([reaction])
            reactants = reaction.reactants
            transition_state = reaction.transition_state
            species = reactants + [transition_state]

            for _, row in gsolv_excel_data.iterrows():
                temperature = row['Temperature (K)']

                for reactant in species:
                    calculated_g = solvation._g(reactant, temperature)
                    expected_g = row[reactant.name]
                    assert np.isclose(
                        calculated_g,
                        expected_g,
                        rtol=0.01, atol=1e-9
                    ), (
                        f"Gsolv mismatch for species '{reactant.name}' "
                        f"(Species: '{reactant}') "
                        f"at {temperature}K: "
                        f"Calculated={calculated_g:.4e}, "
                        f"Excel={expected_g:.4e}"
                    )

    def test_compare_correction_factor(
            self,
            construct_system_1,
            test_data_system_1_gsolv,
    ):
        """Compares calculated correction factors with values from the
        gsolv.xlsx sheet.
        """
        gsolv_excel_data = pd.read_excel(test_data_system_1_gsolv)

        for reaction in construct_system_1:
            solvation = Solvation([reaction])
            reactants = reaction.reactants
            transition_state = reaction.transition_state

            for _, row in gsolv_excel_data.iterrows():
                temperature = row['Temperature (K)']

                expected_reactants_gsolv = sum([
                    row[reactant.name] for reactant in reactants
                ])
                calculated_reactans_gsolv = sum([
                    solvation._g(reactant, temperature)
                    for reactant in reactants
                ])
                assert np.isclose(
                    calculated_reactans_gsolv, expected_reactants_gsolv,
                    rtol=0.01, atol=1e-9
                ), (
                    f"Gibbs mismatch for reaction '{reaction.name}' "
                    f"(Species: '{reaction}') "
                    f"at {temperature}K: "
                    f"Calculated={calculated_reactans_gsolv:.4e}, "
                    f"Excel={expected_reactants_gsolv:.4e}"
                )

                expected_transition_state_gsolv = row[transition_state.name]
                calculated_transition_state_gsolv = solvation._g(
                    transition_state, temperature
                )

                assert np.isclose(
                    calculated_transition_state_gsolv,
                    expected_transition_state_gsolv,
                    rtol=0.01, atol=1e-9
                ), (
                    f"Gibbs mismatch for reaction '{reaction.name}' "
                    f"(Species: '{reaction}') "
                    f"at {temperature}K: "
                    f"Calculated={calculated_transition_state_gsolv:.4e}, "
                    f"Excel={expected_transition_state_gsolv:.4e}"
                )


                expected_delta_g = (
                        expected_transition_state_gsolv -
                        expected_reactants_gsolv
                )
                expected_delta_g_si_units = expected_delta_g * 4184
                expected_correction_factor = np.exp(
                    - expected_delta_g_si_units / (8.3145 * temperature)
                )

                calculated_correction_factor = solvation.correction_factor(
                    reaction, temperature
                )

                assert np.isclose(
                    calculated_correction_factor, expected_correction_factor,
                    rtol=0.01, atol=1e-9
                ), (
                    f"Correction Factor mismatch for reaction "
                    f"'{reaction.name}' "
                    f"(Reaction: '{reaction}') "
                    f"at {temperature}K: "
                    f"Calculated={calculated_correction_factor:.4e}, "
                    f"Excel={expected_correction_factor:.4e}"
                )






