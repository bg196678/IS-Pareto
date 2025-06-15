import pytest
import numpy as np

from insiliciopt.species import (
    Reactant,
    Product,
    TransitionState,
    Reaction,
)
from insiliciopt.kinetics import Kinetics

STOICHIOMETRY = {
    'R1_fwd': {'Substrate': -1, 'Nucleophilic': -1, 'ITS1': 1},
    'R1_rev': {'ITS1': -1, 'Substrate': 1, 'Nucleophilic': 1},
    'R2_fwd': {'Substrate': -1, 'Nucleophilic': -1, 'ITS2': 1},
    'R2_rev': {'ITS2': -1, 'Substrate': 1, 'Nucleophilic': 1},
    'R3_fwd': {'Product2': -1, 'Nucleophilic': -1, 'ITS3': 1},
    'R3_rev': {'ITS3': -1, 'Product2': 1, 'Nucleophilic': 1},
    'R4_fwd': {'Product1': -1, 'Nucleophilic': -1, 'ITS4': 1},
    'R4_rev': {'ITS4': -1, 'Product1': 1, 'Nucleophilic': 1},
    'R5_fwd': {'ITS1': -1, 'Product1': 1, 'Leaving_Group': 1},
    'R5_rev': {'Product1': -1, 'Leaving_Group': -1, 'ITS1': 1},
    'R6_fwd': {'ITS2': -1, 'Product2': 1, 'Leaving_Group': 1},
    'R6_rev': {'Product2': -1, 'Leaving_Group': -1, 'ITS2': 1},
    'R7_fwd': {'ITS3': -1, 'Product3': 1, 'Leaving_Group': 1},
    'R7_rev': {'Product3': -1, 'Leaving_Group': -1, 'ITS3': 1},
    'R8_fwd': {'ITS4': -1, 'Product3': 1, 'Leaving_Group': 1},
    'R8_rev': {'Product3': -1, 'Leaving_Group': -1, 'ITS4': 1},
}

REACTION_TS_MAP = {
    'R1_fwd': 'TS1_fwd',
    'R1_rev': 'TS1_rev',
    'R2_fwd': 'TS2_fwd',
    'R2_rev': 'TS2_rev',
    'R3_fwd': 'TS3_fwd',
    'R3_rev': 'TS3_rev',
    'R4_fwd': 'TS4_fwd',
    'R4_rev': 'TS4_rev',
    'R5_fwd': 'TS12_fwd',
    'R5_rev': 'TS12_rev',
    'R6_fwd': 'TS22_fwd',
    'R6_rev': 'TS22_rev',
    'R7_fwd': 'TS32_fwd',
    'R7_rev': 'TS32_rev',
    'R8_fwd': 'TS42_fwd',
    'R8_rev': 'TS42_rev',
}


class TestKinetics:
    """Test suite for the Kinetics class"""

    @pytest.fixture
    def species_dict(self, test_data_system_1_path):
        """Create dictionary of all species from test data"""
        species = {}

        regular_species = [
            'Substrate', 'Nucleophilic', 'Product1', 'Product2',
            'Product3', 'Leaving_Group', 'ITS1', 'ITS2', 'ITS3', 'ITS4'
        ]

        for name in regular_species:
            fchk_path = test_data_system_1_path / f"{name}.fchk"
            if name.startswith('ITS'):
                species[name] = Product(
                    name=name,
                    mass=100.0,
                    fchk_file_path=fchk_path
                )
            elif name in ['Product1', 'Product2', 'Product3', 'Leaving_Group']:
                species[name] = Product(
                    name=name,
                    mass=100.0,
                    fchk_file_path=fchk_path
                )
            else:
                species[name] = Reactant(
                    name=name,
                    mass=100.0,
                    fchk_file_path=fchk_path
                )

        for ts_name in REACTION_TS_MAP.values():
            fchk_path = test_data_system_1_path / f"{ts_name}.fchk"
            species[ts_name] = TransitionState(
                name=ts_name,
                fchk_file_path=fchk_path
            )

        return species

    @pytest.fixture
    def reactions(self, species_dict):
        """Create list of reactions from stoichiometry and TS mapping"""
        reactions = []

        for rxn_name, stoich in STOICHIOMETRY.items():
            rxn_stoich = {}
            for species_name, coeff in stoich.items():
                rxn_stoich[species_dict[species_name]] = coeff

            ts_name = REACTION_TS_MAP[rxn_name]
            ts = species_dict[ts_name]

            reaction = Reaction(
                name=rxn_name,
                stoichiometry=rxn_stoich,
                transition_state=ts
            )
            reactions.append(reaction)

        return reactions

    @pytest.fixture
    def simple_reaction(self, species_dict):
        """Create a simple reaction for basic tests"""
        stoich = {
            species_dict['Substrate']: -1,
            species_dict['Nucleophilic']: -1,
            species_dict['ITS1']: 1
        }
        return Reaction(
            name='R1_fwd',
            stoichiometry=stoich,
            transition_state=species_dict['TS1_fwd']
        )

    def test_kinetics_initialization(self, simple_reaction):
        """Test basic initialization of Kinetics object"""
        kinetics = Kinetics(reactions=[simple_reaction])

        assert kinetics.reactions == [simple_reaction]
        assert kinetics.tunneling_correction is None
        assert kinetics.gradient_threshold == 4e-3

    def test_kinetics_with_tunneling(self, simple_reaction):
        """Test initialization with different tunneling corrections"""
        kinetics_wigner = Kinetics(
            reactions=[simple_reaction],
            tunneling_correction='wigner'
        )
        assert kinetics_wigner.tunneling_correction == 'wigner'

        kinetics_eckart = Kinetics(
            reactions=[simple_reaction],
            tunneling_correction='eckart'
        )
        assert kinetics_eckart.tunneling_correction == 'eckart'

        kinetics_miller = Kinetics(
            reactions=[simple_reaction],
            tunneling_correction='miller'
        )
        assert kinetics_miller.tunneling_correction == 'miller'

    def test_invalid_tunneling_correction(self, simple_reaction):
        """Test that invalid tunneling correction raises error"""
        with pytest.raises(ValueError, match="Tunneling correction must be"):
            Kinetics(
                reactions=[simple_reaction],
                tunneling_correction='invalid'
            )

    def test_reaction_without_reactants(self, species_dict):
        """Test that reaction without reactants raises error"""
        stoich = {species_dict['Product1']: 1}
        reaction = Reaction(
            name='invalid',
            stoichiometry=stoich,
            transition_state=species_dict['TS1_fwd']
        )

        with pytest.raises(ValueError, match="at least one reactant"):
            Kinetics(reactions=[reaction])

    def test_reaction_without_transition_state(self, species_dict):
        """Test that reaction without transition state raises error"""
        stoich = {
            species_dict['Substrate']: -1,
            species_dict['Product1']: 1
        }
        reaction = Reaction(
            name='invalid',
            stoichiometry=stoich,
            transition_state=None
        )

        with pytest.raises(ValueError, match="one transition state"):
            Kinetics(reactions=[reaction])

    def test_invalid_reaction_type(self):
        """Test that non-Reaction objects raise error"""
        with pytest.raises(TypeError, match="must be a list of"):
            Kinetics(reactions=["not a reaction"])

    @pytest.mark.parametrize("temperature", [298.15, 373.15, 473.15])
    def test_rate_constant_calculation(self, simple_reaction, temperature):
        """Test rate constant calculation at different temperatures"""
        kinetics = Kinetics(reactions=[simple_reaction])

        k = kinetics.k(simple_reaction, temperature)
        assert k > 0
        assert np.isfinite(k)

    def test_rate_constant_temperature_dependence(self, simple_reaction):
        """Test that rate constant increases with temperature"""
        kinetics = Kinetics(reactions=[simple_reaction])

        temperatures = [298.15, 373.15, 473.15]
        rate_constants = [kinetics.k(simple_reaction, T) for T in temperatures]

        for i in range(len(rate_constants) - 1):
            assert rate_constants[i + 1] > rate_constants[i]

    @pytest.mark.parametrize("tunneling", ["wigner", "eckart", "miller"])
    def test_tunneling_effect_on_rate(self, simple_reaction, tunneling):
        """Test that tunneling corrections increase rate constant"""
        kinetics_no_tunnel = Kinetics(reactions=[simple_reaction])
        kinetics_with_tunnel = Kinetics(
            reactions=[simple_reaction],
            tunneling_correction=tunneling
        )

        temperature = 298.15
        k_no_tunnel = kinetics_no_tunnel.k(simple_reaction, temperature)
        k_with_tunnel = kinetics_with_tunnel.k(simple_reaction, temperature)

        assert k_with_tunnel >= k_no_tunnel

    def test_multiple_reactions(self, reactions):
        """Test kinetics with multiple reactions"""
        # Take first 4 reactions for faster test
        test_reactions = reactions[:4]
        kinetics = Kinetics(reactions=test_reactions)

        assert len(kinetics._kinetics_models) == len(test_reactions)

        temperature = 298.15
        for reaction in test_reactions:
            k = kinetics.k(reaction, temperature)
            assert k > 0
            assert np.isfinite(k)

    def test_gradient_threshold_custom(self, simple_reaction):
        """Test custom gradient threshold"""
        custom_threshold = 1e-4
        kinetics = Kinetics(
            reactions=[simple_reaction],
            gradient_threshold=custom_threshold
        )

        assert kinetics.gradient_threshold == custom_threshold

    def test_energy_correction(self, test_data_system_1_path, species_dict):
        """Test species with energy correction"""
        # Create species with energy correction
        substrate_with_energy = Reactant(
            name='Substrate_corrected',
            mass=100.0,
            fchk_file_path=test_data_system_1_path / 'Substrate.fchk',
            energy=-123.456  # Example energy correction
        )

        stoich = {
            substrate_with_energy: -1,
            species_dict['Product1']: 1
        }

        reaction = Reaction(
            name='test_energy',
            stoichiometry=stoich,
            transition_state=species_dict['TS1_fwd']
        )

        kinetics = Kinetics(reactions=[reaction])
        k = kinetics.k(reaction, 298.15)
        assert k > 0

    def test_arrhenius_behavior(self, simple_reaction):
        """Test that rate constants follow Arrhenius-like behavior"""
        kinetics = Kinetics(reactions=[simple_reaction])

        temperatures = np.linspace(250, 400, 10)
        rate_constants = [kinetics.k(simple_reaction, T) for T in temperatures]

        ln_k = np.log(rate_constants)
        inv_T = 1.0 / temperatures
        correlation = np.corrcoef(inv_T, ln_k)[0, 1]

        assert correlation < -0.9

    @pytest.mark.parametrize("rxn_name", ['R1_fwd', 'R5_fwd', 'R8_rev'])
    def test_specific_reactions(self, reactions, rxn_name):
        """Test specific reactions from the system"""
        reaction = next(r for r in reactions if r.name == rxn_name)
        kinetics = Kinetics(reactions=[reaction])

        k = kinetics.k(reaction, 298.15)
        assert k > 0

        k_hot = kinetics.k(reaction, 373.15)
        assert k_hot > k
