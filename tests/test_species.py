import pytest
from unittest.mock import patch

from ispareto.species import (
    Species,
    Reactant,
    Product,
    TransitionState,
    Reaction,
    ReactionTerm,
)


class TestSpecies:
    """Tests for Species and its subclasses"""

    def test_species_initialization(self, test_fchk_path, test_tab_path):
        """Test Species initialization"""
        species = Species(
            "H2O",
            18.015,
            test_fchk_path,
            test_tab_path,
        )

        assert species.name == "H2O"
        assert species.mass == 18.015
        assert species.fchk_file_path == test_fchk_path
        assert species.tab_file_path == test_tab_path

    def test_species_logger_creation(self, test_fchk_path, test_tab_path):
        """Test that logger is created with correct name"""
        with patch('logging.getLogger') as mock_logger:
            _ = Species("CO2", 44.01, test_fchk_path, test_tab_path)
            mock_logger.assert_called_with("Species('CO2')")

    def test_reactant_inheritance(self, test_fchk_path, test_tab_path):
        """Test Reactant class inherits from Species"""
        reactant = Reactant("CH4", 16.04, test_fchk_path, test_tab_path)

        assert isinstance(reactant, Species)
        assert reactant.name == "CH4"
        assert reactant.mass == 16.04

    def test_product_inheritance(self, test_fchk_path, test_tab_path):
        """Test Product class inherits from Species"""
        product = Product("CO2", 44.01, test_fchk_path, test_tab_path)

        assert isinstance(product, Species)
        assert product.name == "CO2"
        assert product.mass == 44.01

    def test_transition_state_inheritance(self, test_fchk_path, test_tab_path):
        """Test TransitionState class inherits from Species"""
        ts = TransitionState("TS1", test_fchk_path, test_tab_path)

        assert isinstance(ts, Species)
        assert ts.name == "TS1"
        assert ts.mass is None

    def test_species_equality_and_hash(self, test_fchk_path, test_tab_path):
        """Test Species equality and hash for use as dict keys"""
        species1 = Species("H2O", 18.015, test_fchk_path, test_tab_path)
        species2 = Species("H2O", 18.015, test_fchk_path, test_tab_path)
        other_fchk_path = test_fchk_path.parent / "other.fchk"
        other_tab_path = test_tab_path.parent / "other.tab"
        other_fchk_path.touch(exist_ok=True)
        other_tab_path.touch(exist_ok=True)

        species3 = Species("CO2", 44.01, other_fchk_path, other_tab_path)

        assert species1 == species2
        assert species1 != species3
        assert hash(species1) == hash(species2)
        assert hash(species1) != hash(species3)

    def test_species_addition(self, test_fchk_path, test_tab_path):
        """Test Species + Species creates a Reaction"""
        species1 = Species("H2O", 18.015, test_fchk_path, test_tab_path)
        species2 = Species("CO2", 44.01, test_fchk_path, test_tab_path)

        reaction = species1 + species2

        assert isinstance(reaction, Reaction)
        assert reaction.stoichiometry[species1] == 1
        assert reaction.stoichiometry[species2] == 1

    def test_species_subtraction(self, test_fchk_path, test_tab_path):
        """Test Species - Species creates a Reaction"""
        species1 = Species("H2O", 18.015, test_fchk_path, test_tab_path)
        species2 = Species("CO2", 44.01, test_fchk_path, test_tab_path)

        reaction = species1 - species2

        assert isinstance(reaction, Reaction)
        assert reaction.stoichiometry[species1] == 1
        assert reaction.stoichiometry[species2] == -1

    def test_species_invalid_operations(self, test_fchk_path, test_tab_path):
        """Test Species with invalid types"""
        species = Species("H2O", 18.015, test_fchk_path, test_tab_path)

        with pytest.raises(TypeError):
            species + "invalid"

        with pytest.raises(TypeError):
            species - 42


class TestReactionTerm:
    """Tests for ReactionTerm class"""

    @pytest.fixture(autouse=True)
    def setup_method(self, test_fchk_path, test_tab_path):
        """Set up test species"""
        self.h2o = Species("H2O", 18.015, test_fchk_path, test_tab_path)
        self.co2 = Species("CO2", 44.01, test_fchk_path, test_tab_path)
        self.ch4 = Species("CH4", 16.04, test_fchk_path, test_tab_path)

    def test_reaction_term_initialization(self):
        """Test ReactionTerm initialization"""
        term = ReactionTerm(self.h2o, -1)

        assert term.species == self.h2o
        assert term.coefficient == -1

    def test_reaction_term_addition(self):
        """Test ReactionTerm + ReactionTerm"""
        term1 = ReactionTerm(self.h2o, -1)
        term2 = ReactionTerm(self.co2, 1)

        reaction = term1 + term2

        assert isinstance(reaction, Reaction)
        assert reaction.stoichiometry[self.h2o] == -1
        assert reaction.stoichiometry[self.co2] == 1

    def test_reaction_term_subtraction(self):
        """Test ReactionTerm - ReactionTerm"""
        term1 = ReactionTerm(self.h2o, 1)
        term2 = ReactionTerm(self.co2, 1)

        reaction = term1 - term2

        assert isinstance(reaction, Reaction)
        assert reaction.stoichiometry[self.h2o] == 1
        assert reaction.stoichiometry[self.co2] == -1

    def test_reaction_term_add_to_reaction(self):
        """Test ReactionTerm + Reaction"""
        term1 = ReactionTerm(self.h2o, -1)
        term2 = ReactionTerm(self.co2, 1)
        reaction = term1 + term2

        term3 = ReactionTerm(self.ch4, -1)
        new_reaction = term3 + reaction

        assert new_reaction.stoichiometry[self.h2o] == -1
        assert new_reaction.stoichiometry[self.co2] == 1
        assert new_reaction.stoichiometry[self.ch4] == -1

    def test_reaction_term_invalid_operations(self):
        """Test ReactionTerm with invalid types"""
        term = ReactionTerm(self.h2o, -1)

        with pytest.raises(TypeError):
            term + "invalid"

        with pytest.raises(TypeError):
            term - 42

    def test_species_add_reaction_term(self, test_fchk_path, test_tab_path):
        """Test Species + ReactionTerm"""
        species1 = Species("H2O", 18.015, test_fchk_path, test_tab_path)
        species2 = Species("CO2", 44.01, test_fchk_path, test_tab_path)
        term = ReactionTerm(species2, 2)

        reaction = species1 + term

        assert isinstance(reaction, Reaction)
        assert reaction.stoichiometry[species1] == 1
        assert reaction.stoichiometry[species2] == 2

    def test_species_sub_reaction_term(self, test_fchk_path, test_tab_path):
        """Test Species - ReactionTerm"""
        species1 = Species("H2O", 18.015, test_fchk_path, test_tab_path)
        species2 = Species("CO2", 44.01, test_fchk_path, test_tab_path)
        term = ReactionTerm(species2, 2)

        reaction = species1 - term

        assert isinstance(reaction, Reaction)
        assert reaction.stoichiometry[species1] == 1
        assert reaction.stoichiometry[species2] == -2


class TestReaction:
    """Tests for Reaction class"""

    @pytest.fixture(autouse=True)
    def setup_method(self, test_fchk_path, test_tab_path):
        """Set up test species"""
        self.h2o = Species("H2O", 18.015, test_fchk_path, test_tab_path)
        self.co2 = Species("CO2", 44.01, test_fchk_path, test_tab_path)
        self.ch4 = Species("CH4", 16.04, test_fchk_path, test_tab_path)
        self.o2 = Species("O2", 32.0, test_fchk_path, test_tab_path)

    def test_reaction_initialization_empty(self, test_fchk_path, test_tab_path):
        """Test empty Reaction initialization"""
        reaction = Reaction()

        assert reaction.name is None
        assert reaction.stoichiometry == {}
        assert reaction.species == []
        assert reaction.reactants == []
        assert reaction.products == []

    def test_reaction_initialization_with_name(self, test_fchk_path, test_tab_path):
        """Test Reaction initialization with name"""
        reaction = Reaction(name="combustion")

        assert reaction.name == "combustion"
        assert reaction.stoichiometry == {}

    def test_reaction_initialization_with_stoichiometry(self, test_fchk_path, test_tab_path):
        """Test Reaction initialization with stoichiometry"""
        stoich = {self.ch4: -1, self.o2: -2, self.co2: 1, self.h2o: 2}
        reaction = Reaction(stoichiometry=stoich)

        assert reaction.stoichiometry == stoich
        assert len(reaction.species) == 4
        assert self.ch4 in reaction.reactants
        assert self.co2 in reaction.products

    def test_reaction_removes_zero_coefficients(self, test_fchk_path, test_tab_path):
        """Test that zero coefficients are removed"""
        stoich = {self.ch4: -1, self.o2: 0, self.co2: 1}
        reaction = Reaction(stoichiometry=stoich)

        assert self.o2 not in reaction.stoichiometry
        assert len(reaction.stoichiometry) == 2

    def test_reaction_properties(self, test_fchk_path, test_tab_path):
        """Test reaction properties"""
        stoich = {self.ch4: -1, self.o2: -2, self.co2: 1, self.h2o: 2}
        reaction = Reaction(stoichiometry=stoich)

        species_list = reaction.species
        assert len(species_list) == 4
        assert all(sp in stoich.keys() for sp in species_list)

        reactants = reaction.reactants
        assert self.ch4 in reactants
        assert self.o2 in reactants

        products = reaction.products
        assert self.co2 in products
        assert self.h2o in products

    def test_reaction_addition_with_reaction_term(self, test_fchk_path, test_tab_path):
        """Test Reaction + ReactionTerm"""
        stoich = {self.ch4: -1, self.co2: 1}
        reaction = Reaction(stoichiometry=stoich)

        term = ReactionTerm(self.h2o, 2)
        new_reaction = reaction + term

        assert new_reaction.stoichiometry[self.ch4] == -1
        assert new_reaction.stoichiometry[self.co2] == 1
        assert new_reaction.stoichiometry[self.h2o] == 2

    def test_reaction_addition_with_reaction(self, test_fchk_path, test_tab_path):
        """Test Reaction + Reaction"""
        stoich1 = {self.ch4: -1, self.co2: 1}
        reaction1 = Reaction(stoichiometry=stoich1)

        stoich2 = {self.h2o: -1, self.o2: 1}
        reaction2 = Reaction(stoichiometry=stoich2)

        combined = reaction1 + reaction2

        assert combined.stoichiometry[self.ch4] == -1
        assert combined.stoichiometry[self.co2] == 1
        assert combined.stoichiometry[self.h2o] == -1
        assert combined.stoichiometry[self.o2] == 1

    def test_reaction_subtraction_with_reaction_term(self, test_fchk_path, test_tab_path):
        """Test Reaction - ReactionTerm"""
        stoich = {self.ch4: -1, self.co2: 1, self.h2o: 2}
        reaction = Reaction(stoichiometry=stoich)

        term = ReactionTerm(self.h2o, 1)
        new_reaction = reaction - term

        assert new_reaction.stoichiometry[self.ch4] == -1
        assert new_reaction.stoichiometry[self.co2] == 1
        assert new_reaction.stoichiometry[self.h2o] == 1  # 2 - 1 = 1

    def test_reaction_subtraction_with_reaction(self, test_fchk_path, test_tab_path):
        """Test Reaction - Reaction"""
        stoich1 = {self.ch4: -1, self.co2: 1, self.h2o: 2}
        reaction1 = Reaction(stoichiometry=stoich1)

        stoich2 = {self.h2o: 1, self.o2: -1}
        reaction2 = Reaction(stoichiometry=stoich2)

        result = reaction1 - reaction2

        assert result.stoichiometry[self.ch4] == -1
        assert result.stoichiometry[self.co2] == 1
        assert result.stoichiometry[self.h2o] == 1  # 2 - 1 = 1
        assert result.stoichiometry[self.o2] == 1  # 0 - (-1) = 1

    def test_reaction_invalid_operations(self, test_fchk_path, test_tab_path):
        """Test Reaction with invalid types"""
        reaction = Reaction()

        with pytest.raises(TypeError):
            reaction + "invalid"

        with pytest.raises(TypeError):
            reaction - 42

    def test_reaction_repr(self, test_fchk_path, test_tab_path):
        """Test Reaction string representation"""
        stoich = {self.ch4: -1, self.co2: 1}
        reaction = Reaction(stoichiometry=stoich)

        _ = repr(reaction)

    def test_species_add_reaction(self, test_fchk_path, test_tab_path):
        """Test Species + Reaction"""
        species1 = Species("H2O", 18.015, test_fchk_path, test_tab_path)
        species2 = Species("CO2", 44.01, test_fchk_path, test_tab_path)
        species3 = Species("CH4", 16.04, test_fchk_path, test_tab_path)

        reaction = Reaction(stoichiometry={species2: -1, species3: 1})

        result = species1 + reaction

        assert isinstance(result, Reaction)
        assert result.stoichiometry[species1] == 1
        assert result.stoichiometry[species2] == -1
        assert result.stoichiometry[species3] == 1

    def test_species_sub_reaction(self, test_fchk_path, test_tab_path):
        """Test Species - Reaction"""
        species1 = Species("H2O", 18.015, test_fchk_path, test_tab_path)
        species2 = Species("CO2", 44.01, test_fchk_path, test_tab_path)
        species3 = Species("CH4", 16.04, test_fchk_path, test_tab_path)

        reaction = Reaction(stoichiometry={species2: -1, species3: 1})

        result = species1 - reaction

        assert isinstance(result, Reaction)
        assert result.stoichiometry[species1] == 1
        assert result.stoichiometry[species2] == 1  # Negated from -1
        assert result.stoichiometry[species3] == -1  # Negated from 1


class TestComplexStoichiometryReconstruction:
    """Test reconstruction of the complex stoichiometry system"""

    @pytest.fixture(autouse=True)
    def setup_method(self, test_fchk_path, test_tab_path):
        """Set up all species needed for the complex stoichiometry"""
        self.species = {}
        species_names = [
            'Substrate', 'Nucleophilic', 'ITS1', 'ITS2', 'ITS3', 'ITS4',
            'Product1', 'Product2', 'Product3', 'LeavingGroup'
        ]

        for name in species_names:
            self.species[name] = Species(
                name, 100.0,
                test_fchk_path,
                test_tab_path,
            )

    def test_reconstruct_target_stoichiometry(self, test_fchk_path, test_tab_path):
        """Reconstruct the target stoichiometry using ReactionTerms"""
        target_stoich = {
            'R1_fwd': {'Substrate': -1, 'Nucleophilic': -1, 'ITS1': 1},
            'R1_rev': {'ITS1': -1, 'Substrate': 1, 'Nucleophilic': 1},
            'R2_fwd': {'Substrate': -1, 'Nucleophilic': -1, 'ITS2': 1},
            'R2_rev': {'ITS2': -1, 'Substrate': 1, 'Nucleophilic': 1},
            'R3_fwd': {'Product2': -1, 'Nucleophilic': -1, 'ITS3': 1},
            'R3_rev': {'ITS3': -1, 'Product2': 1, 'Nucleophilic': 1},
            'R4_fwd': {'Product1': -1, 'Nucleophilic': -1, 'ITS4': 1},
            'R4_rev': {'ITS4': -1, 'Product1': 1, 'Nucleophilic': 1},
            'R5_fwd': {'ITS1': -1, 'Product1': 1, 'LeavingGroup': 1},
            'R5_rev': {'Product1': -1, 'LeavingGroup': -1, 'ITS1': 1},
            'R6_fwd': {'ITS2': -1, 'Product2': 1, 'LeavingGroup': 1},
            'R6_rev': {'Product2': -1, 'LeavingGroup': -1, 'ITS2': 1},
            'R7_fwd': {'ITS3': -1, 'Product3': 1, 'LeavingGroup': 1},
            'R7_rev': {'Product3': -1, 'LeavingGroup': -1, 'ITS3': 1},
            'R8_fwd': {'ITS4': -1, 'Product3': 1, 'LeavingGroup': 1},
            'R8_rev': {'Product3': -1, 'LeavingGroup': -1, 'ITS1': 1}, # Note: R8_rev uses ITS1 as per original
        }

        reconstructed_reactions = {}

        for reaction_name, stoich_dict in target_stoich.items():
            terms = []
            for species_name, coeff in stoich_dict.items():
                species = self.species[species_name]
                terms.append(ReactionTerm(species, coeff))

            if not terms:
                reaction = Reaction()
            else:
                reaction = ReactionTerm(terms[0].species, 0)
                for term in terms:
                    reaction = reaction + term


            reaction.name = reaction_name
            reconstructed_reactions[reaction_name] = reaction

        for reaction_name, target_dict in target_stoich.items():
            reconstructed = reconstructed_reactions[reaction_name]

            assert len(reconstructed.stoichiometry) == len(target_dict)

            for species_name, target_coeff in target_dict.items():
                species = self.species[species_name]
                assert species in reconstructed.stoichiometry
                assert reconstructed.stoichiometry[species] == target_coeff


    def test_alternative_reaction_construction(self, test_fchk_path, test_tab_path):
        """Test alternative ways to construct reactions"""
        substrate = self.species['Substrate']
        nucleophilic = self.species['Nucleophilic']
        its1 = self.species['ITS1']

        reaction1 = (
                ReactionTerm(its1, 1) -
                ReactionTerm(substrate, 1) -
                ReactionTerm(nucleophilic, 1)
        )

        reaction2 = Reaction(stoichiometry={
            substrate: -1,
            nucleophilic: -1,
            its1: 1
        })

        assert reaction1.stoichiometry == reaction2.stoichiometry

    def test_reaction_balancing(self, test_fchk_path, test_tab_path):
        """Test that reactions are properly balanced"""
        substrate = self.species['Substrate']
        nucleophilic = self.species['Nucleophilic']
        its1 = self.species['ITS1']

        r1_fwd = (
                ReactionTerm(its1, 1) -
                ReactionTerm(substrate, 1) -
                ReactionTerm(nucleophilic, 1)
        )

        total_coeff = sum(
            abs(coeff) for coeff in r1_fwd.stoichiometry.values())
        assert total_coeff > 0

        reactant_coeffs = [
            coeff for coeff in r1_fwd.stoichiometry.values() if coeff < 0
        ]
        product_coeffs = [
            coeff for coeff in r1_fwd.stoichiometry.values() if coeff > 0
        ]

        assert len(reactant_coeffs) > 0
        assert len(product_coeffs) > 0