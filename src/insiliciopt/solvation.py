import re
import logging
import numpy as np

from insiliciopt.species import (
    Reaction,
    Species,
)
from insiliciopt.utils import ReactionInput


class Solvation(ReactionInput):

    gas_constant = 8.3145  # J/mol/K
    """Gas Constant"""

    _g_values: dict[Species, dict[float, float]]
    """G Values mapping dict[Species, dict[temperature, gsolv]]"""

    def __init__(
        self,
        reactions: list[Reaction],
    ) -> None:
        super().__init__(reactions)
        self._logger = logging.getLogger(__name__)

        self._g_values = {}
        self._extract_g()

    @staticmethod
    def _parse_cosmo_therm_file(species: Species) -> dict[float, float]:
        temperature_g_solve = {}
        content = species.tab_file_path.read_text()
        lines = content.split("\n")

        temperature = None

        for i, line in enumerate(lines):
            line = line.strip()
            if not line:
                continue

            if "Settings " in line:
                parts = line.split(";")[0]
                parts = parts.split("=")[1]
                parts = parts.strip()
                parts = parts.replace("K", "")
                temperature = float(parts)

            if "Compound" in line:
                parts = lines[i+2].split()
                parts = parts[5]
                gsolv = float(parts)
                temperature_g_solve[temperature] = gsolv
                temperature = None

        return temperature_g_solve

    def _extract_g(self):
        all_species = self.species.copy()
        all_species.update(self.transition_states)
        for species in all_species:
            self._g_values[species] = self._parse_cosmo_therm_file(species)

    def _g(self, species: Species, temperature: float) -> float:
        """Obtain G for one species"""
        return float(
            np.interp(
                temperature,
                list(self._g_values[species].keys()),
                list(self._g_values[species].values()),
            )
        )

    def correction_factor(
            self, reaction: Reaction, temperature: float
    ) -> float:
        """Obtain the delt G for one reaction in SI units"""

        reactants = reaction.reactants
        transition_state = reaction.transition_state

        g_reactants = sum(
            [self._g(reactant, temperature) for reactant in reactants]
        )
        g_transition_state = self._g(transition_state, temperature)

        delta_g = g_transition_state - g_reactants
        delta_g_si_units = delta_g * 4184

        return np.exp(-delta_g_si_units / (self.gas_constant * temperature))
