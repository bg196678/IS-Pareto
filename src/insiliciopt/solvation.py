import logging
import numpy as np
import pandas as pd
from pathlib import Path

from insiliciopt.species import (
    Reaction,
    Species,
)
from insiliciopt.utils import (
    ReactionInput,
    DataOutput,
)


class Solvation(ReactionInput, DataOutput):

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

    def __repr__(self) -> str:
        string = "\n\nSOLVATION:\n"
        string += f"  Num Reactions: {len(self.reactions)}\n"
        string += f"  Num Species: {len(self.species)}\n"
        string += f"  Num Transition States: {len(self.transition_states)}\n"
        string += f"  Gas Constant: {self.gas_constant}\n"
        return string

    def _parse_cosmo_therm_file(self, species: Species) -> dict[float, float]:
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

        self._logger.debug(
            f"Extracted {len(temperature_g_solve)} Temperature G Solvation "
            f"values"
        )

        return temperature_g_solve

    def _extract_g(self):
        all_species = self.species.copy()
        all_species.update(self.transition_states)
        for species in all_species:
            self._g_values[species] = self._parse_cosmo_therm_file(species)

        self._logger.debug(
            f"Extracted G Solvation values for {len(self._g_values)} species"
        )

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

    def _dump_gsolv(self, path: Path) -> None:
        g_values = {
            key.name: value for key, value in self._g_values.items()
        }
        df = pd.DataFrame.from_dict(g_values, orient="index")
        df.to_csv(path / "gsolv.csv")
        self._logger.debug("Dumped Gsolv to csv file.")

    def _dump_correction(self, path: Path) -> None:
        """Dumps the correction factors for each reaction at various
        temperatures to a CSV file.
        """
        all_temps = list(self._g_values.values())[0].keys()

        correction_values = {}
        for reaction in self.reactions:
            correction_values[reaction.name] = {}
            for temp in all_temps:
                correction = self.correction_factor(reaction, temp)
                correction_values[reaction.name][temp] = correction

        df = pd.DataFrame.from_dict(correction_values, orient="index")
        df.to_csv(path / "corrections.csv")
        self._logger.debug(f"Dumped correction factors to csv.")

    def dump(self, path: Path) -> None:
        """Dumps the Gsolv values and correction factors"""
        self._dump_gsolv(path)
        self._dump_correction(path)

