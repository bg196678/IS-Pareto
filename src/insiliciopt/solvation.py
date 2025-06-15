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

    def __init__(
        self,
        reactions: list[Reaction],
    ) -> None:
        super().__init__(reactions)
        self._logger = logging.getLogger(__name__)
        self._reactions = reactions

    def _g(self, species: Species, temperature: float) -> float:
        """Obtain G for one species"""
        # TODO
        raise NotImplementedError


    def correction_factor(
            self, reaction: Reaction, temperature: float
    ) -> float:
        """Obtain the delt G for one reaction in SI units"""

        reactants = reaction.reactants
        transition_state = reaction.transition_state

        G_reactants = sum(
            [self._g(reactant, temperature) for reactant in reactants]
        )
        G_transition_state = self._g(transition_state, temperature)

        delta_G = G_transition_state - G_reactants
        delta_G_SI = delta_G * 4184

        return np.exp(-delta_G_SI / (self.gas_constant * temperature))


