from pathlib import Path
from abc import abstractmethod, ABC

from quantumpareto.species import (
    Reaction,
    Species,
)

_logo = """
 ██████╗ ██╗   ██╗ █████╗ ███╗   ██╗████████╗██╗   ██╗███╗   ███╗
██╔═══██╗██║   ██║██╔══██╗████╗  ██║╚══██╔══╝██║   ██║████╗ ████║
██║   ██║██║   ██║███████║██╔██╗ ██║   ██║   ██║   ██║██╔████╔██║
██║▄▄ ██║██║   ██║██╔══██║██║╚██╗██║   ██║   ██║   ██║██║╚██╔╝██║
╚██████╔╝╚██████╔╝██║  ██║██║ ╚████║   ██║   ╚██████╔╝██║ ╚═╝ ██║
 ╚══▀▀═╝  ╚═════╝ ╚═╝  ╚═╝╚═╝  ╚═══╝   ╚═╝    ╚═════╝ ╚═╝     ╚═╝

██████╗  █████╗ ██████╗ ███████╗████████╗ ██████╗ 
██╔══██╗██╔══██╗██╔══██╗██╔════╝╚══██╔══╝██╔═══██╗
██████╔╝███████║██████╔╝█████╗     ██║   ██║   ██║
██╔═══╝ ██╔══██║██╔══██╗██╔══╝     ██║   ██║   ██║
██║     ██║  ██║██║  ██║███████╗   ██║   ╚██████╔╝
╚═╝     ╚═╝  ╚═╝╚═╝  ╚═╝╚══════╝   ╚═╝    ╚═════╝ 

                Multi-Objective Quantum Optimization
                              Version 0.1.0

"""


class DataOutput(ABC):

    @abstractmethod
    def dump(self, path: Path) -> None:
        """Dumps the data which the class holds."""


class ReactionInput:

    reactions: list[Reaction]
    """Contains the Reactions present in the System"""

    def __init__(self, reactions: list[Reaction]):
        self.reactions = reactions
        self._check_input()

    def _check_input(self):

        for reaction in self.reactions:
            if not isinstance(reaction, Reaction):
                raise TypeError(
                    "Reaction must be a list of quantumpareto.Reactions"
                )

            if not reaction.reactants:
                raise ValueError(
                    "Reaction must have at least one reactant"
                )

            if not reaction.transition_state:
                raise ValueError(
                    "Reaction must have one transition state"
                )

    @property
    def species(self) -> set[Species]:
        species = set()
        for reaction in self.reactions:
            species.update(reaction.species)
        return species

    @property
    def transition_states(self) -> set[Species]:
        transition_states = set()
        for reaction in self.reactions:
            transition_states.add(reaction.transition_state)
        return transition_states