from insiliciopt.species import (
    Reaction,
    Species,
)


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
                    "Reaction must be a list of insiliciopt.Reactions"
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