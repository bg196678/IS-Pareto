import logging

from insiliciopt.species import Reaction


class Kinetics:

    reactions: list[Reaction]
    """Contains the Reactions present in the System"""

    use_tunneling: bool
    """Includes Tunneling into the Kinetics calculation"""

    def __init__(
            self,
            reactions: list[Reaction],
            use_tunneling: bool = True,
    ):
        self._logger: logging.Logger = logging.getLogger(__class__.__name__)

        self.reactions = reactions
