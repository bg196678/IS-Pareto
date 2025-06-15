import logging
from dataclasses import dataclass
from pyomo.opt import SolverFactory
from pyomo.environ import (
    ConcreteModel,
    Var,
    NonNegativeReals,
)
from pyomo.dae.contset import ContinuousSet
from pyomo.dae.diffvar import DerivativeVar

from insiliciopt.species import (
    Reaction,
    Species,
)
from insiliciopt.kinetics import Kinetics
from insiliciopt.solvation import Solvation
from insiliciopt.utils import ReactionInput

@dataclass
class ReactorConditions:
    temperature: float
    """Temperature in Celsius"""
    concentration_reactant_1: float
    """Concentration in #TODO"""
    reactant_1: Species
    """Reactant 1"""
    reactant_2: Species
    """Reactant 2"""
    ratio: float
    """Ratio"""
    time: float
    """Time in Minutes"""
    concentration_reactant_2: float | None = None
    """Concentration Reactant 2 in # TODO"""


class Reactor(ReactionInput):

    kinetics: Kinetics

    def __init__(
            self,
            reactions: list[Reaction],
            kinetics: Kinetics,
            solvation: Solvation,
    ):
        super().__init__(reactions)

        self._logger: logging.Logger = logging.getLogger(__name__)

        self.kinetics = kinetics
        self.solvation = solvation

    @staticmethod
    def _convert_conditions(
            conditions: ReactorConditions
    ) -> ReactorConditions:
        conditions_converted = ReactorConditions(
            temperature=conditions.temperature + 273.15,  # Kelvin
            concentration_reactant_1=conditions.concentration_reactant_1,
            reactant_1=conditions.reactant_1,
            reactant_2=conditions.reactant_2,
            ratio=conditions.ratio,
            time=conditions.time * 60, # Seconds
            concentration_reactant_2=(
                conditions.concentration_reactant_2 * conditions.ratio
            ),
        )
        return conditions_converted

    def _setup_reactor(self, conditions: ReactorConditions) -> ConcreteModel:
        """Setup the pyomo reactor"""
        model = ConcreteModel()
        model.t = ContinuousSet(bounds=(0, conditions.time))
        model.C = Var(self.species, model.t, domain=NonNegativeReals)
        model.dCdt = DerivativeVar(model.C, wrt=model.t)

        return model

    def simulate(self, conditions: ReactorConditions) -> tuple[float, float]:
        """Returns E and STY for given starting conditions"""
        conditions_converted = self._convert_conditions(conditions)
        model = self._setup_reactor(conditions_converted)