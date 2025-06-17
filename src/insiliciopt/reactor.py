import logging
from dataclasses import dataclass

from pyomo.core import BuildAction, TransformationFactory
from pyomo.opt import SolverFactory
from pyomo.environ import (
    ConcreteModel,
    Var,
    NonNegativeReals,
    Constraint,
    value,
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
    concentrations: dict[Species, float]
    """Dictionary of initial concentrations, all other species will be set 
    to 0
    """
    products: list[Species]
    """Products"""
    time: float
    """Time in Minutes"""


class Reactor(ReactionInput):

    conditions: ReactorConditions | None
    """Reactor Conditions"""

    kinetics: Kinetics
    """Kinetics"""

    solvation: Solvation
    """Solvation delta G"""

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
        self.conditions = None

    def __repr__(self):
        string = "\n\nREACTOR:\n"
        string += f"  Num Reactions: {len(self.reactions)}\n"
        string += f"  Num Species: {len(self.species)}\n"
        string += f"  Num Transition States: {len(self.transition_states)}\n"
        return string

    def _convert_conditions(
            self,
            conditions: ReactorConditions
    ) -> ReactorConditions:
        self._logger.debug(
            f"Conditions: {conditions}"
        )

        conditions_converted = ReactorConditions(
            temperature=conditions.temperature + 273.15,  # Kelvin
            concentrations=conditions.concentrations,
            products=conditions.products,
            time=conditions.time * 60, # Seconds
        )
        self.conditions = conditions_converted

        self._logger.debug(
            f"Converted conditions: {conditions_converted}"
        )
        return conditions_converted

    def __rate_rule(
            self,
            model: ConcreteModel,
            reaction: Reaction,
            time: float
    ) -> float:
        k = self.kinetics.k(reaction, self.conditions.temperature)
        correction = self.solvation.correction_factor(
            reaction, self.conditions.temperature
        )
        rate = k * correction
        for reactant in reaction.reactants:
            rate *= model.C[reactant, time]
        return rate

    def __mass_balance(
            self,
            model: ConcreteModel,
            species: Species,
            time: float,
    ) -> bool:
        return model.dCdt[species, time] == sum(
            reaction.coefficient(species) * self.__rate_rule(
                model, reaction, time
            )
            for reaction in self.reactions
        )

    def __init_conditions(self, model) -> None:
        for species in self.species:
            concentration = self.conditions.concentrations.get(species, 0.0)
            model.C[species, 0].fix(concentration)

    def _setup_reactor(self) -> ConcreteModel:
        """Setup the pyomo reactor"""
        model = ConcreteModel()
        model.t = ContinuousSet(bounds=(0, self.conditions.time))
        model.C = Var(self.species, model.t, domain=NonNegativeReals)
        model.dCdt = DerivativeVar(model.C, wrt=model.t)
        model.mass_balance = Constraint(
            self.species, model.t, rule=self.__mass_balance
        )
        model.init = BuildAction(rule=lambda m: self.__init_conditions(m))
        return model

    def _extract_results(self, model: ConcreteModel) -> tuple[float, float]:
        concentration = {
            sp: value(model.C[sp, self.conditions.time])
            for sp in self.species
        }
        mass_product = sum(
            concentration[species] * species.mass
            for species in self.conditions.products
        )
        mass_waste = sum(
            concentration[species] * species.mass
            for species in self.species
            if species not in self.conditions.products
        )

        sty = 3600 * mass_product / self.conditions.time
        e_factor = mass_waste / mass_product

        self._logger.debug(
            f"Result: E={e_factor} - STY={sty}"
        )

        return sty, e_factor

    def simulate(self, conditions: ReactorConditions) -> tuple[float, float]:
        """Returns E and STY for given starting conditions"""
        self._convert_conditions(conditions)
        model = self._setup_reactor()

        TransformationFactory(
            "dae.finite_difference"
        ).apply_to(model, nfe=200, scheme="BACKWARD")
        solver = SolverFactory("ipopt")
        solver.solve(model, tee=False)

        return self._extract_results(model)

