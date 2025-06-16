import logging
import pandas as pd
from pathlib import Path
from dataclasses import dataclass
from abc import abstractmethod, ABC
from summit.domain import (
    ContinuousVariable,
    Domain,
)
from summit.strategies import (
    TSEMO,
    LHS,
)
from summit.utils.dataset import DataSet

from insiliciopt.reactor import (
    Reactor,
    ReactorConditions,
)
from insiliciopt.species import Species


@dataclass
class OptimizationSpecies:
    reactant_1: Species
    """Reactant 1"""
    reactant_2: Species
    """Reactant 2"""
    products: list[Species]
    """List of Products"""

@dataclass
class OptimizationBoundaries:
    """Optimization Boundaries"""
    temperature: tuple[float, float]
    """Temperature bounds in Celsius"""
    concentration_reactant_1: tuple[float, float]
    """Concentration Reactant 1 Boundaries"""
    concentration_ratio: tuple[float, float]
    """Concentration Ration Boundaries"""
    time: tuple[float, float]
    """Time Boundaries in Minutes"""
    STY: tuple[float, float] = (0, 1e8)
    """Space Time Yield Objective"""
    E_factor: tuple[float, float] = (0, 50)
    """Mass product waste factor objective"""

@dataclass
class _OptimizerConditions:
    """Optimizer Conditions"""
    temperature: float
    """Temperature in Celsius"""
    concentration_reactant_1: float
    """Concentration of Reactant 1"""
    concentration_ratio: float
    """Concentration Ratio"""
    time: float
    """Time in Minutes"""

class Optimizer(ABC):

    species: OptimizationSpecies
    """Optimization Species"""

    boundaries: OptimizationBoundaries
    """Optimization Boundaries"""

    reactor: Reactor
    """Reactor Model to optimize"""

    output_directory: Path
    """Stores the data and plots in this directory"""

    results: pd.DataFrame
    """Dataframe holding the results"""

    _results_columns_names: list[str] = [
        "Temperature [Celsius]",
        "Concentration [%]",
        "Concentration Ratio",
        "Time [Minutes]",
        "STY",
        "E",
    ]

    _counter: int = 0
    """Counter of iterations"""

    def __init__(
            self,
            species: OptimizationSpecies,
            boundaries: OptimizationBoundaries,
            reactor: Reactor,
            output_directory: Path,
    ) -> None:
        self._logger: logging.Logger = logging.getLogger(__class__.__name__)

        self.species = species
        self.boundaries = boundaries
        self.reactor = reactor
        self.output_directory = output_directory

        self.results = pd.DataFrame(
            columns=[
                "Temperature",
                "Concentration Reactant 1",
                "Concentration Ratio",
                "Time",
            ]
        )

    def _convert_optimizer_to_reactor_conditions(
            self, conditions: _OptimizerConditions,
    ) -> ReactorConditions:
        """Convert optimizer conditions to reactor conditions"""
        concentrations = {
            self.species.reactant_1: conditions.concentration_reactant_1,
            self.species.reactant_2: (
                conditions.concentration_ratio *
                conditions.concentration_reactant_1
            ),
        }
        rector_conditions = ReactorConditions(
            temperature=conditions.temperature,
            concentrations=concentrations,
            products=self.species.products,
            time=conditions.time,
        )
        return rector_conditions

    def _store_results(self) -> None:
        self.results.to_csv(
            self.output_directory / "results.csv",
            index=False,
        )

    def _add_to_result(
            self,
            E: float,
            STY: float,
            conditions: _OptimizerConditions,
    ) -> None:
        """Adds the results to the Dataframe"""
        results_columns_dataframe = pd.DataFrame(
            [
                {
                    self._results_columns_names[0]: conditions.temperature,
                    self._results_columns_names[1]: (
                        conditions.concentration_reactant_1
                    ),
                    self._results_columns_names[2]: (
                        conditions.concentration_ratio
                    ),
                    self._results_columns_names[3]: conditions.time,
                    self._results_columns_names[4]: STY,
                    self._results_columns_names[5]: E,
                }
            ]
        )
        self.results = pd.concat(
            [self.results, results_columns_dataframe],
            ignore_index=True,
        )

        self._counter += 1

    def _simulate_reactor(
            self, conditions: _OptimizerConditions
    ) -> tuple[float, float]:
        """Runs the reactor model with given conditions"""
        reactor_conditions = self._convert_optimizer_to_reactor_conditions(
            conditions
        )
        return self.reactor.simulate(reactor_conditions)

    @abstractmethod
    def run(self, num_iterations: int) -> None:
        """Runs the optimizer"""


class TSEmoOptimizer(Optimizer):

    domain: Domain
    """Optimization Domain"""

    num_lhs_points: int
    """Number of initial LHS points"""

    _columns: list[str]
    """Domain column Names"""

    def __init__(
            self,
            species: OptimizationSpecies,
            boundaries: OptimizationBoundaries,
            reactor: Reactor,
            output_directory: Path,
            num_lhs_points: int = 20,
    ) -> None:
        super().__init__(species, boundaries, reactor, output_directory)

        self.num_lhs_points = num_lhs_points

    def _create_domain(self):
        domain = Domain()
        domain += ContinuousVariable(
            self._results_columns_names[0],
            bounds=list(self.boundaries.temperature),
            description="Temperature Celsius",
        )
        domain += ContinuousVariable(
            self._results_columns_names[1],
            bounds=list(self.boundaries.concentration_reactant_1),
            description="Concentration Reactant 1",
        )
        domain += ContinuousVariable(
            self._results_columns_names[2],
            bounds=list(self.boundaries.concentration_ratio),
            description="Concentration Ratio",
        )
        domain += ContinuousVariable(
            self._results_columns_names[3],
            bounds=list(self.boundaries.time),
            description="Time in Minutes",
        )
        domain += ContinuousVariable(
            self._results_columns_names[4],
            is_objective=True,
            bounds=list(self.boundaries.STY),
            description="Space Time Yield Objective",
        )
        domain += ContinuousVariable(
            self._results_columns_names[5],
            is_objective=True,
            bounds=list(self.boundaries.E_factor),
            description="Mass product waste factor objective",
        )
        self._columns = [variable.name for variable in domain.variables]
        self.domain = domain

    def _construct_optimization_conditions(
            self, suggestion: DataSet
    ) -> _OptimizerConditions:
        """Converst summit Dataset to Optimizer Conditions"""
        data = suggestion.to_dict()["data"][0]
        temperature = data[self._columns.index("Temperature")]
        concentration = data[self._columns.index("Concentration")]
        ratio = data[self._columns.index("Concentration Ratio")]
        time = data[self._columns.index("Time")]

        optimization_conditions = _OptimizerConditions(
            temperature=temperature,
            concentration_reactant_1=concentration,
            concentration_ratio=ratio,
            time=time,
        )
        return optimization_conditions

    def _run_lhs(self) -> None:
        """Initial LHS sampling"""
        lhs_sampler = LHS(self.domain)
        initial_suggestions = lhs_sampler.suggest_experiments(
            self.num_lhs_points,
        )

        for i in range(self.num_lhs_points):
            suggestion = initial_suggestions.iloc[[i]]
            conditions = self._construct_optimization_conditions(suggestion)
            STY, E = self._simulate_reactor(conditions)
            self._add_to_result(E, STY, conditions)

    def _run_tsemo(self, num_iterations: int) -> None:
        """Optimized TSEMO sampling"""
        tsemo_sampler = TSEMO(self.domain)

        for i in range(num_iterations):

            dataset = DataSet.from_dataframe(self.results)
            suggestion = tsemo_sampler.suggest_experiments(
                1, prev_res=dataset,
            )
            conditions = self._construct_optimization_conditions(suggestion)
            STY, E = self._simulate_reactor(conditions)
            self._add_to_result(E, STY, conditions)

    def run(self, num_iterations: int) -> None:
        """Runs the TSEMO optimizer"""
        self._run_lhs()
        self._run_tsemo(num_iterations=num_iterations)
