import logging
import numpy as np
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
    SOBO,
    ENTMOOT,
)
from summit.utils.dataset import DataSet

from insiliciopt.reactor import (
    Reactor,
    ReactorConditions,
)
from insiliciopt.species import Species
from insiliciopt.utils import _logo
from insiliciopt.visualization import Visualization


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
        "Temperature[Celsius]",
        "Concentration[%]",
        "ConcentrationRatio",
        "Time[Minutes]",
        "STY",
        "E",
    ]

    _counter: int = 0
    """Counter of iterations"""

    _log_level: int = logging.INFO
    """Logging level"""

    _output_data_directory: Path
    """Output Data Storage"""

    _visualization: Visualization
    """Visualization of the Pareto Front"""

    def __init__(
            self,
            species: OptimizationSpecies,
            boundaries: OptimizationBoundaries,
            reactor: Reactor,
            output_directory: Path,
            log_level: int = logging.INFO,
            num_initial_points: int | None = None,
    ) -> None:
        self._logger: logging.Logger = logging.getLogger(__class__.__name__)

        self.species = species
        self.boundaries = boundaries
        self.reactor = reactor
        self.output_directory = output_directory
        self._log_level = log_level

        self.results = pd.DataFrame(
            columns=self._results_columns_names
        )
        self._setup_logger()
        self._log_start()

        # Data
        self._output_data_directory = self.output_directory / "data"
        self._output_data_directory.mkdir(parents=True, exist_ok=True)
        self.reactor.kinetics.dump(self._output_data_directory)
        self.reactor.solvation.dump(self._output_data_directory)

        # Visualization
        self._output_plot_directory = self.output_directory / "plots"
        self._output_plot_directory.mkdir(parents=True, exist_ok=True)
        self._visualization = Visualization(
            self._output_plot_directory,
            num_initial_points=num_initial_points,
        )

    @abstractmethod
    def __repr__(self) -> str:
        """Optimizer"""

    def _setup_logger(self) -> None:
        """Adds a console and file handler"""
        self._logger.handlers = []

        formatter = logging.Formatter(
            '%(asctime)s - %(name)s - %(levelname)s\n'
            '%(message)s'
        )
        self._logger.setLevel(self._log_level)

        console_handler = logging.StreamHandler()
        console_handler.setLevel(self._log_level)
        console_handler.setFormatter(formatter)
        self._logger.addHandler(console_handler)

        log_file = self.output_directory / "insiliciopt.log"
        if log_file.exists():
            log_file.unlink()

        file_handler = logging.FileHandler(log_file)
        file_handler.setLevel(self._log_level)
        file_handler.setFormatter(formatter)
        self._logger.addHandler(file_handler)

    def _log_start(self) -> None:
        """Log initial message with detailed optimization information"""
        string = _logo

        # Species to Optimize
        string += "\nOPTIMIZATION SPECIES:\n"
        string += f"    Reactant 1: {self.species.reactant_1.name}\n"
        string += f"    Reactant 2: {self.species.reactant_2.name}\n"
        string += (
            f"    Products: "
            f"{', '.join([product.name for product in self.species.products])}\n\n"
        )

        # Optimization Boundaries
        string += "\nOPTIMIZATION BOUNDARIES:\n"
        string += f"    Temperature: {self.boundaries.temperature} °C\n"
        string += (
            f"    Concentration Reactant 1: "
                   f"{self.boundaries.concentration_reactant_1} %\n"
        )
        string += f"    Concentration Ratio: {self.boundaries.concentration_ratio}\n"
        string += f"    Time: {self.boundaries.time} minutes\n"
        string += f"    STY: {self.boundaries.STY}\n"
        string += f"    E_factor: {self.boundaries.E_factor}"

        string += f"\n{self.reactor}"
        string += f"{self.reactor.kinetics}"
        string += f"{self.reactor.solvation}"

        string += "\n\nOUTPUT DIRECTORY:\n"
        string += f"    {self.output_directory}\n"

        string += f"{self}\n\n"

        string += "REACTIONS:\n"
        for reaction in self.reactor.reactions:
            string += f"   {reaction}\n"

        string += "\n\nSPECIES:\n"
        for species in self.reactor.species:
            string += f"   {species}\n"

        string += "\n\nTRANSITION STATES:\n"
        for transition_state in self.reactor.transition_states:
            string += f"   {transition_state}\n"

        string += "\n\n===== STARTING OPTIMIZATION =====\n\n"

        self._logger.info(string)

    def _log_iteration(
            self,
            E: float,
            STY: float,
            conditions: _OptimizerConditions,
            iteration_type: str | None = None,
    ) -> None:
        """Log iteration"""
        iteration_str = "=" * 35
        iteration_str += f"\nITERATION {self._counter}\n"

        if iteration_type:
            iteration_str += f"    Type: {iteration_type}\n"

        iteration_str += "    Conditions:\n"
        iteration_str += (
            f"        "
            f"Temperature: {conditions.temperature} °C\n"
        )
        iteration_str += (
            f"        Concentration Reactant 1: "
            f"{conditions.concentration_reactant_1}\n"
        )
        iteration_str += (
            f"        "
            f"Concentration Ratio: {conditions.concentration_ratio:.2f}\n"
        )
        iteration_str += f"        Time: {conditions.time:.2f} minutes\n"

        iteration_str += "    Results:\n"
        iteration_str += f"        STY: {STY:.4e}\n"
        iteration_str += f"        E-factor: {E:.4f}\n"

        self._logger.info(iteration_str)

    def _log_summary(self) -> None:
        """Logs a summary of the optimization"""
        summary = "\n\n===== OPTIMIZATION SUMMARY =====\n\n"
        summary += f"Total iterations: {self._counter}\n\n"

        best_sty_idx = self.results[self._results_columns_names[4]].idxmax()
        best_sty_row = self.results.iloc[best_sty_idx]

        summary += "BEST SPACE TIME YIELD (STY) SOLUTION:\n"
        summary += (
            f"    STY: "
            f"{best_sty_row[self._results_columns_names[4]]}\n"
        )
        summary += (
            f"    "
            f"E-factor: {best_sty_row[self._results_columns_names[5]]}\n"
        )
        summary += (
            f"    "
            f"Temperature: {best_sty_row[self._results_columns_names[0]]} °C\n"
        )
        summary += (
            f"    "
            f"Concentration Reactant 1: "
            f"{best_sty_row[self._results_columns_names[1]]}\n"
        )
        summary += (
            f"    "
            f"Concentration Ratio: "
            f"{best_sty_row[self._results_columns_names[2]]}\n"
        )
        summary += (
            f"    "
            f"Time: "
            f"{best_sty_row[self._results_columns_names[3]]:.2f} minutes\n\n"
        )

        best_e_idx = self.results[self._results_columns_names[5]].idxmin()
        best_e_row = self.results.iloc[best_e_idx]

        summary += "BEST E-FACTOR SOLUTION:\n"
        summary += (
            f"    "
            f"E-factor: {best_e_row[self._results_columns_names[5]]:.4f}\n"
        )
        summary += (
            f"    "
            f"STY: {best_e_row[self._results_columns_names[4]]:.4e}\n"
        )
        summary += (
            f"    "
            f"Temperature: {best_e_row[self._results_columns_names[0]]:.2f} °C\n"
        )
        summary += (
            f"    "
            f"Concentration Reactant 1: "
            f"{best_e_row[self._results_columns_names[1]]:.4f}\n"
        )
        summary += (
            f"    "
            f"Concentration Ratio: "
            f"{best_e_row[self._results_columns_names[2]]:.4f}\n"
        )
        summary += (
            f"    "
            f"Time: "
            f"{best_e_row[self._results_columns_names[3]]:.2f} minutes\n\n"
        )

        summary += "STATISTICS:\n"
        summary += (
            f"    "
            f"STY Range: "
            f"[{self.results[self._results_columns_names[4]].min():.4e}, "
            f"{self.results[self._results_columns_names[4]].max():.4e}]\n"
        )
        summary += (
            f"    "
            f"E-factor Range: "
            f"[{self.results[self._results_columns_names[5]].min():.4f}, "
            f"{self.results[self._results_columns_names[5]].max():.4f}]\n"
        )
        summary += (
            f"    "
            f"Temperature Range: "
            f"[{self.results[self._results_columns_names[0]].min():.2f}, "
            f"{self.results[self._results_columns_names[0]].max():.2f}] °C\n"
        )
        summary += (
            f"    "
            f"Final results stored at: "
            f"{self._output_data_directory / 'insiliciopt.csv'}\n"
        )

        summary += "\n\n============= TERMINATED NORMALLY =============\n\n"

        self._logger.info(summary)

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
        reactor_conditions = ReactorConditions(
            temperature=conditions.temperature,
            concentrations=concentrations,
            products=self.species.products,
            time=conditions.time,
        )
        return reactor_conditions

    def _store_results(self) -> None:
        self.results.to_csv(
            self._output_data_directory / "insiliciopt.csv",
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
        self._store_results()
        self._visualization.plot(
            e=np.array(self.results[self._results_columns_names[-1]].values),
            sty=np.array(self.results[self._results_columns_names[-2]].values),
            iteration=self._counter,
        )

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

    def _end(self) -> None:
        self._visualization.animate()
        self._log_summary()


class SummitOptimizer(Optimizer):

    domain: Domain
    """Optimization Domain"""

    num_initial_points: int
    """Number of initial LHS points"""

    _columns: list[str]
    """Domain column Names"""

    def __init__(
            self,
            species: OptimizationSpecies,
            boundaries: OptimizationBoundaries,
            reactor: Reactor,
            output_directory: Path,
            num_initial_points: int = 20,
    ) -> None:
        self.num_initial_points = num_initial_points

        if self.num_initial_points < 4:
            raise ValueError(
                "The number of LHS points must be at least 3."
            )

        super().__init__(
            species,
            boundaries,
            reactor,
            output_directory,
            num_initial_points=num_initial_points,
        )

        self._create_domain()

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
            maximize=True,
            description="Space Time Yield Objective",
        )
        domain += ContinuousVariable(
            self._results_columns_names[5],
            is_objective=True,
            bounds=list(self.boundaries.E_factor),
            maximize=False,
            description="Mass product waste factor objective",
        )
        self._columns = [variable.name for variable in domain.variables]
        self.domain = domain

    def _construct_optimization_conditions(
            self, suggestion: DataSet
    ) -> _OptimizerConditions:
        """Convert summit Dataset to Optimizer Conditions"""
        data = suggestion.to_dict()["data"][0]
        temperature = data[
            round(self._columns.index(
                self._results_columns_names[0]
            ), 4)
        ]
        concentration = data[
            round(self._columns.index(
                self._results_columns_names[1]
            ), 4)
        ]
        ratio = data[
            round(self._columns.index(
                self._results_columns_names[2]
            ), 4)
        ]
        time = data[
            round(self._columns.index(
                self._results_columns_names[3]
            ), 4)
        ]

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
            self.num_initial_points,
        )

        for i in range(self.num_initial_points):
            suggestion = initial_suggestions.iloc[[i]]
            conditions = self._construct_optimization_conditions(suggestion)
            E, STY = self._simulate_reactor(conditions)
            self._add_to_result(E=E, STY=STY, conditions=conditions)
            self._log_iteration(
                E=E, STY=STY, conditions=conditions, iteration_type="LHS"
            )

    def _run_optimizer(self, num_iterations: int) -> None:
        """Optimized Summit sampling"""
        for i in range(num_iterations):
            dataset = DataSet.from_df(self.results)
            suggestion = self._summit_optimizer_suggest(dataset)
            conditions = self._construct_optimization_conditions(suggestion)
            E, STY = self._simulate_reactor(conditions)
            self._add_to_result(E=E, STY=STY, conditions=conditions)
            self._log_iteration(
                E=E, STY=STY, conditions=conditions, iteration_type="TSEMO"
            )

    @abstractmethod
    def _summit_optimizer_suggest(self, dataset: DataSet) -> DataSet:
        """Summit Optimizer Suggest"""

    def run(self, num_iterations: int) -> None:
        """Runs the Summit optimizer"""
        self._run_lhs()
        self._run_optimizer(num_iterations=num_iterations)
        self._end()

class TSEmoOptimizer(SummitOptimizer):

    _ts_emo_strategy: TSEMO
    """TSEmo Optimizer Strategy"""

    def __init__(
            self,
            species: OptimizationSpecies,
            boundaries: OptimizationBoundaries,
            reactor: Reactor,
            output_directory: Path,
            num_initial_points: int = 20,
    ) -> None:

        super().__init__(
            species,
            boundaries,
            reactor,
            output_directory,
            num_initial_points,
        )

        self._ts_emo_strategy = TSEMO(self.domain)

    def __repr__(self) -> str:
        string = "\n\nOPTIMIZER:\n"
        string += "  Type: TSEMO\n"
        string += f"  Num LHS points: {self.num_initial_points}\n"
        return string

    def _summit_optimizer_suggest(self, dataset: DataSet) -> DataSet:
        return self._ts_emo_strategy.suggest_experiments(
                1, prev_res=dataset,
            )

class SoBoOptimizer(SummitOptimizer):

    _so_bo_optimizer: SOBO
    """SoBo Optimizer Strategy"""

    def __init__(
            self,
            species: OptimizationSpecies,
            boundaries: OptimizationBoundaries,
            reactor: Reactor,
            output_directory: Path,
            num_initial_points: int = 20,
    ) -> None:

        super().__init__(
            species,
            boundaries,
            reactor,
            output_directory,
            num_initial_points,
        )

        self._so_bo_optimizer = SOBO(self.domain)

    def __repr__(self) -> str:
        string = "\n\nOPTIMIZER:\n"
        string += "  Type: SOBO\n"
        string += f"  Num LHS points: {self.num_initial_points}\n"
        return string

    def _summit_optimizer_suggest(self, dataset: DataSet) -> DataSet:
        return self._so_bo_optimizer.suggest_experiments(
                1, prev_res=dataset,
            )

class EntMootOptimizer(SummitOptimizer):

    _ent_moot_optimizer: ENTMOOT
    """EntMoot Optimizer Strategy"""

    def __init__(
            self,
            species: OptimizationSpecies,
            boundaries: OptimizationBoundaries,
            reactor: Reactor,
            output_directory: Path,
            num_initial_points: int = 20,
    ) -> None:
        super().__init__(
            species,
            boundaries,
            reactor,
            output_directory,
            num_initial_points,
        )

        self._ent_moot_optimizer = ENTMOOT(self.domain)

    def __repr__(self) -> str:
        string = "\n\nOPTIMIZER:\n"
        string += "  Type: ENTMOOT\n"
        string += f"  Num LHS points: {self.num_initial_points}\n"
        return string

    def _summit_optimizer_suggest(self, dataset: DataSet) -> DataSet:
        return self._ent_moot_optimizer.suggest_experiments(
            1, prev_res=dataset,
        )
