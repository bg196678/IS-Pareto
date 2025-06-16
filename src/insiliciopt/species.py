import logging
from pathlib import Path

class Species:
    """Abstract class for species"""

    name: str
    """Name Identifier for the Species"""

    mass: float | None
    """Mass of the Species"""

    energy: float | None
    """Energy (high fidelity) correction of the Species (not necessary)"""

    fchk_file_path: Path
    """Gaussian Fchk File Path"""

    tab_file_path: Path | None
    """Cosmotherm Tab File Path"""

    def __init__(
            self,
            name: str,
            mass: float | None,
            fchk_file_path: Path,
            tab_file_path: Path | None = None,
            energy: float | None = None
    ):
        self._logger: logging.Logger = logging.getLogger(
            f"Species('{name}')"
        )

        self.name = name
        self.mass = mass
        self.energy = energy
        self.tab_file_path = tab_file_path
        self.fchk_file_path = fchk_file_path

        if not self.fchk_file_path.exists():
            raise FileNotFoundError(
                f"FCHK File not found: {self.fchk_file_path}"
            )

        if not self.tab_file_path.exists():
            raise FileNotFoundError(
                f"TAB File not found: {self.tab_file_path}"
            )

    def __eq__(self, other):
        if not isinstance(other, Species):
            return False
        return (
            self.name == other.name and
            self.mass == other.mass and
            self.fchk_file_path == other.fchk_file_path
        )

    def __add__(self, other):
        if isinstance(other, Species):
            return Reaction(
                stoichiometry={
                    self: 1,
                    other: 1
                }
            )
        elif isinstance(other, ReactionTerm):
            return ReactionTerm(species=self, coefficient=1) + other
        elif isinstance(other, Reaction):
            result = other.stoichiometry.copy()
            result[self] = result.get(self, 0) + 1
            return Reaction(stoichiometry=result)
        else:
            raise TypeError(f"Cannot add Species to {type(other)}")

    def __sub__(self, other):
        if isinstance(other, Species):
            return Reaction(
                stoichiometry={
                    self: 1,
                    other: -1
                }
            )
        elif isinstance(other, ReactionTerm):
            return ReactionTerm(species=self, coefficient=1) - other
        elif isinstance(other, Reaction):
            result = {k: -v for k, v in other.stoichiometry.items()}
            result[self] = result.get(self, 0) + 1
            return Reaction(stoichiometry=result)
        else:
            raise TypeError(f"Cannot subtract {type(other)} from Species")

    def __hash__(self):
        return hash((self.name, self.mass, str(self.fchk_file_path)))

    def __repr__(self) -> str:
        return f"Species('{self.name}')"


class Reactant(Species):
    """Reactants"""


class Product(Species):
    """Products"""


class TransitionState(Species):
    """Transition states"""

    def __init__(
            self,
            name: str,
            fchk_file_path: Path,
            tab_file_path: Path,
            energy: float | None = None,
    ) -> None:
        super().__init__(
            name=name,
            mass=None, # Transition States do not need mass
            tab_file_path=tab_file_path,
            fchk_file_path=fchk_file_path,
            energy=energy,
        )


class ReactionTerm:
    """Represents a term in a reaction with a species and its stoichiometric
    coefficient
    """

    def __init__(self, species, coefficient):
        self.species = species
        self.coefficient = coefficient

    def __add__(self, other):
        if isinstance(other, ReactionTerm):
            return Reaction(
                stoichiometry=
                {
                    self.species: self.coefficient,
                    other.species: other.coefficient
                }
            )

        elif isinstance(other, Reaction):
            result = other.stoichiometry.copy()
            result[self.species] = result.get(
                self.species, 0
            ) + self.coefficient
            return Reaction(stoichiometry=result)
        else:
            raise TypeError(f"Cannot add ReactionTerm to {type(other)}")

    def __sub__(self, other):
        if isinstance(other, ReactionTerm):
            return Reaction(
                stoichiometry=
                {
                    self.species: self.coefficient,
                    other.species: -other.coefficient
                }
            )

        elif isinstance(other, Reaction):
            result = {k: -v for k, v in other.stoichiometry.items()}
            result[self.species] = result.get(
                self.species,  0
            ) + self.coefficient
            return Reaction(stoichiometry=result)
        else:
            raise TypeError(f"Cannot subtract {type(other)} from ReactionTerm")


class Reaction:
    """Reaction"""

    name: str | None = None
    """Name Identifier for the Reaction"""

    transition_state: TransitionState | None = None

    def __init__(
            self,
            name: str | None = None,
            stoichiometry: dict[Species, int] | None = None,
            transition_state: TransitionState | None = None,
    ):
        """
        Initialize a Reaction with stoichiometry

        stoichiometry: Dictionary mapping species to stoichiometric
        coefficients
        """
        self._logger: logging.Logger = logging.getLogger(
            name or __class__.__name__
        )
        self.name = name
        self.transition_state = transition_state
        self.stoichiometry = stoichiometry or {}
        self.stoichiometry = {
            k: v for k, v in self.stoichiometry.items() if v != 0
        } # no need for 0 coefficients

    @property
    def species(self) -> list[Species]:
        return [
            sp for sp in self.stoichiometry.keys()
        ]

    @property
    def reactants(self) -> list[Species]:
        return [
            sp for sp, coeff in self.stoichiometry.items()
            if coeff < 0
        ]

    @property
    def products(self) -> list[Species]:
        return [
            sp for sp, coeff in self.stoichiometry.items()
            if coeff > 0
        ]

    def coefficient(self, species: Species) -> float:
        return self.stoichiometry.get(species, 0)

    def __add__(self, other):
        if isinstance(other, ReactionTerm):
            result = self.stoichiometry.copy()
            result[other.species] = result.get(
                other.species, 0
            ) + other.coefficient
            return Reaction(stoichiometry=result)

        elif isinstance(other, Reaction):
            result = self.stoichiometry.copy()
            for species, coef in other.stoichiometry.items():
                result[species] = result.get(species, 0) + coef
            return Reaction(stoichiometry=result)

        elif isinstance(other, Species):
            result = self.stoichiometry.copy()
            result[other] = result.get(other, 0) + 1
            return Reaction(stoichiometry=result)

        else:
            raise TypeError(f"Cannot add Reaction to {type(other)}")

    def __sub__(self, other):

        if isinstance(other, ReactionTerm):
            result = self.stoichiometry.copy()
            result[other.species] = result.get(
                other.species, 0
            ) - other.coefficient
            return Reaction(stoichiometry=result)

        elif isinstance(other, Reaction):
            result = self.stoichiometry.copy()
            for species, coef in other.stoichiometry.items():
                result[species] = result.get(species, 0) - coef
            return Reaction(stoichiometry=result)

        elif isinstance(other, Species):
            result = self.stoichiometry.copy()
            result[other] = result.get(other, 0) - 1
            return Reaction(stoichiometry=result)

        else:
            raise TypeError(f"Cannot subtract {type(other)} from Reaction")

    def __repr__(self):
        return f"Reaction({self.stoichiometry})"
