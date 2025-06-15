import logging
from tamkin import (
    PartFun,
    NMA,
    load_molecule_g03fchk,
    ConstrainExt,
    ExtTrans,
    ExtRot,
    KineticModel,
    Wigner,
    Eckart,
    Miller
)

from insiliciopt.species import (
    Reaction,
    Species, TransitionState,
)


class Kinetics:

    reactions: list[Reaction]
    """Contains the Reactions present in the System"""

    tunneling_correction: str | None
    """Includes Tunneling into the Kinetics calculation, types are:
    wigner, eckart or miller
    """

    gradient_threshold: float
    """Gradient Threshold for Tamkin"""

    _kinetics_models: dict[Reaction, KineticModel]
    """Holds the constructed Kinetic Model for every Reaction"""

    def __init__(
            self,
            reactions: list[Reaction],
            tunneling_correction: str | None = None,
            gradient_threshold: float = 4e-3,
    ):
        self._logger: logging.Logger = logging.getLogger(__class__.__name__)

        self.reactions = reactions
        self.tunneling_correction = tunneling_correction
        self.gradient_threshold = gradient_threshold

        self._check_input()
        self._construct_kinetics_model()

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

    def _construct_tunneling_correction(
            self,
            transition_state_partition_function: PartFun,
            reactant_partition_functions: list[PartFun],
            product_partition_functions: list[PartFun],
    ) -> Wigner | Eckart | Miller | None:
        """Constructs the Tunneling correction"""
        tunneling = None

        if self.tunneling_correction not in [
            "wigner", "eckart", "miller", None,
        ]:
            raise ValueError(
                "Tunneling correction must be one of 'wigner', 'eckart', "
                "miller."
            )

        if self.tunneling_correction == "eckart":
            tunneling = Eckart(
                reactant_partition_functions,
                transition_state_partition_function,
                product_partition_functions
            )

        if self.tunneling_correction == "wigner":
            tunneling = Wigner(transition_state_partition_function)

        if self.tunneling_correction == "miller":
            tunneling = Miller(transition_state_partition_function)

        return tunneling

    def _construct_partition_function(self, species: Species) -> PartFun:
        """Constructs the partition function for one species"""
        partition_function = PartFun(
            NMA(
                load_molecule_g03fchk(
                    species.fchk_file_path,
                    energy=species.energy
                ),
                ConstrainExt(
                    gradient_threshold=self.gradient_threshold
                )
            ),
            [ExtTrans(), ExtRot()]
        )
        return partition_function


    def _construct_kinetics_model(self):
        """Constructs the Kinetics Model which can be evaluated on the fly"""

        for reaction in self.reactions:

            reactants_pfs: list[PartFun] = []
            product_pfs: list[PartFun] = []

            reactants = reaction.reactants
            products = reaction.products
            transition_state = reaction.transition_state

            for reactant in reactants:
                reactants_pfs.append(
                    self._construct_partition_function(reactant)
                )

            # Product Partition Functions are only needed for eckart
            # Tunneling correction, otherwise -> unnecessary construction
            if self.tunneling_correction == "eckart":
                for product in products:
                    product_pfs.append(
                        self._construct_partition_function(product)
                    )

            transition_state_pf = self._construct_partition_function(
                transition_state
            )

            kinetic_model = KineticModel(
                pfs_react=reactants_pfs,
                pf_trans=transition_state_pf,
                tunneling=self._construct_tunneling_correction(
                    transition_state_pf,
                    reactants_pfs,
                    product_pfs,
                )
            )
            self._kinetics_models[reaction] = kinetic_model


    def k(self, reaction: Reaction, temperature: float) -> float:
        """Returns the rate constant for given species and temperature in SI
        units
        """
        kinetic_model = self._kinetics_models[reaction]
        k = kinetic_model.rate_constant(temperature)
        k_SI = k / kinetic_model.unit
        return k_SI
