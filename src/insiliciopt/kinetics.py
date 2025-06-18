import logging
import numpy as np
import pandas as pd
from pathlib import Path

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
    Species,
)
from insiliciopt.utils import (
    ReactionInput,
    DataOutput,
)


class Kinetics(ReactionInput, DataOutput):

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
        super().__init__(reactions)
        self._logger: logging.Logger = logging.getLogger(__class__.__name__)

        self.tunneling_correction = tunneling_correction
        self.gradient_threshold = gradient_threshold

        self._kinetics_models = {}

        self._check_input()
        self._construct_kinetics_model()

    def __repr__(self) -> str:
        string = "\n\nKINETICS:\n"
        string += f"  Num Reactions: {len(self.reactions)}\n"
        string += f"  Num Species: {len(self.species)}\n"
        string += f"  Num Transition States: {len(self.transition_states)}\n"
        string += f"  Tunneling correction: {self.tunneling_correction}\n"
        string += f"  Gradient threshold: {self.gradient_threshold}\n"
        return string

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
                "Tunneling correction must be one of 'wigner', 'eckart' or"
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

        self._logger.debug(
            f"Tunneling Setup with {tunneling}."
        )
        return tunneling

    def _construct_partition_function(self, species: Species) -> PartFun:
        """Constructs the partition function for one species"""
        partition_function = PartFun(
            NMA(
                load_molecule_g03fchk(
                    species.fchk_file_path,
                    #energy=species.energy
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

            self._logger.debug(
                f"Constructed Kinetics Model for {Reaction}."
            )

    def k(self, reaction: Reaction, temperature: float) -> float:
        """Returns the rate constant for given species and temperature in SI
        units
        """
        kinetic_model = self._kinetics_models[reaction]
        k = kinetic_model.rate_constant(temperature)
        k_si = k / kinetic_model.unit
        return k_si

    def dump(self, path: Path) -> None:
        """Dumps the kinetics"""
        temperature_range = np.linspace(50, 2500, 2450)
        kinetics = {}
        for reaction in self.reactions:
            rates = []
            for temp in temperature_range:
                k = self.k(reaction, temp)
                rates.append(k)
            kinetics[reaction.name] = {
                'temperature': temperature_range.tolist(),
                'rate_constants': rates
            }
        df_list = []
        for reaction_name, data in kinetics.items():
            temp_df = pd.DataFrame({
                'temperature': data['temperature'],
                f'rate_{reaction_name}': data['rate_constants']
            })
            df_list.append(temp_df)

        result_df = df_list[0]
        for df in df_list[1:]:
            result_df = pd.merge(
                result_df, df, on='temperature', how='outer'
            )

        result_df.to_csv(path / "kinetics.csv", index=False)
        self._logger.debug(f"Kinetics data dumped to {path}")


