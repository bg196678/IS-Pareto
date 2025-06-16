from pathlib import Path

from insiliciopt.species import (
    Species,
    TransitionState,
    Reaction,
)
from insiliciopt.kinetics import Kinetics
from insiliciopt.reactor import Reactor
from insiliciopt.solvation import Solvation
from insiliciopt.optimizer import (
    TSEmoOptimizer,
    OptimizationSpecies,
    OptimizationBoundaries,
)

##############
# Paths
##############
current_dir = Path(__file__).parent
gaussian_dir = current_dir / "data" / "gaussian"
cosmotherm_dir = current_dir / "data" / "cosmotherm"
output_dir = current_dir / "output"

##############
# Species
##############

substrate = Species(
    name="Substrate",
    mass=0.159,
    fchk_file_path=gaussian_dir / "substrate.fchk",
    energy=-635.185999076,
)

nucleophilic = Species(
    name="Nucleophilic",
    mass=0.071,
    fchk_file_path=gaussian_dir / "nucleophilic.fchk",
    energy=-212.576801654,
)

ITS1 = Species(
    name="ITS1",
    mass=0.159 + 0.071,
    fchk_file_path=gaussian_dir / "ITS1.fchk",
    energy=-847.760037351,
)

ITS2 = Species(
    name="ITS2",
    mass=0.159 + 0.071,
    fchk_file_path=gaussian_dir / "ITS2.fchk",
    energy=-847.734323929,
)

# TODO

product_1 = None
product_2 = None
product_3 = None


##################
# Transition State
##################

TS1_fwd = TransitionState(
    name="TS1_fwd",
    fchk_file_path=gaussian_dir / "TS1_fwd.fchk",
    energy=-847.759104357,
)

# TODO


##############
# Reactions
##############

reaction_1_fwd = Reaction(
    name="R1_fwd",
    transition_state=TS1_fwd,
)
 # TODO

reactions = [
    reaction_1_fwd,
]

############
# Kinetics
############

kinetics = Kinetics(
    reactions=reactions,
    tunneling_correction="eckart",
)

############
# Solvation
############

solvation = Solvation(
    reactions=reactions,
)

############
# Reactor
############

reactor = Reactor(
    reactions=reactions,
    kinetics=kinetics,
    solvation=solvation,
)


############
# Optimizer
############
optimzation_species = OptimizationSpecies(
    reactant_1=substrate,
    reactant_2=nucleophilic,
    products=[
        product_1,
        product_2,
        product_3,
    ]
)
optimization_boundaries = OptimizationBoundaries(
    temperature=(60, 140),
    concentration_reactant_1=(100, 500),
    concentration_ratio=(10, 5.0),
    time=(0.5, 2.0),
)
optimizer = TSEmoOptimizer(
    species=optimzation_species,
    boundaries=optimization_boundaries,
    reactor=reactor,
    output_directory=output_dir,
    num_lhs_points=20,
)
