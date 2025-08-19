from pathlib import Path

from ispareto.species import (
    Species,
    TransitionState,
)
from ispareto.kinetics import Kinetics
from ispareto.reactor import Reactor
from ispareto.solvation import Solvation
from ispareto.optimizer import (
    TSEmoOptimizer,
    OptimizationSpecies,
    OptimizationBoundaries,
)
import ispareto.visualization as visualization

##############
# Paths
##############
current_dir = Path(__file__).parent
gaussian_dir = current_dir / "data" / "gaussian"
cosmo_therm_dir = current_dir / "data" / "cosmotherm"
output_dir = current_dir / "output"

##############
# Visualization
##############
visualization.e_factor_bounds = (0.0, 15.0)
visualization.sty_bounds = (0.0, 500)

##############
# Species
##############
hbr = Species(
    name="HBr",
    mass=0.080911,
    fchk_file_path=gaussian_dir / "HBr.fchk",
    tab_file_path=cosmo_therm_dir / "HBr.tab",
    energy=-2574.66039157451,
)

product_1 = Species(
    name="Product1",
    mass=0.12108915 + 0.16997311 - 0.080911,
    fchk_file_path=gaussian_dir / "Product1.fchk",
    tab_file_path=cosmo_therm_dir / "Product1.tab",
    energy=-636.581035187144,
)

product_2 = Species(
    name="Product2",
    mass=0.12108915 + 2*0.16997311 - 2*0.080911,
    fchk_file_path=gaussian_dir / "Product2.fchk",
    tab_file_path=cosmo_therm_dir / "Product2.tab",
    energy=-906.940141884542,
)

reactant_1 = Species(
    name="Reactant1",
    mass=0.12108915,
    fchk_file_path=gaussian_dir / "Reactant1.fchk",
    tab_file_path=cosmo_therm_dir / "Reactant1.tab",
    energy=-366.220519437833,
)

reactant_2 = Species(
    name="Reactant2",
    mass=0.16997311,
    fchk_file_path=gaussian_dir / "Reactant2.fchk",
    tab_file_path=cosmo_therm_dir / "Reactant2.tab",
    energy=-2845.0216955827,
)


##################
# Transition State
##################
ts_1_fwd = TransitionState(
    name="TS1_fwd",
    fchk_file_path=gaussian_dir / "TS1.fchk",
    tab_file_path=cosmo_therm_dir / "TS1.tab",
    energy=-3211.2069180058,
)

ts_1_rev = TransitionState(
    name="TS1_rev",
    fchk_file_path=gaussian_dir / "TS1.fchk",
    tab_file_path=cosmo_therm_dir / "TS1.tab",
    energy=-3211.2069180058,
)

ts_2_fwd = TransitionState(
    name="TS2_fwd",
    fchk_file_path=gaussian_dir / "TS2.fchk",
    tab_file_path=cosmo_therm_dir / "TS2.tab",
    energy=-3481.570039358,
)

ts_2_rev = TransitionState(
    name="TS2_rev",
    fchk_file_path=gaussian_dir / "TS2.fchk",
    tab_file_path=cosmo_therm_dir / "TS2.tab",
    energy=-3481.570039358,
)


##############
# Reactions
##############
# R1_fwd = Reactant 1 + Reactant 2 -> Product 1 + HBr
reaction_1_fwd = product_1 + hbr - reactant_1 - reactant_2
reaction_1_fwd.name = "R1_fwd"
reaction_1_fwd.transition_state = ts_1_fwd

# R2_fwd = Reactant 2 + Product 1 -> Product 2 + HBr
reaction_2_fwd = product_2 + hbr - reactant_2 - product_1
reaction_2_fwd.name = "R2_fwd"
reaction_2_fwd.transition_state = ts_2_fwd

reactions = [
    reaction_1_fwd,
    reaction_2_fwd,
]


############
# Kinetics
############
kinetics = Kinetics(
    reactions=reactions,
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
optimization_species = OptimizationSpecies(
    reactant_1=reactant_1,
    reactant_2=reactant_2,
    products=[product_1,]
)
optimization_boundaries = OptimizationBoundaries(
    temperature=(100, 150),
    concentration_reactant_1=(100, 500),
    concentration_ratio=(1.0, 5.0),
    time=(10.0, 20.0),
)
optimizer = TSEmoOptimizer(
    species=optimization_species,
    boundaries=optimization_boundaries,
    reactor=reactor,
    output_directory=output_dir,
    num_initial_points=5,
)
optimizer.run(num_iterations=110)
