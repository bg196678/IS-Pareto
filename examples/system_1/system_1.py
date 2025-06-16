from pathlib import Path

from insiliciopt.species import (
    Species,
    TransitionState,
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
cosmo_therm_dir = current_dir / "data" / "cosmotherm"
output_dir = current_dir / "output"

##############
# Species
##############
substrate = Species(
    name="Substrate",
    mass=0.159,
    fchk_file_path=gaussian_dir / "Substrate.fchk",
    tab_file_path=cosmo_therm_dir / "Substrate.tab",
    energy=-635.185999076,
)

nucleophilic = Species(
    name="Nucleophilic",
    mass=0.071,
    fchk_file_path=gaussian_dir / "Nucleophilic.fchk",
    tab_file_path=cosmo_therm_dir / "Nucleophilic.tab",
    energy=-212.576801654,
)

its_1 = Species(
    name="ITS1",
    mass=0.159 + 0.071,
    fchk_file_path=gaussian_dir / "ITS1.fchk",
    tab_file_path=cosmo_therm_dir / "ITS1.tab",
    energy=-847.760037351,
)

its_2 = Species(
    name="ITS2",
    mass=0.159 + 0.071,
    fchk_file_path=gaussian_dir / "ITS2.fchk",
    tab_file_path=cosmo_therm_dir / "ITS2.tab",
    energy=-847.734323929,
)

its_3 = Species(
    name="ITS3",
    mass=0.159 + 0.071 + 0.071,
    fchk_file_path=gaussian_dir / "ITS3.fchk",
    # TODO
    energy=-959.9276068,
)

its_4 = Species(
    name="ITS4",
    mass=0.159 + 0.071 + 0.071,
    fchk_file_path=gaussian_dir / "ITS4.fchk",
    # TODO
    energy=-959.881772303,
)

product_1 = Species(
    name="Product1",
    mass=0.21008046,
    fchk_file_path=gaussian_dir / "Product1.fchk",
    tab_file_path=cosmo_therm_dir / "Product1.tab",
    energy=-747.339396317,
)

product_2 = Species(
    name="Product2",
    mass=0.21008046,
    fchk_file_path=gaussian_dir / "Product2.fchk",
    tab_file_path=cosmo_therm_dir / "Product2.tab",
    energy=-747.339639927,
)

product_3 = Species(
    name="Product3",
    mass=0.26114773,
    fchk_file_path=gaussian_dir / "Product3.fchk",
    tab_file_path=cosmo_therm_dir / "Product3.tab",
    energy=-859.487278239,
)

leaving_group = Species(
    name="LeavingGroup",
    mass=0.020,
    fchk_file_path=gaussian_dir / "Leaving_Group.fchk",
    tab_file_path=cosmo_therm_dir / "Leaving_Group.tab",
    energy=-100.453158015,
)


##################
# Transition State
##################
ts_1_fwd = TransitionState(
    name="TS1_fwd",
    fchk_file_path=gaussian_dir / "TS1_fwd.fchk",
    tab_file_path=cosmo_therm_dir / "TS1.tab",
    energy=-847.759104357,
)

ts_1_rev = TransitionState(
    name="TS1_rev",
    fchk_file_path=gaussian_dir / "TS1_rev.fchk",
    tab_file_path=cosmo_therm_dir / "TS1.tab",
    energy=-847.759104357,
)

ts_2_fwd = TransitionState(
    name="TS2_fwd",
    fchk_file_path=gaussian_dir / "TS2_fwd.fchk",
    tab_file_path=cosmo_therm_dir / "TS2.tab",
    energy=-847.748732583,
)

ts_2_rev = TransitionState(
    name="TS2_rev",
    fchk_file_path=gaussian_dir / "TS2_rev.fchk",
    tab_file_path=cosmo_therm_dir / "TS2.tab",
    energy=-847.748732583,
)

ts_3_fwd = TransitionState(
    name="TS3_fwd",
    fchk_file_path=gaussian_dir / "TS3_fwd.fchk",
    tab_file_path=cosmo_therm_dir / "TS3.tab",
    energy=-959.907069926,
)

ts_3_rev = TransitionState(
    name="TS3_rev",
    fchk_file_path=gaussian_dir / "TS3_rev.fchk",
    tab_file_path=cosmo_therm_dir / "TS3.tab",
    energy=-959.907069926,
)

ts_4_fwd = TransitionState(
    name="TS4_fwd",
    fchk_file_path=gaussian_dir / "TS4_fwd.fchk",
    tab_file_path=cosmo_therm_dir / "TS4.tab",
    energy=-959.894825835,
)

ts_4_rev = TransitionState(
    name="TS4_rev",
    fchk_file_path=gaussian_dir / "TS4_rev.fchk",
    tab_file_path=cosmo_therm_dir / "TS4.tab",
    energy=-959.894825835,
)

ts_1_2_fwd = TransitionState(
    name="TS12_fwd",
    fchk_file_path=gaussian_dir / "TS12_fwd.fchk",
    # TODO
    energy=-847.749041202,
)

ts_1_2_rev = TransitionState(
    name="TS12_rev",
    fchk_file_path=gaussian_dir / "TS12_rev.fchk",
    # TODO
    energy=-847.749041202,
)

ts_2_2_fwd = TransitionState(
    name="TS22_fwd",
    fchk_file_path=gaussian_dir / "TS22_fwd.fchk",
    energy=-847.729373816,
    # TODO
)

ts_2_2_rev = TransitionState(
    name="TS22_rev",
    fchk_file_path=gaussian_dir / "TS22_rev.fchk",
    # TODO
    energy=-847.729373816,
)

ts_3_2_fwd = TransitionState(
    name="TS32_fwd",
    fchk_file_path=gaussian_dir / "TS32_fwd.fchk",
    tab_file_path=cosmo_therm_dir / "TS32.tab",
    energy=-959.900414613,
)

ts_3_2_rev = TransitionState(
    name="TS32_2rev",
    fchk_file_path=gaussian_dir / "TS32_rev.fchk",
    tab_file_path=cosmo_therm_dir / "TS32.tab",
    energy=-959.900414613,
)

ts_4_2_fwd = TransitionState(
    name="TS42_fwd",
    fchk_file_path=gaussian_dir / "TS42_fwd.fchk",
    tab_file_path=cosmo_therm_dir / "TS42.tab",
    energy=-959.874575137,
)

ts_4_2_rev = TransitionState(
    name="TS42_rev",
    fchk_file_path=gaussian_dir / "TS42_rev.fchk",
    tab_file_path=cosmo_therm_dir / "TS42.tab",
    energy=-959.874575137,
)


##############
# Reactions
##############
# R1_fwd: Substrate + Nucleophilic -> ITS1
reaction_1_fwd = its_1 - substrate - nucleophilic
reaction_1_fwd.name = "R1_fwd"
reaction_1_fwd.transition_state = ts_1_fwd

# R1_rev: ITS1 -> Substrate + Nucleophilic
reaction_1_rev = substrate + nucleophilic - its_1
reaction_1_rev.name = "R1_rev"
reaction_1_rev.transition_state = ts_1_rev

# R2_fwd: Substrate + Nucleophilic -> ITS2
reaction_2_fwd = its_2 - substrate - nucleophilic
reaction_2_fwd.name = "R2_fwd"
reaction_2_fwd.transition_state = ts_2_fwd

# R2_rev: ITS2 -> Substrate + Nucleophilic
reaction_2_rev = substrate + nucleophilic - its_2
reaction_2_rev.name = "R2_rev"
reaction_2_rev.transition_state = ts_2_rev

# R3_fwd: Product2 + Nucleophilic -> ITS3
reaction_3_fwd = its_3 - product_2 - nucleophilic
reaction_3_fwd.name = "R3_fwd"
reaction_3_fwd.transition_state = ts_3_fwd

# R3_rev: ITS3 -> Product2 + Nucleophilic
reaction_3_rev = product_2 + nucleophilic - its_3
reaction_3_rev.name = "R3_rev"
reaction_3_rev.transition_state = ts_3_rev

# R4_fwd: Product1 + Nucleophilic -> ITS4
reaction_4_fwd = its_4 - product_1 - nucleophilic
reaction_4_fwd.name = "R4_fwd"
reaction_4_fwd.transition_state = ts_4_fwd

# R4_rev: ITS4 -> Product1 + Nucleophilic
reaction_4_rev = product_1 + nucleophilic - its_4
reaction_4_rev.name = "R4_rev"
reaction_4_rev.transition_state = ts_4_rev

# R5_fwd: ITS1 -> Product1 + LeavingGroup
reaction_5_fwd = product_1 + leaving_group - its_1
reaction_5_fwd.name = "R5_fwd"
reaction_5_fwd.transition_state = ts_1_2_fwd

# R5_rev: Product1 + LeavingGroup -> ITS1
reaction_5_rev = its_1 - product_1 - leaving_group
reaction_5_rev.name = "R5_rev"
reaction_5_rev.transition_state = ts_1_2_rev

# R6_fwd: ITS2 -> Product2 + LeavingGroup
reaction_6_fwd = product_2 + leaving_group - its_2
reaction_6_fwd.name = "R6_fwd"
reaction_6_fwd.transition_state = ts_2_2_fwd

# R6_rev: Product2 + LeavingGroup -> ITS2
reaction_6_rev = its_2 - product_2 - leaving_group
reaction_6_rev.name = "R6_rev"
reaction_6_rev.transition_state = ts_2_2_rev

# R7_fwd: ITS3 -> Product3 + LeavingGroup
reaction_7_fwd = product_3 + leaving_group - its_3
reaction_7_fwd.name = "R7_fwd"
reaction_7_fwd.transition_state = ts_3_2_fwd

# R7_rev: Product3 + LeavingGroup -> ITS3
reaction_7_rev = its_3 - product_3 - leaving_group
reaction_7_rev.name = "R7_rev"
reaction_7_rev.transition_state = ts_3_2_rev

# R8_fwd: ITS4 -> Product3 + LeavingGroup
reaction_8_fwd = product_3 + leaving_group - its_4
reaction_8_fwd.name = "R8_fwd"
reaction_8_fwd.transition_state = ts_4_2_fwd

# R8_rev: Product3 + LeavingGroup -> ITS4
reaction_8_rev = its_4 - product_3 - leaving_group
reaction_8_rev.name = "R8_rev"
reaction_8_rev.transition_state = ts_4_2_rev


reactions = [
    reaction_1_fwd, reaction_1_rev,
    reaction_2_fwd, reaction_2_rev,
    reaction_3_fwd, reaction_3_rev,
    reaction_4_fwd, reaction_4_rev,
    reaction_5_fwd, reaction_5_rev,
    reaction_6_fwd, reaction_6_rev,
    reaction_7_fwd, reaction_7_rev,
    reaction_8_fwd, reaction_8_rev
]


############
# Kinetics
############
kinetics = Kinetics(
    reactions=reactions,
    tunneling_correction="wigner",
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
    species=optimization_species,
    boundaries=optimization_boundaries,
    reactor=reactor,
    output_directory=output_dir,
    num_lhs_points=20,
)
optimizer.run(num_iterations=100)
