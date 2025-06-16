from pathlib import Path
from insiliciopt.species import Species, TransitionState

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
    energy=-847.734323929
)

# TODO