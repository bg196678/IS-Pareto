import os
import imageio
import tempfile
import numpy as np
from pathlib import Path
import matplotlib.pyplot as plt


e_factor_bounds: tuple[float, float] = (0, 3.0)
"""Plotting bounds for the E-Factor axis"""

sty_bounds: tuple[float, float] = (0.0, 13000)
"""Plotting bounds for the STY axis"""

plot_style = r"""
# Matplotlib style for scientific plotting
# This is the base style for "SciencePlots"
# see: https://github.com/garrettj403/SciencePlots

# Set color cycle: blue, green, yellow, red, violet, gray
axes.prop_cycle : cycler('color', ['0C5DA5', '00B945', 'FF9500', 'FF2C00', '845B97', '474747', '9e9e9e'])

# Set default figure size
figure.figsize : 3.5, 2.625

# Set x axis
xtick.direction : in
xtick.major.size : 3
xtick.major.width : 0.5
xtick.minor.size : 1.5
xtick.minor.width : 0.5
xtick.minor.visible : True
xtick.top : True

# Set y axis
ytick.direction : in
ytick.major.size : 3
ytick.major.width : 0.5
ytick.minor.size : 1.5
ytick.minor.width : 0.5
ytick.minor.visible : True
ytick.right : True

# Set line widths
axes.linewidth : 0.5
grid.linewidth : 0.5
lines.linewidth : 1.

# Remove legend frame
legend.frameon : False

# Always save as 'tight'
savefig.bbox : tight
savefig.pad_inches : 0.05

# Use serif fonts
# font.serif : Times
font.family : serif
mathtext.fontset : dejavuserif

# Use LaTeX for math formatting
#text.usetex : True
#text.latex.preamble : \usepackage{amsmath} \usepackage{amssymb}
"""

style_file = tempfile.NamedTemporaryFile(
    mode='w', suffix='.mplstyle', delete=False
)
style_file.write(plot_style)
style_file.close()
plt.style.use(style_file.name)
os.unlink(style_file.name)

class Visualization:
    """Visualization of the Development of the Pareto Front with
    iterations
    """

    plot_directory: Path
    """Path to the Plotting directory"""

    num_initial_points: int | None = None
    """Initial Points"""

    def __init__(
            self,
            plot_directory: Path,
            num_initial_points: int | None = None,
    ) -> None:
        self.plot_directory = plot_directory
        self.num_initial_points = num_initial_points

    def _pareto_front(
            self, e: np.ndarray, sty: np.ndarray
    ) -> tuple[np.ndarray, np.ndarray]:
        """Calculate the Pareto front from arrays of E-factor and STY values.
        """
        n = len(e)
        is_pareto = np.ones(n, dtype=bool)

        for i in range(n):
            for j in range(n):
                if i == j:
                    continue
                if (sty[j] >= sty[i] and e[j] <= e[i]) and (
                        sty[j] > sty[i] or e[j] < e[i]):
                    is_pareto[i] = False
                    break

        pareto_indices = np.where(is_pareto)[0]
        sorted_indices = pareto_indices[np.argsort(e[pareto_indices])]
        return e[sorted_indices], sty[sorted_indices]

    def plot(
            self,
            e: np.ndarray,
            sty: np.ndarray,
            iteration: int,
    ) -> None:

        e_pareto, sty_pareto = self._pareto_front(e, sty)

        plt.figure(figsize=(8, 6))

        if self.num_initial_points:
            plt.scatter(
                sty[:self.num_initial_points], e[:self.num_initial_points],
                marker='s', color='black', label=f'Initial Points',
                s=64, zorder=4
            )

        plt.scatter(
            sty[self.num_initial_points:], e[self.num_initial_points:],
            marker='x', color='blue', label=f'Optimized Points',
            alpha=0.7, zorder=3
        )

        plt.scatter(
            sty_pareto, e_pareto,
            facecolors='orange', edgecolors='red', s=100, linewidths=1.5,
            label='Pareto Front', zorder=5, marker='o'
        )

        plt.plot(
            sty_pareto, e_pareto,
            color='red', linewidth=2, linestyle='--', zorder=2,
            label='Pareto Curve'
        )

        plt.xlim(*sty_bounds)
        plt.ylim(*e_factor_bounds)
        plt.xlabel('STY [kg/m³/h]', fontsize=18)
        plt.ylabel('E-Factor', fontsize=18)
        plt.title(
            f'In-silicio closed loop multi-objective optimization - '
            f'Iteration {iteration}',
            fontsize=16
        )
        plt.tick_params(axis='both', which='major', labelsize=15)
        plt.legend(fontsize=15)
        plt.grid(True, linestyle='--', alpha=0.5)
        plt.tight_layout()

        file_path = self.plot_directory / f"pareto_front_iteration_{iteration}.png"
        plt.savefig(str(file_path), dpi=300)
        plt.close()

    def animate(self) -> None:
        frame_paths = []

        for file in self.plot_directory.iterdir():
            if not ".png" in file.name:
                continue

            frame_paths.append(file)

        output_path = self.plot_directory / "pareto_front_animation.gif"

        with imageio.get_writer(
                str(output_path),
                mode='I',
                duration=0.4,
                loop=0
        ) as writer:
            for frame_path in frame_paths:
                image = imageio.imread(frame_path)
                writer.append_data(image)
