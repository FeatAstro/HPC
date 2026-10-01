"""Figure for Step 3b: the electric field of Ohm's law for a smooth plasma, computed by the code
(coloured) against the exact one (black, dashed), and its error against dx.

Run after building:  python3 plots/plot_step3b_ohm.py
"""
import math

import numpy as np

from style import INK, SERIES, plt, save
from common import SmoothPlasma, cell_centre_positions, node_positions, pic, plasma_parameters

LENGTH = 2.0 * math.pi
PARAMETERS = plasma_parameters(electron_temperature=0.5, resistivity=0.1, hyper_resistivity=0.05)
COMPONENTS = ["Ex", "Ey", "Ez"]
CELLS = np.array([32, 64, 128, 256, 512, 1024])
EXPECTED = dict(color=INK, lw=1, ls="--")


def solved_components(plasma, grid):
    E = plasma.solved_E(grid)
    return [E.x, E.y, E.z]


def plot_ohm():
    plasma = SmoothPlasma(LENGTH, PARAMETERS)
    fig, axes = plt.subplots(1, 4, figsize=(18, 4.8))

    coarse_grid = pic.Grid(24, LENGTH / 24)
    points = [cell_centre_positions(coarse_grid), node_positions(coarse_grid), node_positions(coarse_grid)]
    where = ["cell centres", "nodes", "nodes"]
    fine_x = np.linspace(0.0, LENGTH, 400)
    solved = solved_components(plasma, coarse_grid)
    exact = plasma.exact_E_on_the_yee_grid(coarse_grid)
    for k, name in enumerate(COMPONENTS):
        ax = axes[k]
        ax.plot(points[k], solved[k], "o", color=SERIES[k], ms=7, label=f"measured: on the {where[k]}, 24 cells")
        ax.plot(fine_x, plasma.E(fine_x)[k], label="expected: exact Ohm's law", **EXPECTED)
        ax.set_xlabel("x")
        ax.set_ylabel(name)
        ax.set_title(f"{name}\nlargest |measured − expected| = {np.max(np.abs(solved[k] - exact[k])):.3f}")
        ax.legend(loc="upper center", bbox_to_anchor=(0.5, -0.16), fontsize=8)

    ax = axes[3]
    spacings = LENGTH / CELLS
    slopes = []
    for k, name in enumerate(COMPONENTS):
        errors = []
        for nx in CELLS:
            grid = pic.Grid(nx, LENGTH / nx)
            errors.append(np.max(np.abs(solved_components(plasma, grid)[k] - plasma.exact_E_on_the_yee_grid(grid)[k])))
        slopes.append(np.polyfit(np.log(spacings), np.log(errors), 1)[0])
        ax.loglog(spacings, errors, "o-", color=SERIES[k], ms=7, lw=3, label=f"measured: {name} (slope {slopes[k]:.2f})")
        if k == 0:
            first_error = errors[0]
    ax.loglog(spacings, 0.4 * first_error * (spacings / spacings[0]) ** 2, label="expected: ∝ dx²", **EXPECTED)
    ax.set_xlabel("grid spacing dx")
    ax.set_ylabel("largest error over the grid")
    ax.set_title("error against the grid spacing\nexpected slope: 2")
    ax.legend(loc="upper center", bbox_to_anchor=(0.5, -0.16), fontsize=8)

    fig.suptitle("Ohm's law for a smooth plasma: flow, Hall, electron pressure, resistivity and hyper-resistivity "
                 "all non-zero", color=INK)
    save(fig, "step3b_ohm.png")


if __name__ == "__main__":
    plot_ohm()
