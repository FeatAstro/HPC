"""Figure for Step 3a: error of the derivatives and averages between nodes and cell centres
(Yee layout): what the code gives (coloured) against what is expected (black, dashed).

Run after building:  python3 plots/plot_step3a_yee.py
"""
import math

import numpy as np

from style import BLUE, INK, ORANGE, plt, save
from common import cell_centre_positions, node_positions, pic

LENGTH = 2.0 * math.pi
WAVENUMBER = 2.0
CELLS = np.array([16, 32, 64, 128, 256, 512])
EXPECTED = dict(color=INK, lw=1, ls="--")


def wave(x):
    return np.sin(WAVENUMBER * x)


def slope(x):
    return WAVENUMBER * np.cos(WAVENUMBER * x)


def largest_errors(nx):
    grid = pic.Grid(nx, LENGTH / nx)
    nodes, cell_centres = node_positions(grid), cell_centre_positions(grid)
    return {
        "derivative at cell centres": np.max(np.abs(pic.derivative_at_cell_centres(wave(nodes), grid)
                                                    - slope(cell_centres))),
        "derivative at nodes": np.max(np.abs(pic.derivative_at_nodes(wave(cell_centres), grid) - slope(nodes))),
        "average at cell centres": np.max(np.abs(pic.average_at_cell_centres(wave(nodes)) - wave(cell_centres))),
        "average at nodes": np.max(np.abs(pic.average_at_nodes(wave(cell_centres)) - wave(nodes))),
    }


def plot_errors(ax, errors, spacings, operations, expected_error, expected_label, title):
    # the two directions give almost the same error: large markers under small ones keep both visible
    for operation, color, marker_size in zip(operations, (BLUE, ORANGE), (13, 6)):
        measured = np.array([error[operation] for error in errors])
        order = np.polyfit(np.log(spacings), np.log(measured), 1)[0]
        ax.loglog(spacings, measured, "o", color=color, ms=marker_size,
                  label=f"measured: {operation} (slope {order:.2f})")
    ax.loglog(spacings, expected_error, label=expected_label, **EXPECTED)
    ax.set_xlabel("grid spacing dx")
    ax.set_ylabel("largest error over the grid")
    ax.set_title(title)
    ax.legend(loc="upper left", fontsize=8)


def plot_yee():
    spacings = LENGTH / CELLS
    errors = [largest_errors(nx) for nx in CELLS]
    k = WAVENUMBER

    fig, (ax_derivative, ax_average) = plt.subplots(1, 2, figsize=(12, 5))
    plot_errors(ax_derivative, errors, spacings, ("derivative at cell centres", "derivative at nodes"),
                k**3 * spacings**2 / 24.0, "expected: k³ dx² / 24",
                "Derivative of sin(kx), k = 2\nexpected slope: 2")
    plot_errors(ax_average, errors, spacings, ("average at cell centres", "average at nodes"),
                k**2 * spacings**2 / 8.0, "expected: k² dx² / 8",
                "Average of sin(kx), k = 2\nexpected slope: 2")
    save(fig, "step3a_yee.png")


if __name__ == "__main__":
    plot_yee()
