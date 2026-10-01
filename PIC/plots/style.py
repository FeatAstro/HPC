"""Common setup of the plot scripts: imports `pic` and the test helpers, colours, figure folder."""
import os
import sys

import matplotlib

matplotlib.use("Agg")  # write files, no window
import matplotlib.pyplot as plt  # noqa: E402

HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, os.path.join(HERE, "..", "tests"))

FIGURES = os.path.join(HERE, "figures")
os.makedirs(FIGURES, exist_ok=True)

BLUE, ORANGE, AQUA = "#2a78d6", "#eb6834", "#1baf7a"
SERIES = [BLUE, ORANGE, AQUA, "#eda100", "#e87ba4", "#008300"]  # one colour per curve, in this order
INK, INK_SECONDARY, GRID = "#0b0b0b", "#52514e", "#e4e3df"

plt.rcParams.update({
    "figure.dpi": 150,
    "font.size": 10,
    "axes.edgecolor": INK_SECONDARY,
    "axes.labelcolor": INK,
    "axes.titlesize": 11,
    "axes.grid": True,
    "grid.color": GRID,
    "grid.linewidth": 0.6,
    "xtick.color": INK_SECONDARY,
    "ytick.color": INK_SECONDARY,
    "lines.linewidth": 2.0,
    "legend.frameon": False,
})


def save(fig, name, **layout):
    fig.tight_layout(**layout)
    fig.savefig(os.path.join(FIGURES, name))
    plt.close(fig)
    print(f"wrote plots/figures/{name}")
