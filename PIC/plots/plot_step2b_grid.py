"""Figure for Step 2b: order-1 deposit (hat functions) and gather (linear interpolation),
measured from the code and compared with what is expected.

Run after building:  python3 plots/plot_step2b_grid.py
"""
import numpy as np

from style import BLUE, INK, INK_SECONDARY, SERIES, plt, save
from common import pic

NX, DX = 6, 1.0
WEIGHT = 1.0
NODE_VALUES = np.array([0.2, 0.9, 0.5, -0.4, -0.1, 0.6])


def deposited_density(grid, x):
    population = pic.Population("protons", 1.0, 1.0)
    population.x, population.w = [x], [WEIGHT]
    population.vx = population.vy = population.vz = [0.0]
    population.deposit(grid)
    return population.density


def expected_hat(x, node_position, grid):
    """w max(0, 1 - d/dx), d = periodic distance between the particle and the node."""
    distance = np.abs(x - node_position)
    distance = np.minimum(distance, grid.length() - distance)
    return WEIGHT * np.maximum(0.0, 1.0 - distance / grid.dx)


def expected_interpolation(x, grid):
    """Straight line between the values of the two nodes around x (node nx is node 0 again)."""
    node_positions = np.arange(grid.nx + 1) * grid.dx
    return np.interp(x, node_positions, np.append(NODE_VALUES, NODE_VALUES[0]))


def mark_nodes(ax, grid, label_height):
    """Dotted line and name at every node; node 0 appears again at x = L (periodic)."""
    nodes_and_positions = [(node, grid.node_position(node)) for node in range(grid.nx)] + [(0, grid.length())]
    for node, x in nodes_and_positions:
        ax.axvline(x, color=INK_SECONDARY, lw=0.6, ls=":")
        ax.text(x, label_height, f"node {node}", color=INK_SECONDARY, ha="center", fontsize=8)


def plot_deposit(ax, grid, positions):
    measured = np.array([deposited_density(grid, x) for x in positions])
    largest_difference = 0.0
    for node in range(grid.nx):
        expected = expected_hat(positions, grid.node_position(node), grid)
        largest_difference = max(largest_difference, np.max(np.abs(measured[:, node] - expected)))
        ax.plot(positions, measured[:, node], color=SERIES[node], lw=3,
                label="measured: density on each node (one colour per node)" if node == 0 else None)
        ax.plot(positions, expected, color=INK, lw=1, ls="--",
                label="expected: w · max(0, 1 − |x − x_node| / dx)" if node == 0 else None)
    mark_nodes(ax, grid, label_height=1.04)
    ax.plot(positions, measured.sum(axis=1), color=INK_SECONDARY, lw=1.5, ls=":",
            label="measured: sum over the nodes (expected: w = 1)")

    ax.set_xlabel("position x of the particle  (nodes at x = 0, 1, …, 5; periodic, L = 6)")
    ax.set_ylabel("density received by the node")
    ax.set_title("Deposit: one particle of weight w = 1 moved through the domain\n"
                 f"largest |measured − expected| = {largest_difference:.0e}")
    ax.set_ylim(-0.05, 1.5)
    ax.legend(loc="upper center", fontsize=8)


def plot_gather(ax, grid, positions):
    field = pic.VectorField(grid)
    field.x = NODE_VALUES
    measured = np.array([pic.gather(field, grid, x).x for x in positions])
    expected = expected_interpolation(positions, grid)
    largest_difference = np.max(np.abs(measured - expected))

    ax.plot(positions, measured, color=BLUE, lw=3, label="measured: field gathered at the particle position")
    ax.plot(positions, expected, color=INK, lw=1, ls="--",
            label="expected: straight line between the two neighbouring nodes")
    node_positions = np.arange(grid.nx) * grid.dx
    ax.plot(node_positions, NODE_VALUES, "o", color=INK, ms=7, label="field values given on the nodes")
    ax.plot(grid.length(), NODE_VALUES[0], "o", color=INK, ms=7, mfc="white", label="node 0 again (periodic)")
    mark_nodes(ax, grid, label_height=1.02)

    ax.set_xlabel("position x of the particle")
    ax.set_ylabel("field felt by the particle")
    ax.set_title("Gather: field felt by a particle moved through the domain\n"
                 f"largest |measured − expected| = {largest_difference:.0e}")
    ax.set_ylim(-1.15, 1.15)
    ax.legend(loc="lower left", fontsize=8)


def plot_deposit_and_gather():
    grid = pic.Grid(NX, DX)
    positions = np.linspace(0.0, grid.length(), 600, endpoint=False)
    fig, (ax_deposit, ax_gather) = plt.subplots(1, 2, figsize=(13, 5))
    plot_deposit(ax_deposit, grid, positions)
    plot_gather(ax_gather, grid, positions)
    save(fig, "step2b_deposit_gather.png")


if __name__ == "__main__":
    plot_deposit_and_gather()
