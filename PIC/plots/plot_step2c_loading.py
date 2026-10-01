"""Figure for Step 2c: moments of a Maxwellian loaded from user profiles, and their error vs N:
what the code gives (coloured) against what is expected (thin, black, dashed).

Run after building:  python3 plots/plot_step2c_loading.py
"""
import numpy as np

from style import BLUE, INK, ORANGE, plt, save
from common import example_profiles, loader_expectation, pic

GRID = pic.Grid(16, 0.5)  # coarse on purpose, so that the dx^2 floor is visible
PROFILES = example_profiles(GRID.length())
MOMENTS = ["density", "bulk velocity", "temperature"]
SEEDS = range(3)
EXPECTED = dict(color=INK, lw=1, ls="--")


def deposited_moments(particles_per_cell, seed):
    population = pic.Population("protons", 1.0, 1.0)
    pic.load_maxwellian(population, GRID, particles_per_cell, *PROFILES, seed)
    population.deposit(GRID)
    return [population.density, pic.bulk_velocity([population], GRID).x, population.temperature()]


def profiles_on_nodes(x):
    density, bulk_velocity, temperature = PROFILES
    return [np.array([density(xi) for xi in x]), np.array([bulk_velocity(xi).x for xi in x]),
            np.array([temperature(xi) for xi in x])]


def below_the_row(ax):
    """One legend for the three panels of a row, under the middle one, without squeezing the panels."""
    legend = ax.legend(loc="upper center", bbox_to_anchor=(0.5, -0.17), ncol=3, fontsize=9)
    legend.set_in_layout(False)


def plot_loading():
    nodes = np.arange(GRID.nx) * GRID.dx
    fine_x = np.linspace(0.0, GRID.length(), 400)
    expected = loader_expectation(GRID, 1.0, *PROFILES)
    user_on_nodes = profiles_on_nodes(nodes)
    user_curves = profiles_on_nodes(fine_x)

    particles_per_cell = np.array([10, 30, 100, 300, 1000, 3000, 10000, 30000, 100000])
    error_vs_expectation = np.zeros((3, len(particles_per_cell)))
    error_vs_profile = np.zeros((3, len(particles_per_cell)))
    for j, ppc in enumerate(particles_per_cell):
        for seed in SEEDS:
            moments = deposited_moments(ppc, seed)
            for k in range(3):
                error_vs_expectation[k, j] += np.mean((moments[k] - expected[k]) ** 2) / len(SEEDS)
                error_vs_profile[k, j] += np.mean((moments[k] - user_on_nodes[k]) ** 2) / len(SEEDS)
    error_vs_expectation, error_vs_profile = np.sqrt(error_vs_expectation), np.sqrt(error_vs_profile)
    total_particles = particles_per_cell * GRID.nx

    few, many = deposited_moments(10, 0), deposited_moments(1000, 0)
    fig, axes = plt.subplots(2, 3, figsize=(15, 10))
    for k, name in enumerate(MOMENTS):
        ax = axes[0, k]
        ax.plot(nodes, few[k], "o", color=ORANGE, ms=7, label="measured: deposited, 10 particles / cell")
        ax.plot(nodes, many[k], "o", color=BLUE, ms=7, label="measured: deposited, 1000 particles / cell")
        ax.plot(fine_x, user_curves[k], label="expected: the profile given by the user", **EXPECTED)
        ax.set_xlabel("x")
        ax.set_ylabel(name)
        ax.set_title(f"{name}\nlargest difference: "
                     f"{np.max(np.abs(few[k] - user_on_nodes[k])):.2f} (10 / cell), "
                     f"{np.max(np.abs(many[k] - user_on_nodes[k])):.3f} (1000 / cell)")
        if k == 1:
            below_the_row(ax)

        ax = axes[1, k]
        measured_slope = np.polyfit(np.log(total_particles), np.log(error_vs_expectation[k]), 1)[0]
        expected_decrease = error_vs_expectation[k, 0] * (total_particles / total_particles[0]) ** -0.5
        ax.loglog(total_particles, error_vs_expectation[k], "o-", color=BLUE, ms=7, lw=3,
                  label="measured: error against the loader's exact average")
        ax.loglog(total_particles, error_vs_profile[k], "s-", color=ORANGE, ms=6, lw=1.5,
                  label="measured: error against the user profile (stops at the grid error ∝ dx²)")
        ax.loglog(total_particles, expected_decrease, label="expected: particle noise ∝ 1/√N", **EXPECTED)
        ax.set_xlabel("total number of particles N")
        ax.set_ylabel("RMS error over the nodes")
        ax.set_title(f"error on the {name}\nmeasured slope = {measured_slope:.2f}  (expected: −0.50)")
        if k == 1:
            below_the_row(ax)
    fig.suptitle(f"Maxwellian loaded from n(x), u(x), T(x) on {GRID.nx} cells (dx = {GRID.dx}); "
                 f"errors averaged over {len(SEEDS)} seeds", color=INK)
    save(fig, "step2c_loading.png", h_pad=5.0, rect=(0.0, 0.05, 1.0, 1.0))


if __name__ == "__main__":
    plot_loading()
