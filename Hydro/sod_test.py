"""Validation: the Sod shock tube along x, compared with the exact Riemann solution."""
import matplotlib.pyplot as plt
import numpy as np

from exact_riemann import RiemannState, exact_solution
from gudanov2d import Godunov2D, hll_flux, hllc_flux

GAMMA = 1.4
X_MIN, X_MAX = 0.0, 10.0
Y_MIN, Y_MAX = 0.0, 1.0
NY = 4  # the problem is 1-D along x, a few rows are enough
INTERFACE = 5.0
T_END = 2.0
CFL = 0.4
LEFT = RiemannState(density=1.0, velocity=0.0, pressure=1.0, gamma=GAMMA)
RIGHT = RiemannState(density=0.125, velocity=0.0, pressure=0.1, gamma=GAMMA)

SOLVERS = {"HLL": hll_flux, "HLLC": hllc_flux}
PLOT_RESOLUTION = 200


def sod_initial_state(X):
    on_left = X < INTERFACE
    rho = np.where(on_left, LEFT.density, RIGHT.density)
    u = np.where(on_left, LEFT.velocity, RIGHT.velocity)
    v = np.zeros_like(X)
    P = np.where(on_left, LEFT.pressure, RIGHT.pressure)
    return np.array([rho, u, v, P])


def run_sod(nx, riemann_solver):
    """Density, velocity and pressure along the middle row at T_END, with the cell centres x."""
    sim = Godunov2D(nx, NY, X_MIN, X_MAX, Y_MIN, Y_MAX, GAMMA, clouds=[],
                    riemann_solver=riemann_solver, boundary="outflow")
    X, Y = sim.cell_centres()
    sim.set_primitives(sod_initial_state(X))
    sim.run_until(T_END, CFL)

    rho, u, v, P = sim.interior_primitives()
    middle_row = NY // 2
    x = X_MIN + (np.arange(nx) + 0.5) * sim.dx
    return x, rho[middle_row], u[middle_row], P[middle_row]


def density_error(nx, riemann_solver):
    """L1 error of the density against the exact solution."""
    x, rho, u, P = run_sod(nx, riemann_solver)
    exact_rho, exact_u, exact_P = exact_solution(x, INTERFACE, T_END, LEFT, RIGHT, GAMMA)
    return np.mean(np.abs(rho - exact_rho)) * (X_MAX - X_MIN)


def print_errors():
    for solver_name, riemann_solver in SOLVERS.items():
        error = density_error(PLOT_RESOLUTION, riemann_solver)
        print(f"{solver_name:5} nx = {PLOT_RESOLUTION}   L1 density error against the exact solution = {error:.4e}")


def plot_comparison():
    figure, axes = plt.subplots(1, 3, figsize=(15, 4), layout="constrained")
    x_fine = np.linspace(X_MIN, X_MAX, 2000)
    exact_fields = exact_solution(x_fine, INTERFACE, T_END, LEFT, RIGHT, GAMMA)
    labels = ["density", "velocity", "pressure"]
    for ax, label, exact_field in zip(axes, labels, exact_fields):
        ax.plot(x_fine, exact_field, "k-", label="exact")
        ax.set_xlabel("x")
        ax.set_ylabel(label)

    for solver_name, riemann_solver in SOLVERS.items():
        x, rho, u, P = run_sod(PLOT_RESOLUTION, riemann_solver)
        for ax, field in zip(axes, [rho, u, P]):
            ax.plot(x, field, ".", markersize=4, label=solver_name)

    axes[0].legend()
    figure.suptitle(f"Sod shock tube at t = {T_END}, nx = {PLOT_RESOLUTION}")
    figure.savefig("sod_test.png", dpi=130)
    plt.close(figure)
    print("saved sod_test.png")


if __name__ == "__main__":
    plot_comparison()
    print_errors()
