"""Shared by the tests, the plots and the benchmarks: imports the compiled module `pic`
(from build/, or from $PIC_BUILD_DIR) and builds particles and populations."""
import math
import os
import sys

import numpy as np

HERE = os.path.dirname(os.path.abspath(__file__))
BUILD_DIR = os.environ.get("PIC_BUILD_DIR", os.path.join(HERE, "..", "build"))
sys.path.insert(0, os.path.abspath(BUILD_DIR))

import pic  # noqa: E402

NO_FIELD = pic.Vec3(0.0, 0.0, 0.0)


def field_along_z(amplitude):
    return pic.Vec3(0.0, 0.0, amplitude)


def uniform_fields(grid, E, B):
    """Fields on the grid with the same E and B at every node."""
    fields = pic.Fields(grid)
    for component in ("x", "y", "z"):
        setattr(fields.E, component, np.full(grid.nx, getattr(E, component)))
        setattr(fields.B, component, np.full(grid.nx, getattr(B, component)))
    return fields


def cyclotron_period(q, m, B):
    return 2.0 * math.pi * m / (abs(q) * B)


def larmor_radius(v_perp, q, m, B):
    return m * v_perp / (abs(q) * B)


def make_particle(v, q=1.0, m=1.0, x=0.0):
    particle = pic.Particle()
    particle.x = x
    particle.v = pic.Vec3(*v)
    particle.q, particle.m = q, m
    return particle


def thermal_population(n, length, temperature=1.0, q=1.0, m=1.0, seed=0):
    """Uniform positions in [0, length), Maxwellian velocities (each component has variance T/m),
    equal weights 1/n."""
    rng = np.random.default_rng(seed)
    population = pic.Population("protons", q, m)
    population.x = rng.uniform(0.0, length, n)
    thermal_speed = math.sqrt(temperature / m)
    population.vx = rng.normal(0.0, thermal_speed, n)
    population.vy = rng.normal(0.0, thermal_speed, n)
    population.vz = rng.normal(0.0, thermal_speed, n)
    population.w = np.full(n, 1.0 / n)
    return population


def temperature(population):
    """T = (m/3) <|v - <v>|^2>, averages weighted by w."""
    w = population.w
    variance_sum = 0.0
    for v in (population.vx, population.vy, population.vz):
        mean = np.average(v, weights=w)
        variance_sum += np.average((v - mean) ** 2, weights=w)
    return population.m * variance_sum / 3.0


def example_profiles(length):
    """Periodic density, bulk velocity and temperature profiles over [0, length)."""
    k = 2.0 * math.pi / length

    def density(x):
        return 1.0 + 0.5 * math.sin(k * x)

    def bulk_velocity(x):
        return pic.Vec3(0.3 * math.cos(k * x), 0.1, 0.0)

    def temperature(x):
        return 0.5 + 0.2 * math.sin(2.0 * k * x)

    return density, bulk_velocity, temperature


def loader_expectation(grid, m, density, bulk_velocity, temperature):
    """Exact average of the deposited n, u_x and T for load_maxwellian: every cell carries the
    profiles of its centre and gives half of each moment to each of its two nodes."""
    centres = cell_centre_positions(grid)
    n = np.array([density(x) for x in centres])
    u = np.array([[bulk_velocity(x).x, bulk_velocity(x).y, bulk_velocity(x).z] for x in centres])
    T = np.array([temperature(x) for x in centres])
    energy = n * (0.5 * m * np.sum(u**2, axis=1) + 1.5 * T)

    def to_nodes(cell_values):
        return 0.5 * (cell_values + np.roll(cell_values, 1, axis=0))

    node_density = to_nodes(n)
    node_velocity = to_nodes(n[:, None] * u) / node_density[:, None]
    node_energy = to_nodes(energy)
    node_temperature = (2.0 / 3.0) * (node_energy / node_density - 0.5 * m * np.sum(node_velocity**2, axis=1))
    return node_density, node_velocity[:, 0], node_temperature


def node_positions(grid):
    return np.arange(grid.nx) * grid.dx


def cell_centre_positions(grid):
    return (np.arange(grid.nx) + 0.5) * grid.dx


def plasma_parameters(electron_temperature=0.0, resistivity=0.0, hyper_resistivity=0.0):
    parameters = pic.PlasmaParameters()
    parameters.electron_temperature = electron_temperature
    parameters.resistivity = resistivity
    parameters.hyper_resistivity = hyper_resistivity
    return parameters


def vector_field_on_nodes(grid, x, y, z):
    field = pic.VectorField(grid)
    field.x, field.y, field.z = np.broadcast_to(x, grid.nx), np.broadcast_to(y, grid.nx), np.broadcast_to(z, grid.nx)
    return field


def magnetic_field(grid, x, y, z):
    """B with the Yee layout: x is given on the nodes, y and z at the cell centres."""
    B = pic.Fields(grid).B
    B.x, B.y, B.z = np.broadcast_to(x, grid.nx), np.broadcast_to(y, grid.nx), np.broadcast_to(z, grid.nx)
    return B


class SmoothPlasma:
    """A smooth periodic plasma in which every term of Ohm's law is non-zero, with its exact
    electric field: the analytic case the field solver is compared with."""

    def __init__(self, length, parameters):
        self.k = 2.0 * math.pi / length
        self.parameters = parameters

    def density(self, x):
        return 1.0 + 0.3 * np.sin(self.k * x)

    def bulk_velocity(self, x):
        k = self.k
        return 0.2 * np.cos(k * x), 0.1 * np.sin(k * x), -0.15 * np.cos(2.0 * k * x)

    def B(self, x):
        k = self.k
        return 0.8 + 0.0 * x, 0.5 * np.cos(k * x), 0.4 * np.sin(2.0 * k * x)

    def current(self, x):
        """j = curl B = (0, -dBz/dx, dBy/dx)."""
        k = self.k
        return 0.0 * x, -0.8 * k * np.cos(2.0 * k * x), -0.5 * k * np.sin(k * x)

    def E(self, x):
        """E = -u x B + (j x B) / n - T_e grad(n) / n + eta j - nu laplacian(j)."""
        k, p = self.k, self.parameters
        n = self.density(x)
        ux, uy, uz = self.bulk_velocity(x)
        Bx, By, Bz = self.B(x)
        _, jy, jz = self.current(x)
        dn_dx = 0.3 * k * np.cos(k * x)
        laplacian_jy = 3.2 * k**3 * np.cos(2.0 * k * x)
        laplacian_jz = 0.5 * k**3 * np.sin(k * x)
        Ex = -(uy * Bz - uz * By) + (jy * Bz - jz * By) / n - p.electron_temperature * dn_dx / n
        Ey = -(uz * Bx - ux * Bz) + jz * Bx / n + p.resistivity * jy - p.hyper_resistivity * laplacian_jy
        Ez = -(ux * By - uy * Bx) - jy * Bx / n + p.resistivity * jz - p.hyper_resistivity * laplacian_jz
        return Ex, Ey, Ez

    def solved_E(self, grid):
        """E computed by the code from the exact n, u, B given on the Yee grid."""
        nodes, cell_centres = node_positions(grid), cell_centre_positions(grid)
        B = magnetic_field(grid, self.B(nodes)[0], self.B(cell_centres)[1], self.B(cell_centres)[2])
        bulk_velocity = vector_field_on_nodes(grid, *self.bulk_velocity(nodes))
        return pic.ohm(self.density(nodes), bulk_velocity, B, pic.ampere(B, grid), grid, self.parameters)

    def exact_E_on_the_yee_grid(self, grid):
        nodes, cell_centres = node_positions(grid), cell_centre_positions(grid)
        return self.E(cell_centres)[0], self.E(nodes)[1], self.E(nodes)[2]


def ions_at_rest(grid, density=1.0):
    """Uniform density and zero bulk velocity on the nodes: the fixed ion background of Step 3c."""
    return np.full(grid.nx, density), vector_field_on_nodes(grid, 0.0, 0.0, 0.0)


def transverse_wave(grid, mode, amplitude, Bx):
    """Fields with B = (Bx, b cos(kx), b sin(kx)), k = 2 pi mode / L."""
    k = 2.0 * math.pi * mode / grid.length()
    x = cell_centre_positions(grid)
    fields = pic.Fields(grid)
    fields.B = magnetic_field(grid, Bx, amplitude * np.cos(k * x), amplitude * np.sin(k * x))
    return fields, k


def transverse_amplitude(B, grid, mode):
    """Complex amplitude of By + i Bz in the Fourier mode exp(i k x): its modulus is the wave
    amplitude, its phase turns in time at the wave frequency."""
    k = 2.0 * math.pi * mode / grid.length()
    return np.mean((B.y + 1j * B.z) * np.exp(-1j * k * cell_centre_positions(grid)))


def follow_amplitude(fields, grid, parameters, mode, dt, steps_between_samples, n_samples):
    """Advances the fields with the ions at rest; returns the times and the complex amplitude of
    the mode at each sample."""
    density, bulk_velocity = ions_at_rest(grid)
    times, amplitudes = [0.0], [transverse_amplitude(fields.B, grid, mode)]
    for sample in range(1, n_samples + 1):
        pic.advance_fields_with_fixed_ions(fields, density, bulk_velocity, grid, parameters, dt, steps_between_samples)
        times.append(sample * steps_between_samples * dt)
        amplitudes.append(transverse_amplitude(fields.B, grid, mode))
    return np.array(times), np.array(amplitudes)
