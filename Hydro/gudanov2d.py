import matplotlib.pyplot as plt
import numpy as np

NG = 2  # two neighbours on each side of a cell (ghost cells)
DT_MIN = 1e-8  # below this the run is stalling 

# Floors protect divisions and square roots; they are not a positivity-preserving scheme.
DENSITY_FLOOR = 1e-12
PRESSURE_FLOOR = 1e-12
MIN_SPEED_GAP = 1e-12


# ---------- Stage 1: physical model, W = (rho, u, v, P) <-> U = (rho, rho u, rho v, E) 

def sound_speed(gamma, P, rho):
    return np.sqrt(gamma * P / rho)


def total_energy(rho, u, v, P, gamma):
    """E = internal energy P / (gamma - 1) + kinetic energy rho (u^2 + v^2) / 2."""
    return P / (gamma - 1) + 0.5 * rho * (u**2 + v**2)


def prim_to_cons(W, gamma):
    rho, u, v, P = W
    E = total_energy(rho, u, v, P, gamma)
    return np.array([rho, rho * u, rho * v, E])


def cons_to_prim(U, gamma):
    rho = np.maximum(U[0], DENSITY_FLOOR)
    u = U[1] / rho
    v = U[2] / rho
    P = (gamma - 1) * (U[3] - 0.5 * rho * (u**2 + v**2))
    P = np.maximum(P, PRESSURE_FLOOR)
    return np.array([rho, u, v, P])


# ---------- Stage 2: 1-D Riemann problem at a face, HLL and HLLC 
# idir = 0 for an x-face (normal velocity u), idir = 1 for a y-face (normal velocity v).

def normal_velocity(W, idir):
    """Velocity across the face."""
    if idir == 0:
        return W[1]
    return W[2]


def physical_flux(W, gamma, idir):
    """F(U) through an x-face, G(U) through a y-face (slide 14)."""
    rho, u, v, P = W
    E = total_energy(rho, u, v, P, gamma)
    if idir == 0:
        return np.array([rho * u, rho * u**2 + P, rho * u * v, u * (E + P)])
    return np.array([rho * v, rho * u * v, rho * v**2 + P, v * (E + P)])


def outer_wave_speeds(W_L, W_R, gamma, idir):
    """S_L and S_R: estimates of the fastest left- and right-going waves (Davis)."""
    un_L = normal_velocity(W_L, idir)
    un_R = normal_velocity(W_R, idir)
    c_L = sound_speed(gamma, W_L[3], W_L[0])
    c_R = sound_speed(gamma, W_R[3], W_R[0])
    S_L = np.minimum(un_L - c_L, un_R - c_R)
    S_R = np.maximum(un_L + c_L, un_R + c_R)
    return S_L, S_R


def hll_star_flux(F_L, F_R, U_L, U_R, S_L, S_R):
    """Flux of the single averaged state between S_L and S_R."""
    return (S_R * F_L - S_L * F_R + S_L * S_R * (U_R - U_L)) / (S_R - S_L)


def hll_flux(W_L, W_R, gamma, idir):
    S_L, S_R = outer_wave_speeds(W_L, W_R, gamma, idir)
    F_L = physical_flux(W_L, gamma, idir)
    F_R = physical_flux(W_R, gamma, idir)
    U_L = prim_to_cons(W_L, gamma)
    U_R = prim_to_cons(W_R, gamma)
    F_star = hll_star_flux(F_L, F_R, U_L, U_R, S_L, S_R)

    # The face lies left of all waves (S_L >= 0), right of all waves (S_R <= 0), or between.
    F = np.where(S_R <= 0, F_R, F_star)
    F = np.where(S_L >= 0, F_L, F)
    return F


def contact_speed(W_L, W_R, S_L, S_R, idir):
    """S_star: speed of the contact wave in the middle of the star region."""
    rho_L, P_L = W_L[0], W_L[3]
    rho_R, P_R = W_R[0], W_R[3]
    un_L = normal_velocity(W_L, idir)
    un_R = normal_velocity(W_R, idir)
    numerator = P_R - P_L + rho_L * un_L * (S_L - un_L) - rho_R * un_R * (S_R - un_R)
    denominator = rho_L * (S_L - un_L) - rho_R * (S_R - un_R)
    return numerator / denominator


def nonzero_speed_gap(gap):
    """gap, kept at least MIN_SPEED_GAP away from 0 with its sign."""
    too_small = np.abs(gap) < MIN_SPEED_GAP
    sign = np.where(gap < 0, -1.0, 1.0)
    return np.where(too_small, sign * MIN_SPEED_GAP, gap)


def hllc_star_state(W, U, S, S_star, idir):
    """Conserved state between the outer wave S and the contact S_star.

    Every component, energy included, follows from the Rankine-Hugoniot jump across the wave S
    (Toro 10.39); this keeps the energy flux through a contact at rest exactly zero.
    """
    rho, u, v, P = W
    un = normal_velocity(W, idir)
    E = U[3]
    rho_star = rho * (S - un) / nonzero_speed_gap(S - S_star)

    if idir == 0:
        u_star = S_star
        v_star = v
    else:
        u_star = u
        v_star = S_star
    E_star = rho_star * (E / rho + (S_star - un) * (S_star + P / (rho * (S - un))))
    return np.array([rho_star, rho_star * u_star, rho_star * v_star, E_star])


def hllc_flux(W_L, W_R, gamma, idir):
    S_L, S_R = outer_wave_speeds(W_L, W_R, gamma, idir)
    S_star = contact_speed(W_L, W_R, S_L, S_R, idir)
    F_L = physical_flux(W_L, gamma, idir)
    F_R = physical_flux(W_R, gamma, idir)
    U_L = prim_to_cons(W_L, gamma)
    U_R = prim_to_cons(W_R, gamma)

    U_star_L = hllc_star_state(W_L, U_L, S_L, S_star, idir)
    U_star_R = hllc_star_state(W_R, U_R, S_R, S_star, idir)
    F_star_L = F_L + S_L * (U_star_L - U_L)
    F_star_R = F_R + S_R * (U_star_R - U_R)

    F = np.where(S_R >= 0, F_star_R, F_R)
    F = np.where(S_star >= 0, F_star_L, F)
    F = np.where(S_L >= 0, F_L, F)
    return F


# ---------- Stage 3: finite volumes, ghost cells and MUSCL reconstruction

def array_axis(idir):
    """Axis of the (variable, y, x) arrays along which idir runs."""
    if idir == 0:
        return 2
    return 1


def fill_outflow_ghosts(U):
    """The ghost cells copy the nearest domain cell, so waves leave freely."""
    U[:, :, 0] = U[:, :, NG]
    U[:, :, 1] = U[:, :, NG]
    U[:, :, -1] = U[:, :, -1 - NG]
    U[:, :, -2] = U[:, :, -1 - NG]
    U[:, 0, :] = U[:, NG, :]
    U[:, 1, :] = U[:, NG, :]
    U[:, -1, :] = U[:, -1 - NG, :]
    U[:, -2, :] = U[:, -1 - NG, :]
    return U


def fill_periodic_ghosts(U):
    """Gas leaving through one edge comes back through the opposite edge."""
    U[:, :, :NG] = U[:, :, -2 * NG:-NG]
    U[:, :, -NG:] = U[:, :, NG:2 * NG]
    U[:, :NG, :] = U[:, -2 * NG:-NG, :]
    U[:, -NG:, :] = U[:, NG:2 * NG, :]
    return U


def apply_bcs(U, boundary):
    """Fill the NG = 2 ghost layers on all four sides (slide 10)."""
    if boundary == "outflow":
        return fill_outflow_ghosts(U)
    if boundary == "periodic":
        return fill_periodic_ghosts(U)
    raise ValueError(f"unknown boundary {boundary!r}: use 'outflow' or 'periodic'")


def minmod(a, b):
    """The smaller slope if a and b have the same sign, 0 otherwise."""
    smaller = np.where(np.abs(a) < np.abs(b), a, b)
    return np.where(a * b > 0, smaller, 0.0)


def limited_slopes(W, axis):
    """sigma_i = minmod(q_i - q_i-1, q_i+1 - q_i), for every variable and cell."""
    backward_difference = W - np.roll(W, 1, axis=axis)
    forward_difference = np.roll(W, -1, axis=axis) - W
    return minmod(backward_difference, forward_difference)


def face_states(W, idir):
    """Left and right states at the face i+1/2, stored at index i.

    left:  q_i + sigma_i / 2, extrapolated from cell i;
    right: q_i+1 - sigma_i+1 / 2, extrapolated from cell i+1.
    """
    axis = array_axis(idir)
    slopes = limited_slopes(W, axis)
    left_state = W + 0.5 * slopes
    right_state = np.roll(W - 0.5 * slopes, -1, axis=axis)
    return left_state, right_state


# ---------- Stage 4: spatial right-hand side and time integration 

def flux_divergence(F, axis, spacing):
    """(F_i+1/2 - F_i-1/2) / spacing, with F_i+1/2 stored at index i."""
    return (F - np.roll(F, 1, axis=axis)) / spacing


def directional_divergence(W, spacing, gamma, riemann_solver, idir):
    """Reconstruct, solve the Riemann problem at every face, take the divergence."""
    W_left, W_right = face_states(W, idir)
    F = riemann_solver(W_left, W_right, gamma, idir)
    return flux_divergence(F, array_axis(idir), spacing)


def godunov_update(U, dx, dy, gamma, riemann_solver):
    """L(U) = -(F_i+1/2 - F_i-1/2)/dx - (G_j+1/2 - G_j-1/2)/dy (slide 17); U has filled ghosts."""
    W = cons_to_prim(U, gamma)
    divergence_x = directional_divergence(W, dx, gamma, riemann_solver, idir=0)
    divergence_y = directional_divergence(W, dy, gamma, riemann_solver, idir=1)
    return -(divergence_x + divergence_y)


def euler_step(U, dt, dx, dy, gamma, riemann_solver, boundary):
    """U + dt L(U), with the ghost cells refilled for the next evaluation of L."""
    U_next = U + dt * godunov_update(U, dx, dy, gamma, riemann_solver)
    return apply_bcs(U_next, boundary)


def advance(U, dt, dx, dy, gamma, riemann_solver, boundary):
    """SSP-RK2 (slide 18): U1 = U + dt L(U), then U_new = (U + (U1 + dt L(U1))) / 2."""
    U_1 = euler_step(U, dt, dx, dy, gamma, riemann_solver, boundary)
    U_2 = euler_step(U_1, dt, dx, dy, gamma, riemann_solver, boundary)
    U_new = 0.5 * U + 0.5 * U_2
    return apply_bcs(U_new, boundary)


# ---------- Stage 5: the 2-D simulation, initial conditions and outputs

class Cloud:
    """A round cloud of gas (density, pressure, velocity u, v) in the ambient medium.

    Its profile is a Gaussian: it peaks at the centre and falls to 1/e of the peak at radius.
    """

    def __init__(self, x0, y0, radius, density, pressure, u, v=0.0):
        if radius <= 0:
            raise ValueError("radius must be positive")
        if density <= 0 or pressure <= 0:
            raise ValueError("density and pressure must be positive")
        self.x0 = x0
        self.y0 = y0
        self.radius = radius
        self.density = density
        self.pressure = pressure
        self.u = u
        self.v = v

    def profile(self, X, Y):
        """How much of the cloud is in each cell: 1 at the centre, 0 far away."""
        distance = np.sqrt((X - self.x0) ** 2 + (Y - self.y0) ** 2)
        return np.exp(-(distance / self.radius) ** 2)


class Godunov2D:
    def __init__(self, nx, ny, xmin, xmax, ymin, ymax, gamma, clouds,
                 ambient_density=0.01, ambient_pressure=0.01,
                 riemann_solver=hllc_flux, boundary="outflow"):
        """clouds is a list of non-overlapping Cloud, in a uniform ambient medium at rest;
        boundary is "outflow" or "periodic", on all four sides."""
        self._check_inputs(nx, ny, xmin, xmax, ymin, ymax, gamma, boundary)
        self.t = 0.0
        self.clouds = clouds
        self.ambient_density = ambient_density
        self.ambient_pressure = ambient_pressure
        self.dx = (xmax - xmin) / nx
        self.dy = (ymax - ymin) / ny
        self.extent = (xmin, xmax, ymin, ymax)
        self.gamma = gamma
        self.riemann_solver = riemann_solver
        self.boundary = boundary
        self.nx = nx
        self.ny = ny
        X, Y = self.cell_centres()
        self.set_primitives(self.clouds_in_ambient_medium(X, Y))

    def _check_inputs(self, nx, ny, xmin, xmax, ymin, ymax, gamma, boundary):
        if nx < 3 or ny < 3:
            raise ValueError("nx and ny must be at least 3 physical cells")
        if xmax <= xmin or ymax <= ymin:
            raise ValueError("the domain bounds must be increasing: xmin < xmax, ymin < ymax")
        if gamma <= 1:
            raise ValueError("gamma must be larger than 1")
        if boundary not in ("outflow", "periodic"):
            raise ValueError(f"unknown boundary {boundary!r}: use 'outflow' or 'periodic'")

    def cell_centres(self):
        """X and Y of every cell centre, ghost cells included."""
        xmin, xmax, ymin, ymax = self.extent
        x = xmin + (np.arange(-NG, self.nx + NG) + 0.5) * self.dx
        y = ymin + (np.arange(-NG, self.ny + NG) + 0.5) * self.dy
        return np.meshgrid(x, y)

    def set_primitives(self, W):
        """Replace the state by the primitives W = [rho, u, v, P] given on the whole grid."""
        self.U = apply_bcs(prim_to_cons(W, self.gamma), self.boundary)

    def clouds_in_ambient_medium(self, X, Y):
        """Primitives [rho, u, v, P]: each cloud's state inside it, the ambient state elsewhere."""
        rho = np.full_like(X, self.ambient_density)
        P = np.full_like(X, self.ambient_pressure)
        u = np.zeros_like(X)
        v = np.zeros_like(X)
        for cloud in self.clouds:
            inside = cloud.profile(X, Y)
            rho += (cloud.density - self.ambient_density) * inside
            P += (cloud.pressure - self.ambient_pressure) * inside
            u += cloud.u * inside
            v += cloud.v * inside
        return np.array([rho, u, v, P])

    def interior_primitives(self):
        return cons_to_prim(self.U[:, NG:-NG, NG:-NG], self.gamma)

    def compute_dt(self, CFL):
        """CFL times the shortest time a signal takes to cross a cell ."""
        rho, u, v, P = self.interior_primitives()
        c = sound_speed(self.gamma, P, rho)
        dt_x = self.dx / np.max(np.abs(u) + c)
        dt_y = self.dy / np.max(np.abs(v) + c)
        return CFL * min(dt_x, dt_y)

    def step(self, dt):
        """One SSP-RK2 step; keeps the state U and the time t together."""
        self.U = advance(self.U, dt, self.dx, self.dy, self.gamma,
                         self.riemann_solver, self.boundary)
        self.t += dt

    def run_until(self, t_end, CFL):
        """Advance from the current time to exactly t_end."""
        while self.t < t_end:
            dt = self.compute_dt(CFL)
            # Checked before the clipping below, which may legitimately make the last step tiny.
            if dt < DT_MIN:
                raise RuntimeError(f"time step too small: dt = {dt:.3e} at t = {self.t:.4f}")
            if self.t + dt > t_end:
                dt = t_end - self.t
            self.step(dt)

    def total_mass(self):
        rho = self.interior_primitives()[0]
        return np.sum(rho) * self.dx * self.dy

    def fields(self):
        """The plotted fields on the domain (ghost cells removed)."""
        rho, u, v, P = self.interior_primitives()
        return {
            "density": rho,
            "pressure": P,
            "velocity": np.sqrt(u**2 + v**2),
        }

    def evolve(self, output_times, CFL):
        """Run to each output time, keeping the fields and the total mass there."""
        self.output_times = output_times
        self.snapshots = []
        self.masses = []
        for t_end in output_times:
            self.run_until(t_end, CFL)
            self.snapshots.append(self.fields())
            self.masses.append(self.total_mass())

    def field_to_plot(self, snapshot, field_name):
        if field_name in LOG_FIELDS:
            return np.log10(snapshot[field_name])
        return snapshot[field_name]

    def colour_range(self, field_name):
        """Same colour range at every output time, so the panels of a row can be compared."""
        lowest = np.inf
        highest = -np.inf
        for snapshot in self.snapshots:
            values = self.field_to_plot(snapshot, field_name)
            lowest = min(lowest, np.min(values))
            highest = max(highest, np.max(values))
        return lowest, highest

    def plot_field_row(self, axes_row, field_name):
        """Draw the field at each output time, one time per axis; returns the last image."""
        vmin, vmax = self.colour_range(field_name)
        for snapshot, t, ax in zip(self.snapshots, self.output_times, axes_row):
            values = self.field_to_plot(snapshot, field_name)
            image = ax.imshow(values, origin="lower", extent=self.extent,
                              cmap="inferno", vmin=vmin, vmax=vmax)
            ax.set_title(f"t = {t:.1f}")
            ax.set_xlabel("x")
            ax.set_ylabel("y")
        return image

    def plot_fields(self, case_name):
        """One figure: a row per field, a column per output time; saved as a png."""
        figure, axes = plt.subplots(len(FIELD_LABELS), len(self.output_times),
                                    figsize=(15, 7), layout="constrained")
        row = 0
        for field_name, label in FIELD_LABELS.items():
            image = self.plot_field_row(axes[row], field_name)
            figure.colorbar(image, ax=axes[row], label=label)
            row += 1

        figure.suptitle(f"Two-cloud collision: {case_name}")
        file_name = "collision_" + case_name.replace(", ", "_") + ".png"
        figure.savefig(file_name, dpi=130)
        plt.close(figure)
        print(f"saved {file_name}")


NX, NY = 200, 100
X_MIN, X_MAX = 0.0, 10.0
Y_MIN, Y_MAX = 0.0, 5.0
GAMMA = 1.4
CFL = 0.4
OUTPUT_TIMES = [0.0, 1.0, 2.0, 3.0, 4.0, 5.0]

CLOUD_SETUPS = {
    "pressure_balanced": [
        Cloud(x0=2.5, y0=2.5, radius=0.4, density=2.5, pressure=0.01, u=3.0),
        Cloud(x0=7.5, y0=2.5, radius=1.0, density=2.5, pressure=0.01, u=-1.0),
    ],
}
CLOUD_SETUP = "pressure_balanced"
CLOUDS = CLOUD_SETUPS[CLOUD_SETUP]
AMBIENT_DENSITY = 0.01
AMBIENT_PRESSURE = 0.01

SOLVERS = {"HLL": hll_flux, "HLLC": hllc_flux}
BOUNDARIES = ["outflow", "periodic"]

# Density and pressure span two decades (clouds against a 0.01 medium), so they are shown in log10.
FIELD_LABELS = {
    "density": "log10 density",
    "pressure": "log10 pressure",
    "velocity": "velocity magnitude",
}
LOG_FIELDS = ["density", "pressure"]


def print_masses(simulations):
    print("\nTotal mass in the domain")
    header = ""
    for t in OUTPUT_TIMES:
        header += f"{f't = {t}':>10}"
    print(f"{'':36}{header}")
    for case_name, sim in simulations.items():
        row = ""
        for mass in sim.masses:
            row += f"{mass:10.3f}"
        print(f"{case_name:36}{row}")


def main():
    simulations = {}

    case_name = f"{CLOUD_SETUP}, HLL, outflow"
    print(f"running {case_name}")
    sim = Godunov2D(NX, NY, X_MIN, X_MAX, Y_MIN, Y_MAX, GAMMA, CLOUDS,
                    ambient_density=AMBIENT_DENSITY, ambient_pressure=AMBIENT_PRESSURE,
                    riemann_solver=hll_flux, boundary="outflow")
    sim.evolve(OUTPUT_TIMES, CFL)
    sim.plot_fields(case_name)
    simulations[case_name] = sim

    # for solver_name, riemann_solver in SOLVERS.items():
    #     for boundary in BOUNDARIES:
    #         case_name = f"{CLOUD_SETUP}, {solver_name}, {boundary}"
    #         print(f"running {case_name}")
    #         sim = Godunov2D(NX, NY, X_MIN, X_MAX, Y_MIN, Y_MAX, GAMMA, CLOUDS,
    #                         ambient_density=AMBIENT_DENSITY, ambient_pressure=AMBIENT_PRESSURE,
    #                         riemann_solver=riemann_solver, boundary=boundary)
    #         sim.evolve(OUTPUT_TIMES, CFL)
    #         sim.plot_fields(case_name)
    #         simulations[case_name] = sim

    print_masses(simulations)


if __name__ == "__main__":
    main()
