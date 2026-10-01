# Tests — what each one checks and why

All tests are in Python (`tests/`), using `unittest` + `ddt` (`@data` runs the same test for
each set of values). Run them with `python3 -m unittest discover -s tests -v` (the `-v` shows each
case; the measured numbers are printed next to the expected ones).

Units are dimensionless: `q`, `m`, `B` are numbers of order 1, time is measured with the
cyclotron period `T_c = 2π m / (|q| B)`, and the Larmor radius is `r_L = m v⊥ / (|q| B)`.
The code is 1D3V: the position is `x` only, the velocity has 3 components. B is along +z, so the
particle turns in the (x, y) plane: `x` oscillates between −r_L and +r_L and `(v_x, v_y)` turns on a circle.

## Step 1 — Boris pusher (`test_step1_boris.py`)

Three tests, one for each part of the scheme: **B only**, **the numerical error**, **E and B together**.

### Key fact about Boris in a uniform B
Each step turns the velocity by `θ = 2 atan(ω_c dt / 2)` instead of the exact `ω_c dt`, **without
changing its length**. So the orbit is a circle with exactly the right radius and the speed never
drifts; the only error is in the phase (the particle is slightly late), and it is ∝ dt².

### 1. Larmor radius and rotation direction (`test_larmor_radius_and_rotation`)
E = 0, B = 1, run for ions of mass 1, 4, 16, 1836 and an electron (q = −1, m = 1/1836).
Checks that `x` oscillates with amplitude `r_L = m v / (|q| B)`, measured as `(x_max − x_min)/2`
(to 1e-6: the extremes are sampled at quarter periods, where Boris' tiny phase lag leaves ~1e-8),
and that the velocity `(v_x, v_y)` turns clockwise for ions and anticlockwise for electrons
(sign of `v(n) × v(n+1)` along z).
*Why*: it tests the physics (radius ∝ m, as asked by the professor), and the direction catches a sign
error that the radius cannot see (a circle has the same radius both ways).

### 2. Error vs time step (`test_error_is_second_order_in_dt`)
After exactly one period the particle should be back at `x = 0`; `|x|` is the error.
Measured for 16 … 512 steps per period: the log-log slope must be 2.
*Why*: the order of accuracy is the standard proof that a scheme is implemented as designed; a bug
usually drops it to 1 or 0.

### 3. E×B drift (`test_ExB_drift`)
E in the plane, B along z. A particle started at `v_E = E × B / B²` feels no force
(`E + v × B = 0`), so it must keep the velocity `v_E` and move as `x = v_E,x t`, for ions, electrons and heavy
ions alike (drift independent of q and m). Boris keeps this exactly (1e-14).
*Why*: the only test where E and B act together, so it checks that half kick → rotation → half kick
combines them correctly. Physically, E×B is how a magnetised plasma flows.

## Step 2a — a population of N particles (`test_step2a_population.py`)

A `Population` is one ion species (q, m) and its macro-particles, stored as one array per quantity
(`x, vx, vy, vz, w`). `push` does one Boris step for every particle (the same
`boris_velocity_update` as Step 1) and wraps `x` periodically into `[0, L)`. Since Step 2b, E and B
are gathered from the grid at the mid-step position; these tests give the same E and B at every node.

### 1. Energy and temperature of a thermal population (`test_thermal_population_keeps_energy_and_temperature`)
10⁵ particles with a Maxwellian of temperature 0.5 in a pure B field, ~16 cyclotron periods.
The total kinetic energy `Σ w m |v|²/2` changes by less than 1e-12 (measured: 6e-16) and the
temperature `T = (m/3)⟨|v − ⟨v⟩|²⟩` is unchanged to 1e-12.
*Why*: a magnetic field does no work, so a real plasma keeps its energy and temperature. This is the
Step 1 property (|v| exactly conserved) shown for a whole population, through the same diagnostics
(energy, temperature) we will use later.

### 2. Periodic wrap (`test_periodic_wrap`)
Domain L = 1, fast particles with an E×B drift: they cross the boundary up to ~14 times. Every `x`
stays in `[0, L)`, and equals, modulo L, the position of the same particle in a large domain where it
never wraps (to 1e-10, distance measured around the circle so that 0.999999 and 0.000001 are close).
The velocities are untouched by the wrap.
*Why*: the grid will be periodic; a particle leaving at `x = L` must come back at `x = 0`, exactly where
the unwrapped motion says, and never end outside the domain (where the deposit would fail).

## Step 2b — periodic grid, deposit and gather (`test_step2b_grid.py`)

Grid: `nx` nodes at `x_i = i dx`, periodic (node `nx` is node 0). The moments live on the nodes; the
field components live on the nodes or at the cell centres (Yee layout, Step 3a).
Order-1 shape function: a particle at `x` between nodes `i` and `i+1` gives the fraction `1 − r` of
itself to node `i` and `r` to node `i+1`, with `r = x/dx − i`. The same weights are used to deposit
(particles → grid) and to gather (grid → particles). Figure: `plots/plot_step2b_grid.py`.

### 1. Deposit is the hat function (`test_deposit_is_the_hat_function`)
One particle of weight `w` is swept through the whole domain, across the boundary. On every node the
deposited density must be `w · max(0, 1 − d/dx)`, `d` = periodic distance to the node (to 1e-14), and the
flux must be `density × v`. Run for two grids (ddt), one with an awkward `dx = 0.37`.
*Why*: this *is* the order-1 shape function, measured directly. Swapping `r` and `1 − r`, an off-by-one
node index or a broken periodic wrap give a hat that is flipped, shifted or cut at the boundary.

### 2. Gather is linear interpolation (`test_gather_is_linear_interpolation`)
Random values on the nodes, or at the cell centres (Yee layout of Step 3a, where the hat is shifted by
half a cell); a particle swept through the domain must feel the straight line between the two points
around it, across the periodic boundary too (to 1e-14).
*Why*: this is what "order 1" means for the fields a particle feels, including across the boundary.

### 3. The deposit conserves what the particles carry (`test_deposit_conserves_weight_momentum_and_energy`)
Thousands of random particles of two species (protons, and alphas with q = 2, m = 4 and a drift).
Summed over the nodes: `Σ density = Σ w`, `Σ flux = Σ w v`, `Σ kinetic_energy_density = Σ ½ m w |v|²`,
and `Σ total_density × bulk_velocity` = the total flux of both species.
*Why*: tests 1–2 follow one particle; this one checks that with many particles anywhere, nothing is
lost or created, and that the sums over populations (slide 5) are right.

## Step 2c — Maxwellian loading from user profiles (`test_step2c_loading.py`)

`load_maxwellian(population, grid, particles_per_cell, density, bulk_velocity, temperature, seed)`:
in every cell, `particles_per_cell` particles at uniform random positions, all with the profiles of the
cell centre (as in PHARE): weight `w = n / particles_per_cell`, each velocity component
`u + √(T/m) · ξ` with `ξ` a standard normal number. The profiles are Python functions of x.
Knowing only n, u, T, the Maxwellian is an assumption (thermal equilibrium, isotropic T).
Figure: `plots/plot_step2c_loading.py`.

### The error of each moment falls as 1/√N (`test_error_decreases_as_one_over_sqrt_N`, for n, u_x, T)
32 cells, 10 to 10240 particles per cell (N = 320 to 327 680). The deposited density, bulk velocity and
temperature are compared on every node with the **exact average of the loader**
(`loader_expectation`: a node receives half of each of its two neighbouring cells). RMS over the
nodes and over 4 seeds; the log-log slope against N must be −0.5 ± 0.1 (measured −0.49, −0.49, −0.52).
*Why*: a moment computed from N random particles has a statistical error ∝ 1/√N (central limit
theorem). Seeing the −1/2 slope for n, u and T shows that the loader draws the right distribution and
that the deposit turns it into the right moments, with nothing but particle noise left.
*Why against the loader's average and not directly against the user profile*: evaluating the profiles
at cell centres and depositing with the hat function gives a small grid bias ∝ dx² (for the density
`dx² n''/8`), independent of N. Against the user profile the error stops decreasing once the noise falls
below that bias (the floor in the figure), so the slope would depend on the N range chosen.

## Step 3a — Yee layout (`test_step3a_yee.py`)

Each field component lives either on the nodes, `x = i dx`, or at the cell centres, `x = (i + ½) dx`:
`Bx, Ey, Ez` (and the moments n, u) on the nodes; `By, Bz, Ex` at the cell centres. Four operations go
from one set of points to the other: `derivative_at_cell_centres`, `derivative_at_nodes`,
`average_at_cell_centres`, `average_at_nodes`. Figure: `plots/plot_step3a_yee.py`.

### The operations are second order in dx (`test_operations_are_second_order_in_dx`, 4 cases)
`sin(kx)` is given on one set of points, differentiated (or averaged) onto the other set and compared
with the exact value there, for 16 to 512 cells. The log-log slope of the largest error against `dx`
must be 2 ± 0.05 (measured 1.98 to 2.00); the errors follow `k³dx²/24` (derivative) and `k²dx²/8` (average).
*Why*: this is the reason for the Yee layout. `(f[i+1] − f[i])/dx` is only first-order accurate at
point `i`, but second-order accurate half-way between `i` and `i+1`: placing the result there gives the
accuracy of a wide stencil with the two nearest points. A wrong shift (result attributed to the wrong
set of points, or the wrong neighbour) drops the slope to 1. It is the check asked by the professor
("error of the second order approximation in the Yee grid layout").

## Step 3b — Ampère and Ohm (`test_step3b_ohm.py`)

Normalised units (`μ0 = e = m_i = 1`). `ampere(B, grid)` gives `j = ∇ × B` (in 1D: `jx = 0`,
`jy = −∂Bz/∂x`, `jz = ∂By/∂x`, on the nodes). `ohm(density, bulk_velocity, B, j, grid, parameters)` gives
`E = −u × B + (j × B)/n − ∇p_e/n + η j − ν ∇²j` with isothermal electrons `p_e = n T_e`; `Ex` at the cell
centres, `Ey`, `Ez` on the nodes. Figure: `plots/plot_step3b_ohm.py`.

### 1. Uniform plasma gives −u × B (`test_uniform_plasma_gives_minus_u_cross_B`)
Uniform n, u and B (all components non-zero), with non-zero `T_e`, `η`, `ν`: the current is exactly 0 and
`E = −u × B` to 1e-15 on the three components.
*Why*: the one case with an exact answer, and the E×B physics of Step 1 seen from the fields: a plasma
drifting at u feels no electric field in its own frame. A wrong sign or a swapped component in the cross
product shows immediately.

### 2. Every term, against an exact solution (`test_error_is_second_order_in_dx`, for Ex, Ey, Ez)
A smooth periodic plasma (`SmoothPlasma` in `tests/common.py`): varying density, a flow with three
components, a transverse magnetic wave on a constant `Bx`, and `T_e = 0.5`, `η = 0.1`, `ν = 0.05`, so
that the ideal, Hall, pressure, resistive and hyper-resistive terms are all non-zero. The exact E is
written by hand; the code gets the exact n, u, B on the Yee grid and computes j then E. For 32 to 1024
cells the largest error must fall as dx² (measured slopes 1.99, 2.00, 2.00).
*Why*: it exercises every term and every average between nodes and cell centres. A term with the wrong
sign leaves an error that does not decrease; a term taken on the wrong set of points gives slope 1.

## Step 3c — Faraday and the time scheme, ions held fixed (`test_step3c_faraday.py`)

`faraday(B, E, grid, dt)` returns `B − dt ∇ × E` (in 1D: `Bx` constant, `∂By/∂t = ∂Ez/∂x`,
`∂Bz/∂t = −∂Ey/∂x`). `advance_fields_with_fixed_ions` does one step of the iterated Crank–Nicolson scheme
of session 3 slides 6–7: two predictions of B from the mid-step E (average of E at the start and E from
the predicted B), then the correction. No particles: n and u are given and constant.
Figure: `plots/plot_step3c_faraday.py`.

### 1. Whistler frequency (`test_whistler_frequency`, wavenumbers 1, 2, 4)
Ions at rest, `n = 1`, `Bx = 1`, a small wave `By = b cos(kx)`, `Bz = b sin(kx)`, Hall term only
(`η = ν = 0`). The pattern must turn at `ω = k² Bx / n`. The frequency is measured over two periods from
the phase of the Fourier mode, on 128 cells: 0.9998, 3.997, 15.95 for 1, 4, 16 expected (within 0.5 %;
the small deficit is the grid, which sees `k` as `(2/dx) sin(k dx/2)`), and the amplitude is kept to 1e-6.
*Why*: Ampère, the Hall term, Faraday and the time scheme working together, on the wave that sets the
stability limit of the code (`dt` below about `n dx² / (2 B)`). `ω ∝ k²` is specific to whistlers (most
waves have `ω ∝ k`): matching it for three wavelengths cannot happen by accident.

### 2. Resistive decay (`test_resistive_decay`)
`Bx = 0` (no Hall effect on B), `η = 0.05`, `ν = 0.001`: Faraday's law becomes a diffusion equation and a
sine wave must decay as `exp(−(η k² + ν k⁴) t)`. Measured rate 0.21581 for 0.21600 expected (0.1 %),
while the amplitude is divided by 4.5.
*Why*: tests the two dissipation terms and the time scheme on a problem that decays instead of
oscillating; a wrong sign would make the wave grow. Note the time step: the hyper-resistive term limits
`dt` to about `dx⁴ / (8 ν)`, much smaller than the whistler limit on a fine grid.
