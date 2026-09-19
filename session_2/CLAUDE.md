# session_2 — Hybrid-Kinetic PIC (Boris Pusher + Particle-Mesh Coupling)

Source of the physics: `session_2.pdf` (professor's slides, "HPC Master Class").
This is a 2D hybrid-kinetic PIC code: ions are kinetic macro-particles (Boris-pushed),
electrons are a fluid; the mesh carries deposited moments and the Yee-staggered fields.

## Build & test

```sh
cmake -S . -B build
cmake --build build
ctest --test-dir build --output-on-failure
./build/kinetic_fisher   # demo driver: Boris, deposit/gather, hybrid waves
```

## Conventions

- **2D3V, periodic**: particle positions are 2D `(x, y)`; velocities and fields are full
  3-component `Vec3` (`v × B` needs all three). The domain is periodic: `Grid2D::index()` wraps
  indices and `wrap_x/wrap_y` wrap positions, so `nx`, `ny` count nodes = cells and the domain
  length is `nx * dx`.
- **Normalised units**: defaults of `PlasmaParameters` are `mu0 = e = m_i = n0 = 1`, `B` in `B0`,
  velocity in `v_A`, time in `1/Omega_i`. Settings live in structs (`PlasmaParameters`, `IonLoading`,
  `WaveRunSettings`); change a value there, not in the algorithms.
- **Naming**: variables are spelled out (`ion_density`, `electric_field`, `magnetic_field`,
  `current_density`), comments only where the code cannot say it (physics, sign, layout).
- **Shape function**: first-order = bilinear. `interpolation.cpp` implements deposit and gather
  through the same `bilinear_weights()`, which takes the staggering offset of the quantity —
  this symmetry is required to avoid a self-force artifact (slide 6).
- **Weights**: following the slides, `w_p` is already a density contribution, so
  `deposit_moments` computes `n_ij = Σ S·w_p` with *no* division by `dx·dy`; `Σ n_ij = Σ w_p`
  for any spacing. `v_ij` is the weighted mean. Uniform loading uses `w_p = n0 / particles_per_cell`.
  A macro-particle stands for `w_p·dx·dy` real ions (used for energies).
- **Yee layout** (`yee_layout.hpp`, the single place that defines it): `Ex (i+½,j)`, `Ey (i,j+½)`,
  `Ez (i,j)`; `Bx (i,j+½)`, `By (i+½,j)`, `Bz (i+½,j+½)`; `j` like `E`; ion moments on the nodes.
  Curls use forward differences (E → B) and backward differences (B → E/j), see
  `finite_differences.hpp` and `field_solver.cpp`.
- **Hybrid model**: ions are kinetic (Boris-pushed); electrons are a massless, quasi-neutral,
  isothermal fluid that is never stored: `v_e = v_i − j/(ne)`, `E = −v_e×B − ∇P_e/(ne)`,
  `mu0 j = ∇×B`. Ohm's law is evaluated at each E-component location by interpolating `n`,
  `v_i`, `j`, `B` there. The closure is one function (`electron_pressure` in `field_solver.cpp`);
  electron inertia is dropped.
- **Time stepping** (`pic_loop.cpp`): drift half step, then a predictor (moments with `v^n`,
  push a copy to estimate `v^{n+1}`), then a corrector (moments with `(v^n+v^{n+1})/2`, fields,
  final Boris push, drift half step). `B` is advanced with RK4 and frozen moments; `E` is computed
  from the half-step `B`. RK4 is used because a midpoint scheme is weakly unstable for whistlers.
- **Toolchain**: C++20, assert-style test executables run by `ctest` (helpers in `tests/check.hpp`),
  own minimal `Vec3` (Eigen not installed), HDF5 linked in CMake but unused.
- **Boris pusher** (`boris_pusher.cpp`): slide 7 velocity update only; `particle_push.cpp` owns the
  drift half-steps and the E/B gather.

## Scope: what is implemented

- Step 1: single particle, Boris pusher, prescribed analytic fields.
- Step 2: N particles, bilinear deposit and gather.
- Periodic boundaries, Yee-staggered fields, Ampère + generalized Ohm's law + Faraday (RK4),
  the full predictor-corrector PIC loop (`advance_hybrid_step`).
- Validation: uniform-drift equilibrium, parallel ion-cyclotron / whistler waves against the exact
  dispersion relation (`linear_wave.cpp`, ~1% frequency error, energy conserved to ~1e-5).

Known behaviour: a warm, randomly loaded plasma shows particle noise (transverse field grows to a
level that scales as `1/sqrt(particles per cell)`, total energy drifts up accordingly); `div B` stays
at 1e-14. With `electron_temperature > 0` total ion + magnetic energy is not conserved (the isothermal
electrons exchange energy with a heat bath).

## To do

- Boundary conditions other than periodic (reflecting/open) if a problem needs them.
- Higher-order shape functions (slide 5 mentions B-splines) and noise reduction (`delta-f`, filtering).
- More validation: oblique/compressive waves (magnetosonic), Landau-type damping, ion beam instabilities.
- Performance: the Ohm's-law evaluation re-computes bilinear weights per component; OpenMP over
  particles needs a safe deposit (per-thread grids or colouring); HDF5 output.
- Update `cpp_concepts.tex` for the periodic grid, Yee layout, field solver and PIC loop (it still
  describes the old node-centred `Grid2D` field names).

## Documentation

`cpp_concepts.tex` (PDF tracked, aux/log/toc gitignored) explains the physics first, then the
C++ used; rebuild with `latexmk -pdf cpp_concepts.tex` and keep it in sync when the code's
conventions change. The repo root `README.md` describes both sessions.

## Layout

- `src/vec3.hpp` — 3-vector, `axis::x/y/z` component indices.
- `src/particle.hpp` — `Particle` (x, y, v, q, m, w).
- `src/vector_field.hpp` — `ScalarField`/`VectorField` and small arithmetic helpers.
- `src/yee_layout.hpp` — where each component lives (`yee::electric_field`, ...).
- `src/grid.hpp` — periodic `Grid2D` with the moments and fields.
- `src/finite_differences.hpp` — forward/backward differences on nodes.
- `src/interpolation.{hpp,cpp}` — `bilinear_weights`, `deposit_moments`, `gather`, `gather_scalar`.
- `src/boris_pusher.{hpp,cpp}` — velocity update only.
- `src/particle_push.{hpp,cpp}` — `drift_particles`, `accelerate_particles`.
- `src/step1_single_particle.{hpp,cpp}`, `src/step2_pic_loop.{hpp,cpp}` — prescribed-field leapfrog wrappers.
- `src/plasma_parameters.hpp` — physical constants and plasma properties.
- `src/field_solver.{hpp,cpp}` — curls, Ampère, Ohm's law, Faraday (RK4).
- `src/pic_loop.{hpp,cpp}` — `advance_hybrid_step`.
- `src/particle_loading.{hpp,cpp}` — uniform lattice loading of ions.
- `src/diagnostics.{hpp,cpp}` — energies, `max |div B|`.
- `src/linear_wave.{hpp,cpp}` — parallel-wave dispersion relation and run driver.
- `tests/` — `test_boris_pusher`, `test_deposit_gather`, `test_periodic_grid`, `test_yee_grid`,
  `test_field_solver`, `test_hybrid_waves`.
