# session_2 — Boris Pusher + Particle-Mesh Coupling

Source of the physics: `session_2.pdf` (professor's slides, "HPC Master Class").
This is the particle side of a 2D PIC code: ions are macro-particles (Boris-pushed); the mesh carries
the deposited moments and the prescribed (not solved) fields.

The self-consistent hybrid extension (Yee grid, Ampère + Ohm + Faraday, predictor-corrector loop,
wave validation) lives in `../session_2_advanced/`, which is gitignored and not part of the submission.

## Build & test

```sh
cmake -S . -B build
cmake --build build
ctest --test-dir build --output-on-failure
./build/kinetic_fisher   # demo driver: Step 1 and Step 2
```

## Conventions

- **2D3V, periodic**: particle positions are 2D `(x, y)`; velocities and fields are full
  3-component `Vec3` (`v × B` needs all three). The domain is periodic: `Grid2D::index()` wraps
  indices and `wrap_x/wrap_y` wrap positions, so `nx`, `ny` count nodes = cells and the domain
  length is `nx * dx`.
- **Units**: dimensionless, `q = m = 1`, `B` in `B0`, time in `1/omega_c`.
- **Naming**: variables are spelled out (`ion_density`, `electric_field`, `magnetic_field`),
  comments only where the code cannot say it (physics, sign, layout).
- **Shape function**: first-order = bilinear. `interpolation.cpp` implements deposit and gather
  through the same `bilinear_weights()` and `for_each_corner()` — this symmetry is required to avoid
  a self-force artifact (slide 6).
- **Weights**: following the slides, `w_p` is already a density contribution, so
  `deposit_moments` computes `n_ij = Σ S·w_p` with *no* division by `dx·dy`; `Σ n_ij = Σ w_p`
  for any spacing. `v_ij` is the weighted mean.
- **Grid layout**: moments and the prescribed `E`, `B` all live on the nodes.
- **Toolchain**: C++20, assert-style test executables run by `ctest` (helpers in `tests/check.hpp`),
  own minimal `Vec3` (Eigen not installed), HDF5 linked in CMake but unused.
- **Boris pusher** (`boris_pusher.cpp`): slide 7 velocity update only; `particle_push.cpp` owns the
  drift half-steps and the E/B gather.

## Scope: what is implemented

- Step 1: single particle, Boris pusher, prescribed analytic fields.
- Step 2: N particles, bilinear deposit and gather, periodic boundaries, fields prescribed on the grid.

## To do

- Field equations (Ampère, generalized Ohm's law, Faraday) and the self-consistent PIC loop
  (done in `../session_2_advanced/`).
- Higher-order shape functions (slide 5 mentions B-splines).
- Performance: OpenMP over particles needs a safe deposit (per-thread grids or colouring); HDF5 output.

## Documentation

`cpp_concepts.tex` (PDF tracked, aux/log/toc gitignored) explains the physics first, then the
C++ used; rebuild with `latexmk -pdf cpp_concepts.tex` and keep it in sync when the code's
conventions change. `tests_walkthrough.tex` (PDF tracked) derives every test step by step; update it
when a test or its tolerance changes. `code_walkthrough.tex` (PDF tracked) follows `main()` down the
call tree and explains each function physically; update it when the call structure or a demo changes.
The repo root `README.md` describes both sessions.

## Layout

- `src/vec3.hpp` — 3-vector and its operators.
- `src/particle.hpp` — `Particle` (x, y, v, q, m, w).
- `src/vector_field.hpp` — `ScalarField`/`VectorField` aliases.
- `src/grid.hpp` — periodic `Grid2D` with the moments and fields.
- `src/interpolation.{hpp,cpp}` — `bilinear_weights`, `deposit_moments`, `gather`.
- `src/boris_pusher.{hpp,cpp}` — velocity update only.
- `src/particle_push.{hpp,cpp}` — `drift_particles`, `accelerate_particles`.
- `src/step1_single_particle.{hpp,cpp}`, `src/step2_pic_loop.{hpp,cpp}` — leapfrog wrappers.
- `tests/` — `test_boris_pusher`, `test_deposit_gather`, `test_periodic_grid`.
