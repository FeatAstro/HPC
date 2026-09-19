# session_2 — Hybrid-Kinetic PIC (Boris Pusher + Particle-Mesh Coupling)

Source of the physics: `session_2.pdf` (professor's slides, "HPC Master Class").
This is a 2D hybrid-kinetic PIC code: ions are kinetic macro-particles (Boris-pushed),
electrons are a fluid; the mesh carries deposited moments and (eventually) fields.

## Build & test

```sh
cmake -S . -B build
cmake --build build
ctest --test-dir build --output-on-failure
./build/kinetic_fisher   # demo driver, prints Step 1 & Step 2 diagnostics
```

## Conventions

- **2D3V**: particle positions are 2D `(x, y)`; velocities and fields are full
  3-component `Vec3`, since `v × B` needs all three components even though the
  space grid is 2D.
- **Shape function**: first-order = bilinear. `interpolation.cpp` implements
  both deposit (particle → grid moments `n_ij`, `v_ij`) and gather (grid → particle
  field) through the same `bilinear_weights()` function — this symmetry is required
  to avoid a self-force artifact (see slide 6).
- **Boris pusher** (`boris_pusher.cpp`) implements the exact update from slide 7:
  half `E` accel → exact `B` rotation (`t`, `v'`, `s` trick) → half `E` accel.
  It only rotates/accelerates velocity; callers (`step1_single_particle.cpp`,
  `step2_pic_loop.cpp`) own the `r^{n+1/2}` / `r^{n+1}` position half-steps, since
  they also decide where `E`,`B` are sampled from.

## Scope: what's implemented vs. deferred

Implemented:
- Step 1: single particle, Boris pusher, prescribed analytic `E`,`B` fields.
- Step 2: N particles, bilinear deposit of moments onto a node-centered grid,
  bilinear gather of (prescribed/static) fields back to particles.

**Not implemented yet** (deferred to a later project stage): the self-consistent
field solve — Faraday's law for `B`, Ampere's law for `j`, and the generalized
Ohm's law for `E` (slide 1) — and the Yee-staggered grid layout that solve would
need. Right now `Grid2D` is node-centered (not staggered) and `grid.E`/`grid.B`
are just set directly by the caller; nothing evolves them.

## Layout

- `src/vec3.hpp` — minimal 3-vector.
- `src/particle.hpp` — `Particle` (x, y, v, q, m, w).
- `src/boris_pusher.{hpp,cpp}` — velocity update only.
- `src/step1_single_particle.{hpp,cpp}` — leapfrog wrapper for prescribed-field particle.
- `src/grid.hpp` — `Grid2D`.
- `src/interpolation.{hpp,cpp}` — `bilinear_weights`, `deposit_moments`, `gather`.
- `src/step2_pic_loop.{hpp,cpp}` — leapfrog wrapper using grid-gathered fields.
- `tests/test_boris_pusher.cpp` — Step 1 checks (energy conservation, gyroradius, period).
- `tests/test_deposit_gather.cpp` — Step 2 checks (weight conservation, linear-field
  exactness, consistency with Step 1).
