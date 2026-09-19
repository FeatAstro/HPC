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
- **Weights**: following the slides, `w_p` is already a density contribution, so
  `deposit_moments` computes `n_ij = Σ S·w_p` with *no* division by `dx·dy`; `Σ n_ij = Σ w_p`
  for any spacing. `v_ij` is the weighted mean (divided by the summed weight).
- **Grid indexing**: node `(i, j)` is stored at `j * nx + i` (row-major); `nx`, `ny` count
  nodes, not cells. Always go through `Grid2D::index()`.
- **Toolchain choices**: C++20 (`std::numbers::pi`), assert-based test executables run by
  `ctest` (no test framework), own minimal `Vec3` (Eigen is not installed; kept custom on
  purpose), HDF5 is linked in CMake but unused here.
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

## To do

- Self-consistent field solve on a Yee-staggered grid: Faraday, Ampere, generalized Ohm's
  law, then the full PIC loop of slide 6 (moments → fields → gather → push). `Grid2D` is
  not staggered yet.
- Boundary conditions: particles outside the grid are currently clamped to the edge cell
  and extrapolated (weights can be negative) instead of being wrapped or reflected.
- More Step 2 test cases with known analytic answers once the field solve exists.
- Minor: `main.cpp` uses `std::numbers::pi` without `#include <numbers>` (compiles through
  a transitive include; add it explicitly).

## Documentation

`cpp_concepts.tex` (PDF tracked, aux/log/toc gitignored) explains the physics first, then the
C++ used; rebuild with `latexmk -pdf cpp_concepts.tex` and keep it in sync when the code's
conventions change. The repo root `README.md` describes both sessions.

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
