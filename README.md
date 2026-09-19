# HPC Master Class — Projects

Two C++ sessions, each a self-contained CMake project (C++, HDF5 available on the system).

| Folder | Topic |
| --- | --- |
| [`session_1/`](session_1) | Finite-difference convergence study |
| [`session_2/`](session_2) | Boris pusher and particle–mesh coupling (hybrid-kinetic PIC) |

## session_1 — Finite-difference convergence

Approximates the derivative of `sin(x)` at `x = 1` with a first-order (forward) and a
second-order (centred) finite-difference scheme, sweeping the step `h` from `1e-1` down to
`1e-12`. The absolute error against the exact derivative `cos(x)` is written to
`build/results.h5`, and `analysis.ipynb` plots it on a log–log scale to show the truncation-error
slopes and the round-off floor at small `h`.

```sh
cd session_1
cmake -S . -B build && cmake --build build
cd build && ./finite_difference_convergence    # writes results.h5 in the current folder
```

Then open `analysis.ipynb` (it reads `build/results.h5`).

## session_2 — Boris pusher and particle–mesh coupling

Building blocks of a 2D hybrid-kinetic particle-in-cell code (ions as macro-particles, electrons as a
fluid), following the course slides in `session_2/session_2.pdf`:

1. **Boris pusher, one particle**: a charged particle advanced in prescribed electric and magnetic
   fields, tested against the exact gyration in a uniform magnetic field (speed conservation, Larmor
   radius, return to start after one period).
2. **N particles and the mesh**: particles deposit weight and bulk velocity on a 2D grid with a
   first-order (bilinear) shape function, and grid fields are gathered back to the particles with the
   same shape function. Tested for weight conservation (any grid spacing), exactness on linear fields,
   and consistency with step 1.

Not implemented yet: the self-consistent field solve (Faraday, Ampère, generalized Ohm's law) and the
Yee-staggered grid. Fields on the grid are prescribed.

```sh
cd session_2
cmake -S . -B build && cmake --build build
ctest --test-dir build --output-on-failure     # run the tests
./build/kinetic_fisher                         # demo of both steps
```

`session_2/cpp_concepts.tex` explains the physics and the C++ used (build with
`latexmk -pdf cpp_concepts.tex`); `session_2/CLAUDE.md` summarizes the conventions.
