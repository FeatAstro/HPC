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

## session_2 — Hybrid-kinetic particle-in-cell

A 2D hybrid-kinetic particle-in-cell code (ions as macro-particles, electrons as a massless fluid),
following the course slides in `session_2/session_2.pdf`. The domain is periodic, in normalised units
(`mu0 = e = m_i = n0 = B0 = 1`); every physical and numerical setting is a field of one parameter struct
(`PlasmaParameters`, `IonLoading`, `WaveRunSettings`).

1. **Boris pusher, one particle**: tested against the exact gyration in a uniform magnetic field.
2. **N particles and the mesh**: bilinear deposit of density and bulk velocity, bilinear gather of the
   fields, same shape function both ways; weight conservation for any grid spacing.
3. **Periodic boundaries and a Yee grid**: `E` and `j` on cell edges, `B` on cell faces, each component
   interpolated from where it lives.
4. **Field solver**: Ampère's law (no displacement current), generalized Ohm's law with the electron
   fluid velocity `v_e = v_i - j/(ne)` and an isothermal electron pressure, Faraday's law advanced with
   RK4. `div B` stays at round-off.
5. **Full PIC loop** (predictor-corrector): moments, fields, gather, push.
6. **Validation**: a uniform drift stays in equilibrium, and parallel ion-cyclotron and whistler waves
   propagate at the frequency of the exact dispersion relation (about 1% agreement, energy conserved to
   1e-5 of the wave energy).

```sh
cd session_2
cmake -S . -B build && cmake --build build
ctest --test-dir build --output-on-failure     # run the tests
./build/kinetic_fisher                         # demo of the steps above
```

`session_2/cpp_concepts.tex` (PDF included) explains the physics and the C++ of all the steps above
(build with `latexmk -pdf cpp_concepts.tex`); `session_2/CLAUDE.md` summarizes the conventions.
