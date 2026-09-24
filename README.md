# HPC Master Class — Projects

Two C++ sessions (plus an extension of session 2), each a self-contained CMake project (C++, HDF5 available on the system).

| Folder | Topic |
| --- | --- |
| [`session_1/`](session_1) | Finite-difference convergence study |
| [`session_2/`](session_2) | Boris pusher and particle–mesh coupling (PIC, prescribed fields) |
| [`session_2_advanced/`](session_2_advanced) | Extension of session 2: self-consistent hybrid-kinetic PIC (fields solved) |

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

The particle side of a 2D particle-in-cell code, following the course slides in `session_2/session_2.pdf`.
The domain is periodic and the fields are prescribed (not solved yet).

1. **Boris pusher, one particle**: tested against the exact gyration in a uniform magnetic field.
2. **N particles and the mesh**: bilinear deposit of density and bulk velocity, bilinear gather of the
   fields, same shape function both ways; weight conservation for any grid spacing, periodic wrapping.

```sh
cd session_2
cmake -S . -B build && cmake --build build
ctest --test-dir build --output-on-failure     # run the tests
./build/kinetic_fisher                         # demo of the steps above
```

`session_2/cpp_concepts.tex` (PDF included) explains the physics and the C++ of both steps
(build with `latexmk -pdf cpp_concepts.tex`), `tests_walkthrough.pdf` derives every test and
`code_walkthrough.pdf` follows `main()` through each function; `session_2/CLAUDE.md` summarizes the conventions.

## session_2_advanced — Hybrid-kinetic particle-in-cell (extension)

Goes beyond the session: the same Steps 1–2, plus the field equations of the slides, so that the
particles act back on the fields. Ions are kinetic macro-particles, electrons a massless fluid.

- **Yee grid**: `E` and `j` on cell edges, `B` on cell faces, each component interpolated from where it lives.
- **Field solver**: Ampère's law (no displacement current), generalized Ohm's law with
  `v_e = v_i - j/(ne)` and an isothermal electron pressure, Faraday's law advanced with RK4; `div B`
  stays at round-off.
- **Full PIC loop** (predictor-corrector): deposit moments, solve fields, gather, push.
- **Validation**: a uniform drift stays in equilibrium, and parallel ion-cyclotron and whistler waves
  propagate at the frequency of the exact dispersion relation (about 1% agreement, energy conserved to
  1e-5 of the wave energy).

Normalised units (`mu0 = e = m_i = n0 = B0 = 1`); settings live in the structs `PlasmaParameters`,
`IonLoading`, `WaveRunSettings`. Build, test and run as for session_2 (inside `session_2_advanced/`);
its `cpp_concepts.pdf`, `tests_walkthrough.pdf` and `code_walkthrough.pdf` cover the full hybrid code.
