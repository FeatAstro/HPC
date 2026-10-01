# HPC Master Class — Projects

Two C++ projects from the course, each a self-contained CMake project.

| Folder | Topic |
| --- | --- |
| [`FDC/`](FDC) | Finite-difference convergence study |
| [`PIC/`](PIC) | 1D hybrid particle-in-cell code, built and tested step by step |

The course slides and LaTeX notes are kept locally and are not in the repository.

## FDC — Finite-difference convergence

Approximates the derivative of `sin(x)` at `x = 1` with a first-order (forward) and a
second-order (centred) finite-difference scheme, sweeping the step `h` from `1e-1` down to
`1e-12`. The absolute error against the exact derivative `cos(x)` is written to
`build/results.h5`, and `analysis.ipynb` plots it on a log–log scale to show the truncation-error
slopes and the round-off floor at small `h`.

```sh
cd FDC
cmake -S . -B build && cmake --build build
cd build && ./finite_difference_convergence    # writes results.h5 in the current folder
```

Then open `analysis.ipynb` (it reads `build/results.h5`). Needs HDF5.

## PIC — Hybrid particle-in-cell code

A 1D hybrid-kinetic plasma code: the ions are macro-particles, the electrons a massless isothermal
fluid. Positions are 1D and velocities 3D (1D3V), the domain is periodic, and the units are
normalised (`μ0 = e = m_i = 1`). The physics is in C++; it is exposed to Python with pybind11, and
every test and figure is written in Python.

### What is implemented

| Step | Content | Checked by |
| --- | --- | --- |
| 1 | Boris pusher, one particle in uniform fields | Larmor radius against mass, rotation direction, error ∝ dt², E×B drift |
| 2a | A population of particles (one array per quantity), periodic domain | energy and temperature conserved in a magnetic field, periodic wrap |
| 2b | Periodic grid, order-1 deposit of the moments and gather of the fields | hat-function deposit, linear gather, conservation of weight, momentum and energy |
| 2c | Maxwellian loading from user profiles n(x), u(x), T(x) | error of the deposited moments ∝ 1/√N |
| 3a | Yee layout: `Bx, Ey, Ez` on the nodes, `By, Bz, Ex` at the cell centres | derivatives and averages second order in dx |
| 3b | Ampère's law and the generalised Ohm's law (Hall, electron pressure, resistivity, hyper-resistivity) | uniform plasma gives `E = −u × B`; error ∝ dx² against an exact solution |
| 3c | Faraday's law with the iterated Crank–Nicolson scheme, ions held fixed | whistler frequency `ω = k² B / n`, resistive decay rate |

Not done yet: the full loop, where the particles and the fields advance together.

### Build, test, plot

Needs CMake, a C++17 compiler, Python 3 with numpy and matplotlib, and pybind11 and ddt
(on Debian/Ubuntu: `sudo apt install pybind11-dev python3-pybind11 python3-ddt`).

```sh
cd PIC
cmake -S . -B build && cmake --build build
python3 -m unittest discover -s tests -v       # or: ctest --test-dir build --output-on-failure
python3 plots/plot_step1_boris.py              # one script per step, figures go to plots/figures/
```

### Layout

- `src/` — the C++ code: `boris` (pusher), `population` (particles, deposit, push), `grid` (Yee grid,
  shape function), `loading` (Maxwellian loader), `field_solver` (Ampère, Ohm, Faraday).
- `python/bindings.cpp` — the Python module `pic`.
- `tests/` — the tests, with their shared helpers in `common.py`.
- `plots/` — one figure script per step: what the code gives against what is expected.
- [`TESTS.md`](PIC/TESTS.md) — what each test checks and why.
